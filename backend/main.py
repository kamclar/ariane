# ============================================================
# ARIANE - FastAPI application
# Automated ACMG Rule-based Interpretation and Annotation ENgine
# ============================================================
from fastapi import Depends, FastAPI, HTTPException, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.encoders import jsonable_encoder
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, JSONResponse, Response
from pydantic import ValidationError
from pathlib import Path
from datetime import datetime, timezone
import json
import hashlib
import logging
import os
import secrets
import sys
import time
import asyncio
import uuid
from typing import Optional
from fastapi import Header
from backend.admin import router as admin_router
from backend.review_api import router as review_router
from backend.version import ARIANE_VERSION
from backend.api_auth import require_public_api_key
from backend.ui_session import issue_ui_session, require_ui_session
from backend.classification_runtime import (
    ClassificationCacheRepository,
    ClassificationUsageRepository,
    classification_fingerprint,
    public_api_daily_classification_limit,
    resolve_usage_identity,
)

from backend.config import (
    TABLE4_PATH, TABLE9_PATH, ST7_PATH,
    ENIGMA_REFERENCE_TABLES_PATH, ENIGMA_RULE_CATALOG_PATH,
    PS1_PROTEIN_REGISTRY_PATH,
    ST2_SPLICE_EVIDENCE_PATH, EXON_CNV_EVIDENCE_PATH,
    EXON_CNV_EVIDENCE_MANIFEST_PATH,
    GENE_POLICY_MANIFEST_PATH, GENE_POLICY_METADATA_PATH,
)
from backend.gene_policy import (
    active_genes,
    get_gene_policy,
    not_used_criteria,
    validate_policy_source_bindings,
)
from backend.data_validation import validate_required_datasets
from backend.data_health import get_data_issues
from backend.lookups.spliceai import spliceai_runtime_health
from backend.classification_dag import (
    DagNodeExecutionError,
    get_configured_engine_mode,
)
from backend.classification_dag.manual import execute_manual_evidence
from backend.classification_dag.provider_wiring import production_provider_dependencies
from backend.models import (
    VariantRequest, ClassificationResult,
    BatchRequest, BatchResponse, BatchItemResult,
    ManualEvidenceRequest, ManualEvidenceResult,
    ManualEvidenceStatusRequest, ManualEvidenceStatusResponse,
    ManualCriterionStatus, ManualCriterionResult, EvidenceInteractionWarning,
    ClientValidationRequest, VariantNormalizationResponse,
    Ps1ReferenceResolutionRequest, Ps1ReferenceResolutionResponse,
)
from backend.public_api import (
    PUBLIC_API_VERSION,
    PublicApiError,
    PublicApiErrorResponse,
    build_metadata,
    create_public_api_router,
)
from backend.services import (
    ClassificationCommand,
    EvidenceOrchestrationService,
    EvidenceExecutionError,
    RequiredEvidenceUnavailableError,
    VariantPreparationError,
    execute_variant_classification,
    resolve_ps1_reference,
)
from backend.modules.hgvs_provider import load_panel_provider
from backend.modules.manual_evidence import (
    manual_criterion_statuses,
    manual_criteria_for_gene,
    resource_links_for_gene,
)
from backend.modules.ps1_splice_evidence import list_splice_ps1_candidate_discovery

validate_required_datasets({
    "table4": TABLE4_PATH,
    "table9": TABLE9_PATH,
    "enigma_rule_catalog": ENIGMA_RULE_CATALOG_PATH,
    "enigma_reference_tables": ENIGMA_REFERENCE_TABLES_PATH,
    "st7": ST7_PATH,
    "ps1_protein_registry": PS1_PROTEIN_REGISTRY_PATH,
    "st2_splice_evidence": ST2_SPLICE_EVIDENCE_PATH,
    "exon_cnv_evidence": EXON_CNV_EVIDENCE_PATH,
    "exon_cnv_evidence_manifest": EXON_CNV_EVIDENCE_MANIFEST_PATH,
    "gene_policy_manifest": GENE_POLICY_MANIFEST_PATH,
    "gene_policy_metadata": GENE_POLICY_METADATA_PATH,
})
validate_policy_source_bindings()

# Validate once at process startup. The engine cannot change underneath a
# running classification process.
CLASSIFIER_ENGINE_MODE = get_configured_engine_mode()

# Initialize local sources before serving requests so /api/health reports
# degraded caches even before the first classification.
from backend.population_frequency import PopulationFrequencyService  # noqa: E402
from backend.lookups import coordinates as _coordinate_data_source  # noqa: E402,F401
from backend.lookups import bayesdel as _bayesdel_data_source  # noqa: E402,F401
from backend.lookups import spliceai as _spliceai_data_source  # noqa: E402
from backend.modules.pp4_bp5 import load_pp4_bp5_snapshot  # noqa: E402
from backend.modules.residues import initialize_residue_data  # noqa: E402
from backend.modules.hgvs_engine import validate_hgvs_engine  # noqa: E402
from backend.modules.enigma_rules import (  # noqa: E402
    get_decision_tree,
    public_catalog,
    search_reference_table,
    search_table9,
    validate_rule_catalog,
)
from backend.cache_registry import clear_runtime_caches  # noqa: E402

_spliceai_data_source._load_api_cache()
load_pp4_bp5_snapshot()
initialize_residue_data()
validate_hgvs_engine()
validate_rule_catalog()

# Population snapshots are loaded explicitly and owned by the application.
# Provider dependencies retain the service method, so an administrative reload
# atomically replaces the repository used by subsequent requests.
POPULATION_FREQUENCY_SERVICE = PopulationFrequencyService.load_default()
CLASSIFICATION_ORCHESTRATION = EvidenceOrchestrationService(
    engine_mode=CLASSIFIER_ENGINE_MODE,
    provider_dependencies=production_provider_dependencies(
        population_frequency_lookup=POPULATION_FREQUENCY_SERVICE.get_frequencies,
    ),
)
CLASSIFICATION_CACHE: ClassificationCacheRepository | None = None
CLASSIFICATION_USAGE: ClassificationUsageRepository | None = None
PUBLIC_API_QUOTA: ClassificationUsageRepository | None = None
PUBLIC_API_DAILY_CLASSIFICATION_LIMIT = public_api_daily_classification_limit()

# ── App setup ──────────────────────────────────────────────────────────────
app = FastAPI(
    title="ARIANE",
    description="Automated ACMG Rule-based Interpretation and Annotation ENgine",
    version=ARIANE_VERSION,
)
app.include_router(admin_router)
app.include_router(review_router)

AUDIT_LOGGER = logging.getLogger("ariane.audit")
AUDIT_LOGGER.setLevel(logging.INFO)
AUDIT_LOGGER.propagate = False
if not AUDIT_LOGGER.handlers:
    audit_handler = logging.StreamHandler(sys.stdout)
    audit_handler.setFormatter(logging.Formatter("%(message)s"))
    AUDIT_LOGGER.addHandler(audit_handler)

AUDIT_LOG_PATH = Path(os.getenv("ARIANE_AUDIT_LOG", "/var/log/ariane/audit.jsonl"))
try:
    if AUDIT_LOG_PATH.parent.is_dir():
        audit_file_handler = logging.FileHandler(AUDIT_LOG_PATH, encoding="utf-8")
        audit_file_handler.setFormatter(logging.Formatter("%(message)s"))
        AUDIT_LOGGER.addHandler(audit_file_handler)
except OSError:
    AUDIT_LOGGER.exception("Failed to open the audit log file")

try:
    CLASSIFICATION_CACHE = ClassificationCacheRepository()
except Exception:
    logging.getLogger("ariane.classification_cache").exception(
        "Classification cache could not be initialized"
    )
try:
    CLASSIFICATION_USAGE = ClassificationUsageRepository()
    PUBLIC_API_QUOTA = CLASSIFICATION_USAGE
except Exception:
    logging.getLogger("ariane.classification_usage").exception(
        "Classification usage storage could not be initialized"
    )

def _request_context(request: Request) -> dict:
    return {
        "request_id": getattr(request.state, "request_id", ""),
        "source_ip": request.client.host if request.client else "unknown",
        "method": request.method,
        "path": request.url.path,
        "user_agent": request.headers.get("user-agent", "")[:300],
    }


def _audit(request: Request, event: str, level: str = "info", **fields) -> None:
    record = {
        "log_type": "ariane_audit",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "event": event,
        **_request_context(request),
        **fields,
    }
    message = json.dumps(record, ensure_ascii=True, separators=(",", ":"))
    getattr(AUDIT_LOGGER, level)(message)


def _request_id(value: str | None) -> str:
    """Accept a compact log-safe client ID or generate a server ID."""
    if value:
        candidate = value.strip()
        if (
            1 <= len(candidate) <= 128
            and all(char.isalnum() or char in "-._:" for char in candidate)
        ):
            return candidate
    return uuid.uuid4().hex


@app.middleware("http")
async def audit_request(request: Request, call_next):
    request.state.request_id = _request_id(request.headers.get("x-request-id"))
    started = time.monotonic()
    try:
        response = await call_next(request)
    except Exception as exc:
        _audit(
            request,
            "request_exception",
            level="exception",
            duration_ms=round((time.monotonic() - started) * 1000, 1),
            error_type=type(exc).__name__,
            error=str(exc)[:2000],
        )
        raise
    response.headers["X-Request-ID"] = request.state.request_id
    if request.url.path.startswith("/api/v1/"):
        response.headers["X-ARIANE-API-Version"] = PUBLIC_API_VERSION
        response.headers["X-ARIANE-Version"] = ARIANE_VERSION
    log_completion = (
        request.url.path.startswith("/admin/")
        or request.url.path.startswith("/ui-api/")
        or request.url.path.startswith("/api/") and request.url.path != "/api/health"
    )
    if log_completion:
        _audit(
            request,
            "request_completed",
            level="warning" if response.status_code >= 400 else "info",
            status_code=response.status_code,
            duration_ms=round((time.monotonic() - started) * 1000, 1),
        )
    return response


@app.exception_handler(RequestValidationError)
async def validation_error_handler(request: Request, exc: RequestValidationError):
    _audit(
        request,
        "validation_error",
        level="warning",
        input=jsonable_encoder(exc.body),
        errors=jsonable_encoder(exc.errors()),
    )
    if request.url.path.startswith("/api/v1/"):
        messages = []
        for error in exc.errors():
            location = ".".join(str(value) for value in error.get("loc", ()))
            message = str(error.get("msg", "Invalid request"))
            if message.startswith("Value error, "):
                message = message[len("Value error, "):]
            rendered = f"{location}: {message}" if location else message
            if rendered not in messages:
                messages.append(rendered)
        response = PublicApiErrorResponse(
            metadata=build_metadata(
                request_id=request.state.request_id,
                engine=CLASSIFIER_ENGINE_MODE.value,
            ),
            error=PublicApiError(
                code="invalid_request",
                message="; ".join(messages) or "Invalid request",
                retryable=False,
            ),
        )
        return JSONResponse(
            status_code=422,
            content=response.model_dump(mode="json", exclude_none=True),
        )
    return JSONResponse(
        status_code=422,
        content={"detail": jsonable_encoder(exc.errors())},
    )


@app.exception_handler(HTTPException)
async def http_error_handler(request: Request, exc: HTTPException):
    if not request.url.path.startswith("/api/v1/"):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": jsonable_encoder(exc.detail)},
            headers=exc.headers,
        )
    code_by_status = {
        401: "authentication_required",
        404: "not_found",
        422: "variant_not_classifiable",
        429: "rate_limited",
        503: "evidence_unavailable",
    }
    error_headers = exc.headers or {}
    if isinstance(exc.detail, dict):
        error_code = str(exc.detail.get("code") or code_by_status.get(
            exc.status_code, "request_failed"
        ))
        error_message = str(exc.detail.get("message") or "Request failed")
        retryable = bool(exc.detail.get("retryable", False))
    else:
        error_code = error_headers.get(
            "X-ARIANE-Error-Code",
            code_by_status.get(exc.status_code, "request_failed"),
        )
        error_message = str(exc.detail)
        retryable_header = error_headers.get("X-ARIANE-Retryable")
        retryable = (
            retryable_header.lower() == "true"
            if retryable_header is not None
            else exc.status_code in {429, 502, 503, 504}
        )
    response = PublicApiErrorResponse(
        metadata=build_metadata(
            request_id=request.state.request_id,
            engine=CLASSIFIER_ENGINE_MODE.value,
        ),
        error=PublicApiError(
            code=error_code,
            message=error_message,
            retryable=retryable,
        ),
    )
    return JSONResponse(
        status_code=exc.status_code,
        content=response.model_dump(mode="json", exclude_none=True),
        headers=exc.headers,
    )


FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"
INDEX_HTML_TEMPLATE = (FRONTEND_DIR / "index.html").read_text(encoding="utf-8")


def _frontend_asset_version() -> str:
    """Bind browser asset URLs to the release and their actual contents."""
    static_dir = FRONTEND_DIR / "static"
    asset_paths = [
        static_dir / "css" / "style.css",
        static_dir / "images" / "ariane-icon-v1.svg",
        *sorted((static_dir / "js").glob("*.js")),
    ]
    digest = hashlib.sha256()
    for path in asset_paths:
        digest.update(path.relative_to(static_dir).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
    return f"{ARIANE_VERSION}-{digest.hexdigest()[:12]}"


FRONTEND_ASSET_VERSION = _frontend_asset_version()
app.mount("/static", StaticFiles(directory=FRONTEND_DIR / "static"), name="static")


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    response = HTMLResponse(
        INDEX_HTML_TEMPLATE.replace(
            "__ARIANE_ASSET_VERSION__",
            FRONTEND_ASSET_VERSION,
        ),
        headers={
            "Cache-Control": "no-store, max-age=0",
            "Pragma": "no-cache",
        },
    )
    issue_ui_session(request, response)
    return response


@app.get("/api/health")
async def health():
    issues = get_data_issues()
    panel = load_panel_provider()
    spliceai = spliceai_runtime_health()
    return {
        "status": (
            "degraded"
            if issues or (spliceai["local"] and spliceai["status"] != "ok")
            else "ok"
        ),
        "version": ARIANE_VERSION,
        "classification_engine": CLASSIFIER_ENGINE_MODE.value,
        "data": {
            "table4": TABLE4_PATH.exists(),
            "table9": TABLE9_PATH.exists(),
            "enigma_rule_catalog": ENIGMA_RULE_CATALOG_PATH.exists(),
            "st7":    ST7_PATH.exists(),
            "ps1_protein_registry": PS1_PROTEIN_REGISTRY_PATH.exists(),
            "st2_splice_evidence": ST2_SPLICE_EVIDENCE_PATH.exists(),
            "exon_cnv_evidence": EXON_CNV_EVIDENCE_PATH.exists(),
            "reference_bundle": panel.provenance.get("reference_bundle", ""),
            "normalization_engine": panel.provenance.get("normalization_engine", ""),
            "spliceai": spliceai,
        },
        "data_issues": issues,
    }


@app.get(
    "/ui-api/resources",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.get("/api/resources", dependencies=[Depends(require_public_api_key)])
async def resources(gene: Optional[str] = None):
    return {
        "version": ARIANE_VERSION,
        "build_revision": os.getenv("ARIANE_BUILD_REVISION", "").strip(),
        "issue_tracker_url": os.getenv("ARIANE_ISSUE_TRACKER_URL", "").strip(),
        "manual_criteria": manual_criteria_for_gene(gene) if gene else {},
        "genes": [
            {
                "symbol": symbol,
                "reference_transcript": get_gene_policy(symbol)["gene_config"]["reference_transcript"],
                "reference_protein": get_gene_policy(symbol)["gene_config"]["reference_protein"],
                "policy_id": get_gene_policy(symbol)["policy"]["runtime_policy_id"],
                "policy_name": get_gene_policy(symbol)["policy"]["name"],
                "policy_version": get_gene_policy(symbol)["policy"]["version"],
                "policy_source_url": get_gene_policy(symbol)["policy"]["source_url"],
                "not_used_criteria": not_used_criteria(symbol),
            }
            for symbol in active_genes()
        ],
        "links": resource_links_for_gene(gene),
        "splice_ps1_candidates": list_splice_ps1_candidate_discovery(),
    }


@app.get(
    "/ui-api/rules",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.get("/api/rules", dependencies=[Depends(require_public_api_key)])
async def enigma_rules_catalog():
    """Return the validated source and rule index without local file paths."""
    return public_catalog()


@app.get(
    "/ui-api/rules/trees/{tree_id}",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.get(
    "/api/rules/trees/{tree_id}",
    dependencies=[Depends(require_public_api_key)],
)
async def enigma_decision_tree(tree_id: str):
    tree = get_decision_tree(tree_id)
    if tree is None:
        raise HTTPException(status_code=404, detail="Decision tree not found")
    return tree


@app.get(
    "/ui-api/rules/tables/table9",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.get(
    "/api/rules/tables/table9",
    dependencies=[Depends(require_public_api_key)],
)
async def enigma_table9_records(
    gene: Optional[str] = Query(default=None, max_length=40),
    query: str = Query(default="", max_length=200),
    code: Optional[str] = Query(default=None, max_length=20),
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=25, ge=1, le=100),
):
    if gene is not None:
        gene = gene.strip().upper()
        if gene not in set(active_genes()):
            raise HTTPException(
                status_code=422,
                detail=f"Gene must be one of: {', '.join(active_genes())}",
            )
    return search_table9(
        gene=gene,
        query=query,
        code=code,
        page=page,
        page_size=page_size,
    )


@app.get(
    "/ui-api/rules/tables/{table_id}",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.get(
    "/api/rules/tables/{table_id}",
    dependencies=[Depends(require_public_api_key)],
)
async def enigma_reference_table_records(
    table_id: str,
    section: Optional[str] = Query(default=None, max_length=80),
    query: str = Query(default="", max_length=200),
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=25, ge=1, le=100),
):
    payload = search_reference_table(
        table_id,
        section_id=section,
        query=query,
        page=page,
        page_size=page_size,
    )
    if payload is None:
        raise HTTPException(status_code=404, detail="ENIGMA table or section not found")
    return payload


@app.post(
    "/ui-api/audit/client-validation",
    status_code=204,
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.post(
    "/api/audit/client-validation",
    status_code=204,
    dependencies=[Depends(require_public_api_key)],
)
async def client_validation_error(
    req: ClientValidationRequest,
    request: Request,
):
    _audit(
        request,
        "client_validation_error",
        level="warning",
        form=req.form,
        input=req.input,
        error=req.error,
    )
    return Response(status_code=204)


@app.post(
    "/ui-api/manual-evidence/evaluate",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.post(
    "/api/manual-evidence/evaluate",
    dependencies=[Depends(require_public_api_key)],
)
async def evaluate_manual_evidence_endpoint(
    req: ManualEvidenceRequest,
    request: Request,
) -> ManualEvidenceResult:
    try:
        execution = execute_manual_evidence(
            [criterion.model_dump() for criterion in req.base_criteria],
            [criterion.model_dump() for criterion in req.manual_criteria],
            req.variant_context.model_dump() if req.variant_context else None,
        )
        result = execution.result
    except DagNodeExecutionError as exc:
        cause = exc.__cause__
        if isinstance(cause, ValueError):
            detail = str(cause)
            status_code = 422
        else:
            detail = (
                f"Manual evidence evaluation could not complete at internal "
                f"step {exc.node_id}. No adjusted classification was returned."
            )
            status_code = 503
        _audit(
            request,
            "manual_evidence_error",
            level="warning",
            input=req.model_dump(mode="json"),
            error=str(exc),
            trace=[entry.as_dict() for entry in exc.trace],
        )
        raise HTTPException(status_code=status_code, detail=detail) from exc

    response = ManualEvidenceResult(
        predicted_class=result["predicted_class"],
        predicted_label=result["predicted_label"],
        total_points=result["total_points"],
        classification_note=result["classification_note"],
        manual_criteria=[
            ManualCriterionResult(**criterion)
            for criterion in result["manual_criteria"]
        ],
        evidence_interactions=[
            EvidenceInteractionWarning(**warning)
            for warning in result["evidence_interactions"]
        ],
        assessor=req.assessor,
        assessed_at=req.assessed_at,
    )
    _audit(
        request,
        "manual_evidence_completed",
        input=req.model_dump(mode="json"),
        result={
            "predicted_class": response.predicted_class,
            "predicted_label": response.predicted_label,
            "total_points": response.total_points,
            "evidence_interactions": [
                warning.model_dump(mode="json")
                for warning in response.evidence_interactions
            ],
        },
    )
    return response


@app.post(
    "/ui-api/manual-evidence/status",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.post(
    "/api/manual-evidence/status",
    dependencies=[Depends(require_public_api_key)],
)
async def manual_evidence_status_endpoint(
    req: ManualEvidenceStatusRequest,
) -> ManualEvidenceStatusResponse:
    statuses = manual_criterion_statuses(
        [criterion.model_dump() for criterion in req.base_criteria],
        [criterion.model_dump() for criterion in req.manual_criteria],
        req.variant_context.model_dump() if req.variant_context else None,
    )
    return ManualEvidenceStatusResponse(
        criteria=[ManualCriterionStatus(**status) for status in statuses]
    )


@app.post(
    "/ui-api/manual-evidence/resolve-ps1-reference",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.post(
    "/api/manual-evidence/resolve-ps1-reference",
    dependencies=[Depends(require_public_api_key)],
)
async def resolve_ps1_reference_endpoint(
    req: Ps1ReferenceResolutionRequest,
    request: Request,
) -> Ps1ReferenceResolutionResponse:
    try:
        result = await resolve_ps1_reference(
            req.gene,
            req.assessed_c_notation,
            req.reference_c_notation,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        _audit(
            request,
            "ps1_reference_resolution_error",
            level="exception",
            input=req.model_dump(mode="json"),
            error_type=type(exc).__name__,
            error=str(exc)[:2000],
        )
        raise HTTPException(
            status_code=503,
            detail="PS1 reference facts could not be resolved; no criterion was added.",
        ) from exc
    response = Ps1ReferenceResolutionResponse(**result)
    _audit(
        request,
        "ps1_reference_resolved",
        input=req.model_dump(mode="json"),
        result=response.model_dump(mode="json"),
    )
    return response


# Semaphore limits concurrent external API calls during batch processing
BATCH_CONCURRENCY = 3
BATCH_SEMAPHORE = asyncio.Semaphore(BATCH_CONCURRENCY)


def _batch_variant_label(index: int, raw_item: object) -> str:
    if isinstance(raw_item, dict):
        gene = str(raw_item.get("gene", "")).strip().upper()
        notation = str(raw_item.get("c_notation", "")).strip()
        label = " ".join(value for value in (gene, notation) if value)
        if label:
            return label
    return f"Batch item {index + 1}"


def _batch_validation_message(exc: ValidationError) -> str:
    messages: list[str] = []
    for error in exc.errors(include_url=False):
        location = ".".join(str(value) for value in error.get("loc", ()))
        message = str(error.get("msg", "Invalid variant input"))
        if message.startswith("Value error, "):
            message = message[len("Value error, "):]
        rendered = f"{location}: {message}" if location else message
        if rendered not in messages:
            messages.append(rendered)
    return "; ".join(messages) or "Invalid variant input"


def _execution_error(exc: Exception) -> tuple[str, str, bool]:
    if isinstance(exc, HTTPException):
        message = str(exc.detail)
        headers = exc.headers or {}
        explicit_code = headers.get("X-ARIANE-Error-Code")
        retryable_header = headers.get("X-ARIANE-Retryable")
        explicit_retryable = (
            retryable_header.lower() == "true"
            if retryable_header is not None
            else None
        )
        if exc.status_code == 422:
            return explicit_code or "variant_not_classifiable", message, False
        if exc.status_code == 503:
            return (
                explicit_code or "evidence_unavailable",
                message,
                True if explicit_retryable is None else explicit_retryable,
            )
        return (
            explicit_code or "request_failed",
            message,
            (
                exc.status_code in {429, 502, 504}
                if explicit_retryable is None
                else explicit_retryable
            ),
        )
    return "internal_error", str(exc) or type(exc).__name__, False


def _policy_metadata(gene: str) -> dict[str, str]:
    configured = get_gene_policy(gene)
    return {
        "gene": gene,
        "policy_id": configured["policy"]["runtime_policy_id"],
        "policy_version": configured["policy"]["version"],
        "reference_transcript": configured["gene_config"]["reference_transcript"],
    }


async def _classify_one(
    gene: str,
    c_notation: str,
    p_notation: str = "",
    dup_type: str = "Unknown",
) -> ClassificationResult:
    """Compatibility facade for single and batch API handlers."""
    try:
        return await execute_variant_classification(
            ClassificationCommand(
                gene=gene,
                c_notation=c_notation,
                p_notation=p_notation,
                dup_type=dup_type,
            ),
            engine_mode=CLASSIFIER_ENGINE_MODE,
            orchestration=CLASSIFICATION_ORCHESTRATION,
        )
    except VariantPreparationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except RequiredEvidenceUnavailableError as exc:
        headers = {
            "X-ARIANE-Error-Code": exc.code,
            "X-ARIANE-Retryable": str(exc.retryable).lower(),
        }
        if exc.retryable:
            headers["Retry-After"] = "5"
        raise HTTPException(
            status_code=503,
            detail=str(exc),
            headers=headers,
        ) from exc
    except EvidenceExecutionError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc


async def _classify_one_cached(
    gene: str,
    c_notation: str,
    p_notation: str = "",
    dup_type: str = "Unknown",
    reference_transcript: str = "",
) -> tuple[ClassificationResult, str, str]:
    """Return a current, fingerprint-bound result and its cache status."""
    fingerprint = classification_fingerprint(gene, CLASSIFIER_ENGINE_MODE.value)
    cache_read_status = "disabled" if CLASSIFICATION_CACHE is None else "miss"
    try:
        cached = CLASSIFICATION_CACHE.get(
            gene=gene,
            transcript=reference_transcript,
            c_notation=c_notation,
            p_notation=p_notation,
            dup_type=dup_type,
            fingerprint=fingerprint,
        ) if CLASSIFICATION_CACHE is not None else None
        if cached is not None:
            cache_read_status = cached.status
    except Exception:
        logging.getLogger("ariane.classification_cache").exception(
            "Classification cache read failed for %s:%s", gene, c_notation
        )
        cached = None
        cache_read_status = "read_error"
    if cached is not None and cached.result is not None:
        return cached.result, cached.status, fingerprint

    response = await _classify_one(gene, c_notation, p_notation, dup_type)
    cache_status = cache_read_status
    try:
        if CLASSIFICATION_CACHE is not None:
            CLASSIFICATION_CACHE.put(
                response,
                transcript=reference_transcript,
                dup_type=dup_type,
                fingerprint=fingerprint,
            )
    except Exception:
        logging.getLogger("ariane.classification_cache").exception(
            "Classification cache write failed for %s:%s", gene, c_notation
        )
        cache_status = "write_error"
    return response, cache_status, fingerprint


def _record_classification_usage(**values) -> None:
    """Keep statistics auxiliary so storage failure cannot change a result."""
    if CLASSIFICATION_USAGE is None:
        return
    try:
        CLASSIFICATION_USAGE.record(**values)
    except Exception:
        logging.getLogger("ariane.classification_usage").exception(
            "Classification usage event could not be stored"
        )


def _reserve_public_api_classifications(
    request: Request,
    response: Response,
    units: int,
) -> None:
    """Reserve per-key work before classification; browser traffic is separate."""
    api_key_id = str(getattr(request.state, "api_key_id", "") or "").strip()
    if not api_key_id:
        return
    if PUBLIC_API_QUOTA is None:
        raise HTTPException(
            status_code=503,
            detail={
                "code": "api_quota_unavailable",
                "message": "Public API quota storage is unavailable",
                "retryable": True,
            },
            headers={"Retry-After": "60"},
        )
    try:
        reservation = PUBLIC_API_QUOTA.reserve_public_api_classifications(
            api_key_id=api_key_id,
            units=units,
            limit=PUBLIC_API_DAILY_CLASSIFICATION_LIMIT,
        )
    except Exception as exc:
        logging.getLogger("ariane.api_quota").exception(
            "Public API quota reservation failed for key %s", api_key_id
        )
        raise HTTPException(
            status_code=503,
            detail={
                "code": "api_quota_unavailable",
                "message": "Public API quota storage is unavailable",
                "retryable": True,
            },
            headers={"Retry-After": "60"},
        ) from exc
    reset_epoch = str(int(reservation.reset_at.timestamp()))
    quota_headers = {
        "X-RateLimit-Limit": str(reservation.limit),
        "X-RateLimit-Remaining": str(reservation.remaining),
        "X-RateLimit-Reset": reset_epoch,
    }
    if not reservation.allowed:
        retry_after = max(
            1,
            int(
                reservation.reset_at.timestamp()
                - datetime.now(timezone.utc).timestamp()
            ),
        )
        raise HTTPException(
            status_code=429,
            detail={
                "code": "daily_classification_quota_exceeded",
                "message": (
                    "The API key has reached its daily classification quota. "
                    "The quota resets at 00:00 UTC."
                ),
                "retryable": True,
            },
            headers={**quota_headers, "Retry-After": str(retry_after)},
        )
    for name, value in quota_headers.items():
        response.headers[name] = value


@app.post(
    "/ui-api/normalize",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.post("/api/normalize", dependencies=[Depends(require_public_api_key)])
async def normalize_variant(
    req: VariantRequest, request: Request
) -> VariantNormalizationResponse:
    response = VariantNormalizationResponse(
        gene=req.gene,
        submitted_notation=req.submitted_notation,
        c_notation=req.c_notation,
        p_notation=req.p_notation,
        reference_transcript=req.reference_transcript,
        normalization_source=req.normalization_source,
        consequence_status=req.consequence_status,
        normalization_provenance=req.normalization_provenance,
        protein_consequence_explanation=req.protein_consequence_explanation,
        assembly=req.assembly,
    )
    _audit(
        request,
        "variant_normalized",
        input={"gene": req.gene, "notation": req.submitted_notation, "assembly": req.assembly},
        result=response.model_dump(mode="json"),
    )
    return response


@app.post(
    "/ui-api/classify",
    dependencies=[Depends(require_ui_session)],
    include_in_schema=False,
)
@app.post("/api/classify", dependencies=[Depends(require_public_api_key)])
async def classify_variant(
    req: VariantRequest,
    request: Request,
    http_response: Response,
) -> ClassificationResult:
    _reserve_public_api_classifications(request, http_response, 1)
    input_data = req.model_dump(mode="json")
    identity = resolve_usage_identity(request)
    identity.set_cookie(http_response)
    policy = get_gene_policy(req.gene)["policy"]
    started = time.monotonic()
    cache_status = "not_attempted"
    fingerprint = classification_fingerprint(req.gene, CLASSIFIER_ENGINE_MODE.value)
    try:
        response, cache_status, fingerprint = await _classify_one_cached(
            req.gene,
            req.c_notation,
            req.p_notation or "",
            req.dup_type,
            req.reference_transcript,
        )
        response.reference_transcript = req.reference_transcript
        response.submitted_notation = req.submitted_notation
        response.normalization_source = req.normalization_source
        response.consequence_status = req.consequence_status
        response.normalization_provenance = req.normalization_provenance
        response.protein_consequence_explanation = req.protein_consequence_explanation
    except Exception as exc:
        _record_classification_usage(
            actor_id=identity.actor_id,
            actor_type=identity.actor_type,
            request_id=request.state.request_id,
            request_mode="single",
            gene=req.gene,
            transcript=req.reference_transcript,
            c_notation=req.c_notation,
            status="error",
            cache_status=cache_status,
            predicted_class=None,
            total_points=None,
            duration_ms=(time.monotonic() - started) * 1000,
            policy_id=policy["runtime_policy_id"],
            policy_version=policy["version"],
            classifier_fingerprint=fingerprint,
        )
        _audit(
            request,
            "classification_error",
            level="exception",
            input=input_data,
            error_type=type(exc).__name__,
            error=str(exc)[:2000],
        )
        raise
    _record_classification_usage(
        actor_id=identity.actor_id,
        actor_type=identity.actor_type,
        request_id=request.state.request_id,
        request_mode="single",
        gene=req.gene,
        transcript=req.reference_transcript,
        c_notation=req.c_notation,
        status="completed",
        cache_status=cache_status,
        predicted_class=response.predicted_class,
        total_points=response.total_points,
        duration_ms=(time.monotonic() - started) * 1000,
        policy_id=policy["runtime_policy_id"],
        policy_version=policy["version"],
        classifier_fingerprint=fingerprint,
    )
    duration_ms = (time.monotonic() - started) * 1000
    request.state.classification_execution = {
        "policy": _policy_metadata(req.gene),
        "classifier_fingerprint": fingerprint,
        "cache_status": cache_status,
        "duration_ms": duration_ms,
    }
    http_response.headers["X-ARIANE-Cache-Status"] = cache_status
    http_response.headers["X-ARIANE-Classifier-Fingerprint"] = fingerprint
    _audit(
        request,
        "classification_completed",
        input=input_data,
        result={
            "predicted_class": response.predicted_class,
            "predicted_label": response.predicted_label,
            "total_points": response.total_points,
            "evidence_direction": response.evidence_direction,
            "mixed_evidence": response.mixed_evidence,
            "pathogenic_points": response.pathogenic_points,
            "benign_points": response.benign_points,
            "classification_cache": cache_status,
            "evidence_interactions": [
                warning.model_dump(mode="json")
                for warning in response.evidence_interactions
            ],
            "spliceai_audit": response.spliceai_audit.model_dump(mode="json")
            if response.spliceai_audit else None,
        },
    )
    return response


@app.post(
    "/api/classify/batch",
    dependencies=[Depends(require_public_api_key)],
)
async def classify_batch(
    req: BatchRequest,
    request: Request,
    http_response: Response,
) -> BatchResponse:
    """
    Classify multiple variants. Up to 200 per request.
    Results preserve input order. Per-variant errors are reported inline.
    Concurrency is limited to avoid overwhelming external APIs.
    """
    batch_started = time.monotonic()
    identity = resolve_usage_identity(request)
    identity.set_cookie(http_response)

    async def _one(idx: int, item: VariantRequest):
        started = time.monotonic()
        cache_status = "not_attempted"
        fingerprint = classification_fingerprint(item.gene, CLASSIFIER_ENGINE_MODE.value)
        async with BATCH_SEMAPHORE:
            try:
                res, cache_status, fingerprint = await _classify_one_cached(
                    item.gene,
                    item.c_notation,
                    item.p_notation or "",
                    item.dup_type,
                    item.reference_transcript,
                )
                res.reference_transcript = item.reference_transcript
                res.submitted_notation = item.submitted_notation
                res.normalization_source = item.normalization_source
                res.consequence_status = item.consequence_status
                res.normalization_provenance = item.normalization_provenance
                res.protein_consequence_explanation = item.protein_consequence_explanation
                output = BatchItemResult(
                    index=idx, status="ok",
                    variant=f"{item.gene} {item.c_notation}",
                    result=res,
                )
            except Exception as exc:
                error_code, error_message, retryable = _execution_error(exc)
                output = BatchItemResult(
                    index=idx, status="error",
                    variant=f"{item.gene} {item.c_notation}",
                    error=error_message,
                    error_code=error_code,
                    error_retryable=retryable,
                )
            policy = get_gene_policy(item.gene)["policy"]
            duration_ms = (time.monotonic() - started) * 1000
            output.cache_status = cache_status
            output.classifier_fingerprint = fingerprint
            output.policy_id = policy["runtime_policy_id"]
            output.policy_version = policy["version"]
            output.duration_ms = duration_ms
            return (
                item,
                output,
                cache_status,
                fingerprint,
                duration_ms,
            )

    valid_items: list[tuple[int, VariantRequest]] = []
    validation_errors: list[BatchItemResult] = []
    for idx, raw_item in enumerate(req.variants):
        try:
            item = VariantRequest.model_validate(raw_item)
        except ValidationError as exc:
            error = _batch_validation_message(exc)
            validation_errors.append(BatchItemResult(
                index=idx,
                status="error",
                variant=_batch_variant_label(idx, raw_item),
                error=error,
                error_code="invalid_variant",
                error_retryable=False,
            ))
            if isinstance(raw_item, dict):
                safe_input = {
                    "gene": str(raw_item.get("gene", ""))[:40],
                    "c_notation": str(raw_item.get("c_notation", ""))[:200],
                }
            else:
                safe_input = {"input_type": type(raw_item).__name__}
            _audit(
                request,
                "batch_item_validation_error",
                level="warning",
                item_index=idx,
                input=safe_input,
                error=error,
            )
        else:
            valid_items.append((idx, item))

    if valid_items:
        _reserve_public_api_classifications(
            request,
            http_response,
            len(valid_items),
        )

    executions = await asyncio.gather(*[
        _one(idx, item) for idx, item in valid_items
    ])
    items = sorted(
        validation_errors + [execution[1] for execution in executions],
        key=lambda value: value.index,
    )
    success = sum(1 for r in items if r.status == "ok")
    for input_item, output_item, cache_status, fingerprint, duration_ms in executions:
        policy = get_gene_policy(input_item.gene)["policy"]
        _record_classification_usage(
            actor_id=identity.actor_id,
            actor_type=identity.actor_type,
            request_id=request.state.request_id,
            request_mode="batch",
            gene=input_item.gene,
            transcript=input_item.reference_transcript,
            c_notation=input_item.c_notation,
            status="completed" if output_item.status == "ok" else "error",
            cache_status=cache_status,
            predicted_class=(
                output_item.result.predicted_class if output_item.result else None
            ),
            total_points=(output_item.result.total_points if output_item.result else None),
            duration_ms=duration_ms,
            policy_id=policy["runtime_policy_id"],
            policy_version=policy["version"],
            classifier_fingerprint=fingerprint,
        )
        _audit(
            request,
            "batch_item_completed" if output_item.status == "ok" else "batch_item_error",
            level="info" if output_item.status == "ok" else "warning",
            item_index=output_item.index,
            input=input_item.model_dump(mode="json"),
            result={
                "predicted_class": output_item.result.predicted_class,
                "predicted_label": output_item.result.predicted_label,
                "total_points": output_item.result.total_points,
                "classification_cache": cache_status,
            } if output_item.result else None,
            error=output_item.error,
        )
    response = BatchResponse(
        total=len(items),
        success_count=success,
        error_count=len(items) - success,
        results=list(items),
    )
    request.state.batch_duration_ms = (time.monotonic() - batch_started) * 1000
    return response


app.include_router(create_public_api_router(
    classify_single=classify_variant,
    classify_batch=classify_batch,
    engine=CLASSIFIER_ENGINE_MODE.value,
    batch_concurrency=BATCH_CONCURRENCY,
    per_key_classifications_per_utc_day=PUBLIC_API_DAILY_CLASSIFICATION_LIMIT,
))


@app.post("/api/clear-cache")
async def clear_cache(
    request: Request,
    x_ariane_admin_token: Optional[str] = Header(default=None),
):
    admin_token = os.getenv("ARIANE_ADMIN_TOKEN", "")
    if not admin_token:
        raise HTTPException(status_code=503, detail="Administrative API is disabled")
    if not x_ariane_admin_token or not secrets.compare_digest(
        x_ariane_admin_token, admin_token
    ):
        raise HTTPException(status_code=403, detail="Invalid administrative token")
    clear_runtime_caches()

    POPULATION_FREQUENCY_SERVICE.reload()

    _audit(request, "cache_cleared")
    return {"status": "ok", "message": "All caches cleared"}
