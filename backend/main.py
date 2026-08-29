# ============================================================
# ARIANE - FastAPI application
# Automated ACMG Rule-based Interpretation and Annotation ENgine
# ============================================================
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.encoders import jsonable_encoder
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse, Response
from pathlib import Path
from datetime import datetime, timezone
import json
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
from backend.version import ARIANE_VERSION

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
from backend.classification_dag import (
    DagNodeExecutionError,
    execute_manual_evidence,
    get_configured_engine_mode,
)
from backend.classification_dag.provider_wiring import production_provider_dependencies
from backend.models import (
    VariantRequest, ClassificationResult,
    BatchRequest, BatchResponse, BatchItemResult,
    ManualEvidenceRequest, ManualEvidenceResult,
    ManualCriterionResult, EvidenceInteractionWarning,
    ClientValidationRequest, VariantNormalizationResponse,
    Ps1ReferenceResolutionRequest, Ps1ReferenceResolutionResponse,
)
from backend.services import (
    ClassificationCommand,
    EvidenceOrchestrationService,
    EvidenceExecutionError,
    VariantPreparationError,
    execute_variant_classification,
    resolve_ps1_reference,
)
from backend.modules.hgvs_provider import load_panel_provider
from backend.modules.manual_evidence import (
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
from backend.lookups.indels import load_indel_snapshot  # noqa: E402
from backend.lookups.precomputed import validate_classification_snapshot  # noqa: E402
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
validate_classification_snapshot()
load_indel_snapshot()
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

# ── App setup ──────────────────────────────────────────────────────────────
app = FastAPI(
    title="ARIANE",
    description="Automated ACMG Rule-based Interpretation and Annotation ENgine",
    version=ARIANE_VERSION,
)
app.include_router(admin_router)

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


@app.middleware("http")
async def audit_request(request: Request, call_next):
    request.state.request_id = request.headers.get("x-request-id") or uuid.uuid4().hex
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
    log_completion = (
        request.url.path.startswith("/admin/")
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
    return JSONResponse(
        status_code=422,
        content={"detail": jsonable_encoder(exc.errors())},
    )


FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"
app.mount("/static", StaticFiles(directory=FRONTEND_DIR / "static"), name="static")


@app.get("/", response_class=HTMLResponse)
async def index():
    return FileResponse(FRONTEND_DIR / "index.html")


@app.get("/api/health")
async def health():
    issues = get_data_issues()
    panel = load_panel_provider()
    return {
        "status": "degraded" if issues else "ok",
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
        },
        "data_issues": issues,
    }


@app.get("/api/resources")
async def resources(gene: Optional[str] = None):
    return {
        "version": ARIANE_VERSION,
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


@app.get("/api/rules")
async def enigma_rules_catalog():
    """Public, versioned source and rule index without local file paths."""
    return public_catalog()


@app.get("/api/rules/trees/{tree_id}")
async def enigma_decision_tree(tree_id: str):
    tree = get_decision_tree(tree_id)
    if tree is None:
        raise HTTPException(status_code=404, detail="Decision tree not found")
    return tree


@app.get("/api/rules/tables/table9")
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


@app.get("/api/rules/tables/{table_id}")
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


@app.post("/api/audit/client-validation", status_code=204)
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


@app.post("/api/manual-evidence/evaluate")
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


@app.post("/api/manual-evidence/resolve-ps1-reference")
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
BATCH_SEMAPHORE = asyncio.Semaphore(3)

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
    except EvidenceExecutionError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc


@app.post("/api/normalize")
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


@app.post("/api/classify")
async def classify_variant(req: VariantRequest, request: Request) -> ClassificationResult:
    input_data = req.model_dump(mode="json")
    try:
        response = await _classify_one(
            req.gene, req.c_notation, req.p_notation or "", req.dup_type
        )
        response.reference_transcript = req.reference_transcript
        response.submitted_notation = req.submitted_notation
        response.normalization_source = req.normalization_source
        response.consequence_status = req.consequence_status
        response.normalization_provenance = req.normalization_provenance
        response.protein_consequence_explanation = req.protein_consequence_explanation
    except Exception as exc:
        _audit(
            request,
            "classification_error",
            level="exception",
            input=input_data,
            error_type=type(exc).__name__,
            error=str(exc)[:2000],
        )
        raise
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
            "evidence_interactions": [
                warning.model_dump(mode="json")
                for warning in response.evidence_interactions
            ],
            "spliceai_audit": response.spliceai_audit.model_dump(mode="json")
            if response.spliceai_audit else None,
        },
    )
    return response


@app.post("/api/classify/batch")
async def classify_batch(req: BatchRequest, request: Request) -> BatchResponse:
    """
    Classify multiple variants. Up to 200 per request.
    Results preserve input order. Per-variant errors are reported inline.
    Concurrency is limited to avoid overwhelming external APIs.
    """
    async def _one(idx: int, item: VariantRequest) -> BatchItemResult:
        async with BATCH_SEMAPHORE:
            try:
                res = await _classify_one(
                    item.gene, item.c_notation, item.p_notation or "", item.dup_type
                )
                res.reference_transcript = item.reference_transcript
                res.submitted_notation = item.submitted_notation
                res.normalization_source = item.normalization_source
                res.consequence_status = item.consequence_status
                res.normalization_provenance = item.normalization_provenance
                res.protein_consequence_explanation = item.protein_consequence_explanation
                return BatchItemResult(
                    index=idx, status="ok",
                    variant=f"{item.gene} {item.c_notation}",
                    result=res,
                )
            except Exception as exc:
                return BatchItemResult(
                    index=idx, status="error",
                    variant=f"{item.gene} {item.c_notation}",
                    error=str(exc),
                )

    items = await asyncio.gather(*[_one(i, v) for i, v in enumerate(req.variants)])
    items = sorted(items, key=lambda r: r.index)
    success = sum(1 for r in items if r.status == "ok")
    for input_item, output_item in zip(req.variants, items):
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
            } if output_item.result else None,
            error=output_item.error,
        )
    return BatchResponse(
        total=len(items),
        success_count=success,
        error_count=len(items) - success,
        results=list(items),
    )


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
