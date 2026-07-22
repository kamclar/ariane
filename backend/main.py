# ============================================================
# ARIANE - FastAPI application
# Automated ACMG Rule-based Interpretation and Annotation ENgine
# ============================================================
from fastapi import FastAPI, HTTPException, Request
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

from backend.config import (
    DATA_DIR, TABLE4_PATH, TABLE9_PATH, ST7_PATH,
    ALLOWED_GENES, TRANSCRIPTS,
)
from backend.data_validation import validate_required_datasets
from backend.data_health import get_data_issues, get_user_warnings
from backend.lookup_execution import lookup_or_unavailable
from backend.models import (
    VariantRequest, ClassificationResult, CriterionResult,
    ExternalComparison, ExternalSubmitter, CLASS_LABELS,
    BatchRequest, BatchResponse, BatchItemResult,
    AlphaMissenseResult, VusExplanation,
    RnaReviewRecommendation,
    ManualEvidenceRequest, ManualEvidenceResult,
    ManualCriterionResult, ClientValidationRequest,
)

validate_required_datasets({"table4": TABLE4_PATH, "table9": TABLE9_PATH, "st7": ST7_PATH})

# Initialize local sources before serving requests so /api/health reports
# degraded caches even before the first classification.
from backend.modules import frequency as _frequency_data_source  # noqa: E402,F401
from backend.lookups import coordinates as _coordinate_data_source  # noqa: E402,F401
from backend.lookups import bayesdel as _bayesdel_data_source  # noqa: E402,F401
from backend.lookups import spliceai as _spliceai_data_source  # noqa: E402
from backend.modules.residues import initialize_residue_data  # noqa: E402

_spliceai_data_source._load_precomputed_cache()
_spliceai_data_source._load_api_cache()
initialize_residue_data()

# ── App setup ──────────────────────────────────────────────────────────────
app = FastAPI(
    title="ARIANE",
    description="Automated ACMG Rule-based Interpretation and Annotation ENgine for BRCA1/2",
    version="1.8.0",
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
    return {
        "status": "degraded" if issues else "ok",
        "version": "1.8.0",
        "data": {
            "table4": TABLE4_PATH.exists(),
            "table9": TABLE9_PATH.exists(),
            "st7":    ST7_PATH.exists(),
        },
        "data_issues": issues,
    }


@app.get("/api/resources")
async def resources():
    from backend.modules.manual_evidence import MANUAL_CRITERIA, RESOURCE_LINKS
    from backend.modules.splice_ps1_reference import (
        load_splice_ps1_reference_candidates,
    )

    return {
        "manual_criteria": MANUAL_CRITERIA,
        "links": RESOURCE_LINKS,
        "splice_ps1_reference_candidates": load_splice_ps1_reference_candidates(),
    }


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
    from backend.modules.manual_evidence import evaluate_manual_evidence

    try:
        result = evaluate_manual_evidence(
            [criterion.model_dump() for criterion in req.base_criteria],
            [criterion.model_dump() for criterion in req.manual_criteria],
        )
    except ValueError as exc:
        _audit(
            request,
            "manual_evidence_error",
            level="warning",
            input=req.model_dump(mode="json"),
            error=str(exc),
        )
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    response = ManualEvidenceResult(
        predicted_class=result["predicted_class"],
        predicted_label=result["predicted_label"],
        total_points=result["total_points"],
        classification_note=result["classification_note"],
        manual_criteria=[
            ManualCriterionResult(**criterion)
            for criterion in result["manual_criteria"]
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
        },
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
    """Core classification logic shared by single and batch endpoints."""
    from backend.lookups.coordinates import resolve_variant, get_grch37, get_grch38
    from backend.modules.variant_type import infer_variant_type
    from backend.lookups.spliceai import get_spliceai_score
    from backend.lookups.bayesdel import get_bayesdel_and_alphamissense
    from backend.lookups.clinvar import clinvar_lookup, clinvar_review_stars
    from backend.lookups.clingen import clingen_erepo_lookup
    from backend.modules.frequency import get_gnomad_frequencies
    from backend.modules.table9 import table9_lookup_ps3_bs3
    from backend.modules.pp4_bp5 import evaluate_pp4_bp5
    from backend.modules.ps1 import evaluate_ps1
    from backend.modules.residues import check_important_residue
    from backend.modules.classifier import evaluate_variant as _evaluate
    from backend.modules.external import external_comparison
    from backend.lookups.precomputed import lookup_classification_snapshot
    from backend.modules.reference_validation import validate_reference_allele

    # Reject a wrong stated reference before coordinates, external lookups, or
    # evidence evaluation. This is deliberately fail-closed.
    try:
        validate_reference_allele(gene, c_notation)
    except (ValueError, RuntimeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    # The versioned coding-SNV snapshot carries the reference-transcript
    # protein consequence. Use it when the optional p. field was omitted and
    # reject contradictory user input instead of silently choosing one.
    snapshot = lookup_classification_snapshot(gene, c_notation)
    snapshot_p = ""
    if snapshot:
        snapshot_p = str(snapshot.get("record", {}).get("p_notation") or "")
    if not p_notation and snapshot_p:
        p_notation = snapshot_p
    elif p_notation and snapshot_p and p_notation != snapshot_p:
        raise HTTPException(
            status_code=422,
            detail=(
                f"Protein consequence mismatch for {gene} {c_notation}: "
                f"the configured reference transcript gives {snapshot_p}, not {p_notation}."
            ),
        )

    # Step 1: resolve coordinates
    resolved = {}
    lookup_diagnostics = []
    try:
        rv = resolve_variant(gene, c_notation)
        if rv:
            resolved[f"{gene}:{c_notation}"] = rv
            if rv.status != "ok" or rv.source == "Mutalyzer":
                lookup_diagnostics.extend(
                    f"Coordinate resolver: {warning}" for warning in rv.warnings
                )
    except Exception as exc:
        message = f"Coordinate lookup failed: {type(exc).__name__}: {exc}"
        logging.getLogger(__name__).exception(message)
        lookup_diagnostics.append(message)
    grch37 = get_grch37(resolved, gene, c_notation)
    grch38 = get_grch38(resolved, gene, c_notation)

    # Step 2: variant type
    variant_type = infer_variant_type(c_notation, p_notation)

    # Step 3: parallel external lookups
    # get_bayesdel_and_alphamissense returns (bayesdel_score, alphamissense_dict)
    # in a single myvariant.info call - no extra API overhead for AlphaMissense.
    spliceai_score, (bayesdel_score, alphamissense), cv, er = await asyncio.gather(
        lookup_or_unavailable(get_spliceai_score, None, "SpliceAI", lookup_diagnostics, gene, c_notation),
        lookup_or_unavailable(get_bayesdel_and_alphamissense, (None, None), "MyVariant/BayesDel", lookup_diagnostics, gene, c_notation),
        lookup_or_unavailable(
            clinvar_lookup,
            {"status": "api_timeout", "error": "ClinVar lookup timed out"},
            "ClinVar", lookup_diagnostics,
            gene,
            c_notation,
        ),
        lookup_or_unavailable(
            clingen_erepo_lookup,
            {"status": "api_timeout", "error": "ClinGen ERepo lookup timed out"},
            "ClinGen ERepo", lookup_diagnostics,
            gene,
            c_notation,
        ),
    )

    from backend.lookups.spliceai import SPLICEAI_STATUS_CACHE
    from backend.lookups.bayesdel import BAYESDEL_STATUS_CACHE
    variant_key = f"{gene}:{c_notation}"
    splice_status = SPLICEAI_STATUS_CACHE.get(variant_key, {})
    if splice_status.get("status") not in {None, "ok"}:
        lookup_diagnostics.append(
            f"SpliceAI unavailable: status={splice_status.get('status')}; "
            f"{splice_status.get('reason', 'no reason reported')}"
        )
    bayesdel_status = BAYESDEL_STATUS_CACHE.get(variant_key, {})
    if bayesdel_status.get("status") in {"api_error", "no_grch37_coords"}:
        lookup_diagnostics.append(
            f"MyVariant/BayesDel unavailable: status={bayesdel_status.get('status')}; "
            f"{bayesdel_status.get('reason', 'no reason reported')}"
        )

    # Step 4: fast local lookups
    gnomad_data = None
    if grch37 or grch38:
        gnomad_data = get_gnomad_frequencies(grch37=grch37, grch38=grch38)

    table9_result  = table9_lookup_ps3_bs3(gene, c_notation)
    pp4_bp5_result = evaluate_pp4_bp5(gene, c_notation)
    ps1_result     = evaluate_ps1(
        gene, c_notation, p_notation,
        variant_type=variant_type,
        spliceai_score=spliceai_score,
    )
    residue_info = check_important_residue(gene, p_notation)

    # Step 5: evaluate
    result = _evaluate(
        gene=gene, variant_type=variant_type,
        p_notation=p_notation, c_notation=c_notation,
        spliceai_score=spliceai_score, bayesdel_score=bayesdel_score,
        gnomad_data=gnomad_data, table9_result=table9_result,
        pp4_bp5_result=pp4_bp5_result, ps1_result=ps1_result,
        residue_info=residue_info, dup_type=dup_type,
    )
    for diagnostic in lookup_diagnostics:
        result["warnings"].append(f"External evidence unavailable: {diagnostic}")
    for warning in get_user_warnings():
        if warning not in result["warnings"]:
            result["warnings"].append(warning)
    if cv.get("status") == "ambiguous":
        result["warnings"].append(
            "ClinVar lookup was ambiguous; no external ClinVar record was selected. "
            f"Candidate IDs: {', '.join(cv.get('candidate_ids', [])) or 'not reported'}."
        )
    elif cv.get("status") not in {"ok", "not_found"}:
        result["warnings"].append(
            f"ClinVar evidence unavailable: status={cv.get('status')}; {cv.get('error', '')}".rstrip("; ")
        )
    if er.get("status") not in {"ok", "not_found"}:
        result["warnings"].append(
            f"ClinGen ERepo evidence unavailable: status={er.get('status')}; {er.get('error', '')}".rstrip("; ")
        )

    # Step 6: external comparison
    ext = external_comparison(gene, c_notation, result["predicted_class"], cv, er)

    # Step 7: narrative summary
    from backend.modules.vus_explanation import explain_vus
    vus_explanation = explain_vus(result)
    from backend.modules.narrative import generate_narrative
    narrative = generate_narrative(
        gene=gene,
        c_notation=c_notation,
        p_notation=p_notation,
        variant_type=variant_type,
        result=result,
        spliceai_score=spliceai_score,
        bayesdel_score=bayesdel_score,
        alphamissense=alphamissense,
    )

    # Step 8: build response model
    criteria = [
        CriterionResult(
            name=name,
            applies=crit.get("applies", True),
            strength=crit.get("strength"),
            points=crit.get("points", 0),
            reason=crit.get("reason", ""),
        )
        for name, crit in result["criteria"].items()
    ]

    ext_model = None
    if cv.get("status") == "ok":
        submitters = [
            ExternalSubmitter(
                scv=s.get("scv", ""),
                org=s.get("org", ""),
                classification=s.get("class") or "",
                date_eval=s.get("date_eval", ""),
                is_enigma_ep=s.get("is_enigma_ep", False),
                review_status=s.get("review", ""),
                curated_status=(
                    "ClinGen/ENIGMA curated submitter"
                    if s.get("is_enigma_ep", False)
                    else ""
                ),
                comment=s.get("comment", "")[:200],
            )
            for s in cv.get("submissions", [])
        ]
        ext_model = ExternalComparison(
            clinvar_classification=cv.get("aggregate", {}).get("classification", ""),
            clinvar_review_status=cv.get("aggregate", {}).get("review_status", ""),
            clinvar_review_stars=clinvar_review_stars(
                cv.get("aggregate", {}).get("review_status", "")
            ),
            clinvar_n_submitters=cv.get("aggregate", {}).get("n_submitters", 0),
            clinvar_has_conflict=cv.get("has_conflict", False),
            clinvar_submitters=submitters,
            enigma_ep_class=ext.get("enigma_class", ""),
            enigma_ep_source=ext.get("enigma_source", ""),
        )

    return ClassificationResult(
        variant=result["variant"],
        gene=gene,
        c_notation=c_notation,
        p_notation=p_notation,
        predicted_class=result["predicted_class"],
        predicted_label=CLASS_LABELS.get(result["predicted_class"], ""),
        total_points=result["total_points"],
        criteria=criteria,
        warnings=result["warnings"],
        external=ext_model,
        has_functional_evidence=result.get("has_functional_evidence", False),
        classification_note=result.get("classification_note", ""),
        narrative=narrative,
        alphamissense=AlphaMissenseResult(
            am_score=alphamissense.get("am_score") if alphamissense else None,
            am_class=alphamissense.get("am_class", "") if alphamissense else "",
        ) if alphamissense else None,
        vus_explanation=VusExplanation(**vus_explanation) if vus_explanation else None,
        rna_review=RnaReviewRecommendation(**result["rna_review"])
        if result.get("rna_review") else None,
        splice_ps1_review=RnaReviewRecommendation(**result["splice_ps1_review"])
        if result.get("splice_ps1_review") else None,
        initiation_review=RnaReviewRecommendation(**result["initiation_review"])
        if result.get("initiation_review") else None,
    )


@app.post("/api/classify")
async def classify_variant(req: VariantRequest, request: Request) -> ClassificationResult:
    input_data = req.model_dump(mode="json")
    try:
        response = await _classify_one(
            req.gene, req.c_notation, req.p_notation or "", req.dup_type
        )
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
    from backend.lookups.spliceai import SPLICEAI_CACHE, SPLICEAI_STATUS_CACHE
    from backend.lookups.bayesdel import BAYESDEL_CACHE
    from backend.lookups.clinvar import CLINVAR_CACHE
    from backend.lookups.clingen import EREPO_CACHE
    from backend.modules.frequency import GNOMAD_CACHE, load_gnomad_local_cache, load_gnomad_coverage_cache

    SPLICEAI_CACHE.clear()
    SPLICEAI_STATUS_CACHE.clear()
    BAYESDEL_CACHE.clear()
    CLINVAR_CACHE.clear()
    EREPO_CACHE.clear()

    load_gnomad_local_cache()
    load_gnomad_coverage_cache()

    _audit(request, "cache_cleared")
    return {"status": "ok", "message": "All caches cleared"}
