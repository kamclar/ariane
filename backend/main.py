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
    PS1_PROTEIN_REGISTRY_PATH, PS1_SPLICE_REFERENCE_PATH,
    ST2_SPLICE_EVIDENCE_PATH,
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
    SpliceAIAudit,
    RnaReviewRecommendation, ProteinPs1ReviewRecommendation,
    ManualEvidenceRequest, ManualEvidenceResult,
    ManualCriterionResult, EvidenceInteractionWarning,
    ClientValidationRequest, VariantNormalizationResponse,
)

validate_required_datasets({
    "table4": TABLE4_PATH,
    "table9": TABLE9_PATH,
    "st7": ST7_PATH,
    "ps1_protein_registry": PS1_PROTEIN_REGISTRY_PATH,
    "ps1_splice_reference": PS1_SPLICE_REFERENCE_PATH,
    "st2_splice_evidence": ST2_SPLICE_EVIDENCE_PATH,
})

# Initialize local sources before serving requests so /api/health reports
# degraded caches even before the first classification.
from backend.modules import frequency as _frequency_data_source  # noqa: E402,F401
from backend.lookups import coordinates as _coordinate_data_source  # noqa: E402,F401
from backend.lookups import bayesdel as _bayesdel_data_source  # noqa: E402,F401
from backend.lookups import spliceai as _spliceai_data_source  # noqa: E402
from backend.lookups.indels import load_indel_snapshot  # noqa: E402
from backend.lookups.precomputed import validate_classification_snapshot  # noqa: E402
from backend.modules.pp4_bp5 import load_pp4_bp5_snapshot  # noqa: E402
from backend.modules.residues import initialize_residue_data  # noqa: E402
from backend.modules.hgvs_engine import validate_hgvs_engine  # noqa: E402

_spliceai_data_source._load_precomputed_cache()
_spliceai_data_source._load_api_cache()
validate_classification_snapshot()
load_indel_snapshot()
load_pp4_bp5_snapshot()
initialize_residue_data()
validate_hgvs_engine()

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
    from backend.modules.hgvs_provider import load_panel_provider
    panel = load_panel_provider()
    return {
        "status": "degraded" if issues else "ok",
        "version": "1.8.0",
        "data": {
            "table4": TABLE4_PATH.exists(),
            "table9": TABLE9_PATH.exists(),
            "st7":    ST7_PATH.exists(),
            "ps1_protein_registry": PS1_PROTEIN_REGISTRY_PATH.exists(),
            "ps1_splice_reference": PS1_SPLICE_REFERENCE_PATH.exists(),
            "st2_splice_evidence": ST2_SPLICE_EVIDENCE_PATH.exists(),
            "reference_bundle": panel.provenance.get("reference_bundle", ""),
            "normalization_engine": panel.provenance.get("normalization_engine", ""),
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
    from backend.modules.ps1 import evaluate_ps1, select_vua_spliceai_for_ps1
    from backend.modules.residues import check_important_residue
    from backend.modules.classifier import evaluate_variant as _evaluate
    from backend.modules.external import external_comparison
    from backend.modules.variant_input import normalize_variant_input

    # Every entry point, including internal/batch calls, uses the same local
    # reference-transcript normalizer before any evidence or external lookup.
    try:
        normalized_input = normalize_variant_input(
            gene, c_notation, p_notation=p_notation or None
        )
    except (ValueError, RuntimeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    c_notation = normalized_input.c_notation
    p_notation = normalized_input.p_notation

    # Determine the variant type before planning external lookups. Exon-level
    # CNVs with uncertain breakpoints cannot be represented by one genomic
    # allele, so SNV/small-variant coordinate services are not applicable.
    variant_type = infer_variant_type(c_notation, p_notation)
    is_exon_cnv = variant_type.lower() in {"exon_deletion", "exon_duplication"}

    # Step 1: resolve coordinates where the HGVS description has exact bounds.
    resolved = {}
    lookup_diagnostics = []
    if not is_exon_cnv:
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

    # Step 2: parallel external lookups
    # get_bayesdel_and_alphamissense returns (bayesdel_score, alphamissense_dict)
    # in a single myvariant.info call - no extra API overhead for AlphaMissense.
    external_tasks = [
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
    ]
    if is_exon_cnv:
        cv, er = await asyncio.gather(*external_tasks)
        spliceai_score = None
        bayesdel_score, alphamissense = None, None
    else:
        spliceai_score, (bayesdel_score, alphamissense), cv, er = await asyncio.gather(
            lookup_or_unavailable(
                get_spliceai_score, None, "SpliceAI", lookup_diagnostics,
                gene, c_notation,
            ),
            lookup_or_unavailable(
                get_bayesdel_and_alphamissense, (None, None),
                "MyVariant/BayesDel", lookup_diagnostics, gene, c_notation,
            ),
            *external_tasks,
        )

    from backend.lookups.spliceai import SPLICEAI_STATUS_CACHE
    from backend.lookups.bayesdel import BAYESDEL_STATUS_CACHE
    variant_key = f"{gene}:{c_notation}"
    splice_status = {} if is_exon_cnv else SPLICEAI_STATUS_CACHE.get(variant_key, {})
    if splice_status.get("status") not in {None, "ok"}:
        lookup_diagnostics.append(
            f"SpliceAI unavailable: status={splice_status.get('status')}; "
            f"{splice_status.get('reason', 'no reason reported')}"
        )
    bayesdel_status = {} if is_exon_cnv else BAYESDEL_STATUS_CACHE.get(variant_key, {})
    if bayesdel_status.get("status") in {"api_error", "no_grch37_coords"}:
        lookup_diagnostics.append(
            f"MyVariant/BayesDel unavailable: status={bayesdel_status.get('status')}; "
            f"{bayesdel_status.get('reason', 'no reason reported')}"
        )

    # Step 4: fast local lookups
    gnomad_data = None
    if grch37 or grch38:
        gnomad_data = get_gnomad_frequencies(
            gene=gene,
            grch37=grch37,
            grch38=grch38,
        )

    table9_result  = table9_lookup_ps3_bs3(gene, c_notation)
    from backend.modules.ps1_splice_evidence import evaluate_defined_splice_sources
    ps1_vua_splice_evidence = evaluate_defined_splice_sources(
        gene, c_notation, table9_result
    )
    ps1_spliceai_score, ps1_spliceai_source = select_vua_spliceai_for_ps1(
        spliceai_score, table9_result
    )
    pp4_bp5_result = evaluate_pp4_bp5(gene, c_notation)
    ps1_result     = evaluate_ps1(
        gene, c_notation, p_notation,
        variant_type=variant_type,
        spliceai_score=ps1_spliceai_score,
        vua_spliceai_source=ps1_spliceai_source,
        vua_splice_evidence_status=ps1_vua_splice_evidence["status"],
        vua_splice_sources_checked=ps1_vua_splice_evidence["sources_checked"],
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
    # Detailed provider responses belong in the server log, not in the clinical
    # result. The public result receives one actionable summary per source.
    for diagnostic in lookup_diagnostics:
        logging.getLogger(__name__).warning("External lookup diagnostic: %s", diagnostic)
    if is_exon_cnv:
        result["warnings"].append(
            "Coordinate-dependent evidence was not evaluated: this exon-level "
            "copy-number variant has uncertain genomic breakpoints."
        )
    elif not grch37 and not grch38:
        result["warnings"].append(
            "Coordinate-dependent evidence was not evaluated because genomic "
            "coordinates could not be resolved."
        )
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
            "ClinVar comparison is temporarily unavailable."
        )
    if er.get("status") not in {"ok", "not_found"}:
        result["warnings"].append(
            "ClinGen ERepo comparison is temporarily unavailable."
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
        reference_transcript=TRANSCRIPTS.get(gene, ""),
        normalization_source=normalized_input.normalization_source,
        consequence_status=normalized_input.consequence_status,
        normalization_provenance=normalized_input.normalization_provenance or {},
        protein_consequence_explanation=normalized_input.protein_consequence_explanation,
        predicted_class=result["predicted_class"],
        predicted_label=CLASS_LABELS.get(result["predicted_class"], ""),
        total_points=result["total_points"],
        criteria=criteria,
        warnings=result["warnings"],
        external=ext_model,
        has_functional_evidence=result.get("has_functional_evidence", False),
        classification_note=result.get("classification_note", ""),
        evidence_direction=result.get("evidence_direction", "none"),
        mixed_evidence=result.get("mixed_evidence", False),
        pathogenic_points=result.get("pathogenic_points", 0),
        benign_points=result.get("benign_points", 0),
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
        protein_ps1_review=ProteinPs1ReviewRecommendation(**result["protein_ps1_review"])
        if result.get("protein_ps1_review") else None,
        initiation_review=RnaReviewRecommendation(**result["initiation_review"])
        if result.get("initiation_review") else None,
        spliceai_audit=SpliceAIAudit(**{
            field: splice_status.get(field)
            for field in SpliceAIAudit.model_fields
            if splice_status.get(field) is not None
        }) if splice_status else None,
        population_frequency_audit=(
            gnomad_data.get("population_frequency_audit", {})
            if gnomad_data else {}
        ),
        evidence_interactions=[
            EvidenceInteractionWarning(**warning)
            for warning in result.get("evidence_interactions", [])
        ],
    )


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
