# ============================================================
# ARIANE - FastAPI application
# Automated ACMG Rule-based Interpretation and Annotation ENgine
# ============================================================
from fastapi import FastAPI, HTTPException, Request
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, FileResponse
from pathlib import Path
import json
import os
import secrets
import time
import asyncio
from typing import Optional
from fastapi import Header

from backend.config import (
    DATA_DIR, TABLE4_PATH, TABLE9_PATH, ST7_PATH,
    ALLOWED_GENES, TRANSCRIPTS,
)
from backend.models import (
    VariantRequest, ClassificationResult, CriterionResult,
    ExternalComparison, ExternalSubmitter, CLASS_LABELS,
    BatchRequest, BatchResponse, BatchItemResult,
    AlphaMissenseResult, VusExplanation,
    RnaReviewRecommendation,
    ManualEvidenceRequest, ManualEvidenceResult,
    ManualCriterionResult,
)

# ── App setup ──────────────────────────────────────────────────────────────
app = FastAPI(
    title="ARIANE",
    description="Automated ACMG Rule-based Interpretation and Annotation ENgine for BRCA1/2",
    version="1.8.0",
)

FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"
app.mount("/static", StaticFiles(directory=FRONTEND_DIR / "static"), name="static")


@app.get("/", response_class=HTMLResponse)
async def index():
    return FileResponse(FRONTEND_DIR / "index.html")


@app.get("/api/health")
async def health():
    return {
        "status": "ok",
        "version": "1.8.0",
        "data": {
            "table4": TABLE4_PATH.exists(),
            "table9": TABLE9_PATH.exists(),
            "st7":    ST7_PATH.exists(),
        }
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


@app.post("/api/manual-evidence/evaluate")
async def evaluate_manual_evidence_endpoint(
    req: ManualEvidenceRequest,
) -> ManualEvidenceResult:
    from backend.modules.manual_evidence import evaluate_manual_evidence

    try:
        result = evaluate_manual_evidence(
            [criterion.model_dump() for criterion in req.base_criteria],
            [criterion.model_dump() for criterion in req.manual_criteria],
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    return ManualEvidenceResult(
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


# Semaphore limits concurrent external API calls during batch processing
BATCH_SEMAPHORE = asyncio.Semaphore(3)
EXTERNAL_LOOKUP_TIMEOUT = 12


async def _lookup_or_default(func, default, *args):
    try:
        return await asyncio.wait_for(
            asyncio.to_thread(func, *args),
            timeout=EXTERNAL_LOOKUP_TIMEOUT,
        )
    except Exception:
        return default


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

    # Step 1: resolve coordinates
    resolved = {}
    try:
        rv = resolve_variant(gene, c_notation)
        if rv:
            resolved[f"{gene}:{c_notation}"] = rv
    except Exception:
        pass
    grch37 = get_grch37(resolved, gene, c_notation)
    grch38 = get_grch38(resolved, gene, c_notation)

    # Step 2: variant type
    variant_type = infer_variant_type(c_notation, p_notation)

    # Step 3: parallel external lookups
    # get_bayesdel_and_alphamissense returns (bayesdel_score, alphamissense_dict)
    # in a single myvariant.info call - no extra API overhead for AlphaMissense.
    spliceai_score, (bayesdel_score, alphamissense), cv, er = await asyncio.gather(
        _lookup_or_default(get_spliceai_score, None, gene, c_notation),
        _lookup_or_default(get_bayesdel_and_alphamissense, (None, None), gene, c_notation),
        _lookup_or_default(
            clinvar_lookup,
            {"status": "api_timeout", "error": "ClinVar lookup timed out"},
            gene,
            c_notation,
        ),
        _lookup_or_default(
            clingen_erepo_lookup,
            {"status": "api_timeout", "error": "ClinGen ERepo lookup timed out"},
            gene,
            c_notation,
        ),
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
async def classify_variant(req: VariantRequest) -> ClassificationResult:
    return await _classify_one(
        req.gene, req.c_notation, req.p_notation or "", req.dup_type
    )


@app.post("/api/classify/batch")
async def classify_batch(req: BatchRequest) -> BatchResponse:
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
    return BatchResponse(
        total=len(items),
        success_count=success,
        error_count=len(items) - success,
        results=list(items),
    )


@app.post("/api/clear-cache")
async def clear_cache(x_ariane_admin_token: Optional[str] = Header(default=None)):
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

    return {"status": "ok", "message": "All caches cleared"}
