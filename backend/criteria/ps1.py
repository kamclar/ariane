"""ENIGMA v1.2 protein-level PS1 evaluation.

ST7 P/LP records are accepted as the reference classification basis. The
versioned registry records the independent protein-mechanism and splice checks
that remain mandatory before a reference can be used for automatic scoring.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional

from backend.policy.criterion_points import points_for_strength
from backend.policy.gene import (
    spliceai_thresholds,
    vcep_specification,
)
from backend.reference_data.ps1_registry_validation import (
    NO_DAMAGING_SPLICE_STATUSES,
    extract_missense_aa_change as _extract_aa_change,
)
from backend.reference_data.ps1_protein import (
    Ps1ReferenceMatches,
    matching_ps1_references,
    public_ps1_reference,
)


def select_vua_spliceai_for_ps1(
    service_or_cache_score: Optional[float],
    table9_result: Optional[Dict[str, Any]],
) -> tuple[Optional[float], str]:
    """Return the configured SpliceAI result used for protein PS1.

    Table 9 is used separately for curated PS3/BS3 and published splice
    evidence. Its recorded prediction does not replace the configured result.
    The second argument remains in the public helper signature for callers that
    already provide the reviewed Table 9 record.
    """
    return service_or_cache_score, "configured SpliceAI source"
def _deduplicate_candidates(records: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    deduplicated: Dict[tuple, Dict[str, Any]] = {}
    for record in records:
        key = (record.get("gene"), record.get("c_notation"))
        current = deduplicated.get(key)
        if current is None or record.get("reference_status") == "approved":
            deduplicated[key] = record
    return list(deduplicated.values())
def evaluate_ps1(
    gene: str,
    c_notation: str,
    p_notation: str,
    variant_type: str,
    spliceai_score: Optional[float] = None,
    vua_spliceai_source: str = "configured SpliceAI source",
    vua_splice_evidence_status: str = "not_assessed",
    vua_splice_sources_checked: Optional[List[str]] = None,
    reference_spliceai_scores: Optional[Dict[str, Optional[float]]] = None,
    reference_matches: Ps1ReferenceMatches | None = None,
) -> Dict[str, Any]:
    """Apply PS1 only after every recorded and runtime condition is satisfied."""
    result: Dict[str, Any] = {
        "applies": False,
        "code": "PS1",
        "strength": None,
        "points": 0,
        "reason": "",
        "reference_variant": None,
        "reference_status": "not_found",
        "application_status": "not_applied",
        "review_required": False,
        "blocking_reasons": [],
        "candidates": [],
        "vua_splice_evidence_status": vua_splice_evidence_status,
        "vua_splice_sources_checked": vua_splice_sources_checked or [],
        "vua_spliceai_score": spliceai_score,
        "vua_spliceai_source": vua_spliceai_source,
        "reference_spliceai_scores": reference_spliceai_scores or {},
        "source_url": vcep_specification(gene)["url"],
        "assessed_variant": {
            "gene": gene,
            "c_notation": c_notation,
            "p_notation": p_notation,
        },
    }
    splice_low = spliceai_thresholds(gene)["bp4"]

    if (variant_type or "").lower() != "missense":
        result["reason"] = "Protein-level PS1 applies only to missense variants"
        return result
    aa_key = _extract_aa_change(p_notation)
    if not aa_key:
        result["reason"] = f"Could not extract a normalized missense change from {p_notation}"
        return result

    matches = reference_matches or matching_ps1_references(
        gene,
        c_notation,
        p_notation,
    )
    st7_matches = list(matches.candidates)
    approved_records = list(matches.approved)
    candidates = _deduplicate_candidates(
        [*st7_matches, *(public_ps1_reference(item) for item in approved_records)]
    )
    result["candidates"] = candidates
    if not candidates:
        result["reason"] = f"No P/LP reference with the same missense change {aa_key} was found"
        return result

    review_candidates = [
        item for item in candidates if item.get("reference_status") == "review_required"
    ]
    excluded_candidates = [
        item for item in candidates if item.get("reference_status") == "excluded"
    ]
    result["reference_status"] = (
        "approved"
        if approved_records
        else "review_required"
        if review_candidates
        else "excluded"
    )

    if spliceai_score is None:
        result["blocking_reasons"].append(
            "SpliceAI is unavailable for the variant under assessment"
        )
    elif spliceai_score > splice_low:
        result["reason"] = (
            f"Protein-level PS1 not applicable: assessed-variant SpliceAI "
            f"{spliceai_score:.3f} > {splice_low}"
        )
        result["reference_status"] = "ineligible_for_this_application"
        result["application_status"] = "not_applicable"
        return result

    if vua_splice_evidence_status == "abnormal":
        result["reason"] = (
            "Protein-level PS1 not applicable: a confirmed splice effect is recorded "
            "for the variant under assessment"
        )
        result["reference_status"] = "ineligible_for_this_application"
        result["application_status"] = "not_applicable"
        return result
    if vua_splice_evidence_status == "conflicting":
        result["blocking_reasons"].append(
            "Confirmed splice evidence for the variant under assessment is conflicting"
        )
    elif vua_splice_evidence_status not in NO_DAMAGING_SPLICE_STATUSES:
        result["blocking_reasons"].append(
            "The defined confirmed RNA/splice evidence sources have not been completely checked"
        )

    if not approved_records and excluded_candidates and not review_candidates:
        result["application_status"] = "reference_ineligible"
        result["reason"] = "; ".join(
            item.get("status_reason") or "Matching ST7 reference is ineligible"
            for item in excluded_candidates
        )
        return result

    runtime_scores = reference_spliceai_scores or {}
    if not approved_records and review_candidates:
        review_unavailable = [
            str(item["c_notation"])
            for item in review_candidates
            if runtime_scores.get(str(item["c_notation"])) is None
        ]
        review_high = [
            (str(item["c_notation"]), float(runtime_scores[str(item["c_notation"])]))
            for item in review_candidates
            if runtime_scores.get(str(item["c_notation"])) is not None
            and float(runtime_scores[str(item["c_notation"])]) > splice_low
        ]
        if len(review_high) == len(review_candidates):
            result["reference_status"] = "ineligible_for_this_application"
            result["application_status"] = "not_applicable"
            result["reason"] = (
                "Protein-level PS1 not applicable: every matching ST7 reference "
                f"has SpliceAI above {splice_low} ("
                + ", ".join(
                    f"{reference_c}={score:.3f}"
                    for reference_c, score in review_high
                )
                + ")"
            )
            return result
        result["blocking_reasons"].append(
            "A matching ST7 reference was found, but its recorded protein or "
            "splice evidence remains unresolved"
        )
        if review_unavailable:
            result["blocking_reasons"].append(
                "SpliceAI is unavailable for matching reference variant(s): "
                + ", ".join(review_unavailable)
            )

    runtime_approved_records: List[Dict[str, Any]] = []
    unavailable_references: List[str] = []
    high_score_references: List[tuple[str, float]] = []
    for reference in approved_records:
        reference_c = str(reference["c_notation"])
        reference_score = runtime_scores.get(reference_c)
        if reference_score is None:
            unavailable_references.append(reference_c)
        elif reference_score > splice_low:
            high_score_references.append((reference_c, reference_score))
        else:
            runtime_approved_records.append(reference)

    if approved_records and unavailable_references:
        result["blocking_reasons"].append(
            "SpliceAI is unavailable for matching reference variant(s): "
            + ", ".join(unavailable_references)
        )

    if result["blocking_reasons"]:
        result["review_required"] = True
        result["application_status"] = "manual_review_required"
        result["reason"] = "; ".join(result["blocking_reasons"])
        return result

    if approved_records and not runtime_approved_records:
        result["reference_status"] = "ineligible_for_this_application"
        result["application_status"] = "not_applicable"
        result["reason"] = (
            "Protein-level PS1 not applicable: matching reference variant "
            f"SpliceAI exceeds {splice_low} ("
            + ", ".join(
                f"{reference_c}={score:.3f}"
                for reference_c, score in high_score_references
            )
            + ")"
        )
        return result

    best = max(
        runtime_approved_records,
        key=lambda item: 5 if item["classification"] == "Pathogenic" else 4,
    )
    best_reference_score = runtime_scores[str(best["c_notation"])]
    strength = "Strong" if best["classification"] == "Pathogenic" else "Moderate"
    result.update(
        {
            "applies": True,
            "strength": strength,
            "points": points_for_strength(strength),
            "reference_variant": public_ps1_reference(best),
            "reference_status": "approved",
            "application_status": "auto_applied",
            "reason": (
                f"Same normalized missense change {aa_key} as approved "
                f"{best['classification']} reference {best['c_notation']} "
                f"{best['p_notation']} from "
                f"{best.get('candidate_source', 'the curated PS1 reference registry')} "
                f"(classification source: {best['classification_source']}); both variants "
                "meet the recorded ENIGMA protein-level PS1 splice conditions "
                f"(assessed-variant SpliceAI {spliceai_score:.3f}, reference-variant "
                f"SpliceAI {best_reference_score:.3f}, source {vua_spliceai_source})"
            ),
        }
    )
    return result
