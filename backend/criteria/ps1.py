"""ENIGMA v1.2 protein-level PS1 evaluation.

ST7 P/LP records are accepted as the reference classification basis. The
versioned registry records the independent protein-mechanism and splice checks
that remain mandatory before a reference can be used for automatic scoring.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from backend.policy.gene import (
    spliceai_thresholds,
    vcep_specification,
)
from backend.reference_data.ps1_registry_validation import (
    NO_DAMAGING_SPLICE_STATUSES,
    compute_approval_basis_checksum as compute_approval_basis_checksum,
    extract_missense_aa_change as _extract_aa_change,
    validate_ps1_reference_registry,
)


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
ST7_PATH = DATA_DIR / "st7_reference_set.json"
PS1_REGISTRY_PATH = DATA_DIR / "ps1_protein_reference_registry.json"

PS1_POINTS = {"Strong": 4, "Moderate": 2}
_ST7_LOOKUP: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
_APPROVED_LOOKUP: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
_LOADED = False


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




def _load_references() -> None:
    global _LOADED, _ST7_LOOKUP, _APPROVED_LOOKUP
    if _LOADED:
        return
    if not ST7_PATH.is_file():
        raise RuntimeError(f"Required ENIGMA ST7 dataset is missing: {ST7_PATH}")
    if not PS1_REGISTRY_PATH.is_file():
        raise RuntimeError(f"Required PS1 protein reference registry is missing: {PS1_REGISTRY_PATH}")

    st7 = json.loads(ST7_PATH.read_text(encoding="utf-8"))
    registry = json.loads(PS1_REGISTRY_PATH.read_text(encoding="utf-8"))
    validate_ps1_reference_registry(registry)

    expected_st7 = {
        (record.get("gene"), record.get("c_notation")): record
        for record in st7.get("variants", [])
        if record.get("iarc_class") in {4, 5}
        and _extract_aa_change(record.get("p_notation") or "")
    }
    registry_records = registry.get("references", [])
    registry_by_variant = {
        (record.get("gene"), record.get("c_notation")): record
        for record in registry_records
    }
    st7_registry = {
        key: record
        for key, record in registry_by_variant.items()
        if "enigma_st7_v1_2_reference_set" in record.get("source_memberships", [])
    }
    if set(st7_registry) != set(expected_st7):
        raise RuntimeError(
            "PS1 protein reference registry does not contain the complete ST7 P/LP missense set"
        )
    for key, source in expected_st7.items():
        record = st7_registry[key]
        expected_class = "Pathogenic" if source["iarc_class"] == 5 else "Likely Pathogenic"
        if (
            record.get("p_notation") != source.get("p_notation")
            or record.get("st7_source_classification") != expected_class
            or record.get("st7_source_classification_source")
            != (source.get("source") or "ENIGMA ST7 v1.2")
        ):
            raise RuntimeError(f"PS1 registry record does not match ST7 for {key[0]}:{key[1]}")

    st7_lookup: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    for record in registry_records:
        aa_key = _extract_aa_change(record["p_notation"])
        st7_lookup.setdefault(record["gene"], {}).setdefault(aa_key, []).append(
            _public_registry_candidate(record)
        )

    approved_lookup: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    for record in registry_records:
        if record.get("status") != "eligible":
            continue
        aa_key = _extract_aa_change(record["p_notation"])
        approved_lookup.setdefault(record["gene"], {}).setdefault(aa_key, []).append(record)

    _ST7_LOOKUP = st7_lookup
    _APPROVED_LOOKUP = approved_lookup
    _LOADED = True


def reset_ps1_reference_cache_for_tests() -> None:
    global _LOADED, _ST7_LOOKUP, _APPROVED_LOOKUP
    _LOADED = False
    _ST7_LOOKUP = {}
    _APPROVED_LOOKUP = {}


def discover_ps1_reference_variants(
    gene: str,
    c_notation: str,
    p_notation: str,
    variant_type: str,
) -> List[str]:
    """Return candidate reference c. variants that require runtime SpliceAI.

    Discovery is local and deterministic. It does not apply PS1 and it does not
    call an external service. The caller can therefore schedule reference and
    assessed-variant SpliceAI lookups together.
    """
    if (variant_type or "").lower() != "missense":
        return []
    aa_key = _extract_aa_change(p_notation)
    if not aa_key:
        return []
    _load_references()
    return sorted(
        {
            str(record["c_notation"])
            for record in _ST7_LOOKUP.get(gene, {}).get(aa_key, [])
            if record.get("reference_status") != "excluded"
            if record.get("c_notation") != c_notation
        }
    )


def _deduplicate_candidates(records: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    deduplicated: Dict[tuple, Dict[str, Any]] = {}
    for record in records:
        key = (record.get("gene"), record.get("c_notation"))
        current = deduplicated.get(key)
        if current is None or record.get("reference_status") == "approved":
            deduplicated[key] = record
    return list(deduplicated.values())


def _public_registry_candidate(record: Dict[str, Any]) -> Dict[str, Any]:
    registry_status = record.get("status", "review_required")
    return {
        "key": f"REGISTRY|{record['reference_id']}",
        "reference_id": record["reference_id"],
        "gene": record["gene"],
        "transcript": record["transcript"],
        "c_notation": record["c_notation"],
        "p_notation": record["p_notation"],
        "classification": record["classification"],
        "iarc_class": 5 if record["classification"] == "Pathogenic" else 4,
        "classification_basis": record["classification_verification"],
        "classification_source": record["classification_source"],
        "reference_status": "approved" if registry_status == "eligible" else registry_status,
        "status_reason": record.get("status_reason", ""),
        "source_dataset": record.get("candidate_source", "Curated PS1 protein reference registry"),
        "reference_splice_evidence_status": (
            record.get("reference_splice_evidence", {}).get("confirmed_status", "not_assessed")
        ),
        "reference_splice_sources_checked": (
            record.get("reference_splice_evidence", {}).get("sources_checked", [])
        ),
        "classification_ps1_dependency_used": (
            record.get("classification_ps1_dependency", {}).get("used", "unknown")
        ),
    }


def lookup_ps1_reference_variant(gene: str, c_notation: str) -> Optional[Dict[str, Any]]:
    """Return one exact registry reference for guided manual PS1 review."""
    _load_references()
    for records_by_change in _ST7_LOOKUP.get(gene, {}).values():
        for record in records_by_change:
            if record.get("c_notation") == c_notation:
                return dict(record)
    return None


def _public_approved_candidate(record: Dict[str, Any]) -> Dict[str, Any]:
    return _public_registry_candidate(record)


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
) -> Dict[str, Any]:
    """Apply PS1 only after every recorded and runtime condition is satisfied."""
    _load_references()
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

    st7_matches = [
        item for item in _ST7_LOOKUP.get(gene, {}).get(aa_key, [])
        if item["c_notation"] != c_notation
    ]
    approved_records = [
        item for item in _APPROVED_LOOKUP.get(gene, {}).get(aa_key, [])
        if item["c_notation"] != c_notation
    ]
    candidates = _deduplicate_candidates(
        [*st7_matches, *(_public_approved_candidate(item) for item in approved_records)]
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
            "points": PS1_POINTS[strength],
            "reference_variant": _public_approved_candidate(best),
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
