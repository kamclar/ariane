"""Automatic PP4/BP5 lookup from the validated local clinical-LR snapshot."""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Optional
from backend.gene_policy import (
    clinical_lr_thresholds,
    policy_name,
    policy_version,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_pp4_clinical_lr_snapshot.index.json"
METADATA_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_pp4_clinical_lr_snapshot.metadata.json"
SOURCE_MANIFEST_PATH = REPOSITORY_ROOT / "data" / "sources" / "enigma" / "clinical_lr_sources.manifest.json"
INDEL_SNAPSHOT_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_normalized_indel_snapshot.index.json"
INDEL_METADATA_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_normalized_indel_snapshot.metadata.json"

PP4_POINTS = {"Very Strong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1}
BP5_POINTS = {"Very Strong": -8, "Strong": -4, "Moderate": -2, "Supporting": -1}

_SNAPSHOT: dict[str, dict[str, Any]] | None = None
_ALIASES: dict[str, str] | None = None

_SOURCE_LABEL_PATTERN = re.compile(
    r"^(PP4|BP5)\s*-\s*(Pathogenic|Benign)\s*-\s*"
    r"(Supporting|Moderate|Strong|Very strong)$",
    re.IGNORECASE,
)
_CANONICAL_STRENGTHS = {
    "supporting": "Supporting",
    "moderate": "Moderate",
    "strong": "Strong",
    "very strong": "Very Strong",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def lr_to_pp4_strength(gene: str, lr: float) -> Optional[str]:
    thresholds = clinical_lr_thresholds(gene)["pp4"]
    if lr >= thresholds["very_strong_min_inclusive"]:
        return "Very Strong"
    if lr >= thresholds["strong_min_inclusive"]:
        return "Strong"
    if lr >= thresholds["moderate_min_inclusive"]:
        return "Moderate"
    if lr >= thresholds["supporting_min_inclusive"]:
        return "Supporting"
    return None


def lr_to_bp5_strength(gene: str, lr: float) -> Optional[str]:
    thresholds = clinical_lr_thresholds(gene)["bp5"]
    if lr <= thresholds["very_strong_max_inclusive"]:
        return "Very Strong"
    if lr <= thresholds["strong_max_inclusive"]:
        return "Strong"
    if lr <= thresholds["moderate_max_inclusive"]:
        return "Moderate"
    if lr <= thresholds["supporting_max_inclusive"]:
        return "Supporting"
    return None


def _parse_source_acmg_label(label: str) -> tuple[Optional[str], Optional[str]]:
    """Parse the pinned track label or reject an unsupported source schema."""
    value = str(label or "").strip()
    if value.casefold() == "not informative":
        return None, None
    match = _SOURCE_LABEL_PATTERN.fullmatch(value)
    if not match:
        raise RuntimeError(f"Unsupported PP4/BP5 source ACMG label: {value!r}")
    code, direction, strength = match.groups()
    expected_direction = "Pathogenic" if code.upper() == "PP4" else "Benign"
    if direction.casefold() != expected_direction.casefold():
        raise RuntimeError(
            f"PP4/BP5 source ACMG label has an invalid direction: {value!r}"
        )
    return code.upper(), _CANONICAL_STRENGTHS[strength.casefold()]


def _threshold_rule(
    gene: str,
    code: Optional[str],
    strength: Optional[str],
) -> dict[str, Any]:
    thresholds = clinical_lr_thresholds(gene)
    if code == "PP4" and strength:
        key = {
            "Supporting": "supporting_min_inclusive",
            "Moderate": "moderate_min_inclusive",
            "Strong": "strong_min_inclusive",
            "Very Strong": "very_strong_min_inclusive",
        }[strength]
        boundary = thresholds["pp4"][key]
        return {
            "operator": ">=",
            "boundary": boundary,
            "text": f"{code} {strength} applies at combined LR >= {boundary:g}",
        }
    if code == "BP5" and strength:
        key = {
            "Supporting": "supporting_max_inclusive",
            "Moderate": "moderate_max_inclusive",
            "Strong": "strong_max_inclusive",
            "Very Strong": "very_strong_max_inclusive",
        }[strength]
        boundary = thresholds["bp5"][key]
        return {
            "operator": "<=",
            "boundary": boundary,
            "text": f"{code} {strength} applies at combined LR <= {boundary:g}",
        }
    return {
        "operator": "between",
        "boundary": None,
        "text": (
            "PP4/BP5 is not informative when combined LR is greater than 0.48 "
            "and less than 2.08"
        ),
    }


def _threshold_comparison(
    gene: str,
    lr: float,
    calculated_code: Optional[str],
    calculated_strength: Optional[str],
    source_label: str,
) -> dict[str, Any]:
    source_code, source_strength = _parse_source_acmg_label(source_label)
    status = (
        "match"
        if (source_code, source_strength) == (calculated_code, calculated_strength)
        else "different"
    )
    vcep_label = (
        f"{calculated_code} {calculated_strength}"
        if calculated_code and calculated_strength
        else "Not informative"
    )
    threshold = _threshold_rule(gene, calculated_code, calculated_strength)
    policy_label = f"{policy_name(gene)} v{policy_version(gene)}"
    reason = ""
    if status == "different":
        reason = (
            f"The source track labels this record {source_label}, while {policy_label} "
            f"maps combined LR={lr:.6g} to {vcep_label}. {threshold['text']}. "
            "ARIANE uses the VCEP threshold result for classification and retains "
            "the source label for audit."
        )
    return {
        "status": status,
        "source_label": source_label,
        "source_code": source_code,
        "source_strength": source_strength,
        "vcep_policy": policy_label,
        "vcep_label": vcep_label,
        "vcep_code": calculated_code,
        "vcep_strength": calculated_strength,
        "combined_lr": lr,
        "threshold_operator": threshold["operator"],
        "threshold_value": threshold["boundary"],
        "threshold_rule": threshold["text"],
        "reason": reason,
    }


def load_pp4_bp5_snapshot() -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    """Load and validate the snapshot. Missing or corrupted data is fatal."""
    global _SNAPSHOT, _ALIASES
    if _SNAPSHOT is not None and _ALIASES is not None:
        return _SNAPSHOT, _ALIASES
    if not SNAPSHOT_PATH.is_file() or not METADATA_PATH.is_file():
        raise RuntimeError("PP4/BP5 clinical LR snapshot or its metadata is missing")

    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    if metadata.get("status") != "validated_derived_snapshot":
        raise RuntimeError("PP4/BP5 clinical LR snapshot is not validated")
    if metadata.get("index_sha256") != _sha256(SNAPSHOT_PATH):
        raise RuntimeError("PP4/BP5 clinical LR snapshot checksum does not match metadata")
    if not SOURCE_MANIFEST_PATH.is_file():
        raise RuntimeError("PP4/BP5 clinical LR source manifest is missing")
    if metadata.get("source_manifest_sha256") != _sha256(SOURCE_MANIFEST_PATH):
        raise RuntimeError("PP4/BP5 clinical LR source manifest checksum mismatch")
    source_manifest = json.loads(SOURCE_MANIFEST_PATH.read_text(encoding="utf-8"))
    if source_manifest.get("schema_version") != 3:
        raise RuntimeError("PP4/BP5 clinical LR source manifest schema is unsupported")
    datasets = source_manifest.get("datasets")
    if not isinstance(datasets, dict) or len(datasets) != 1:
        raise RuntimeError("PP4/BP5 clinical LR source manifest must pin one active dataset")
    if source_manifest.get("update_policy", {}).get("automatic_release_activation") is not False:
        raise RuntimeError("PP4/BP5 clinical LR source updates must not activate automatically")
    normalization = metadata.get("normalization")
    if not isinstance(normalization, dict) or not normalization.get("provenance"):
        raise RuntimeError("PP4/BP5 clinical LR snapshot normalization provenance is missing")
    indel_dependency = normalization.get("normalized_indel_dependency")
    if not isinstance(indel_dependency, dict):
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel dependency is missing")
    if not INDEL_SNAPSHOT_PATH.is_file() or not INDEL_METADATA_PATH.is_file():
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel dependency is unavailable")
    if indel_dependency.get("index_sha256") != _sha256(INDEL_SNAPSHOT_PATH):
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel dependency checksum mismatch")
    if indel_dependency.get("metadata_sha256") != _sha256(INDEL_METADATA_PATH):
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel metadata checksum mismatch")

    snapshot = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    if metadata.get("records") != len(snapshot):
        raise RuntimeError("PP4/BP5 clinical LR snapshot record count does not match metadata")
    for record in snapshot.values():
        _parse_source_acmg_label(record.get("source_acmg_label", ""))
        lr = record.get("combined_lr")
        lr_status = record.get("likelihood_ratio_status")
        expected_lr_status = (
            "unavailable_conflict"
            if lr is None
            else "source_reported_zero"
            if lr == 0
            else "available"
        )
        if lr_status != expected_lr_status:
            raise RuntimeError(
                "PP4/BP5 clinical LR snapshot has an inconsistent likelihood-ratio status"
            )
    lr_status_counts = Counter(
        record["likelihood_ratio_status"] for record in snapshot.values()
    )
    if metadata.get("likelihood_ratio_statuses") != dict(sorted(lr_status_counts.items())):
        raise RuntimeError("PP4/BP5 likelihood-ratio status counts do not match metadata")
    conflict_records = [
        record
        for record in snapshot.values()
        if record.get("overlap_status") == "conflicting_normalized_source_rows"
    ]
    conflict_source_row_count = sum(
        len(record.get("conflicting_source_rows", []))
        for record in conflict_records
    )
    expected_conflict_counts = {
        "variant_count": len(conflict_records),
        "source_row_count": conflict_source_row_count,
        "excess_source_row_count": conflict_source_row_count - len(conflict_records),
    }
    if metadata.get("normalization_conflicts") != expected_conflict_counts:
        raise RuntimeError("PP4/BP5 normalization conflict counts do not match metadata")

    aliases: dict[str, str] = {key: key for key in snapshot}
    ambiguous: set[str] = set()
    for canonical_key, record in snapshot.items():
        for notation in record.get("input_c_notations", []):
            alias = f"{record['gene']}:{notation}"
            previous = aliases.get(alias)
            if previous and previous != canonical_key:
                ambiguous.add(alias)
            else:
                aliases[alias] = canonical_key
    for alias in ambiguous:
        aliases.pop(alias, None)
    if ambiguous:
        raise RuntimeError(f"PP4/BP5 clinical LR snapshot has ambiguous aliases: {len(ambiguous)}")

    _SNAPSHOT, _ALIASES = snapshot, aliases
    return snapshot, aliases


def evaluate_pp4_bp5(gene: str, c_notation: str) -> Dict:
    snapshot, aliases = load_pp4_bp5_snapshot()
    query_key = f"{gene}:{c_notation}"
    canonical_key = query_key if query_key in snapshot else aliases.get(query_key)
    entry = snapshot.get(canonical_key) if canonical_key else None
    result = {
        "applies": False, "code": None, "strength": None, "points": 0,
        "reason": "", "likelihood_ratio": None,
        "likelihood_ratio_status": "not_found",
        "source_components": [],
        "source_bundle_ids": [],
        "source_bundle_count": 0,
        "independent_source_group_count": 0,
        "candidate_likelihood_ratio": None,
        "overlap_status": "not_assessed",
        "double_counting_risk": False,
        "source_reported_overlap_caveat": False,
        "data_release": "",
        "overlap_assessment_note": "",
        "overlap_assessment_sources": [],
        "automatic_combination_allowed": False,
        "application_status": "not_found",
        "clinical_evidence_types": [],
        "distinct_clinical_evidence_type_count": 0,
        "likelihood_ratio_contribution_count": 0,
        "independent_evidence_contribution_count": 0,
        "single_strong_likely_benign_eligible": False,
        "single_strong_likely_benign_basis": "",
        "threshold_comparison": {"status": "not_available"},
    }
    if entry is None:
        result["reason"] = (
            "No informative variant-specific combined clinical LR is available; "
            "PP4/BP5 is not applied under ENIGMA v1.2"
        )
        return result

    result.update({
        "source_bundle_ids": entry.get("source_bundle_ids", []),
        "source_bundle_count": entry.get("source_bundle_count", 0),
        "independent_source_group_count": entry.get(
            "independent_source_group_count", 0
        ),
        "candidate_likelihood_ratio": entry.get("candidate_combined_lr"),
        "overlap_status": entry.get("overlap_status", "unknown"),
        "double_counting_risk": entry.get("double_counting_risk", True),
        "source_reported_overlap_caveat": entry.get(
            "source_reported_overlap_caveat", False
        ),
        "data_release": entry.get("source", {}).get("track_release", ""),
        "overlap_assessment_note": entry.get("assessment_note", ""),
        "overlap_assessment_sources": entry.get("assessment_sources", []),
        "automatic_combination_allowed": entry.get(
            "automatic_combination_allowed", False
        ),
        "application_status": entry.get(
            "automatic_application_status", "review_required"
        ),
        "likelihood_ratio_status": entry["likelihood_ratio_status"],
        "threshold_comparison": {
            "status": "not_assessed",
            "source_label": entry.get("source_acmg_label", ""),
        },
    })
    lr = entry.get("combined_lr")
    result["likelihood_ratio"] = lr
    result["source_components"] = entry.get("source_components", [])
    clinical_data = [
        item
        for component in result["source_components"]
        for item in component.get("clinical_data", [])
        if isinstance(item, dict) and item.get("lr") is not None
    ]
    clinical_evidence_types = sorted({
        str(item.get("data_type") or "").strip()
        for item in clinical_data
        if str(item.get("data_type") or "").strip()
    })
    result["clinical_evidence_types"] = clinical_evidence_types
    result["distinct_clinical_evidence_type_count"] = len(clinical_evidence_types)
    result["likelihood_ratio_contribution_count"] = len(clinical_data)
    result["independent_evidence_contribution_count"] = result[
        "independent_source_group_count"
    ]
    if not result["automatic_combination_allowed"]:
        result["reason"] = result["overlap_assessment_note"] or (
            "The published clinical LR record requires expert source review, so "
            "ARIANE did not apply PP4/BP5 automatically."
        )
        return result
    if lr is None:
        raise RuntimeError(
            "Clinical LR record permits automatic combination but has no combined LR"
        )
    pp4_strength = lr_to_pp4_strength(gene, lr)
    bp5_strength = lr_to_bp5_strength(gene, lr)
    code = "PP4" if pp4_strength else "BP5" if bp5_strength else None
    strength = pp4_strength or bp5_strength
    result["threshold_comparison"] = _threshold_comparison(
        gene,
        lr,
        code,
        strength,
        entry["source_acmg_label"],
    )
    if not code:
        result["reason"] = f"Combined clinical LR={lr:.6g} is not informative for PP4 or BP5"
        return result

    result.update({
        "applies": True,
        "code": code,
        "strength": strength,
        "points": (PP4_POINTS if code == "PP4" else BP5_POINTS)[strength],
        "application_status": "applied",
    })
    if (
        code == "BP5"
        and strength == "Strong"
        and len(clinical_data) >= 2
        and len(clinical_evidence_types) >= 2
    ):
        result["single_strong_likely_benign_eligible"] = True
        result["single_strong_likely_benign_basis"] = (
            "Multiple recorded clinical evidence types and likelihood-ratio "
            "contributions support the BP5 Strong code"
        )
    pmids = sorted({
        str(pmid)
        for component in result["source_components"]
        for pmid in (
            component.get("pmids")
            or ([component.get("pmid")] if component.get("pmid") else [])
        )
        if str(pmid).strip()
    })
    result["reason"] = (
        f"ENIGMA v1.2 combined clinical evidence from the UCSC ENIGMA "
        f"data release {result['data_release']}: "
        f"combined LR={lr:.6g}; {code} {strength}; "
        f"PMID {', '.join(pmids)}"
    )
    if result["threshold_comparison"]["status"] == "different":
        result["reason"] += ". " + result["threshold_comparison"]["reason"]
    return result
