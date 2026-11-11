"""Validated protein-PS1 reference records."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from typing import Any, Mapping

from backend.reference_data.paths import PS1_PROTEIN_REGISTRY_PATH, ST7_PATH
from backend.reference_data.ps1_registry_validation import (
    extract_missense_aa_change,
    validate_ps1_reference_registry,
)


@dataclass(frozen=True)
class Ps1ReferenceMatches:
    candidates: tuple[dict[str, Any], ...]
    approved: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class _Ps1ReferenceIndex:
    candidates: Mapping[str, Mapping[str, tuple[dict[str, Any], ...]]]
    approved: Mapping[str, Mapping[str, tuple[dict[str, Any], ...]]]


def public_ps1_reference(record: Mapping[str, Any]) -> dict[str, Any]:
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
        "reference_splice_evidence_status": record.get("reference_splice_evidence", {}).get(
            "confirmed_status", "not_assessed"
        ),
        "reference_splice_sources_checked": record.get("reference_splice_evidence", {}).get(
            "sources_checked", []
        ),
        "classification_ps1_dependency_used": record.get("classification_ps1_dependency", {}).get(
            "used", "unknown"
        ),
    }


def _group_by_change(
    records: list[dict[str, Any]],
    *,
    public: bool,
) -> dict[str, dict[str, tuple[dict[str, Any], ...]]]:
    grouped: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for record in records:
        aa_change = extract_missense_aa_change(str(record.get("p_notation") or ""))
        if not aa_change:
            raise RuntimeError("PS1 reference registry contains an invalid missense record")
        value = public_ps1_reference(record) if public else record
        grouped.setdefault(str(record["gene"]), {}).setdefault(aa_change, []).append(value)
    return {
        gene: {change: tuple(values) for change, values in changes.items()}
        for gene, changes in grouped.items()
    }


@lru_cache(maxsize=1)
def _reference_index() -> _Ps1ReferenceIndex:
    if not ST7_PATH.is_file():
        raise RuntimeError(f"Required ENIGMA ST7 dataset is missing: {ST7_PATH}")
    if not PS1_PROTEIN_REGISTRY_PATH.is_file():
        raise RuntimeError(
            f"Required PS1 protein reference registry is missing: {PS1_PROTEIN_REGISTRY_PATH}"
        )

    st7 = json.loads(ST7_PATH.read_text(encoding="utf-8"))
    registry = json.loads(PS1_PROTEIN_REGISTRY_PATH.read_text(encoding="utf-8"))
    validate_ps1_reference_registry(registry)

    expected_st7 = {
        (record.get("gene"), record.get("c_notation")): record
        for record in st7.get("variants", [])
        if record.get("iarc_class") in {4, 5}
        and extract_missense_aa_change(str(record.get("p_notation") or ""))
    }
    records = list(registry.get("references", []))
    registry_by_variant = {
        (record.get("gene"), record.get("c_notation")): record for record in records
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

    approved = [record for record in records if record.get("status") == "eligible"]
    return _Ps1ReferenceIndex(
        candidates=_group_by_change(records, public=True),
        approved=_group_by_change(approved, public=False),
    )


def matching_ps1_references(
    gene: str,
    c_notation: str,
    p_notation: str,
) -> Ps1ReferenceMatches:
    aa_change = extract_missense_aa_change(p_notation)
    if not aa_change:
        return Ps1ReferenceMatches((), ())
    index = _reference_index()
    return Ps1ReferenceMatches(
        candidates=tuple(
            record
            for record in index.candidates.get(gene, {}).get(aa_change, ())
            if record.get("c_notation") != c_notation
        ),
        approved=tuple(
            record
            for record in index.approved.get(gene, {}).get(aa_change, ())
            if record.get("c_notation") != c_notation
        ),
    )


def discover_ps1_reference_variants(
    gene: str,
    c_notation: str,
    p_notation: str,
    variant_type: str,
) -> list[str]:
    if (variant_type or "").lower() != "missense":
        return []
    matches = matching_ps1_references(gene, c_notation, p_notation)
    return sorted(
        {
            str(record["c_notation"])
            for record in matches.candidates
            if record.get("reference_status") != "excluded"
        }
    )


def lookup_ps1_reference_variant(
    gene: str,
    c_notation: str,
) -> dict[str, Any] | None:
    for records_by_change in _reference_index().candidates.get(gene, {}).values():
        for record in records_by_change:
            if record.get("c_notation") == c_notation:
                return dict(record)
    return None
