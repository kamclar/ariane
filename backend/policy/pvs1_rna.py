"""ENIGMA PVS1 RNA strength policy."""

from __future__ import annotations


_UNQUANTIFIED_PATIENT_RNA_STRENGTH = {
    "Very Strong": "Strong",
    "Strong": "Moderate",
    "Moderate": "Moderate",
    "Supporting": "Supporting",
}


def unquantified_patient_rna_strength(
    baseline_strength: str | None,
) -> str | None:
    return _UNQUANTIFIED_PATIENT_RNA_STRENGTH.get(baseline_strength or "")
