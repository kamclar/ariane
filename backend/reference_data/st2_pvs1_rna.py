"""ST2/ST3 evidence prepared for PVS1 RNA evaluation."""

from __future__ import annotations

import re
from typing import Any

from backend.policy.pvs1_rna import unquantified_patient_rna_strength
from backend.reference_data.ps1_splice_evidence import (
    get_st2_source_metadata,
    get_st2_splice_record,
)
from backend.reference_data.table4 import (
    Table4LookupStatus,
    table4_deletion_exons,
    table4_lookup_deletion,
)


APPENDIX_URL = (
    "https://cspec.genome.network/cspec/File/id/9e6119dc-90b9-42b5-a3b7-1a2eb28b1b12/data"
)

_PATIENT_UNQUANTIFIED_LOF = (
    "patient not allele-specific; aberrant transcripts consistent with loss of function"
)
_DAMAGING_RNA_RESULT = "aberrant transcripts consistent with loss of function"


def normalized_st2_text(value: Any) -> str:
    return " ".join(str(value or "").split()).strip().lower()


def is_damaging_st2_assay(value: Any) -> bool:
    return _DAMAGING_RNA_RESULT in normalized_st2_text(value)


def is_unquantified_patient_st2_assay(value: Any) -> bool:
    return normalized_st2_text(value) == _PATIENT_UNQUANTIFIED_LOF


def reported_whole_exon_deletion(value: Any) -> int | None:
    match = re.fullmatch(r"exon\s+(\d+)\s+deletion", normalized_st2_text(value))
    return int(match.group(1)) if match else None


def table4_exon_for_reported_number(gene: str, exon_number: int) -> str | None:
    matches = []
    for exon in table4_deletion_exons(gene):
        match = re.fullmatch(r"E(\d+)(?:\((\d+)\))?", exon)
        if not match:
            continue
        if int(match.group(2) or match.group(1)) == exon_number:
            matches.append(exon)
    return matches[0] if len(matches) == 1 else None


def lookup_st2_pvs1_rna_evidence(gene: str, c_notation: str) -> dict[str, Any]:
    provenance = get_st2_source_metadata()
    record = get_st2_splice_record(gene, c_notation)
    base: dict[str, Any] = {
        "status": "not_in_source",
        "record": record,
        "table4_exon": None,
        "table4_baseline": None,
        "source_id": "enigma-st2-st3-v1.2",
        "source_version": str(provenance.get("version") or ""),
        "source_checksum": str(provenance.get("source_file_sha256") or ""),
        "source_url": str(provenance.get("source_url") or APPENDIX_URL),
        "reason": ("The exact variant has no record in ENIGMA Supplementary Table 2 v1.2."),
    }
    if record is None:
        return base
    assay_category = record.get("splicing_assay_result_category")
    if not is_damaging_st2_assay(assay_category):
        return {
            **base,
            "status": "not_applicable",
            "reason": (
                "The exact ST2 record does not report an mRNA result in the "
                "damaging loss-of-function category."
            ),
        }
    if not is_unquantified_patient_st2_assay(assay_category):
        return {
            **base,
            "status": "review_required",
            "reason": (
                "The exact ST2 record contains damaging RNA evidence but uses "
                "an assay category that requires a different Appendix E review."
            ),
        }
    if not record.get("st3_references"):
        return {
            **base,
            "status": "review_required",
            "reason": "The ST2 result has no linked supporting ST3 source row.",
        }
    reported_exon = reported_whole_exon_deletion(record.get("result"))
    if reported_exon is None:
        return {
            **base,
            "status": "review_required",
            "reason": (
                "The ST2 transcript result is complex or partial and cannot be "
                "mapped unambiguously to one Table 4 deletion row."
            ),
        }
    table4_exon = table4_exon_for_reported_number(gene, reported_exon)
    if table4_exon is None:
        return {
            **base,
            "status": "review_required",
            "reason": (
                f"Reported exon {reported_exon} does not map uniquely to a Table 4 deletion row."
            ),
        }
    baseline = table4_lookup_deletion(gene, table4_exon)
    if (
        baseline.get("status") != Table4LookupStatus.MATCHED
        or unquantified_patient_rna_strength(baseline.get("pvs1_strength")) is None
    ):
        return {
            **base,
            "status": "review_required",
            "table4_exon": table4_exon,
            "table4_baseline": baseline,
            "reason": (
                f"The {gene} {table4_exon} transcript deletion has no applicable "
                "baseline PVS1 weight in Table 4."
            ),
        }
    return {
        **base,
        "status": "eligible",
        "table4_exon": table4_exon,
        "table4_baseline": baseline,
        "assay_interpretation": ("enigma_st2_curated_apparent_near_complete_loss_of_function"),
        "reason": (
            "The exact checksum-bound ST2 row contains ENIGMA-curated patient "
            "mRNA evidence, a loss-of-function transcript category, linked ST3 "
            "sources and an unambiguous Table 4 consequence."
        ),
    }
