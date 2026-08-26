"""Automatic ENIGMA PVS1 (RNA) from structured, official RNA evidence.

This module contains no variant allowlist. It applies the qualitative patient
mRNA branch of ENIGMA v1.2 Appendix E Table 9 to exact records from
Supplementary Table 2, then obtains the baseline loss-of-function weight from
Specifications Table 4.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Optional

from backend.modules.ps1_splice_evidence import get_st2_splice_record
from backend.modules.table4 import TABLE4_DATA, table4_lookup_deletion
from backend.gene_policy import policy_name, policy_version, vcep_specification


APPENDIX_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "9e6119dc-90b9-42b5-a3b7-1a2eb28b1b12/data"
)

_PATIENT_UNQUANTIFIED_LOF = (
    "patient not allele-specific; aberrant transcripts consistent with loss of function"
)
_RNA_DOWNWEIGHT = {
    "Very Strong": ("Strong", 4),
    "Strong": ("Moderate", 2),
    "Moderate": ("Supporting", 1),
    "Supporting": ("Supporting", 1),
}


def _normalized_text(value: Any) -> str:
    return " ".join(str(value or "").split()).strip().lower()


def _reported_whole_exon_deletion(value: Any) -> Optional[int]:
    """Parse only an unambiguous whole-exon deletion transcript result."""
    match = re.fullmatch(r"exon\s+(\d+)\s+deletion", _normalized_text(value))
    return int(match.group(1)) if match else None


def _table4_exon_for_reported_number(gene: str, exon_number: int) -> Optional[str]:
    """Map clinical exon numbering in ST2 to the corresponding Table 4 row."""
    matches = []
    for exon in TABLE4_DATA.get("deletion_rules", {}).get(gene, {}):
        match = re.fullmatch(r"E(\d+)(?:\((\d+)\))?", exon)
        if not match:
            continue
        reported_number = int(match.group(2) or match.group(1))
        if reported_number == exon_number:
            matches.append(exon)
    return matches[0] if len(matches) == 1 else None


def evaluate_pvs1_rna(gene: str, c_notation: str) -> Dict[str, Any]:
    """Apply the ENIGMA qualitative patient-mRNA PVS1 (RNA) branch.

    Automatic scoring is intentionally limited to information explicitly
    represented by the official datasets:

    * an exact ST2 variant row;
    * patient mRNA without allele-specific quantitation;
    * an abnormal transcript categorised by ENIGMA as loss-of-function;
    * one unambiguous whole-exon deletion outcome;
    * a matching Table 4 deletion row with an applicable baseline PVS1 weight.

    Other RNA records remain unscored because ST2 does not retain the numerical
    transcript proportions needed by the quantitative Appendix E branches.
    """
    specification = vcep_specification(gene)
    result: Dict[str, Any] = {
        "applies": False,
        "code": "PVS1_RNA",
        "strength": None,
        "points": 0,
        "reason": "",
        "source": specification["url"],
        "source_record": None,
        "table4_exon": None,
        "appendix_branch": None,
    }

    record = get_st2_splice_record(gene, c_notation)
    if record is None:
        result["reason"] = (
            "PVS1 (RNA) was not applied because the exact variant has no record "
            "in ENIGMA Supplementary Table 2 v1.2."
        )
        return result

    result["source_record"] = {
        "dataset": "ENIGMA Supplementary Table 2 v1.2",
        "source_row": record.get("source_row"),
        "assay_category": record.get("splicing_assay_result_category"),
        "transcript_result": record.get("result"),
    }

    assay_category = _normalized_text(record.get("splicing_assay_result_category"))
    if assay_category != _PATIENT_UNQUANTIFIED_LOF:
        result["reason"] = (
            "PVS1 (RNA) was not applied automatically because the exact ENIGMA "
            "ST2 assay category does not enter the qualitative patient-mRNA "
            "loss-of-function branch. Quantitative or curator assessment is required."
        )
        return result

    reported_exon = _reported_whole_exon_deletion(record.get("result"))
    if reported_exon is None:
        result["reason"] = (
            "PVS1 (RNA) was not applied automatically because ENIGMA ST2 reports "
            "a complex or partial transcript consequence that cannot be mapped "
            "unambiguously to one Table 4 deletion row."
        )
        return result

    table4_exon = _table4_exon_for_reported_number(gene, reported_exon)
    result["table4_exon"] = table4_exon
    if table4_exon is None:
        result["reason"] = (
            f"PVS1 (RNA) was not applied because reported exon {reported_exon} "
            "does not map uniquely to an ENIGMA Table 4 deletion row."
        )
        return result

    baseline = table4_lookup_deletion(gene, table4_exon)
    baseline_strength = baseline.get("pvs1_strength")
    if not baseline.get("found") or baseline_strength not in _RNA_DOWNWEIGHT:
        result["reason"] = (
            f"PVS1 (RNA) was not applied because the {gene} {table4_exon} "
            "transcript deletion has no applicable baseline PVS1 weight in Table 4."
        )
        return result

    strength, points = _RNA_DOWNWEIGHT[baseline_strength]
    result.update({
        "applies": True,
        "strength": strength,
        "points": points,
        "appendix_branch": "patient_mrna_without_allele_specific_quantitation",
        "reason": (
            f"ENIGMA Supplementary Table 2 row {record.get('source_row')} reports "
            f"patient mRNA without allele-specific quantitation and a loss-of-function "
            f"transcript ({record.get('result')}). Table 4 assigns baseline "
            f"{baseline.get('pvs1_code')} to {gene} {table4_exon}; Appendix E Table 9 "
            f"downweights the qualitative apparent near-complete patient-mRNA branch "
            f"to PVS1 {strength} (RNA)."
        ),
        "source": APPENDIX_URL,
        "decision_path": {
            "tree_id": "figure-1b",
            "tree_version": "ENIGMA VCEP 1.2.0",
            "branch_id": "other-nucleotide-position",
            "criterion": "PVS1_RNA",
            "outcome": "applied",
            "outcome_node": "rna-other-aberrant",
            "steps": [
                {
                    "node_id": "rna-other-quality",
                    "question": "Review assay design, wild-type aberration and transcripts",
                    "result": "reviewed",
                    "observed": (
                        f"ENIGMA Supplementary Table 2 row {record.get('source_row')}; "
                        f"{record.get('splicing_assay_result_category')}"
                    ),
                },
                {
                    "node_id": "rna-other-result",
                    "question": "Observed mRNA result?",
                    "result": "aberrant",
                    "observed": str(record.get("result") or "aberrant transcript"),
                },
            ],
            "sources": [
                {
                    "source_id": "enigma-v1.2-specifications",
                    "label": (
                        f"{policy_name(gene)} v{policy_version(gene)} Specifications"
                    ),
                    "url": specification["url"],
                    "location": "Figure 1B",
                    "figure_url": "/static/enigma/figure-1b-rna.jpg",
                },
                {
                    "source_id": "enigma-v1.2-appendix",
                    "label": "ENIGMA BRCA1/2 VCEP Appendix v1.2",
                    "url": APPENDIX_URL,
                    "location": "Appendix E Table 9",
                },
            ],
        },
    })
    return result
