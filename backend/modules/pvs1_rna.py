"""Prepare ENIGMA PVS1 (RNA) review from structured official RNA evidence.

Supplementary Table 2 can identify an exact RNA-evidence candidate and Table 4
can provide its loss-of-function context. An unquantified ST2 result does not,
however, identify the Appendix E Table 9 branch or criterion strength. Those
records therefore prefill expert review and never create automatic points.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Optional

from backend.modules.ps1_splice_evidence import get_st2_splice_record
from backend.modules.table4 import TABLE4_DATA, table4_lookup_deletion
from backend.gene_policy import reference_transcript, vcep_specification


APPENDIX_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "9e6119dc-90b9-42b5-a3b7-1a2eb28b1b12/data"
)

_PATIENT_UNQUANTIFIED_LOF = (
    "patient not allele-specific; aberrant transcripts consistent with loss of function"
)
_DAMAGING_RNA_RESULT = "aberrant transcripts consistent with loss of function"


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


def _st3_citation(reference: Dict[str, Any]) -> str:
    author = str(reference.get("author") or "Unknown author")
    year = str(reference.get("year") or "unknown year")
    pmid = str(reference.get("pmid") or "").strip()
    citation = f"{author} {year}"
    if pmid:
        citation += f", PMID {pmid}"
    return citation


def _st3_reference_url(reference: Dict[str, Any]) -> Optional[str]:
    pmid = str(reference.get("pmid") or "").strip()
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None


def evaluate_pvs1_rna(gene: str, c_notation: str) -> Dict[str, Any]:
    """Return a non-scoring PVS1 (RNA) review candidate from exact ST2 data."""
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
        "application_status": "not_applied",
        "review_required": False,
        "manual_review_prefill": {},
    }

    record = get_st2_splice_record(gene, c_notation)
    if record is None:
        result["reason"] = (
            "PVS1 (RNA) was not applied because the exact variant has no record "
            "in ENIGMA Supplementary Table 2 v1.2."
        )
        return result

    st3_references = list(record.get("st3_references") or [])
    result["source_record"] = {
        "dataset": "ENIGMA Supplementary Table 2 v1.2",
        "source_row": record.get("source_row"),
        "gene": record.get("gene"),
        "c_notation": record.get("c_notation"),
        "p_notation": record.get("p_notation"),
        "position_category": record.get("included_in_analysis"),
        "assay_category": record.get("splicing_assay_result_category"),
        "assay_summary": record.get("variant_assay_summary"),
        "transcript_result": record.get("result"),
        "evidence_mechanism": "rna_splicing",
        "supporting_st3_references": st3_references,
    }

    assay_category = _normalized_text(record.get("splicing_assay_result_category"))
    if _DAMAGING_RNA_RESULT not in assay_category:
        result["reason"] = (
            "PVS1 (RNA) was not applied automatically because the exact ENIGMA "
            "ST2 record does not report an mRNA result in the damaging "
            "loss-of-function category. The record must be considered under the "
            "appropriate RNA branch and cannot create PVS1 (RNA)."
        )
        return result

    is_unquantified = assay_category == _PATIENT_UNQUANTIFIED_LOF
    review_summary = (
        f"ENIGMA Supplementary Table 2 row {record.get('source_row')}: "
        f"{record.get('result')}. Assay category: "
        f"{record.get('splicing_assay_result_category')}. "
        f"Supporting ST3 sources: {', '.join(_st3_citation(ref) for ref in st3_references)}. "
        "The curator must determine the applicable Appendix E RNA weighting "
        "branch and document the supported strength."
    )
    result.update({
        "application_status": "review_required",
        "review_required": True,
        "appendix_branch": (
            "unquantified_patient_mrna_requires_consensus_review"
            if is_unquantified
            else "structured_rna_result_requires_curator_weighting"
        ),
        "source": APPENDIX_URL,
        "manual_review_prefill": {
            "assay_scope": "mrna_only",
            "rna_conclusion": "damaging",
            "transcript_accession": reference_transcript(gene),
            "transcript_result_summary": review_summary,
            "source_citation": "ENIGMA Supplementary Tables 2 and 3 v1.2",
            "source_references": [
                url for url in (_st3_reference_url(ref) for ref in st3_references)
                if url
            ],
        },
    })

    reported_exon = _reported_whole_exon_deletion(record.get("result"))
    if reported_exon is None:
        result["reason"] = (
            "PVS1 (RNA) was not applied automatically because ENIGMA ST2 reports "
            "a complex or partial transcript consequence that cannot be mapped "
            "unambiguously to one Table 4 deletion row. The exact RNA record has "
            "been prepared for expert review without a code weight."
        )
        return result

    table4_exon = _table4_exon_for_reported_number(gene, reported_exon)
    result["table4_exon"] = table4_exon
    if table4_exon is None:
        result["reason"] = (
            f"PVS1 (RNA) was not applied because reported exon {reported_exon} "
            "does not map uniquely to an ENIGMA Table 4 deletion row. The exact "
            "RNA record has been prepared for expert review without a code weight."
        )
        return result

    baseline = table4_lookup_deletion(gene, table4_exon)
    baseline_strength = baseline.get("pvs1_strength")
    if not baseline.get("found") or not baseline_strength:
        result["reason"] = (
            f"PVS1 (RNA) was not applied because the {gene} {table4_exon} "
            "transcript deletion has no applicable baseline PVS1 weight in Table 4. "
            "The exact RNA record has been prepared for expert review without a "
            "code weight."
        )
        return result

    result["manual_review_prefill"]["table4_context"] = (
        f"{gene} {table4_exon}: {baseline.get('pvs1_code')}"
    )
    if is_unquantified:
        result["reason"] = (
            f"ENIGMA Supplementary Table 2 row {record.get('source_row')} reports "
            f"patient mRNA without allele-specific quantitation and an aberrant "
            f"transcript consistent with loss of function ({record.get('result')}). "
            f"Table 4 supplies the {baseline.get('pvs1_code')} loss-of-function "
            f"context for {gene} {table4_exon}, but ST2 does not establish whether "
            "the unquantified result is apparent near-complete or incomplete. "
            "Appendix E Table 9 requires consensus curator judgement, so PVS1 "
            "(RNA) was not applied automatically."
        )
    else:
        result["reason"] = (
            f"ENIGMA Supplementary Table 2 row {record.get('source_row')} reports "
            f"damaging RNA evidence ({record.get('result')}). Table 4 supplies "
            f"the {baseline.get('pvs1_code')} loss-of-function context for {gene} "
            f"{table4_exon}. The ST2 record does not itself provide a final PVS1 "
            "(RNA) weight, so the result has been prepared for expert review and "
            "no points were applied automatically."
        )
    return result
