"""Evaluate ENIGMA PVS1 (RNA) from structured, versioned RNA evidence.

Exact VCEP assertions retain their published strength. The qualitative
patient-mRNA branch uses the official ST2 curator coding together with the
Table 4 loss-of-function context and the Appendix E RNA weighting matrix.
Complex or insufficiently coded transcript results remain manual-review only.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Optional

from backend.reference_data.ps1_splice_evidence import (
    get_st2_source_metadata,
    get_st2_splice_record,
)
from backend.reference_data.table4 import (
    TABLE4_DATA,
    parse_pvs1_code_strength,
    table4_lookup_deletion,
)
from backend.policy.gene import (
    policy_name,
    policy_version,
    reference_transcript,
    vcep_specification,
)


APPENDIX_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "9e6119dc-90b9-42b5-a3b7-1a2eb28b1b12/data"
)

_PATIENT_UNQUANTIFIED_LOF = (
    "patient not allele-specific; aberrant transcripts consistent with loss of function"
)
_DAMAGING_RNA_RESULT = "aberrant transcripts consistent with loss of function"

# Appendix E, assay results with patient mRNA without allele-specific
# quantitation, apparent (near) complete splicing column. ``Very Strong`` is
# the parsed strength of the unsuffixed Table 4 code ``PVS1``.
_UNQUANTIFIED_PATIENT_RNA_STRENGTH = {
    "Very Strong": ("Strong", 4),
    "Strong": ("Moderate", 2),
    "Moderate": ("Moderate", 2),
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


def lookup_st2_pvs1_rna_evidence(gene: str, c_notation: str) -> Dict[str, Any]:
    """Return exact ST2 facts eligible for the qualitative Appendix E branch.

    This lookup does not assign a criterion. It verifies that the official ST2
    curator coding, supporting ST3 rows and Table 4 transcript consequence are
    sufficiently explicit for the rule layer to derive the RNA strength.
    """
    provenance = get_st2_source_metadata()
    record = get_st2_splice_record(gene, c_notation)
    base: Dict[str, Any] = {
        "status": "not_in_source",
        "record": record,
        "table4_exon": None,
        "table4_baseline": None,
        "source_id": "enigma-st2-st3-v1.2",
        "source_version": str(provenance.get("version") or ""),
        "source_checksum": str(provenance.get("source_file_sha256") or ""),
        "source_url": str(provenance.get("source_url") or APPENDIX_URL),
        "reason": (
            "The exact variant has no record in ENIGMA Supplementary Table 2 v1.2."
        ),
    }
    if record is None:
        return base

    assay_category = _normalized_text(record.get("splicing_assay_result_category"))
    if _DAMAGING_RNA_RESULT not in assay_category:
        return {
            **base,
            "status": "not_applicable",
            "reason": (
                "The exact ST2 record does not report an mRNA result in the "
                "damaging loss-of-function category."
            ),
        }
    if assay_category != _PATIENT_UNQUANTIFIED_LOF:
        return {
            **base,
            "status": "review_required",
            "reason": (
                "The exact ST2 record contains damaging RNA evidence but uses "
                "an assay category that requires a different Appendix E review."
            ),
        }

    references = list(record.get("st3_references") or [])
    if not references:
        return {
            **base,
            "status": "review_required",
            "reason": "The ST2 result has no linked supporting ST3 source row.",
        }

    reported_exon = _reported_whole_exon_deletion(record.get("result"))
    if reported_exon is None:
        return {
            **base,
            "status": "review_required",
            "reason": (
                "The ST2 transcript result is complex or partial and cannot be "
                "mapped unambiguously to one Table 4 deletion row."
            ),
        }
    table4_exon = _table4_exon_for_reported_number(gene, reported_exon)
    if table4_exon is None:
        return {
            **base,
            "status": "review_required",
            "reason": (
                f"Reported exon {reported_exon} does not map uniquely to a "
                "Table 4 deletion row."
            ),
        }
    baseline = table4_lookup_deletion(gene, table4_exon)
    if (
        not baseline.get("found")
        or baseline.get("pvs1_strength") not in _UNQUANTIFIED_PATIENT_RNA_STRENGTH
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
        "assay_interpretation": (
            "enigma_st2_curated_apparent_near_complete_loss_of_function"
        ),
        "reason": (
            "The exact checksum-bound ST2 row contains ENIGMA-curated patient "
            "mRNA evidence, a loss-of-function transcript category, linked ST3 "
            "sources and an unambiguous Table 4 consequence."
        ),
    }


def _evaluate_approved_erepo_record(
    gene: str,
    c_notation: str,
    registry_result: Optional[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Convert one exact approved registry record into an automatic decision."""
    if not registry_result or registry_result.get("status") == "not_in_registry":
        return None
    specification = vcep_specification(gene)
    record = registry_result.get("record")
    problems = []
    if registry_result.get("status") != "eligible":
        problems.append(f"registry status is {registry_result.get('status') or 'unknown'}")
    if not isinstance(record, dict):
        problems.append("the exact assertion record is missing")
        record = {}
    if record.get("gene") != gene or record.get("c_notation") != c_notation:
        problems.append("the assertion does not match the assessed variant exactly")
    if record.get("reference_transcript") != reference_transcript(gene):
        problems.append("the assertion does not use the configured reference transcript")
    if record.get("guideline_id") != specification["id"]:
        problems.append("the assertion uses a different VCEP guideline")
    if record.get("guideline_version") != specification["version"]:
        problems.append("the assertion does not use the active VCEP specification version")
    if record.get("code") != "PVS1_RNA" or record.get("evidence_mechanism") != "rna_splicing":
        problems.append("the assertion does not explicitly identify PVS1 RNA evidence")
    strength, points, requires_rna = parse_pvs1_code_strength(
        str(record.get("source_code") or "")
    )
    if not requires_rna or strength is None:
        problems.append("the assertion has no explicit applicable RNA strength")
    if record.get("strength") != strength or record.get("points") != points:
        problems.append("the recorded strength and points are inconsistent")
    if not record.get("assertion_uuid") or not record.get("published_date"):
        problems.append("the assertion has no published identity and date")

    source = str(record.get("assertion_url") or registry_result.get("source_url") or "")
    source_record = {
        **record,
        "dataset": registry_result.get("registry_id"),
        "dataset_version": registry_result.get("registry_version"),
        "dataset_sha256": registry_result.get("registry_sha256"),
        "coverage_status": registry_result.get("coverage_status"),
    }
    if problems:
        return {
            "applies": False,
            "code": "PVS1_RNA",
            "strength": None,
            "points": 0,
            "reason": (
                "PVS1 (RNA) was not applied automatically because the curated "
                "ERepo candidate did not satisfy every runtime guard: "
                + "; ".join(problems)
                + "."
            ),
            "source": source or specification["url"],
            "source_record": source_record,
            "table4_exon": None,
            "appendix_branch": "curated_erepo_record_requires_review",
            "application_status": "review_required",
            "review_required": True,
            "manual_review_prefill": {
                "assay_scope": "mrna_only",
                "rna_conclusion": "damaging",
                "transcript_accession": reference_transcript(gene),
                "transcript_result_summary": str(record.get("evidence_summary") or ""),
                "source_citation": "ClinGen ERepo ENIGMA BRCA1/2 VCEP assertion",
                "source_references": [source] if source else [],
            },
        }

    return {
        "applies": True,
        "code": "PVS1_RNA",
        "strength": strength,
        "points": points,
        "reason": (
            f"ClinGen ERepo assertion {record['assertion_uuid']} reports "
            f"{record['source_code']} for this exact {reference_transcript(gene)} "
            f"variant under ENIGMA {specification['id']} v{specification['version']}. "
            f"{record['evidence_summary']}"
        ),
        "source": source,
        "source_record": source_record,
        "table4_exon": None,
        "appendix_branch": "published_vcep_pvs1_rna_assertion",
        "application_status": "applied_from_approved_registry",
        "review_required": False,
        "manual_review_prefill": {},
    }


def evaluate_pvs1_rna(
    gene: str,
    c_notation: str,
    erepo_registry_result: Optional[Dict[str, Any]] = None,
    st2_evidence_result: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Apply approved ERepo or eligible ST2 RNA evidence and prepare review."""
    curated_result = _evaluate_approved_erepo_record(
        gene,
        c_notation,
        erepo_registry_result,
    )
    if curated_result is not None:
        return curated_result

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

    st2_evidence = (
        st2_evidence_result
        if st2_evidence_result is not None
        else lookup_st2_pvs1_rna_evidence(gene, c_notation)
    )
    record = st2_evidence.get("record")
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
        "source_version": st2_evidence.get("source_version"),
        "source_checksum": st2_evidence.get("source_checksum"),
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
        expected_status = "eligible"
        expected_exon = st2_evidence.get("table4_exon")
        expected_baseline = st2_evidence.get("table4_baseline") or {}
        runtime_problems = []
        if st2_evidence.get("status") != expected_status:
            runtime_problems.append(
                f"ST2 evidence status is {st2_evidence.get('status') or 'unknown'}"
            )
        if record.get("gene") != gene or record.get("c_notation") != c_notation:
            runtime_problems.append("the ST2 evidence does not match the variant")
        if not st3_references:
            runtime_problems.append("no supporting ST3 source is linked")
        if expected_exon != table4_exon:
            runtime_problems.append("the recorded Table 4 exon does not match")
        if expected_baseline.get("pvs1_code") != baseline.get("pvs1_code"):
            runtime_problems.append("the recorded Table 4 baseline does not match")
        if baseline_strength not in _UNQUANTIFIED_PATIENT_RNA_STRENGTH:
            runtime_problems.append("the Table 4 baseline strength is unsupported")
        if runtime_problems:
            result["reason"] = (
                "PVS1 (RNA) was not applied automatically because the curated "
                "ST2 candidate failed runtime validation: "
                + "; ".join(runtime_problems)
                + ". The evidence remains available for manual review."
            )
            return result

        strength, points = _UNQUANTIFIED_PATIENT_RNA_STRENGTH[baseline_strength]
        result.update({
            "applies": True,
            "strength": strength,
            "points": points,
            "application_status": "applied_from_curated_st2",
            "review_required": False,
            "appendix_branch": (
                "patient_mrna_without_allele_specific_quantitation_"
                "apparent_near_complete"
            ),
            "reason": (
                f"ENIGMA Supplementary Table 2 row {record.get('source_row')} "
                "codes patient mRNA without allele-specific quantitation as an "
                "aberrant transcript consistent with loss of function "
                f"({record.get('result')}) and links "
                f"{len(st3_references)} supporting ST3 source(s). Table 4 assigns "
                f"baseline {baseline.get('pvs1_code')} to {gene} {table4_exon}; "
                f"the Appendix E qualitative patient-mRNA matrix yields PVS1 "
                f"{strength} (RNA)."
            ),
            "source": APPENDIX_URL,
            "manual_review_prefill": {},
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
                        "question": "Review assay design, source and transcript result",
                        "result": "curated_in_st2",
                        "observed": (
                            f"ENIGMA ST2 row {record.get('source_row')}; "
                            f"{record.get('splicing_assay_result_category')}"
                        ),
                    },
                    {
                        "node_id": "rna-other-result",
                        "question": "Observed mRNA result?",
                        "result": "aberrant_loss_of_function",
                        "observed": str(record.get("result") or ""),
                    },
                    {
                        "node_id": "rna-other-weight",
                        "question": "Appendix E qualitative RNA weight?",
                        "result": strength,
                        "observed": (
                            f"Table 4 baseline {baseline.get('pvs1_code')}"
                        ),
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
