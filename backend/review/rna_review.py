"""Informational recommendation for review or generation of RNA evidence."""

import re
from collections.abc import Iterable
from typing import Any, Dict, List, Optional
from backend.policy.gene import (
    reference_transcript,
    spliceai_thresholds,
    vcep_specification,
)
from backend.reference_data.enigma_rules import lookup_st16_table4_rna_references
from backend.review.recommendation import empty_review_recommendation

SPLICE_RELEVANT_TYPES = {
    "splice_site",
    "intronic",
    "synonymous",
    "silent",
    "missense",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_delins",
    "delins",
}


_PMID_CONTEXT = re.compile(r"PMIDs?|pubmed\.ncbi\.nlm\.nih\.gov/", re.I)
_PMID_NUMBER = re.compile(r"\b\d{6,9}\b")


def _pmids(value: Any) -> set[str]:
    """Extract explicitly labelled PubMed identifiers from nested audit values."""
    values: Iterable[Any]
    if isinstance(value, dict):
        values = value.values()
    elif isinstance(value, (list, tuple, set)):
        values = value
    else:
        text = str(value or "")
        return set(_PMID_NUMBER.findall(text)) if _PMID_CONTEXT.search(text) else set()
    return set().union(*(_pmids(item) for item in values)) if values else set()


def _table9_functional_splice_overlap(criteria: Dict) -> list[str]:
    """Return PMIDs used in both Table 9 assay and splice-context fields.

    This is a provenance-overlap signal, not proof that the same experimental
    observation was counted twice. Table 9 may describe a combined mRNA and
    protein assay or a paper containing more than one experiment.
    """
    overlaps: set[str] = set()
    for code in ("PS3", "BS3"):
        criterion = criteria.get(code) or {}
        audit = criterion.get("table9_audit") or {}
        functional_pmids = _pmids(audit.get("assay_results") or [])
        splice_pmids = _pmids(audit.get("splice_result_published") or "")
        overlaps.update(functional_pmids & splice_pmids)
    return sorted(overlaps, key=int)


def _table9_published_splice_prefill(gene: str, criteria: Dict) -> Dict:
    """Prepare independent published splice context for manual RNA review.

    Table 9 is an automatic PS3/BS3 source, not an automatic PVS1 RNA source.
    Its published splice-result field can nevertheless save a curator from
    re-entering an already structured transcript result and its references.
    PMIDs also used for the functional assay are deliberately excluded so that
    the prefill does not present overlapping evidence as an independent source.
    """
    summaries: list[str] = []
    independent_pmids: set[str] = set()

    for code in ("PS3", "BS3"):
        criterion = criteria.get(code) or {}
        if not criterion.get("applies"):
            continue
        audit = criterion.get("table9_audit") or {}
        splice_flag = str(
            audit.get("predicted_or_observed_splicing") or ""
        ).strip().upper()
        summary = str(audit.get("splice_result_published") or "").strip()
        if not splice_flag.startswith("Y") or not summary:
            continue

        splice_pmids = _pmids(summary)
        functional_pmids = _pmids(audit.get("assay_results") or [])
        independent = splice_pmids - functional_pmids
        if not independent:
            continue

        if summary not in summaries:
            summaries.append(summary)
        independent_pmids.update(independent)

    if not summaries:
        return {}

    return {
        "transcript_accession": reference_transcript(gene),
        "transcript_result_summary": (
            "ENIGMA Specifications Table 9 published splice context: "
            + " | ".join(summaries)
        ),
        "source_citation": (
            "ENIGMA Specifications Table 9 v1.2 (published splice context)"
        ),
        "source_references": [
            f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            for pmid in sorted(independent_pmids, key=int)
        ],
    }


def _merge_review_prefill(existing: Dict, additional: Dict) -> Dict:
    """Merge non-scoring review context without replacing curated ST2 data."""
    if not existing:
        return dict(additional)
    merged = dict(existing)

    additional_summary = str(additional.get("transcript_result_summary") or "")
    existing_summary = str(merged.get("transcript_result_summary") or "")
    if additional_summary and additional_summary not in existing_summary:
        merged["transcript_result_summary"] = (
            f"{existing_summary}\n\n{additional_summary}"
            if existing_summary
            else additional_summary
        )

    if not merged.get("transcript_accession"):
        merged["transcript_accession"] = additional.get("transcript_accession")
    if not merged.get("source_citation"):
        merged["source_citation"] = additional.get("source_citation")

    references = list(merged.get("source_references") or [])
    for reference in additional.get("source_references") or []:
        if reference not in references:
            references.append(reference)
    merged["source_references"] = references
    return merged


def _st16_table4_rna_prefill(
    gene: str,
    c_notation: str,
    pvs1_result: Dict,
) -> Dict:
    """Prepare exact Table 4 RNA publication context without assigning points."""
    if not pvs1_result.get("requires_rna") or not c_notation:
        return {}
    record = lookup_st16_table4_rna_references(gene, c_notation)
    if record is None:
        return {}
    pmids = list(record["pmids"])
    return {
        "transcript_accession": reference_transcript(gene),
        "transcript_result_summary": (
            "ENIGMA Table 4 identifies this exact allele as PVS1 (RNA). "
            "Supplementary Table 16 links the publication"
            + ("s" if len(pmids) != 1 else "")
            + " below. Review the reported transcript result under Appendix E; "
            "ARIANE has not inferred an RNA strength."
        ),
        "source_citation": (
            f"ENIGMA BRCA1/2 VCEP v1.2 Supplementary Table 16, "
            f"row {record['source_row']}"
        ),
        "source_references": [
            f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" for pmid in pmids
        ],
    }


def evaluate_rna_review(
    gene: str,
    variant_type: str,
    spliceai_score: Optional[float],
    c_notation: str = "",
    pvs1_result: Optional[Dict] = None,
    pvs1_rna_result: Optional[Dict] = None,
    criteria: Optional[Dict] = None,
) -> Dict:
    """Identify cases in which RNA evidence may clarify a splice effect."""
    variant_type = (variant_type or "").lower()
    pvs1_result = pvs1_result or {}
    pvs1_rna_result = pvs1_rna_result or {}
    criteria = criteria or {}
    reasons: List[str] = []
    potential_branches: List[str] = []
    priority = "none"
    manual_review_prefill: Dict = {}
    splice_high = spliceai_thresholds(gene)["pp3"]
    source_url = vcep_specification(gene)["url"]

    if criteria.get("PVS1_RNA", {}).get("applies"):
        return empty_review_recommendation(gene)

    if pvs1_rna_result.get("review_required"):
        priority = "high"
        reasons.append(pvs1_rna_result.get("reason") or (
            "An exact official RNA-evidence record requires curator review."
        ))
        potential_branches.append("PVS1 (RNA)")
        manual_review_prefill = dict(
            pvs1_rna_result.get("manual_review_prefill") or {}
        )

    if pvs1_result.get("requires_rna"):
        priority = "high"
        reasons.append(
            "ENIGMA Table 4 routes this variant to an RNA-dependent PVS1 assessment."
        )
        potential_branches.append("PVS1 (RNA)")

    st16_prefill = _st16_table4_rna_prefill(gene, c_notation, pvs1_result)
    if st16_prefill:
        manual_review_prefill = _merge_review_prefill(
            manual_review_prefill,
            st16_prefill,
        )
        reasons.append(
            "ENIGMA Supplementary Table 16 links publication references to this "
            "exact Table 4 PVS1 (RNA) allele. ARIANE has copied the references "
            "into the manual RNA form without inferring a strength or adding points."
        )

    pvs1_reason = (pvs1_result.get("reason") or "").lower()
    if variant_type == "splice_site" and not pvs1_result.get("applies"):
        priority = "high"
        if "not found in table 4" in pvs1_reason:
            reasons.append(
                "The canonical splice-site variant was not resolved by the local "
                "Table 4 lookup, so PVS1 must not be inferred from position alone."
            )
        elif not pvs1_result.get("requires_rna"):
            reasons.append(
                "The splice-site variant has no automatically applicable PVS1 "
                "result and needs review of its actual transcript effect."
            )
        if "PVS1 (RNA)" not in potential_branches:
            potential_branches.append("PVS1 (RNA)")
        potential_branches.append("BP7 (RNA)")

    if (
        variant_type in SPLICE_RELEVANT_TYPES
        and spliceai_score is not None
        and spliceai_score >= splice_high
    ):
        if priority == "none":
            priority = "medium"
        reasons.append(
            f"Reference-transcript SpliceAI is {spliceai_score:.3f}, indicating a "
            "predicted splice effect that is not direct experimental evidence."
        )
        if "PVS1 (RNA)" not in potential_branches:
            potential_branches.append("PVS1 (RNA)")

    has_functional_evidence = any(
        criteria.get(code, {}).get("applies") for code in ("PS3", "BS3")
    )
    if (
        has_functional_evidence
        and spliceai_score is not None
        and spliceai_score >= splice_high
    ):
        if priority == "none":
            priority = "medium"
        reasons.append(
            "A Table 9 functional result and a predicted splice effect are both "
            "present. Review whether the assay measured splicing and whether it "
            "could detect nonsense-mediated decay."
        )

    table9_prefill = _table9_published_splice_prefill(gene, criteria)
    if table9_prefill and reasons:
        manual_review_prefill = _merge_review_prefill(
            manual_review_prefill,
            table9_prefill,
        )
        reasons.append(
            "ENIGMA Specifications Table 9 contains a published splice result "
            "with references that are distinct from its PS3/BS3 functional-assay "
            "reference. ARIANE has copied this context into the manual RNA form "
            "without inferring a PVS1 RNA strength or adding points."
        )

    publication_overlap = _table9_functional_splice_overlap(criteria)
    if publication_overlap and (
        pvs1_rna_result.get("review_required") or pvs1_result.get("requires_rna")
    ):
        reasons.append(
            "Table 9 cites PMID "
            + ", ".join(publication_overlap)
            + " both as functional-assay evidence and in its published splice "
            "context. ARIANE counts the Table 9 functional result only once. "
            "Before adding PVS1 (RNA), use independent RNA evidence or document "
            "that the publication contains distinct non-overlapping assay data."
        )

    if not reasons:
        return empty_review_recommendation(gene)

    if "BP7 (RNA)" not in potential_branches:
        potential_branches.append("BP7 (RNA)")

    what_to_test = [
        "Confirm the effect on the configured reference transcript and identify all abnormal transcript products.",
        "Quantify the proportion of normal and abnormal transcript where the assay permits.",
        "Document tissue or cell type, assay method, transcript accession, and whether nonsense-mediated decay could be detected.",
        "Determine whether the abnormal transcript is in-frame or out-of-frame and whether functional transcript remains.",
    ]
    if publication_overlap:
        what_to_test.append(
            "Compare the publications supporting PS3/BS3 and PVS1 (RNA). Record "
            "which assay and observations support each code and exclude reused evidence."
        )

    return {
        "recommended": True,
        "priority": priority,
        "title": "RNA evidence review recommended",
        "summary": (
            "RNA evidence may help determine the actual transcript consequence. "
            "This recommendation is an ARIANE review aid, not an ACMG/AMP or "
            "ENIGMA evidence criterion, and it adds no points."
        ),
        "reasons": reasons,
        "what_to_test": what_to_test,
        "potential_branches": potential_branches,
        "limitations": (
            "A negative RNA result is not automatically benign. Interpretation "
            "depends on assay sensitivity, relevant tissue expression, transcript "
            "coverage, quantification, and the ability to detect transcripts "
            "subject to nonsense-mediated decay."
        ),
        "source_url": source_url,
        "is_evidence_criterion": False,
        "manual_review_prefill": manual_review_prefill,
    }
