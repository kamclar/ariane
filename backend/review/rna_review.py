"""Informational recommendation for review or generation of RNA evidence."""

import re
from typing import Dict, List, Optional
from backend.policy.gene import spliceai_thresholds, vcep_specification

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


def _pmids(value) -> set[str]:
    """Extract explicitly labelled PubMed identifiers from nested audit values."""
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


def evaluate_rna_review(
    gene: str,
    variant_type: str,
    spliceai_score: Optional[float],
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
        return {
            "recommended": False,
            "priority": "none",
            "title": "",
            "summary": "",
            "reasons": [],
            "what_to_test": [],
            "potential_branches": [],
            "limitations": "",
            "source_url": source_url,
            "is_evidence_criterion": False,
            "manual_review_prefill": {},
        }

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
        return {
            "recommended": False,
            "priority": "none",
            "title": "",
            "summary": "",
            "reasons": [],
            "what_to_test": [],
            "potential_branches": [],
            "limitations": "",
            "source_url": source_url,
            "is_evidence_criterion": False,
            "manual_review_prefill": {},
        }

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
