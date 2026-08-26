"""Informational recommendation for splice PS1 review."""

from typing import Dict, List, Optional
from backend.gene_policy import spliceai_thresholds, vcep_specification
SPLICE_PS1_RULE_SOURCE = (
    "Rule and weighting source: ENIGMA BRCA1/2 VCEP v1.2 PS1 and "
    "Appendix J Table 17. ARIANE has no active curated splice-PS1 "
    "reference registry."
)

SPLICE_PS1_RELEVANT_TYPES = {
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


def evaluate_splice_ps1_review(
    gene: str,
    variant_type: str,
    spliceai_score: Optional[float],
    ps1_result: Optional[Dict] = None,
) -> Dict:
    """Identify variants that may need non-automated splice PS1 review."""
    variant_type = (variant_type or "").lower()
    ps1_result = ps1_result or {}
    reasons: List[str] = []
    splice_high = spliceai_thresholds(gene)["pp3"]

    if ps1_result.get("applies"):
        return _empty(gene)

    if variant_type == "splice_site":
        reasons.append(
            "This is a splice-site variant. PS1 may require review against known "
            "P/LP variants with the same experimentally or confidently predicted "
            "splicing consequence."
        )

    if (
        variant_type in SPLICE_PS1_RELEVANT_TYPES
        and spliceai_score is not None
        and spliceai_score >= splice_high
    ):
        reasons.append(
            f"Reference-transcript SpliceAI is {spliceai_score:.3f}. Current "
            "automatic PS1 covers only protein-level same-amino-acid change and "
            "does not evaluate same-splice-effect PS1."
        )

    if not reasons:
        return _empty(gene)

    return {
        "recommended": True,
        "priority": "medium",
        "title": "Splice PS1 review candidate",
        "summary": (
            "This variant may need review for a same-splicing-impact PS1 branch. "
            "ARIANE does not currently score splice PS1 automatically and this "
            "notice adds no points."
        ),
        "reasons": reasons,
        "what_to_test": [
            "Find any known Pathogenic/Likely Pathogenic BRCA1/2 variant at the same splice event or splice site.",
            "Compare whether the assessed variant and reference variant cause the same exon skipping, cryptic splice-site use, or other transcript consequence.",
            "Check whether RNA evidence or ENIGMA/ClinGen expert-panel evidence documents that same splicing consequence.",
            "Confirm the reference variant source, classification strength, transcript accession, and whether the assessed variant prediction is at least comparable.",
        ],
        "potential_branches": ["PS1 (splice)"],
        "limitations": (
            "SpliceAI alone is not enough for PS1. The review needs a specific "
            "known P/LP reference variant and evidence that the splice consequence "
            "is the same."
        ),
        "reference_source": SPLICE_PS1_RULE_SOURCE,
        "source_url": vcep_specification(gene)["url"],
        "is_evidence_criterion": False,
    }


def _empty(gene: str) -> Dict:
    return {
        "recommended": False,
        "priority": "none",
        "title": "",
        "summary": "",
        "reasons": [],
        "what_to_test": [],
        "potential_branches": [],
        "limitations": "",
        "reference_source": "",
        "source_url": vcep_specification(gene)["url"],
        "is_evidence_criterion": False,
    }
