"""Informational recommendation for review or generation of RNA evidence."""

from typing import Dict, List, Optional


ENIGMA_CSPEC_URL = "https://cspec.genome.network/cspec/ui/svi/doc/GN097"

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


def evaluate_rna_review(
    variant_type: str,
    spliceai_score: Optional[float],
    pvs1_result: Optional[Dict] = None,
    criteria: Optional[Dict] = None,
) -> Dict:
    """Identify cases in which RNA evidence may clarify a splice effect."""
    variant_type = (variant_type or "").lower()
    pvs1_result = pvs1_result or {}
    criteria = criteria or {}
    reasons: List[str] = []
    potential_branches: List[str] = []
    priority = "none"

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
        and spliceai_score >= 0.20
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
        and spliceai_score >= 0.20
    ):
        if priority == "none":
            priority = "medium"
        reasons.append(
            "A Table 9 functional result and a predicted splice effect are both "
            "present. Review whether the assay measured splicing and whether it "
            "could detect nonsense-mediated decay."
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
            "source_url": ENIGMA_CSPEC_URL,
            "is_evidence_criterion": False,
        }

    if "BP7 (RNA)" not in potential_branches:
        potential_branches.append("BP7 (RNA)")

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
        "what_to_test": [
            "Confirm the effect on the reference BRCA transcript and identify all abnormal transcript products.",
            "Quantify the proportion of normal and abnormal transcript where the assay permits.",
            "Document tissue or cell type, assay method, transcript accession, and whether nonsense-mediated decay could be detected.",
            "Determine whether the abnormal transcript is in-frame or out-of-frame and whether functional transcript remains.",
        ],
        "potential_branches": potential_branches,
        "limitations": (
            "A negative RNA result is not automatically benign. Interpretation "
            "depends on assay sensitivity, relevant tissue expression, transcript "
            "coverage, quantification, and the ability to detect transcripts "
            "subject to nonsense-mediated decay."
        ),
        "source_url": ENIGMA_CSPEC_URL,
        "is_evidence_criterion": False,
    }
