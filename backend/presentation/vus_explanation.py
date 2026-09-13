"""VUS explanation layer.

This module describes why a class 3 result remains VUS and what kind of
evidence would be useful next. It is informational only and never changes
criteria, points, or the final classification.
"""

from __future__ import annotations

from typing import Dict, Optional


PATH_CODES = {"PVS1", "PM5_PTC", "PM2_Supporting", "PP3", "PP4", "PS1", "PS3"}
BENIGN_CODES = {"BA1", "BS1_Strong", "BS1_Supporting", "BS3", "BP1", "BP4", "BP5", "BP7"}


def _code_set(criteria: Dict) -> set[str]:
    return {
        name
        for name, criterion in criteria.items()
        if criterion.get("applies", True) and criterion.get("points", 0) != 0
    }


def _has_pathogenic(criteria: Dict) -> bool:
    return any(criterion.get("points", 0) > 0 for criterion in criteria.values())


def _has_benign(criteria: Dict) -> bool:
    return any(criterion.get("points", 0) < 0 for criterion in criteria.values())


def _base(
    category: str,
    tier: str,
    title: str,
    summary: str,
    what_to_check: str,
    review_priority: str,
) -> dict:
    return {
        "category": category,
        "tier": tier,
        "title": title,
        "summary": summary,
        "what_to_check": what_to_check,
        "review_priority": review_priority,
    }


def explain_vus(result: Dict) -> Optional[dict]:
    """Return a VUS explanation dict for class 3 results."""
    if result.get("predicted_class") != 3:
        return None

    criteria = result.get("criteria", {})
    codes = _code_set(criteria)
    has_pathogenic = _has_pathogenic(criteria)
    has_benign = _has_benign(criteria)

    if "PVS1" in codes and "BS3" in codes:
        return _base(
            "conflicting_pvs1_bs3",
            "C",
            "VUS with conflicting PVS1 and BS3 evidence",
            "PVS1 provides very strong pathogenic evidence, but calibrated functional evidence supports a benign effect.",
            "Review the PVS1 context and the BS3 functional evidence before changing the classification.",
            "high",
        )

    if "PP3" in codes and "BS3" in codes:
        return _base(
            "conflicting_pp3_bs3",
            "C",
            "VUS with conflicting computational and functional evidence",
            "A computational pathogenic signal is present, but functional evidence supports a benign effect.",
            "Review functional assay details and consider RNA or independent clinical evidence.",
            "high",
        )

    if "PP3" in codes and "PS3" in codes:
        return _base(
            "ps3_pp3_one_step_short",
            "A",
            "VUS: functional and computational evidence point pathogenic",
            "PS3 and PP3 both point toward pathogenicity, but the ENIGMA combination is still incomplete.",
            "Look for RNA confirmation, PP1, PM3, PS4, PP4, curated PS1, or other independent pathogenic evidence.",
            "high",
        )

    if "PVS1" in codes:
        return _base(
            "pvs1_needs_support",
            "B",
            "VUS with PVS1: why this notice is shown",
            "This notice does not mean that the classification or the PVS1 assignment is incorrect. It is shown because PVS1 is strong pathogenic evidence, but PVS1 alone does not satisfy an ENIGMA combination for Likely Pathogenic. With the currently accepted evidence, VUS is therefore the expected result. High review priority means that the result may be close to a classification threshold, not that an error was detected.",
            "Verify the assigned PVS1 strength and review whether any independent, already-supported pathogenic or benign evidence is available. Do not add evidence solely to move the variant out of VUS.",
            "high",
        )

    if "PS3" in codes:
        return _base(
            "strong_pathogenic_one_step_short",
            "B",
            "VUS: strong pathogenic evidence is one step short",
            "Strong pathogenic evidence is present, but another independent evidence item is needed for likely pathogenic.",
            "Look for PP1, PM3, PS4, PP4, curated PS1, or another accepted supporting pathogenic criterion.",
            "high",
        )

    if {"BP4", "BP7", "PM2_Supporting"}.issubset(codes):
        return _base(
            "bp4_bp7_pm2_benign_leaning",
            "D",
            "Benign-leaning VUS: BP4/BP7 plus PM2",
            "Benign splice prediction is present, but PM2 background keeps the result in VUS.",
            "RNA evidence, functional benign evidence, BS2, BS4, or other benign support would be useful.",
            "medium",
        )

    if has_pathogenic and has_benign:
        return _base(
            "mixed_benign_pathogenic_evidence",
            "C",
            "VUS with mixed benign and pathogenic evidence",
            "The automated evidence contains both benign and pathogenic directions, so expert review is required.",
            "Adjudicate the conflicting criteria before any upgrade or downgrade.",
            "high",
        )

    if codes == {"PM2_Supporting"}:
        return _base(
            "pm2_only",
            "E",
            "VUS: PM2 is the only evidence",
            "Absence or rarity in population data is not enough to classify the variant.",
            "Do not prioritize unless another independent evidence source appears.",
            "low",
        )

    if "PP3" in codes:
        return _base(
            "computational_evidence_not_enough",
            "B",
            "VUS: computational evidence is not enough",
            "A computational pathogenic signal is present, but this does not resolve the variant by itself.",
            "Look for functional, RNA, segregation, case-control, trans, or curated clinical evidence.",
            "medium",
        )

    if not codes:
        return _base(
            "no_automated_evidence",
            "E",
            "VUS: no automated Module 1 criteria applied",
            "No automated Module 1 evidence was strong enough to move the classification.",
            "External clinical, functional, segregation, case-control, or curated evidence would be needed.",
            "low",
        )

    return _base(
        "vus_unresolved",
        "E",
        "VUS: available automated evidence is insufficient",
        "The applied automated criteria do not meet an ENIGMA benign or pathogenic combination.",
        "Manual review and non-automated evidence may be needed.",
        "medium",
    )
