"""Pure ENIGMA BRCA v1.2 evidence-combination policy.

This module has no evidence lookup and no variant-specific branches.  Both the
native DAG and the temporary regression oracle can call it during migration.
"""

from __future__ import annotations

from typing import Mapping
from backend.gene_policy import mixed_evidence_point_thresholds


def classify_by_points(
    points: int, has_ba1: bool = False, *, gene: str | None = None
) -> tuple[int, str, str]:
    if has_ba1:
        return (1, "Benign", "BA1 stand-alone benign")
    thresholds = mixed_evidence_point_thresholds(gene)
    if points >= thresholds["pathogenic_min_inclusive"]:
        return (5, "Pathogenic", "")
    if points >= thresholds["likely_pathogenic_min_inclusive"]:
        return (4, "Likely Pathogenic", "")
    if points >= thresholds["vus_min_inclusive"]:
        return (3, "VUS", "")
    if points >= thresholds["likely_benign_min_inclusive"]:
        return (2, "Likely Benign", "")
    return (1, "Benign", "")


def criterion_strength(criterion: Mapping) -> str:
    strength = (criterion.get("strength") or "").lower().replace("-", " ")
    points = abs(criterion.get("points", 0))
    if points >= 8 or "very strong" in strength:
        return "very_strong"
    if points >= 4 or "strong" in strength:
        return "strong"
    if points >= 2 or "moderate" in strength:
        return "moderate"
    return "supporting"


def _classify_pathogenic_combination(criteria: Mapping) -> tuple[int, str, str] | None:
    counts = {"very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0}
    for criterion in criteria.values():
        if criterion.get("points", 0) > 0:
            counts[criterion_strength(criterion)] += 1
    vs, strong, moderate, supporting = (
        counts["very_strong"], counts["strong"],
        counts["moderate"], counts["supporting"],
    )
    if (
        vs >= 2
        or (vs >= 1 and (strong >= 1 or moderate >= 1 or supporting >= 2))
        or strong >= 3
        or (strong >= 2 and (moderate >= 1 or supporting >= 2))
        or (strong >= 1 and (
            moderate >= 3
            or (moderate >= 2 and supporting >= 2)
            or (moderate >= 1 and supporting >= 4)
        ))
    ):
        return (5, "Pathogenic", "")
    if (
        (vs >= 1 and (moderate >= 1 or supporting >= 1))
        or strong >= 2
        or (strong >= 1 and (moderate >= 1 or supporting >= 2))
        or moderate >= 3
        or (moderate >= 2 and supporting >= 2)
        or (moderate >= 1 and supporting >= 4)
    ):
        return (4, "Likely Pathogenic", "")
    return None


def _classify_benign_combination(criteria: Mapping) -> tuple[int, str, str] | None:
    counts = {"very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0}
    has_eligible_single_strong = False
    for criterion in criteria.values():
        if criterion.get("points", 0) >= 0:
            continue
        strength = criterion_strength(criterion)
        counts[strength] += 1
        if (
            strength == "strong"
            and criterion.get("single_strong_likely_benign_eligible") is True
        ):
            has_eligible_single_strong = True
    very_strong, strong, moderate, supporting = (
        counts["very_strong"], counts["strong"],
        counts["moderate"], counts["supporting"],
    )
    if (
        very_strong >= 2
        or (very_strong >= 1 and (strong >= 1 or moderate >= 1 or supporting >= 1))
        or strong >= 2
        or (strong >= 1 and (
            moderate >= 2
            or (moderate >= 1 and supporting >= 1)
            or supporting >= 3
        ))
    ):
        return (1, "Benign", "")
    if (
        very_strong >= 1
        or has_eligible_single_strong
        or ((very_strong + strong) >= 1 and (moderate >= 1 or supporting >= 1))
        or moderate >= 2
        or (moderate >= 1 and supporting >= 1)
        or supporting >= 2
    ):
        return (2, "Likely Benign", "")
    return None


def classify_by_enigma_combination(
    criteria: Mapping, points: int, *, gene: str | None = None
) -> tuple[int, str, str]:
    ba1 = criteria.get("BA1")
    if ba1 and ba1.get("applies", True):
        return (
            1,
            "Benign",
            "BA1 is stand-alone benign evidence. Other criteria remain visible for audit "
            "but do not replace the BA1 classification.",
        )

    has_pathogenic = any(item.get("points", 0) > 0 for item in criteria.values())
    has_benign = any(item.get("points", 0) < 0 for item in criteria.values())
    if has_pathogenic and has_benign:
        cls, label, _ = classify_by_points(points, gene=gene)
        return (
            cls,
            label,
            "Mixed pathogenic and benign evidence. ENIGMA VCEP v1.2 second "
            "classification approach applies the Tavtigian 2020 point system; "
            "expert review is required.",
        )
    if has_pathogenic:
        result = _classify_pathogenic_combination(criteria)
        if result:
            return result
        return (
            3,
            "VUS",
            "Pathogenic evidence does not yet meet an ENIGMA VCEP v1.2 Table 3 "
            "combination for Likely Pathogenic. For example, PVS1 Very Strong "
            "requires at least one additional Supporting criterion.",
        )
    if has_benign:
        result = _classify_benign_combination(criteria)
        if result:
            return result
        strong = [
            criterion
            for criterion in criteria.values()
            if criterion.get("points", 0) < 0
            and criterion_strength(criterion) == "strong"
        ]
        if len(strong) == 1 and len(criteria) == 1:
            return (
                3,
                "VUS",
                "A single Strong benign criterion reaches Likely Benign under "
                "ENIGMA VCEP v1.2 Table 3 only when multiple evidence "
                "contributions are documented for that criterion. This record "
                "does not meet that combination requirement.",
            )
        return (
            3,
            "VUS",
            "Benign evidence does not yet meet an ENIGMA VCEP v1.2 Table 3 "
            "combination for Likely Benign.",
        )
    return (3, "VUS", "")


def verify_acmg_combination(criteria: Mapping, points: int, predicted_class: int) -> str | None:
    path_vs = path_s = path_m = path_p = ben_s = ben_p = 0
    for criterion in criteria.values():
        criterion_points = criterion.get("points", 0)
        strength = (criterion.get("strength") or "").lower()
        if criterion_points > 0:
            if criterion_points >= 8 or "very strong" in strength:
                path_vs += 1
            elif criterion_points >= 4 or "strong" in strength:
                path_s += 1
            elif criterion_points >= 2 or "moderate" in strength:
                path_m += 1
            else:
                path_p += 1
        elif criterion_points < 0:
            if criterion_points <= -4 or "strong" in strength:
                ben_s += 1
            else:
                ben_p += 1
    if predicted_class >= 4 and ben_s > 0 and path_s == 0 and path_vs == 0:
        return (
            "Classification note: pathogenic criteria are all moderate/supporting "
            "but benign strong evidence also present - manual review recommended."
        )
    return None
