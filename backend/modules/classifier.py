# ============================================================
# ARIANE classifier - main evaluation orchestrator
#
# Evidence hierarchy (ENIGMA VCEP v1.2):
#   1. BA1 check - stand-alone benign, if met → class 1, stop
#   2. Table 9 - PS3/BS3 functional evidence; flags has_functional_evidence
#   3. Table 4 - PVS1/PM5 structural rules
#   4. gnomAD - BS1, PM2
#   5. ST7 - PP4/BP5 multifactorial likelihood
#   6. ST7 - PS1 same amino acid change as known P/LP
#   7. SpliceAI/BayesDel - PP3 always; BP4/BP7 only if no functional evidence
#   8. BP1 - domain check
#   9. Classification from adapted ACMG/AMP combinations
#  10. Tavtigian 2020 points only for contradictory evidence
# ============================================================
from typing import Optional, Dict, List, Tuple


def classify_by_points(points: int, has_ba1: bool = False) -> tuple:
    if has_ba1:
        return (1, "Benign", "BA1 stand-alone benign")
    if points >= 10:
        return (5, "Pathogenic", "")
    elif points >= 6:
        return (4, "Likely Pathogenic", "")
    elif points >= -1:
        return (3, "VUS", "")
    elif points >= -6:
        return (2, "Likely Benign", "")
    else:
        return (1, "Benign", "")


def _criterion_strength(criterion: Dict) -> str:
    strength = (criterion.get("strength") or "").lower().replace("-", " ")
    points = abs(criterion.get("points", 0))
    if points >= 8 or "very strong" in strength:
        return "very_strong"
    if points >= 4 or "strong" in strength:
        return "strong"
    if points >= 2 or "moderate" in strength:
        return "moderate"
    return "supporting"


def _classify_pathogenic_combination(criteria: Dict) -> Optional[tuple]:
    counts = {"very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0}
    for criterion in criteria.values():
        if criterion.get("points", 0) > 0:
            counts[_criterion_strength(criterion)] += 1

    vs = counts["very_strong"]
    strong = counts["strong"]
    moderate = counts["moderate"]
    supporting = counts["supporting"]

    if (
        (vs >= 1 and (strong >= 1 or moderate >= 1 or supporting >= 2))
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
        (vs >= 1 and supporting >= 1)
        or (strong >= 1 and (moderate >= 1 or supporting >= 2))
        or moderate >= 3
        or (moderate >= 2 and supporting >= 2)
        or (moderate >= 1 and supporting >= 4)
    ):
        return (4, "Likely Pathogenic", "")

    return None


def _classify_benign_combination(criteria: Dict) -> Optional[tuple]:
    counts = {"very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0}
    composite_strong_codes = {"BP1", "BS4", "BP5"}
    has_composite_strong = False

    for name, criterion in criteria.items():
        if criterion.get("points", 0) >= 0:
            continue
        strength = _criterion_strength(criterion)
        counts[strength] += 1
        if strength == "strong" and name in composite_strong_codes:
            has_composite_strong = True

    strong = counts["very_strong"] + counts["strong"]
    moderate = counts["moderate"]
    supporting = counts["supporting"]

    if strong >= 2:
        return (1, "Benign", "")

    if (
        has_composite_strong
        or (strong >= 1 and (moderate >= 1 or supporting >= 1))
        or (moderate >= 1 and supporting >= 1)
        or supporting >= 2
    ):
        return (2, "Likely Benign", "")

    return None


def classify_by_enigma_combination(criteria: Dict, points: int) -> tuple:
    """
    Apply ENIGMA VCEP v1.2 Table 3 combinations by default.
    Tavtigian points are reserved for contradictory benign/pathogenic evidence.
    """
    has_pathogenic = any(criterion.get("points", 0) > 0 for criterion in criteria.values())
    has_benign = any(criterion.get("points", 0) < 0 for criterion in criteria.values())

    if has_pathogenic and has_benign:
        cls, label, _ = classify_by_points(points)
        return (
            cls,
            label,
            "Contradictory benign and pathogenic evidence - classification uses the "
            "Tavtigian 2020 point system and requires expert review.",
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

    return (3, "VUS", "")


def verify_acmg_combination(criteria: Dict, points: int, predicted_class: int) -> Optional[str]:
    """
    Verify point-based classification against ACMG combination rules.
    Returns warning if inconsistent.
    """
    path_vs, path_s, path_m, path_p = 0, 0, 0, 0
    ben_s, ben_p = 0, 0

    for name, crit in criteria.items():
        pts = crit.get("points", 0)
        strength = (crit.get("strength") or "").lower()
        if pts > 0:
            if pts >= 8 or "very strong" in strength:
                path_vs += 1
            elif pts >= 4 or "strong" in strength:
                path_s += 1
            elif pts >= 2 or "moderate" in strength:
                path_m += 1
            else:
                path_p += 1
        elif pts < 0:
            if pts <= -4 or "strong" in strength:
                ben_s += 1
            else:
                ben_p += 1

    if predicted_class >= 4 and ben_s > 0 and path_s == 0 and path_vs == 0:
        return (
            "Classification note: pathogenic criteria are all moderate/supporting "
            "but benign strong evidence also present - manual review recommended."
        )
    return None


def evaluate_variant(
    gene: str,
    variant_type: str,
    p_notation: str,
    c_notation: str,
    spliceai_score: Optional[float] = None,
    bayesdel_score: Optional[float] = None,
    gnomad_data: Optional[Dict] = None,
    table9_result: Optional[Dict] = None,
    pp4_bp5_result: Optional[Dict] = None,
    ps1_result: Optional[Dict] = None,
    residue_info: Optional[Dict] = None,
    dup_type: str = "Unknown",
) -> Dict:
    from backend.modules.pvs1 import evaluate_pvs1
    from backend.modules.bp1 import evaluate_bp1
    from backend.modules.pp3_bp4 import evaluate_pp3_bp4
    from backend.modules.bp7 import evaluate_bp7
    from backend.modules.frequency import evaluate_frequency_criteria
    from backend.modules.rna_review import evaluate_rna_review
    from backend.modules.splice_ps1_review import evaluate_splice_ps1_review
    from backend.modules.initiation_review import evaluate_initiation_review
    from backend.modules.utils import get_amino_acid_position, is_in_functional_domain

    results = {
        "variant": f"{gene} {c_notation} {p_notation}",
        "gene": gene,
        "c_notation": c_notation,
        "p_notation": p_notation,
        "criteria": {},
        "total_points": 0,
        "warnings": [],
        "has_functional_evidence": False,
        "classification_note": "",
        "residue_info": residue_info,
    }
    results["warnings"].append(
        "FIRST PASS - automatable ENIGMA VCEP v1.2 rules only. "
        "The following criteria are NOT automated and require expert review: "
        "PS4 (case-control data), PM3 (Fanconi anemia / trans variants), "
        "PP1 (co-segregation), BS2 (healthy carriers), BS4 (segregation absence). "
        "This automated result must not replace a full expert variant classification."
    )

    # ── Residue info (informational only) ──────────────────────────────
    if residue_info and residue_info.get("is_important_residue"):
        results["warnings"].append(residue_info["message"])

    # ── Step 1: Frequency (BA1 check) ──────────────────────────────────
    if gnomad_data:
        freq_criteria = evaluate_frequency_criteria(gnomad_data, variant_type)
        for crit_name, crit_data in freq_criteria.items():
            if crit_data.get("applies"):
                results["criteria"][crit_name] = crit_data
                results["total_points"] += crit_data["points"]
            elif crit_name == "PM2" and not crit_data.get("applies"):
                results["warnings"].append(crit_data["reason"])

        gnomad_info = freq_criteria.get("_gnomad_info")
        if gnomad_info:
            results["warnings"].append(gnomad_info["reason"])

        if "BA1" in results["criteria"]:
            cls, label, note = classify_by_points(0, has_ba1=True)
            results["predicted_class"] = cls
            results["predicted_label"] = label
            results["classification_note"] = note
            return results

    # ── Step 2: Table 9 - PS3/BS3 ─────────────────────────────────────
    if table9_result and table9_result.get("applies"):
        results["criteria"][table9_result["code"]] = {
            "applies": True,
            "strength": table9_result["strength"],
            "points": table9_result["points"],
            "reason": table9_result["reason"],
        }
        results["total_points"] += table9_result["points"]
        results["has_functional_evidence"] = True

    # ── Step 3: Table 4 - PVS1/PM5 ────────────────────────────────────
    pvs1 = evaluate_pvs1(
        gene, variant_type, p_notation,
        c_notation=c_notation,
        spliceai_score=spliceai_score,
        dup_type=dup_type,
    )
    if pvs1["applies"]:
        results["criteria"]["PVS1"] = pvs1
        results["total_points"] += pvs1["points"]
    elif pvs1.get("requires_rna") or variant_type.lower() in [
        "initiation_codon", "exon_deletion", "exon_duplication"
    ]:
        results["warnings"].append(pvs1["reason"])

    if pvs1.get("pm5_code") and pvs1.get("pm5_strength"):
        results["criteria"]["PM5_PTC"] = {
            "applies": True,
            "strength": pvs1["pm5_strength"],
            "points": pvs1["pm5_points"],
            "reason": f"Table 4: {pvs1['pm5_code']} for PTC in {pvs1.get('exon', 'unknown exon')}",
        }
        results["total_points"] += pvs1["pm5_points"]

    # ── Step 4: ST7 - PP4/BP5 ─────────────────────────────────────────
    if pp4_bp5_result and pp4_bp5_result.get("applies"):
        code = pp4_bp5_result["code"]
        results["criteria"][code] = {
            "applies": True,
            "strength": pp4_bp5_result["strength"],
            "points": pp4_bp5_result["points"],
            "reason": pp4_bp5_result["reason"],
        }
        results["total_points"] += pp4_bp5_result["points"]

    # ── Step 5: ST7 - PS1 ─────────────────────────────────────────────
    if ps1_result and ps1_result.get("applies"):
        results["criteria"]["PS1"] = {
            "applies": True,
            "strength": ps1_result["strength"],
            "points": ps1_result["points"],
            "reason": ps1_result["reason"],
        }
        results["total_points"] += ps1_result["points"]

    # ── Step 6: Computational predictions ──────────────────────────────
    suppress_benign_comp = results["has_functional_evidence"]
    pp3_bp4 = evaluate_pp3_bp4(
        gene, variant_type, p_notation,
        bayesdel_score=bayesdel_score,
        spliceai_score=spliceai_score,
    )
    for crit_name, crit_data in pp3_bp4.items():
        if crit_data.get("applies"):
            if crit_name == "PP3":
                results["criteria"][crit_name] = crit_data
                results["total_points"] += crit_data["points"]
            elif not suppress_benign_comp:
                results["criteria"][crit_name] = crit_data
                results["total_points"] += crit_data["points"]
            else:
                reason = "functional evidence in Table 9"
                results["warnings"].append(
                    f"{crit_name} not applied - {reason} overrides benign computational prediction."
                )

    # BP7: synonymous - only when no suppression
    if variant_type.lower() in ["synonymous", "silent", "intronic"] and not suppress_benign_comp:
        aa_pos = get_amino_acid_position(p_notation)
        in_domain = False
        if aa_pos:
            in_domain, _ = is_in_functional_domain(gene, aa_pos)
        bp4_met = "BP4" in results["criteria"] and results["criteria"]["BP4"].get("applies", False)

        bp7 = evaluate_bp7(
            variant_type,
            spliceai_score=spliceai_score,
            in_domain=in_domain,
            bp4_met=bp4_met,
            c_notation=c_notation,
        )
        if bp7["applies"]:
            results["criteria"]["BP7"] = bp7
            results["total_points"] += bp7["points"]

    # ── Step 7: BP1 ────────────────────────────────────────────────────
    bp1 = evaluate_bp1(gene, variant_type, p_notation, spliceai_score=spliceai_score)
    if bp1["applies"]:
        results["criteria"]["BP1"] = bp1
        results["total_points"] += bp1["points"]

    # ── Step 8: Warnings ───────────────────────────────────────────────
    if spliceai_score is None:
        results["warnings"].append(
            f"SpliceAI not available for {gene} {c_notation} - "
            "benign criteria BP1/BP4/BP7 require confirmed low score"
        )
    if bayesdel_score is None:
        results["warnings"].append(
            f"BayesDel_noAF not available for {gene} {c_notation}"
        )

    # ── Step 9: Classification ─────────────────────────────────────────
    cls, label, note = classify_by_enigma_combination(
        results["criteria"], results["total_points"]
    )
    results["predicted_class"] = cls
    results["predicted_label"] = label
    results["classification_note"] = note

    # ── Step 9b: first-pass note ───────────────────────────────────────
    if not results["classification_note"]:
        results["classification_note"] = (
            "First pass - automatable ENIGMA VCEP v1.2 rules only. "
            "Non-automated criteria (PS4, PM3, PP1, BS2, BS4) may affect final classification."
        )

    # ── Step 10: ACMG verification ─────────────────────────────────────
    acmg_note = verify_acmg_combination(
        results["criteria"], results["total_points"], results["predicted_class"]
    )
    if acmg_note:
        results["warnings"].append(acmg_note)

    results["rna_review"] = evaluate_rna_review(
        variant_type=variant_type,
        spliceai_score=spliceai_score,
        pvs1_result=pvs1,
        criteria=results["criteria"],
    )
    results["splice_ps1_review"] = evaluate_splice_ps1_review(
        variant_type=variant_type,
        spliceai_score=spliceai_score,
        ps1_result=ps1_result,
    )
    results["initiation_review"] = evaluate_initiation_review(
        variant_type=variant_type,
    )

    return results
