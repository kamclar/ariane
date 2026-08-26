# ============================================================
# ARIANE classifier - main evaluation orchestrator
#
# Evidence hierarchy (ENIGMA VCEP v1.2):
#   1. BA1 check - stand-alone benign, if met → class 1, stop
#   2. Table 9 - PS3/BS3 calibrated functional evidence
#   3. Table 4 - PVS1/PM5 structural rules; ST2 + Appendix E - PVS1 RNA
#   4. gnomAD - BS1, PM2
#   5. Local clinical-LR snapshot - PP4/BP5
#   6. Approved protein-PS1 registry; ST7 supplies review candidates only
#   7. SpliceAI/BayesDel - PP3/BP4/BP7 per variant-type decision tree
#   8. BP1 - domain check
#   9. Classification from adapted ACMG/AMP combinations
#  10. Tavtigian 2020 points only for contradictory evidence
# ============================================================
from typing import Optional, Dict, List, Tuple

from backend.modules.evidence_interactions import (
    apply_automatic_rna_interactions,
    automatic_functional_interactions,
    clinical_functional_risk_interactions,
    pvs1_prediction_deduplication,
)
from backend.modules.pvs1_rna import evaluate_pvs1_rna
from backend.modules.spliceai_policy import compare_table9_spliceai
from backend.gene_policy import policy_name, policy_version, vcep_specification


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


def _classify_benign_combination(criteria: Dict) -> Optional[tuple]:
    counts = {"very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0}
    has_eligible_single_strong = False

    for criterion in criteria.values():
        if criterion.get("points", 0) >= 0:
            continue
        strength = _criterion_strength(criterion)
        counts[strength] += 1
        if (
            strength == "strong"
            and criterion.get("single_strong_likely_benign_eligible") is True
        ):
            has_eligible_single_strong = True

    very_strong = counts["very_strong"]
    strong = counts["strong"]
    moderate = counts["moderate"]
    supporting = counts["supporting"]

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


def classify_by_enigma_combination(criteria: Dict, points: int) -> tuple:
    """
    Apply ENIGMA VCEP v1.2 Table 3 combinations by default.
    Tavtigian points are reserved for contradictory benign/pathogenic evidence.
    """
    ba1 = criteria.get("BA1")
    if ba1 and ba1.get("applies", True):
        return (
            1,
            "Benign",
            "BA1 is stand-alone benign evidence. Other criteria remain visible for audit "
            "but do not replace the BA1 classification.",
        )

    has_pathogenic = any(criterion.get("points", 0) > 0 for criterion in criteria.values())
    has_benign = any(criterion.get("points", 0) < 0 for criterion in criteria.values())

    if has_pathogenic and has_benign:
        cls, label, _ = classify_by_points(points)
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
            and _criterion_strength(criterion) == "strong"
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
    exon_cnv_result: Optional[Dict] = None,
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
    from backend.modules.protein_ps1_review import evaluate_protein_ps1_review
    from backend.modules.initiation_review import evaluate_initiation_review
    from backend.modules.utils import get_amino_acid_position, is_in_functional_domain

    results = {
        "variant": f"{gene} {c_notation} {p_notation}",
        "gene": gene,
        "c_notation": c_notation,
        "p_notation": p_notation,
        "criteria": {},
        "excluded_criteria": {},
        "total_points": 0,
        "warnings": [],
        "has_functional_evidence": False,
        "classification_note": "",
        "evidence_direction": "none",
        "mixed_evidence": False,
        "pathogenic_points": 0,
        "benign_points": 0,
        "evidence_interactions": [],
        "residue_info": residue_info,
    }
    results["warnings"].append(
        "FIRST PASS - automatable ENIGMA VCEP v1.2 rules only. "
        "The following criteria are NOT automated and require expert review: "
        "PS4 (case-control data), PM3 (Fanconi anemia / trans variants), "
        "PP1 (co-segregation), BS2 (healthy carriers), BS4 (segregation absence). "
        "This automated result must not replace a full expert variant classification."
    )

    # Table 9 is authoritative for its PS3/BS3 recommendation and records the
    # SpliceAI context used in that functional review.  Figure 1A and PS1 use
    # the configured ENIGMA-compatible SpliceAI source; Table 9 never replaces
    # an unavailable or different prediction.
    effective_spliceai_score = spliceai_score
    _, table9_spliceai_warnings = compare_table9_spliceai(
        gene, spliceai_score, table9_result
    )
    results["warnings"].extend(table9_spliceai_warnings)

    # ── Residue info (informational only) ──────────────────────────────
    if residue_info and residue_info.get("is_important_residue"):
        results["warnings"].append(residue_info["message"])

    # ── Step 1: Frequency (BA1 check) ──────────────────────────────────
    if gnomad_data:
        freq_criteria = evaluate_frequency_criteria(
            gnomad_data,
            variant_type,
            gene=gene,
            c_notation=c_notation,
        )
        for crit_name, crit_data in freq_criteria.items():
            if crit_name.startswith("_"):
                continue
            if crit_data.get("applies"):
                results["criteria"][crit_name] = crit_data
                results["total_points"] += crit_data["points"]
            elif crit_name == "PM2" and not crit_data.get("applies"):
                results["warnings"].append(crit_data["reason"])

        results["excluded_criteria"].update(
            freq_criteria.get("_excluded_criteria", {})
        )

        gnomad_info = freq_criteria.get("_gnomad_info")
        if gnomad_info:
            results["warnings"].append(gnomad_info["reason"])

        if "BA1" in results["criteria"]:
            cls, label, note = classify_by_points(0, has_ba1=True)
            results["predicted_class"] = cls
            results["predicted_label"] = label
            results["classification_note"] = note
            return results

    # Exon deletions use the generic ENIGMA Appendix G population path.  This
    # input contains no variant-specific criterion assignments.
    if exon_cnv_result and exon_cnv_result.get("found"):
        for criterion in exon_cnv_result.get("criteria", []):
            code = criterion["code"]
            results["criteria"][code] = {
                "applies": True,
                "strength": criterion["strength"],
                "points": criterion["points"],
                "reason": criterion["reason"],
                "source": criterion["source"],
            }
            results["total_points"] += criterion["points"]

    # ── Step 2: Table 9 - PS3/BS3 ─────────────────────────────────────
    if table9_result and table9_result.get("applies"):
        code = table9_result["code"]
        functional_branch = (
            "intronic-silent"
            if variant_type.lower() in {"intronic", "synonymous", "silent"}
            else "exonic-missense-inframe"
        )
        if functional_branch == "intronic-silent":
            functional_steps = [
                {
                    "node_id": "func-rna-assay",
                    "question": "Does the assay measure effects via both mRNA and protein?",
                    "result": "yes",
                    "observed": table9_result["reason"],
                }
            ]
            functional_outcome = "func-rna-code"
        else:
            splice_flag = str(
                table9_result.get("predicted_or_observed_splicing") or ""
            ).strip().upper()
            splicing_present = splice_flag not in {"", "N", "NO"}
            functional_steps = [
                {
                    "node_id": "func-protein-splice",
                    "question": "Is splicing predicted or observed?",
                    "result": "yes" if splicing_present else "no",
                    "observed": (
                        f"Table 9 predicted/observed splicing: {splice_flag or 'not reported'}"
                    ),
                },
                {
                    "node_id": (
                        "func-protein-combined" if splicing_present else "func-protein-only"
                    ),
                    "question": (
                        "Assay measures both mRNA and protein effects"
                        if splicing_present
                        else "Assay measures protein-only effect"
                    ),
                    "result": "eligible",
                    "observed": table9_result["reason"],
                },
            ]
            functional_outcome = "func-protein-code"
        results["criteria"][code] = {
            "applies": True,
            "strength": table9_result["strength"],
            "points": table9_result["points"],
            "reason": table9_result["reason"],
            "decision_path": {
                "tree_id": "figure-1c",
                "tree_version": "ENIGMA VCEP 1.2.0",
                "branch_id": functional_branch,
                "criterion": code,
                "outcome": "applied",
                "outcome_node": functional_outcome,
                "steps": functional_steps,
                "sources": [
                    {
                        "source_id": "enigma-v1.2-specifications",
                        "label": (
                            f"{policy_name(gene)} v{policy_version(gene)} Specifications"
                        ),
                        "url": vcep_specification(gene)["url"],
                        "location": "Figure 1C",
                        "figure_url": "/static/enigma/figure-1c-functional.jpeg",
                    },
                    {
                        "source_id": "enigma-v1.2-table9",
                        "label": "ENIGMA Specifications Table 9 v1.2",
                        "url": "https://cspec.genome.network/cspec/File/id/0a35d6a8-5050-44b6-8a9d-babe8cdc06b2/data",
                        "location": f"{gene}:{c_notation}",
                    },
                ],
            },
        }
        results["total_points"] += table9_result["points"]
        results["has_functional_evidence"] = True

    # ── Step 3: Table 4 - PVS1/PM5 ────────────────────────────────────
    pvs1 = evaluate_pvs1(
        gene, variant_type, p_notation,
        c_notation=c_notation,
        spliceai_score=effective_spliceai_score,
        dup_type=dup_type,
    )
    pvs1_rna = evaluate_pvs1_rna(gene, c_notation)
    if pvs1["applies"]:
        results["criteria"]["PVS1"] = pvs1
        results["total_points"] += pvs1["points"]
    elif pvs1_rna.get("applies"):
        results["criteria"]["PVS1_RNA"] = pvs1_rna
        results["total_points"] += pvs1_rna["points"]
        results["has_functional_evidence"] = True
    elif pvs1.get("requires_rna") or variant_type.lower() in [
        "nonsense", "frameshift", "splice_site", "initiation_codon",
        "exon_deletion", "exon_duplication"
    ]:
        results["warnings"].append(pvs1["reason"])
        if "N/A" in str(pvs1.get("pvs1_code") or ""):
            results["excluded_criteria"]["PVS1"] = {
                "applies": False,
                "strength": "N/A",
                "points": 0,
                "reason": pvs1["reason"],
                "source": pvs1.get("source", ""),
            }

    if (
        not pvs1.get("applies")
        and pvs1_rna.get("source_record")
        and not pvs1_rna.get("applies")
    ):
        results["warnings"].append(pvs1_rna["reason"])

    if pvs1.get("pm5_code") and pvs1.get("pm5_strength"):
        results["criteria"]["PM5_PTC"] = {
            "applies": True,
            "strength": pvs1["pm5_strength"],
            "points": pvs1["pm5_points"],
            "reason": (
                f"Table 4: {pvs1['pm5_code']} for PTC in "
                f"{pvs1.get('pm5_exon') or 'unknown exon'}"
            ),
        }
        results["total_points"] += pvs1["pm5_points"]

    # ── Step 4: local clinical-LR snapshot - PP4/BP5 ───────────────────
    if pp4_bp5_result and pp4_bp5_result.get("applies"):
        code = pp4_bp5_result["code"]
        results["criteria"][code] = dict(pp4_bp5_result)
        results["total_points"] += pp4_bp5_result["points"]

    # ── Step 5: approved protein PS1 references ───────────────────────
    if ps1_result and ps1_result.get("applies"):
        results["criteria"]["PS1"] = {
            "applies": True,
            "strength": ps1_result["strength"],
            "points": ps1_result["points"],
            "reason": ps1_result["reason"],
        }
        results["total_points"] += ps1_result["points"]

    # ── Step 6: Computational predictions ──────────────────────────────
    pp3_bp4 = evaluate_pp3_bp4(
        gene, variant_type, p_notation,
        bayesdel_score=bayesdel_score,
        spliceai_score=effective_spliceai_score,
        c_notation=c_notation,
    )
    for crit_name, crit_data in pp3_bp4.items():
        if crit_data.get("applies"):
            if crit_name == "PP3" and pvs1.get("applies"):
                results["warnings"].append(
                    "PP3 not applied because PVS1 is met; ENIGMA does not stack predictive PP3 with PVS1."
                )
                results["evidence_interactions"].append(
                    pvs1_prediction_deduplication()
                )
            else:
                results["criteria"][crit_name] = crit_data
                results["total_points"] += crit_data["points"]

    # BP7 is independent of calibrated protein functional evidence.
    if variant_type.lower() in ["synonymous", "silent", "intronic"]:
        aa_pos = get_amino_acid_position(p_notation)
        in_domain = False
        if aa_pos:
            in_domain, _ = is_in_functional_domain(gene, aa_pos)
        bp4_met = "BP4" in results["criteria"] and results["criteria"]["BP4"].get("applies", False)

        bp7 = evaluate_bp7(
            variant_type,
            spliceai_score=effective_spliceai_score,
            in_domain=in_domain,
            bp4_met=bp4_met,
            c_notation=c_notation,
            gene=gene,
        )
        if bp7["applies"]:
            results["criteria"]["BP7"] = bp7
            results["total_points"] += bp7["points"]

    # ── Step 7: BP1 ────────────────────────────────────────────────────
    bp1 = evaluate_bp1(
        gene, variant_type, p_notation, spliceai_score=effective_spliceai_score
    )
    if bp1["applies"]:
        results["criteria"]["BP1"] = bp1
        results["total_points"] += bp1["points"]

    results["evidence_interactions"].extend(
        apply_automatic_rna_interactions(results["criteria"])
    )
    results["evidence_interactions"].extend(
        automatic_functional_interactions(results["criteria"])
    )
    results["evidence_interactions"].extend(
        clinical_functional_risk_interactions(results["criteria"])
    )

    # ── Step 8: Warnings ───────────────────────────────────────────────
    figure1a_types = {
        "missense", "inframe_deletion", "inframe_insertion",
        "inframe_delins", "delins", "synonymous", "silent", "intronic",
    }
    if (
        effective_spliceai_score is None
        and variant_type.lower() in figure1a_types
    ):
        results["warnings"].append(
            f"Figure 1A bioinformatic result unavailable for {gene} {c_notation}: "
            "SpliceAI is unavailable. Missing data was not treated as an ENIGMA "
            "prediction band; PP3, BP4, BP1 and BP7 were not applied."
        )
    bayesdel_applicable_types = {
        "missense",
        "inframe_deletion", "inframe_insertion", "inframe_delins", "delins",
    }
    if bayesdel_score is None and variant_type.lower() in bayesdel_applicable_types:
        results["warnings"].append(
            f"BayesDel_noAF not available for {gene} {c_notation}"
        )

    # ── Step 9: Classification ─────────────────────────────────────────
    # Interaction rules may remove weaker evidence after it was initially
    # evaluated, so the final score is always recomputed from retained codes.
    results["total_points"] = sum(
        criterion.get("points", 0)
        for criterion in results["criteria"].values()
    )
    results["pathogenic_points"] = sum(
        max(criterion.get("points", 0), 0)
        for criterion in results["criteria"].values()
    )
    results["benign_points"] = sum(
        min(criterion.get("points", 0), 0)
        for criterion in results["criteria"].values()
    )
    results["mixed_evidence"] = bool(
        results["pathogenic_points"] > 0 and results["benign_points"] < 0
    )
    results["evidence_direction"] = (
        "mixed" if results["mixed_evidence"]
        else "pathogenic" if results["pathogenic_points"] > 0
        else "benign" if results["benign_points"] < 0
        else "none"
    )
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
        gene=gene,
        variant_type=variant_type,
        spliceai_score=spliceai_score,
        pvs1_result=pvs1_rna if pvs1_rna.get("applies") else pvs1,
        criteria=results["criteria"],
    )
    results["splice_ps1_review"] = evaluate_splice_ps1_review(
        gene=gene,
        variant_type=variant_type,
        spliceai_score=spliceai_score,
        ps1_result=ps1_result,
    )
    results["protein_ps1_review"] = evaluate_protein_ps1_review(ps1_result, gene=gene)
    results["initiation_review"] = evaluate_initiation_review(
        gene=gene,
        variant_type=variant_type,
    )

    return results
