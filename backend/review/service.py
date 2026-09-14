"""Build the amended working classification from reviewed manual evidence."""

from __future__ import annotations

from typing import Any, Dict, List, Mapping

from backend.criteria.bp7_rna import evaluate_bp7_rna_variant_context
from backend.criteria.evidence_interactions import (
    apply_rna_interactions,
    automatic_functional_interactions,
    clinical_functional_risk_interactions,
)
from backend.policy.classification import classify_by_enigma_combination
from backend.policy.gene import resolve_policy_gene, rule_is_applicable
from backend.presentation.criterion_order import criterion_sort_key
from backend.review.definitions import MANUAL_CRITERIA, STRENGTH_POINTS, STRUCTURED_CURATED_CODES
from backend.review.strength import (
    _clinical_lr_record_is_complete,
    evaluate_bs4_likelihood_ratio,
    suggest_strength,
)

def evaluate_manual_evidence(
    base_criteria: List[Dict[str, Any]],
    manual_criteria: List[Dict[str, Any]],
    variant_context: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    policy_gene = resolve_policy_gene(
        str((variant_context or {}).get("gene") or "").strip().upper() or None
    )
    combined = {
        criterion["name"]: {
            "applies": criterion.get("applies", True),
            "strength": criterion.get("strength"),
            "points": criterion.get("points", 0),
            "reason": criterion.get("reason", ""),
            "single_strong_likely_benign_eligible": criterion.get(
                "single_strong_likely_benign_eligible", False
            ),
            "single_strong_likely_benign_basis": criterion.get(
                "single_strong_likely_benign_basis", ""
            ),
            "independent_evidence_contribution_count": criterion.get(
                "independent_evidence_contribution_count", 0
            ),
        }
        for criterion in base_criteria
        if criterion.get("applies", True)
    }
    results = []

    enabled_manual = [item for item in manual_criteria if item.get("enabled")]
    bp7_context_criteria = list(base_criteria)
    for item in enabled_manual:
        if item.get("code") != "BS3":
            continue
        manual_bs3_strength = suggest_strength(
            "BS3",
            item.get("evidence", {}),
            variant_context=variant_context,
            base_criteria=base_criteria,
        )
        if manual_bs3_strength:
            bp7_context_criteria.append({
                "name": "BS3",
                "applies": True,
                "strength": manual_bs3_strength,
                "decision_path": {
                    "sources": [{
                        "source_id": "manual-enigma-vcep-functional-review",
                    }],
                },
            })
    clinical_lr_is_used = any(code in combined for code in {"PP4", "BP5"}) or any(
        item.get("code") in {"PP4", "BP5"} for item in enabled_manual
    )
    if clinical_lr_is_used:
        for item in enabled_manual:
            if item.get("code") not in {"PP1", "PS4"}:
                continue
            evidence = item.get("evidence", {})
            rationale = str(evidence.get("independence_rationale") or "").strip()
            if evidence.get("independent_from_pp4_bp5") is not True or not rationale:
                raise ValueError(
                    f"{item['code']} cannot be combined with PP4/BP5 until the reviewer "
                    "confirms independent observations and records an independence rationale"
                )

    for item in manual_criteria:
        code = item["code"]
        if item.get("enabled") and not rule_is_applicable(policy_gene, code):
            raise ValueError(
                f"{code} is not applicable under the configured VCEP policy for {policy_gene}"
            )
        definition = MANUAL_CRITERIA[code]
        if item.get("override_strength") not in {None, ""}:
            raise ValueError(
                "Manual strength overrides are not permitted. "
                f"ARIANE derives criterion strength from the configured VCEP "
                f"policy for {policy_gene}."
            )
        suggested = suggest_strength(
            code,
            item.get("evidence", {}),
            variant_context=variant_context,
            base_criteria=(
                bp7_context_criteria if code == "BP7_RNA" else base_criteria
            ),
        )
        evidence = item.get("evidence", {})
        if (
            code in {"PP4", "BP5"}
            and item.get("enabled")
            and not _clinical_lr_record_is_complete(evidence)
        ):
            raise ValueError(
                f"{code} requires a clinical LR value and scale, recorded source, "
                "and clinical data summary"
            )
        if code == "BP7_RNA" and item.get("enabled") and not suggested:
            context_result = evaluate_bp7_rna_variant_context(
                variant_context, bp7_context_criteria
            )
            if not context_result["eligible"]:
                raise ValueError(context_result["reason"])
        if code in STRUCTURED_CURATED_CODES - {"PP4", "BP5"} and item.get("enabled") and not suggested:
            raise ValueError(
                f"{code} requires a complete structured curated evidence record"
            )
        selected = suggested
        applies = bool(item.get("enabled") and selected)
        points = STRENGTH_POINTS.get(selected, 0)
        if definition["direction"] == "benign":
            points *= -1
        reason = (
            f"Reviewer-curated evidence meets {code} {suggested} requirements"
            if suggested
            else f"Reviewer-curated evidence does not meet the structured {code} requirements"
        )
        single_strong_eligible = False
        single_strong_basis = ""
        contribution_count = 0
        if code == "BS4" and applies:
            _, has_multiple_lrs, contribution_count = evaluate_bs4_likelihood_ratio(
                evidence
            )
            single_strong_eligible = bool(
                selected == "Strong" and has_multiple_lrs
            )
            if single_strong_eligible:
                single_strong_basis = (
                    "Multiple independently identified segregation likelihood ratios "
                    "contribute to BS4 Strong"
                )
                reason += "; multiple independent LR components satisfy the ENIGMA Table 3 single-Strong condition"
            elif selected == "Strong":
                reason += "; BS4 Strong is valid, but the ENIGMA Table 3 single-Strong condition is not documented"
        if code == "BP5" and applies:
            evidence_types = {
                str(item).strip().lower()
                for item in evidence.get("clinical_evidence_types", [])
                if str(item).strip()
            }
            contribution_count = len(evidence_types)
            single_strong_eligible = bool(
                selected == "Strong"
                and contribution_count >= 2
                and evidence.get("independence_review_confirmed") is True
            )
            if single_strong_eligible:
                single_strong_basis = (
                    "At least two independently reviewed clinical evidence types "
                    "contribute to BP5 Strong"
                )
                reason += "; multiple clinical evidence types satisfy the ENIGMA Table 3 single-Strong condition"
            elif selected == "Strong":
                reason += "; BP5 Strong is valid, but the ENIGMA Table 3 single-Strong condition is not documented"
        results.append(
            {
                "code": code,
                "applies": applies,
                "suggested_strength": suggested,
                "selected_strength": selected,
                "points": points if applies else 0,
                "reason": reason,
                "threshold_note": definition["threshold"],
                "overridden": False,
                "notes": item.get("notes", ""),
                "references": item.get("references", []),
                "single_strong_likely_benign_eligible": single_strong_eligible,
                "single_strong_likely_benign_basis": single_strong_basis,
                "independent_evidence_contribution_count": contribution_count,
            }
        )
        if applies:
            if code == "PVS1_INIT":
                combined.pop("PP3", None)
            output_code = "PS1" if code == "PS1_PROTEIN" else code
            if output_code in combined:
                raise ValueError(
                    f"{output_code} is already present in the automated result and cannot be counted twice"
                )
            combined[output_code] = {
                "applies": True,
                "strength": selected,
                "points": points,
                "reason": reason,
                "single_strong_likely_benign_eligible": single_strong_eligible,
                "single_strong_likely_benign_basis": single_strong_basis,
                "independent_evidence_contribution_count": contribution_count,
            }

    applied_manual_codes = {
        result["code"] for result in results if result["applies"]
    }
    evidence_interactions = apply_rna_interactions(
        combined, applied_manual_codes
    )
    evidence_interactions.extend(automatic_functional_interactions(combined))
    evidence_interactions.extend(clinical_functional_risk_interactions(combined))
    total_points = sum(c.get("points", 0) for c in combined.values())
    predicted_class, label, note = classify_by_enigma_combination(
        combined, total_points, gene=policy_gene
    )
    return {
        "predicted_class": predicted_class,
        "predicted_label": label,
        "total_points": total_points,
        "classification_note": note,
        "manual_criteria": sorted(
            results, key=lambda item: criterion_sort_key(item["code"])
        ),
        "evidence_interactions": evidence_interactions,
    }
