"""Validate completeness and readiness of manual evidence forms."""

from __future__ import annotations

from typing import Any, Dict, List, Mapping

from backend.criteria.bp7_rna import evaluate_bp7_rna_variant_context
from backend.policy.gene import resolve_policy_gene, rule_is_applicable
from backend.review.definitions import STRUCTURED_CURATED_CODES
from backend.review.strength import _clinical_lr_record_is_complete, suggest_strength

def manual_criterion_statuses(
    base_criteria: List[Dict[str, Any]],
    manual_criteria: List[Dict[str, Any]],
    variant_context: Mapping[str, Any] | None = None,
) -> List[Dict[str, Any]]:
    """Return backend-derived form readiness without assigning a criterion."""
    policy_gene = resolve_policy_gene(
        str((variant_context or {}).get("gene") or "").strip().upper() or None
    )
    enabled = [item for item in manual_criteria if item.get("enabled")]
    bp7_context_criteria = list(base_criteria)
    for item in enabled:
        if item.get("code") != "BS3":
            continue
        strength = suggest_strength(
            "BS3",
            item.get("evidence", {}),
            variant_context=variant_context,
            base_criteria=base_criteria,
        )
        if strength:
            bp7_context_criteria.append({
                "name": "BS3",
                "applies": True,
                "strength": strength,
                "decision_path": {
                    "sources": [{"source_id": "manual-enigma-vcep-functional-review"}],
                },
            })

    clinical_lr_is_used = any(
        criterion.get("name") in {"PP4", "BP5"}
        and criterion.get("applies", True)
        for criterion in base_criteria
    ) or any(item.get("code") in {"PP4", "BP5"} for item in enabled)

    statuses = []
    for item in manual_criteria:
        code = item["code"]
        if not item.get("enabled"):
            statuses.append({
                "code": code,
                "status": "not_started",
                "message": "Select this criterion to review its evidence.",
                "suggested_strength": None,
            })
            continue

        evidence = item.get("evidence", {})
        issues = []
        if not rule_is_applicable(policy_gene, code):
            issues.append(
                f"{code} is not applicable under the configured VCEP policy for {policy_gene}."
            )
        if not str(item.get("notes") or "").strip():
            issues.append("Add evidence notes.")
        references = item.get("references", [])
        if not any(str(reference).strip() for reference in references):
            issues.append("Add at least one evidence reference.")

        if (
            clinical_lr_is_used
            and code in {"PP1", "PS4"}
            and (
                evidence.get("independent_from_pp4_bp5") is not True
                or not str(evidence.get("independence_rationale") or "").strip()
            )
        ):
            issues.append(
                "Confirm independence from PP4/BP5 and record the rationale."
            )

        context_criteria = (
            bp7_context_criteria if code == "BP7_RNA" else base_criteria
        )
        suggested = suggest_strength(
            code,
            evidence,
            variant_context=variant_context,
            base_criteria=context_criteria,
        )
        if code in {"PP4", "BP5"} and not _clinical_lr_record_is_complete(evidence):
            issues.append(
                "Add the clinical LR value and scale, source, and clinical data summary."
            )
        elif code == "BP7_RNA" and not suggested:
            context_result = evaluate_bp7_rna_variant_context(
                variant_context, bp7_context_criteria
            )
            if not context_result["eligible"]:
                issues.append(context_result["reason"])
            else:
                issues.append("Complete the structured BP7 RNA evidence record.")
        elif code in STRUCTURED_CURATED_CODES and not suggested:
            issues.append(f"Complete the structured {code} evidence record.")
        elif not suggested:
            issues.append(
                "The entered evidence does not yet meet the configured ENIGMA threshold and stipulations."
            )

        output_code = "PS1" if code == "PS1_PROTEIN" else code
        if any(
            criterion.get("name") == output_code and criterion.get("applies", True)
            for criterion in base_criteria
        ):
            issues.append(
                f"{output_code} is already present in the automated result and cannot be counted twice."
            )

        if issues:
            statuses.append({
                "code": code,
                "status": "incomplete",
                "message": " ".join(dict.fromkeys(issues)),
                "suggested_strength": suggested,
            })
        else:
            statuses.append({
                "code": code,
                "status": "ready",
                "message": (
                    f"Required evidence is complete. The backend derives {code} "
                    f"{suggested}."
                ),
                "suggested_strength": suggested,
            })
    return statuses

