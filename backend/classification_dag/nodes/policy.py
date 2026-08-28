"""Evidence interaction and final ENIGMA combination DAG nodes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from backend.classification_dag.nodes.support import (
    FIRST_PASS_NOTE,
    FIRST_PASS_WARNING,
    RetainedEvidence,
    criteria_dict,
    excluded_by_policy,
)
from backend.classification_dag.types import NodeResult
from backend.gene_policy import rule_is_applicable


@dataclass(frozen=True)
class EvidenceInteractionNode:
    id: str = "policy.evidence_interactions"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {
            "classification_inputs",
            "splice_context",
            "frequency_family",
            "exon_cnv_family",
            "functional_family",
            "pvs1_family",
            "clinical_lr_family",
            "protein_ps1_family",
            "bioinformatic_family",
        }
    )
    provides: frozenset[str] = frozenset({"retained_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        from backend.modules.evidence_interactions import (
            apply_automatic_rna_interactions,
            automatic_functional_interactions,
            clinical_functional_risk_interactions,
        )

        ci = inputs["classification_inputs"]
        splice = inputs["splice_context"]
        frequency = inputs["frequency_family"]
        families = tuple(
            inputs[key]
            for key in (
                "frequency_family",
                "exon_cnv_family",
                "functional_family",
                "pvs1_family",
                "clinical_lr_family",
                "protein_ps1_family",
                "bioinformatic_family",
            )
        )
        criteria = {}
        excluded = {}
        not_applicable = {}
        for family in families:
            for item in family.criteria:
                if rule_is_applicable(ci.gene, item.code):
                    criteria[item.code] = item
                else:
                    excluded[item.code] = excluded_by_policy(item, ci.gene)
            for item in family.excluded_criteria:
                excluded[item.code] = item
            for item in family.not_applicable_criteria:
                not_applicable[item.code] = item

        for code in set(criteria) | set(excluded):
            not_applicable.pop(code, None)

        warnings: list[str] = [FIRST_PASS_WARNING]
        warnings.extend(splice.warnings)
        if ci.residue_info and ci.residue_info.get("is_important_residue"):
            warnings.append(ci.residue_info["message"])
        warnings.extend(frequency.warnings)
        warnings.extend(inputs["clinical_lr_family"].warnings)

        if "BA1" in criteria:
            retained = RetainedEvidence(
                criteria=(criteria["BA1"],),
                excluded_criteria=tuple(excluded.values()),
                not_applicable_criteria=tuple(not_applicable.values()),
                warnings=tuple(warnings),
                evidence_interactions=(),
                has_functional_evidence=False,
                metadata={"ba1_terminal": True},
            )
            return NodeResult.succeeded(
                {"retained_evidence": retained},
                provenance={"terminal": "BA1"},
            )

        for family in families[1:]:
            if family.family_id != "rule.clinical_lr":
                warnings.extend(family.warnings)
        interactions: list[Mapping[str, Any]] = list(
            inputs["bioinformatic_family"].evidence_interactions
        )
        public = criteria_dict(tuple(criteria.values()))
        interactions.extend(apply_automatic_rna_interactions(public))
        interactions.extend(automatic_functional_interactions(public))
        interactions.extend(clinical_functional_risk_interactions(public))
        retained_codes = set(public)
        retained_decisions = tuple(
            item for code, item in criteria.items() if code in retained_codes
        )
        retained = RetainedEvidence(
            criteria=retained_decisions,
            excluded_criteria=tuple(excluded.values()),
            not_applicable_criteria=tuple(not_applicable.values()),
            warnings=tuple(warnings),
            evidence_interactions=tuple(interactions),
            has_functional_evidence=any(
                family.has_functional_evidence for family in families
            ),
            metadata={
                "ba1_terminal": False,
                "pvs1": inputs["pvs1_family"].metadata["pvs1"],
                "pvs1_rna": inputs["pvs1_family"].metadata["pvs1_rna"],
                "ps1": inputs["protein_ps1_family"].metadata["ps1"],
            },
        )
        return NodeResult.succeeded(
            {"retained_evidence": retained},
            provenance={"retained": [item.code for item in retained_decisions]},
        )


@dataclass(frozen=True)
class ClassificationPolicyNode:
    id: str = "policy.enigma_combination"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {"normalized_variant", "classification_inputs", "retained_evidence"}
    )
    provides: frozenset[str] = frozenset({"classification_without_reviews"})

    def evaluate(self, context, inputs) -> NodeResult:
        from backend.classification_dag.policy import (
            classify_by_enigma_combination,
            classify_by_points,
            verify_acmg_combination,
        )

        variant = inputs["normalized_variant"]
        ci = inputs["classification_inputs"]
        retained = inputs["retained_evidence"]
        criteria = criteria_dict(retained.criteria)
        excluded = criteria_dict(retained.excluded_criteria)
        not_applicable = criteria_dict(retained.not_applicable_criteria)
        result: dict[str, Any] = {
            "variant": f"{variant.gene} {variant.c_notation} {variant.p_notation}",
            "gene": variant.gene,
            "c_notation": variant.c_notation,
            "p_notation": variant.p_notation,
            "criteria": criteria,
            "excluded_criteria": excluded,
            "not_applicable_criteria": not_applicable,
            "total_points": sum(item.points for item in retained.criteria),
            "warnings": list(retained.warnings),
            "has_functional_evidence": retained.has_functional_evidence,
            "classification_note": "",
            "evidence_direction": "none",
            "mixed_evidence": False,
            "pathogenic_points": 0,
            "benign_points": 0,
            "evidence_interactions": [
                dict(item) for item in retained.evidence_interactions
            ],
            "residue_info": ci.residue_info,
        }
        if retained.metadata.get("ba1_terminal"):
            cls, label, note = classify_by_points(
                0,
                has_ba1=True,
                gene=variant.gene,
            )
            result["predicted_class"] = cls
            result["predicted_label"] = label
            result["classification_note"] = note
            return NodeResult.succeeded(
                {"classification_without_reviews": result},
                provenance={"method": "BA1-stand-alone"},
            )

        result["pathogenic_points"] = sum(
            max(item.points, 0) for item in retained.criteria
        )
        result["benign_points"] = sum(
            min(item.points, 0) for item in retained.criteria
        )
        result["mixed_evidence"] = bool(
            result["pathogenic_points"] > 0 and result["benign_points"] < 0
        )
        result["evidence_direction"] = (
            "mixed"
            if result["mixed_evidence"]
            else "pathogenic"
            if result["pathogenic_points"] > 0
            else "benign"
            if result["benign_points"] < 0
            else "none"
        )
        cls, label, note = classify_by_enigma_combination(
            criteria,
            result["total_points"],
            gene=variant.gene,
        )
        result["predicted_class"] = cls
        result["predicted_label"] = label
        result["classification_note"] = note or FIRST_PASS_NOTE
        acmg_note = verify_acmg_combination(
            criteria,
            result["total_points"],
            cls,
        )
        if acmg_note:
            result["warnings"].append(acmg_note)
        return NodeResult.succeeded(
            {"classification_without_reviews": result},
            provenance={
                "method": (
                    "Tavtigian-2020"
                    if result["mixed_evidence"]
                    else "ENIGMA-Table-3"
                )
            },
        )
