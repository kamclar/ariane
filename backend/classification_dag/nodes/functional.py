"""Calibrated functional-evidence DAG rule nodes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from backend.classification_dag.domain import CriterionDecision, CriterionFamilyResult, NormalizedVariant
from backend.classification_dag.nodes.support import decision
from backend.classification_dag.types import NodeResult
from backend.gene_policy import policy_name, policy_version, vcep_specification


def functional_decision_path(
    variant: NormalizedVariant,
    table9: Mapping[str, Any],
) -> dict[str, Any]:
    specification = vcep_specification(variant.gene)
    branch = (
        "intronic-silent"
        if variant.variant_type.lower() in {"intronic", "synonymous", "silent"}
        else "exonic-missense-inframe"
    )
    if branch == "intronic-silent":
        steps = [
            {
                "node_id": "func-rna-assay",
                "question": "Does the assay measure effects via both mRNA and protein?",
                "result": "yes",
                "observed": table9["reason"],
            }
        ]
        outcome = "func-rna-code"
    else:
        flag = str(table9.get("predicted_or_observed_splicing") or "").strip().upper()
        present = flag not in {"", "N", "NO"}
        steps = [
            {
                "node_id": "func-protein-splice",
                "question": "Is splicing predicted or observed?",
                "result": "yes" if present else "no",
                "observed": f"Table 9 predicted/observed splicing: {flag or 'not reported'}",
            },
            {
                "node_id": "func-protein-combined" if present else "func-protein-only",
                "question": (
                    "Assay measures both mRNA and protein effects"
                    if present
                    else "Assay measures protein-only effect"
                ),
                "result": "eligible",
                "observed": table9["reason"],
            },
        ]
        outcome = "func-protein-code"
    return {
        "tree_id": "figure-1c",
        "tree_version": "ENIGMA VCEP 1.2.0",
        "branch_id": branch,
        "criterion": table9["code"],
        "outcome": "applied",
        "outcome_node": outcome,
        "steps": steps,
        "sources": [
            {
                "source_id": "enigma-v1.2-specifications",
                "label": (
                    f"{policy_name(variant.gene)} "
                    f"v{policy_version(variant.gene)} Specifications"
                ),
                "url": specification["url"],
                "location": "Figure 1C",
                "figure_url": "/static/enigma/figure-1c-functional.jpeg",
            },
            {
                "source_id": "enigma-v1.2-table9",
                "label": "ENIGMA Specifications Table 9 v1.2",
                "url": "https://cspec.genome.network/cspec/File/id/0a35d6a8-5050-44b6-8a9d-babe8cdc06b2/data",
                "location": f"{variant.gene}:{variant.c_notation}",
            },
        ],
    }


@dataclass(frozen=True)
class FunctionalCriteriaNode:
    id: str = "rule.functional.table9"
    version: str = "1"
    requires: frozenset[str] = frozenset({"normalized_variant", "table9_evidence"})
    provides: frozenset[str] = frozenset({"functional_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = inputs["normalized_variant"]
        evidence = inputs["table9_evidence"]
        table9 = evidence.value or {}
        decisions: tuple[CriterionDecision, ...] = ()
        if table9.get("applies"):
            payload = {
                "applies": True,
                "strength": table9["strength"],
                "points": table9["points"],
                "reason": table9["reason"],
                "decision_path": functional_decision_path(variant, table9),
            }
            decisions = (
                decision(
                    table9["code"],
                    payload,
                    gene=variant.gene,
                    family_id=self.id,
                    evidence_item_ids=("enigma_table9",),
                ),
            )
        family = CriterionFamilyResult(
            self.id,
            criteria=decisions,
            has_functional_evidence=bool(decisions),
        )
        return NodeResult.succeeded({"functional_family": family})
