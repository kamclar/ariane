"""Clinical likelihood-ratio and protein PS1 DAG rule nodes."""

from __future__ import annotations

from dataclasses import dataclass

from backend.classification_dag.domain import CriterionFamilyResult
from backend.classification_dag.nodes.support import bundle_value, decision
from backend.classification_dag.types import NodeResult


@dataclass(frozen=True)
class ClinicalLrCriteriaNode:
    id: str = "rule.clinical_lr"
    version: str = "1"
    requires: frozenset[str] = frozenset({"normalized_variant", "evidence_bundle"})
    provides: frozenset[str] = frozenset({"clinical_lr_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = inputs["normalized_variant"]
        value = bundle_value(inputs["evidence_bundle"], "clinical_lr") or {}
        decisions = ()
        warnings = ()
        if value.get("applies"):
            decisions = (
                decision(
                    value["code"],
                    dict(value),
                    gene=variant.gene,
                    family_id=self.id,
                    evidence_item_ids=("clinical_lr",),
                ),
            )
        elif value.get("application_status") == "review_required":
            warnings = (value.get("reason") or "Clinical LR review is required.",)
        return NodeResult.succeeded(
            {
                "clinical_lr_family": CriterionFamilyResult(
                    self.id,
                    criteria=decisions,
                    warnings=warnings,
                    metadata={"clinical_lr": value},
                )
            }
        )


@dataclass(frozen=True)
class ProteinPs1CriteriaNode:
    id: str = "rule.protein_ps1"
    version: str = "1"
    requires: frozenset[str] = frozenset({"normalized_variant", "evidence_bundle"})
    provides: frozenset[str] = frozenset({"protein_ps1_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = inputs["normalized_variant"]
        value = bundle_value(inputs["evidence_bundle"], "protein_ps1") or {}
        decisions = ()
        if value.get("applies"):
            payload = {
                "applies": True,
                "strength": value["strength"],
                "points": value["points"],
                "reason": value["reason"],
            }
            decisions = (
                decision(
                    "PS1",
                    payload,
                    gene=variant.gene,
                    family_id=self.id,
                    evidence_item_ids=("protein_ps1",),
                ),
            )
        return NodeResult.succeeded(
            {
                "protein_ps1_family": CriterionFamilyResult(
                    self.id,
                    criteria=decisions,
                    metadata={"ps1": value},
                )
            }
        )
