"""PVS1, RNA-PVS1 and PTC-PM5 DAG rule nodes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from backend.domain.classification import (
    ClassificationInputs,
    CriterionDecision,
    CriterionDecisionStatus,
    CriterionFamilyResult,
    EvidenceItem,
    NormalizedVariant,
)
from backend.domain.pvs1_rna import Pvs1RnaSelection
from backend.classification_dag.nodes.support import decision
from backend.classification_dag.types import NodeResult
from backend.criteria.pvs1 import evaluate_pvs1
from backend.criteria.pvs1_rna import (
    lookup_st2_pvs1_rna_evidence,
    resolve_pvs1_rna,
)


def _evidence_value(
    inputs: Mapping[str, Any],
    key: str,
) -> dict[str, Any] | None:
    item = inputs[key]
    if not isinstance(item, EvidenceItem):
        raise TypeError(f"{key} must be an EvidenceItem")
    return item.value if isinstance(item.value, dict) else None


@dataclass(frozen=True)
class Pvs1RnaEvidenceResolutionNode:
    """Resolve provider evidence once before dependent rule nodes run."""

    id: str = "context.pvs1_rna"
    version: str = "1"
    requires: frozenset[str] = frozenset({
        "normalized_variant",
        "erepo_pvs1_rna_evidence",
        "table4_pvs1_rna_evidence",
        "st2_pvs1_rna_evidence",
    })
    provides: frozenset[str] = frozenset({"pvs1_rna_selection"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = inputs["normalized_variant"]
        if not isinstance(variant, NormalizedVariant):
            raise TypeError("normalized_variant must be a NormalizedVariant")
        selection = resolve_pvs1_rna(
            variant.gene,
            variant.c_notation,
            erepo_registry_result=_evidence_value(
                inputs,
                "erepo_pvs1_rna_evidence",
            ),
            table4_evidence_result=_evidence_value(
                inputs,
                "table4_pvs1_rna_evidence",
            ),
            st2_evidence_result=_evidence_value(
                inputs,
                "st2_pvs1_rna_evidence",
            ),
        )
        return NodeResult.succeeded(
            {"pvs1_rna_selection": selection},
            provenance={
                "source": selection.source,
                "replaces_prediction": selection.replaces_prediction,
            },
        )


@dataclass(frozen=True)
class Pvs1RnaInputResolutionNode:
    """Resolve RNA evidence supplied directly to the synchronous DAG."""

    id: str = "context.pvs1_rna"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_inputs"})
    provides: frozenset[str] = frozenset({"pvs1_rna_selection"})

    def evaluate(self, context, inputs) -> NodeResult:
        ci = inputs["classification_inputs"]
        if not isinstance(ci, ClassificationInputs):
            raise TypeError("classification_inputs must be ClassificationInputs")
        st2_result = (
            ci.st2_pvs1_rna_result
            if ci.st2_pvs1_rna_result is not None
            else lookup_st2_pvs1_rna_evidence(ci.gene, ci.c_notation)
        )
        selection = resolve_pvs1_rna(
            ci.gene,
            ci.c_notation,
            erepo_registry_result=ci.erepo_pvs1_rna_result,
            table4_evidence_result=ci.table4_pvs1_rna_result,
            st2_evidence_result=st2_result,
        )
        return NodeResult.succeeded(
            {"pvs1_rna_selection": selection},
            provenance={
                "source": selection.source,
                "replaces_prediction": selection.replaces_prediction,
            },
        )


@dataclass(frozen=True)
class Pvs1CriteriaNode:
    id: str = "rule.pvs1_pm5"
    version: str = "1"
    requires: frozenset[str] = frozenset({
        "classification_inputs",
        "splice_context",
        "pvs1_rna_selection",
    })
    provides: frozenset[str] = frozenset({"pvs1_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        ci = inputs["classification_inputs"]
        splice = inputs["splice_context"]
        rna_selection = inputs["pvs1_rna_selection"]
        if not isinstance(rna_selection, Pvs1RnaSelection):
            raise TypeError("pvs1_rna_selection must be a Pvs1RnaSelection")
        pvs1 = evaluate_pvs1(
            ci.gene,
            ci.variant_type,
            ci.p_notation,
            c_notation=ci.c_notation,
            spliceai_score=splice.effective_score,
            dup_type=ci.dup_type,
        )
        pvs1_rna = rna_selection.decision
        decisions: list[CriterionDecision] = []
        excluded: list[CriterionDecision] = []
        not_applicable: list[CriterionDecision] = []
        warnings: list[str] = []
        if pvs1_rna.get("applies"):
            evidence_item_ids = (
                (rna_selection.evidence_item_id,)
                if rna_selection.evidence_item_id
                else ()
            )
            decisions.append(
                decision(
                    "PVS1_RNA",
                    pvs1_rna,
                    gene=ci.gene,
                    family_id=self.id,
                    evidence_item_ids=evidence_item_ids,
                )
            )
        elif pvs1["applies"]:
            decisions.append(
                decision(
                    "PVS1",
                    pvs1,
                    gene=ci.gene,
                    family_id=self.id,
                    evidence_item_ids=("spliceai",),
                )
            )
        elif pvs1.get("requires_rna") or ci.variant_type.lower() in {
            "nonsense",
            "frameshift",
            "ptc_delins",
            "splice_site",
            "initiation_codon",
            "exon_deletion",
            "exon_duplication",
        }:
            warnings.append(pvs1["reason"])
            if "N/A" in str(pvs1.get("pvs1_code") or ""):
                not_applicable.append(
                    decision(
                        "PVS1",
                        {
                            "applies": False,
                            "strength": "N/A",
                            "points": 0,
                            "reason": pvs1["reason"],
                            "source": pvs1.get("source", ""),
                        },
                        gene=ci.gene,
                        family_id=self.id,
                        status=CriterionDecisionStatus.NOT_APPLICABLE,
                    )
                )
        if (
            not pvs1.get("applies")
            and pvs1_rna.get("source_record")
            and not pvs1_rna.get("applies")
            and pvs1_rna.get("review_required")
        ):
            warnings.append(pvs1_rna["reason"])
        if pvs1.get("pm5_code") and pvs1.get("pm5_strength"):
            decisions.append(
                decision(
                    "PM5_PTC",
                    {
                        "applies": True,
                        "strength": pvs1["pm5_strength"],
                        "points": pvs1["pm5_points"],
                        "reason": (
                            f"Table 4: {pvs1['pm5_code']} for PTC in "
                            f"{pvs1.get('pm5_exon') or 'unknown exon'}"
                        ),
                    },
                    gene=ci.gene,
                    family_id=self.id,
                )
            )
        family = CriterionFamilyResult(
            self.id,
            criteria=tuple(decisions),
            excluded_criteria=tuple(excluded),
            not_applicable_criteria=tuple(not_applicable),
            warnings=tuple(warnings),
            has_functional_evidence=bool(pvs1_rna.get("applies")),
            metadata={"pvs1": pvs1, "pvs1_rna": pvs1_rna},
        )
        return NodeResult.succeeded(
            {"pvs1_family": family},
            provenance={"criteria": [item.code for item in decisions]},
        )
