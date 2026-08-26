"""PVS1, RNA-PVS1 and PTC-PM5 DAG rule nodes."""

from __future__ import annotations

from dataclasses import dataclass

from backend.classification_dag.domain import CriterionDecision, CriterionDecisionStatus, CriterionFamilyResult
from backend.classification_dag.nodes.support import decision
from backend.classification_dag.types import NodeResult


@dataclass(frozen=True)
class Pvs1CriteriaNode:
    id: str = "rule.pvs1_pm5"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_inputs", "splice_context"})
    provides: frozenset[str] = frozenset({"pvs1_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        from backend.modules.pvs1 import evaluate_pvs1
        from backend.modules.pvs1_rna import evaluate_pvs1_rna

        ci = inputs["classification_inputs"]
        splice = inputs["splice_context"]
        pvs1 = evaluate_pvs1(
            ci.gene,
            ci.variant_type,
            ci.p_notation,
            c_notation=ci.c_notation,
            spliceai_score=splice.effective_score,
            dup_type=ci.dup_type,
        )
        pvs1_rna = evaluate_pvs1_rna(ci.gene, ci.c_notation)
        decisions: list[CriterionDecision] = []
        excluded: list[CriterionDecision] = []
        warnings: list[str] = []
        functional = False
        if pvs1["applies"]:
            decisions.append(
                decision(
                    "PVS1",
                    pvs1,
                    gene=ci.gene,
                    family_id=self.id,
                    evidence_item_ids=("spliceai",),
                )
            )
        elif pvs1_rna.get("applies"):
            decisions.append(
                decision(
                    "PVS1_RNA",
                    pvs1_rna,
                    gene=ci.gene,
                    family_id=self.id,
                )
            )
            functional = True
        elif pvs1.get("requires_rna") or ci.variant_type.lower() in {
            "nonsense",
            "frameshift",
            "splice_site",
            "initiation_codon",
            "exon_deletion",
            "exon_duplication",
        }:
            warnings.append(pvs1["reason"])
            if "N/A" in str(pvs1.get("pvs1_code") or ""):
                excluded.append(
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
                        status=CriterionDecisionStatus.EXCLUDED,
                    )
                )
        if (
            not pvs1.get("applies")
            and pvs1_rna.get("source_record")
            and not pvs1_rna.get("applies")
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
            warnings=tuple(warnings),
            has_functional_evidence=functional,
            metadata={"pvs1": pvs1, "pvs1_rna": pvs1_rna},
        )
        return NodeResult.succeeded(
            {"pvs1_family": family},
            provenance={"criteria": [item.code for item in decisions]},
        )
