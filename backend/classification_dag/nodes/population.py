"""Population-frequency and exon-CNV DAG rule nodes."""

from __future__ import annotations

from dataclasses import dataclass

from backend.classification_dag.domain import (
    CriterionDecision,
    CriterionDecisionStatus,
    CriterionFamilyResult,
    EvidenceBundle,
)
from backend.classification_dag.nodes.support import bundle_value, decision
from backend.classification_dag.types import NodeResult


@dataclass(frozen=True)
class FrequencyCriteriaNode:
    id: str = "rule.population_frequency"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {"classification_inputs", "normalized_variant", "evidence_bundle"}
    )
    provides: frozenset[str] = frozenset({"frequency_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        from backend.modules.frequency import (
            evaluate_frequency_criteria,
            pm2_not_applicable_decision,
        )

        ci = inputs["classification_inputs"]
        evidence = inputs["evidence_bundle"]
        if not isinstance(evidence, EvidenceBundle):
            raise TypeError("evidence_bundle must be an EvidenceBundle")
        criteria: list[CriterionDecision] = []
        excluded: list[CriterionDecision] = []
        not_applicable: list[CriterionDecision] = []
        warnings: list[str] = []
        pm2_not_applicable = pm2_not_applicable_decision(
            ci.variant_type,
            gene=ci.gene,
            c_notation=ci.c_notation,
        )
        if pm2_not_applicable:
            not_applicable.append(
                decision(
                    "PM2",
                    pm2_not_applicable,
                    gene=ci.gene,
                    family_id=self.id,
                    status=CriterionDecisionStatus.NOT_APPLICABLE,
                )
            )
        gnomad_data = bundle_value(evidence, "gnomad")
        if gnomad_data:
            evaluated = evaluate_frequency_criteria(
                gnomad_data,
                ci.variant_type,
                gene=ci.gene,
                c_notation=ci.c_notation,
            )
            for code, value in evaluated.items():
                if code.startswith("_"):
                    continue
                if value.get("applies"):
                    criteria.append(
                        decision(
                            code,
                            value,
                            gene=ci.gene,
                            family_id=self.id,
                            evidence_item_ids=("gnomad",),
                        )
                    )
                elif code == "PM2":
                    if (
                        "not applicable" in str(value.get("reason") or "").lower()
                        and not pm2_not_applicable
                    ):
                        not_applicable.append(
                            decision(
                                code,
                                value,
                                gene=ci.gene,
                                family_id=self.id,
                                status=CriterionDecisionStatus.NOT_APPLICABLE,
                                evidence_item_ids=("gnomad",),
                            )
                        )
                    else:
                        warnings.append(value["reason"])
            for code, value in evaluated.get("_excluded_criteria", {}).items():
                excluded.append(
                    decision(
                        code,
                        value,
                        gene=ci.gene,
                        family_id=self.id,
                        status=CriterionDecisionStatus.EXCLUDED,
                        evidence_item_ids=("gnomad",),
                    )
                )
            info = evaluated.get("_gnomad_info")
            if info:
                warnings.append(info["reason"])
        family = CriterionFamilyResult(
            family_id=self.id,
            criteria=tuple(criteria),
            excluded_criteria=tuple(excluded),
            not_applicable_criteria=tuple(not_applicable),
            warnings=tuple(warnings),
        )
        return NodeResult.succeeded(
            {"frequency_family": family},
            provenance={"criteria": [item.code for item in criteria]},
        )


@dataclass(frozen=True)
class ExonCnvCriteriaNode:
    id: str = "rule.exon_cnv.population"
    version: str = "1"
    requires: frozenset[str] = frozenset({"normalized_variant", "evidence_bundle"})
    provides: frozenset[str] = frozenset({"exon_cnv_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = inputs["normalized_variant"]
        value = bundle_value(inputs["evidence_bundle"], "exon_cnv") or {}
        decisions = (
            tuple(
                decision(
                    item["code"],
                    {
                        "applies": True,
                        "strength": item["strength"],
                        "points": item["points"],
                        "reason": item["reason"],
                        "source": item["source"],
                    },
                    gene=variant.gene,
                    family_id=self.id,
                    evidence_item_ids=("exon_cnv",),
                )
                for item in value.get("criteria", [])
            )
            if value.get("found")
            else ()
        )
        return NodeResult.succeeded(
            {
                "exon_cnv_family": CriterionFamilyResult(
                    self.id,
                    criteria=decisions,
                )
            }
        )
