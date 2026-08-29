"""Figure 1A bioinformatic DAG rule nodes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from backend.classification_dag.domain import CriterionDecision, CriterionFamilyResult
from backend.classification_dag.nodes.support import bundle_value, decision
from backend.classification_dag.types import NodeResult
from backend.modules.bp1 import evaluate_bp1
from backend.modules.bp7 import evaluate_bp7
from backend.modules.evidence_interactions import pvs1_prediction_deduplication
from backend.modules.pp3_bp4 import evaluate_pp3_bp4
from backend.modules.utils import get_amino_acid_position, is_in_functional_domain


@dataclass(frozen=True)
class BioinformaticCriteriaNode:
    id: str = "rule.bioinformatic.figure1a"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {"classification_inputs", "evidence_bundle", "splice_context", "pvs1_family"}
    )
    provides: frozenset[str] = frozenset({"bioinformatic_family"})

    def evaluate(self, context, inputs) -> NodeResult:
        ci = inputs["classification_inputs"]
        bundle = inputs["evidence_bundle"]
        bayesdel_score = bundle_value(bundle, "bayesdel")
        splice = inputs["splice_context"]
        pvs1 = inputs["pvs1_family"].metadata["pvs1"]
        decisions: list[CriterionDecision] = []
        warnings: list[str] = []
        interactions: list[Mapping[str, Any]] = []
        pp3_bp4 = evaluate_pp3_bp4(
            ci.gene,
            ci.variant_type,
            ci.p_notation,
            bayesdel_score=bayesdel_score,
            spliceai_score=splice.effective_score,
            c_notation=ci.c_notation,
        )
        for code, value in pp3_bp4.items():
            if not value.get("applies"):
                continue
            if code == "PP3" and pvs1.get("applies"):
                warnings.append(
                    "PP3 not applied because PVS1 is met; ENIGMA does not stack "
                    "predictive PP3 with PVS1."
                )
                interactions.append(pvs1_prediction_deduplication())
            else:
                decisions.append(
                    decision(
                        code,
                        value,
                        gene=ci.gene,
                        family_id=self.id,
                        evidence_item_ids=("spliceai", "bayesdel"),
                    )
                )
        if ci.variant_type.lower() in {"synonymous", "silent", "intronic"}:
            aa_pos = get_amino_acid_position(ci.p_notation)
            in_domain = False
            if aa_pos:
                in_domain, _ = is_in_functional_domain(ci.gene, aa_pos)
            bp4_met = any(item.code == "BP4" for item in decisions)
            bp7 = evaluate_bp7(
                ci.variant_type,
                spliceai_score=splice.effective_score,
                in_domain=in_domain,
                bp4_met=bp4_met,
                c_notation=ci.c_notation,
                gene=ci.gene,
            )
            if bp7["applies"]:
                decisions.append(
                    decision(
                        "BP7",
                        bp7,
                        gene=ci.gene,
                        family_id=self.id,
                        evidence_item_ids=("spliceai",),
                    )
                )
        bp1 = evaluate_bp1(
            ci.gene,
            ci.variant_type,
            ci.p_notation,
            spliceai_score=splice.effective_score,
        )
        if bp1["applies"]:
            decisions.append(
                decision(
                    "BP1",
                    bp1,
                    gene=ci.gene,
                    family_id=self.id,
                    evidence_item_ids=("spliceai",),
                )
            )
        figure1a_types = {
            "missense",
            "inframe_deletion",
            "inframe_insertion",
            "inframe_delins",
            "delins",
            "synonymous",
            "silent",
            "intronic",
        }
        figure1a_unavailable = (
            splice.effective_score is None
            and ci.variant_type.lower() in figure1a_types
        )
        if figure1a_unavailable:
            warnings.append(
                f"Figure 1A bioinformatic result unavailable for {ci.gene} "
                f"{ci.c_notation}: SpliceAI is unavailable. Missing data was not "
                "treated as an ENIGMA prediction band; PP3, BP4, BP1 and BP7 "
                "were not applied."
            )
        bayesdel_types = {
            "missense",
            "inframe_deletion",
            "inframe_insertion",
            "inframe_delins",
            "delins",
        }
        if bayesdel_score is None and ci.variant_type.lower() in bayesdel_types:
            warnings.append(
                f"BayesDel_noAF not available for {ci.gene} {ci.c_notation}"
            )
        family = CriterionFamilyResult(
            self.id,
            criteria=tuple(decisions),
            warnings=tuple(warnings),
            evidence_interactions=tuple(interactions),
            metadata={
                "evaluation_status": (
                    "unavailable" if figure1a_unavailable else "evaluated"
                ),
                "spliceai_available": splice.effective_score is not None,
            },
        )
        return NodeResult.succeeded(
            {"bioinformatic_family": family},
            provenance={"criteria": [item.code for item in decisions]},
        )
