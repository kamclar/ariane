"""Context-building nodes shared by several evidence rule families."""

from __future__ import annotations

from dataclasses import dataclass

from backend.classification_dag.domain import EvidenceItem
from backend.classification_dag.nodes.support import SpliceContext, bundle_value
from backend.classification_dag.types import NodeResult
from backend.modules.spliceai_policy import compare_table9_spliceai


@dataclass(frozen=True)
class SpliceContextNode:
    id: str = "context.spliceai.provenance"
    version: str = "2"
    requires: frozenset[str] = frozenset(
        {"classification_inputs", "evidence_bundle", "table9_evidence"}
    )
    provides: frozenset[str] = frozenset({"splice_context"})

    def evaluate(self, context, inputs) -> NodeResult:
        classification_inputs = inputs["classification_inputs"]
        bundle = inputs["evidence_bundle"]
        table9_evidence = inputs["table9_evidence"]
        if not isinstance(table9_evidence, EvidenceItem):
            raise TypeError("table9_evidence must be an EvidenceItem")
        table9 = table9_evidence.value or {}
        if (
            classification_inputs.table9_result is not None
            and dict(classification_inputs.table9_result) != table9
        ):
            raise ValueError(
                "Pre-DAG Table 9 evidence conflicts with the authoritative "
                "DAG provider result"
            )
        service_score = bundle_value(bundle, "spliceai")
        table9_score, warnings = compare_table9_spliceai(
            classification_inputs.gene, service_score, table9
        )
        value = SpliceContext(
            service_score,
            service_score,
            table9_score,
            tuple(warnings),
        )
        return NodeResult.succeeded(
            {"splice_context": value},
            provenance={
                "decision_score_source": "configured_spliceai",
                "table9_score_available_for_audit": table9_score is not None,
            },
        )
