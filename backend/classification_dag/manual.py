"""Fail-closed DAG wrapper for review-adjusted classification."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence
import uuid

from backend.classification_dag.engine import DagDefinition, DagExecutor
from backend.classification_dag.types import DagExecutionContext, DagTraceEntry, NodeResult
from backend.gene_policy import resolve_policy_identity
from backend.modules.manual_evidence import evaluate_manual_evidence


@dataclass(frozen=True)
class ManualEvidenceInputs:
    base_criteria: tuple[Mapping[str, Any], ...]
    manual_criteria: tuple[Mapping[str, Any], ...]
    variant_context: Mapping[str, Any] | None

    @classmethod
    def create(
        cls,
        base_criteria: Sequence[Mapping[str, Any]],
        manual_criteria: Sequence[Mapping[str, Any]],
        variant_context: Mapping[str, Any] | None = None,
    ) -> "ManualEvidenceInputs":
        return cls(
            tuple(dict(item) for item in base_criteria),
            tuple(dict(item) for item in manual_criteria),
            dict(variant_context) if variant_context is not None else None,
        )


@dataclass(frozen=True)
class ManualEvidenceExecution:
    result: Mapping[str, Any]
    graph_id: str
    graph_version: str
    trace: tuple[DagTraceEntry, ...]


@dataclass(frozen=True)
class ManualInputContractNode:
    id: str = "contract.manual_evidence_inputs"
    version: str = "2"
    requires: frozenset[str] = frozenset({"manual_inputs"})
    provides: frozenset[str] = frozenset({"validated_manual_inputs"})

    def evaluate(self, context, inputs) -> NodeResult:
        value = inputs["manual_inputs"]
        if not isinstance(value, ManualEvidenceInputs):
            raise TypeError("manual_inputs must be ManualEvidenceInputs")
        for item in value.base_criteria:
            if not isinstance(item.get("name"), str):
                raise ValueError("Every base criterion requires a name")
        for item in value.manual_criteria:
            if not isinstance(item.get("code"), str):
                raise ValueError("Every manual criterion requires a code")
        return NodeResult.succeeded({"validated_manual_inputs": value})


@dataclass(frozen=True)
class ManualEvidenceEvaluationNode:
    id: str = "rule.manual_evidence"
    version: str = "2"
    requires: frozenset[str] = frozenset({"validated_manual_inputs"})
    provides: frozenset[str] = frozenset({"manual_classification_result"})

    def evaluate(self, context, inputs) -> NodeResult:
        source = inputs["validated_manual_inputs"]
        result = evaluate_manual_evidence(
            [dict(item) for item in source.base_criteria],
            [dict(item) for item in source.manual_criteria],
            source.variant_context,
        )
        applied = [
            item["code"]
            for item in result["manual_criteria"]
            if item.get("applies")
        ]
        return NodeResult.succeeded(
            {"manual_classification_result": result},
            provenance={
                "applied": applied,
                "method": "ENIGMA-Table-3-or-mixed-points",
            },
        )


def build_manual_evidence_graph() -> DagDefinition:
    return DagDefinition(
        id="ariane.vcep.manual-evidence",
        version="3.0.0-single-rule-owner",
        seed_keys={"manual_inputs"},
        nodes=(ManualInputContractNode(), ManualEvidenceEvaluationNode()),
    )


def execute_manual_evidence(
    base_criteria: Sequence[Mapping[str, Any]],
    manual_criteria: Sequence[Mapping[str, Any]],
    variant_context: Mapping[str, Any] | None = None,
) -> ManualEvidenceExecution:
    gene = str((variant_context or {}).get("gene") or "").strip().upper() or None
    selected_policy_id, selected_policy_version = resolve_policy_identity(gene)
    graph = build_manual_evidence_graph()
    run = DagExecutor(graph).run(
        {
            "manual_inputs": ManualEvidenceInputs.create(
                base_criteria, manual_criteria, variant_context
            )
        },
        context=DagExecutionContext(
            run_id=uuid.uuid4().hex,
            variant_key="manual-review",
            policy_id=selected_policy_id,
            metadata={
                "workflow": "review-adjusted-classification",
                "policy_version": selected_policy_version,
                "gene": gene,
            },
        ),
    )
    result = run.values.get("manual_classification_result")
    if not isinstance(result, dict):
        raise RuntimeError("Manual evidence DAG did not produce a result")
    return ManualEvidenceExecution(result, run.graph_id, run.graph_version, run.trace)
