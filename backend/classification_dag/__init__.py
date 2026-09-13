"""Typed, fail-closed execution primitives for ARIANE classification."""

from backend.classification_dag.engine import (
    DagDefinition,
    DagDefinitionError,
    DagExecutor,
    DagNodeExecutionError,
)
from backend.domain.classification import (
    CLASSIFICATION_ENGINE_ID,
    ClassificationInputs,
    CriterionDecision,
    CriterionDecisionStatus,
    CriterionFamilyResult,
    EvidenceBundle,
    EvidenceItem,
    EvidenceStatus,
    NormalizedVariant,
    VariantAssertion,
)
from backend.classification_dag.runtime import (
    ClassificationExecution,
    execute_classification,
    execute_classification_request,
)
from backend.classification_dag.providers import (
    ClassificationRequest,
    ProviderDependencies,
    Table9EvidenceNode,
)
from backend.classification_dag.types import (
    DagExecutionContext,
    DagRun,
    DagTraceEntry,
    NodeResult,
    NodeStatus,
)

__all__ = [
    "ClassificationExecution",
    "CLASSIFICATION_ENGINE_ID",
    "ClassificationInputs",
    "ClassificationRequest",
    "CriterionDecision",
    "CriterionDecisionStatus",
    "CriterionFamilyResult",
    "DagDefinition",
    "DagDefinitionError",
    "DagExecutionContext",
    "DagExecutor",
    "DagNodeExecutionError",
    "DagRun",
    "DagTraceEntry",
    "EvidenceBundle",
    "EvidenceItem",
    "EvidenceStatus",
    "NormalizedVariant",
    "NodeResult",
    "NodeStatus",
    "Table9EvidenceNode",
    "ProviderDependencies",
    "VariantAssertion",
    "execute_classification",
    "execute_classification_request",
]
