"""Typed, fail-closed execution primitives for ARIANE classification."""

from backend.classification_dag.engine import (
    DagDefinition,
    DagDefinitionError,
    DagExecutor,
    DagNodeExecutionError,
)
from backend.classification_dag.domain import (
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
    ClassifierEngineMode,
    execute_classification,
    execute_classification_request,
    get_configured_engine_mode,
)
from backend.classification_dag.providers import (
    ClassificationRequest,
    ProviderDependencies,
    Table9EvidenceNode,
)
from backend.classification_dag.manual import (
    ManualEvidenceExecution,
    ManualEvidenceInputs,
    execute_manual_evidence,
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
    "ClassificationInputs",
    "ClassificationRequest",
    "ClassifierEngineMode",
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
    "ManualEvidenceExecution",
    "ManualEvidenceInputs",
    "VariantAssertion",
    "execute_classification",
    "execute_classification_request",
    "execute_manual_evidence",
    "get_configured_engine_mode",
]
