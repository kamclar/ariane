"""Data contracts shared by classification DAG nodes and the executor."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Awaitable, Mapping, Protocol


class NodeStatus(str, Enum):
    """Explicit outcome of one DAG node.

    Only ``succeeded`` publishes outputs. Other states remain distinct in the
    trace so unavailable evidence can never be mistaken for evidence that was
    evaluated and did not meet a criterion.
    """

    SUCCEEDED = "succeeded"
    NOT_APPLICABLE = "not_applicable"
    UNAVAILABLE = "unavailable"
    AMBIGUOUS = "ambiguous"
    REVIEW_REQUIRED = "review_required"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass(frozen=True)
class NodeResult:
    status: NodeStatus
    outputs: Mapping[str, Any] = field(default_factory=dict)
    reason: str = ""
    provenance: Mapping[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    @classmethod
    def succeeded(
        cls,
        outputs: Mapping[str, Any],
        *,
        provenance: Mapping[str, Any] | None = None,
        warnings: tuple[str, ...] = (),
    ) -> "NodeResult":
        return cls(
            status=NodeStatus.SUCCEEDED,
            outputs=dict(outputs),
            provenance=dict(provenance or {}),
            warnings=warnings,
        )

    @classmethod
    def unavailable(cls, reason: str, **provenance: Any) -> "NodeResult":
        return cls(
            status=NodeStatus.UNAVAILABLE,
            reason=reason,
            provenance=provenance,
        )


@dataclass(frozen=True)
class DagExecutionContext:
    run_id: str
    variant_key: str
    policy_id: str
    metadata: Mapping[str, Any] = field(default_factory=dict)


class DagNode(Protocol):
    id: str
    version: str
    requires: frozenset[str]
    provides: frozenset[str]

    def evaluate(
        self,
        context: DagExecutionContext,
        inputs: Mapping[str, Any],
    ) -> NodeResult | Awaitable[NodeResult]:
        """Evaluate one node without mutating upstream values."""


@dataclass(frozen=True)
class DagTraceEntry:
    node_id: str
    node_version: str
    status: NodeStatus
    requires: tuple[str, ...]
    provides: tuple[str, ...]
    started_at: str
    finished_at: str
    duration_ms: float
    reason: str = ""
    provenance: Mapping[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "node_version": self.node_version,
            "status": self.status.value,
            "requires": list(self.requires),
            "provides": list(self.provides),
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "duration_ms": self.duration_ms,
            "reason": self.reason,
            "provenance": dict(self.provenance),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class DagRun:
    graph_id: str
    graph_version: str
    values: Mapping[str, Any]
    trace: tuple[DagTraceEntry, ...]
