"""Small deterministic DAG executor with strict node contracts."""

from __future__ import annotations

import asyncio
from dataclasses import dataclass
from datetime import datetime, timezone
import inspect
from time import perf_counter
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from backend.classification_dag.types import (
    DagExecutionContext,
    DagNode,
    DagRun,
    DagTraceEntry,
    NodeResult,
    NodeStatus,
)


class DagDefinitionError(ValueError):
    """The graph is incomplete, ambiguous, or cyclic."""


class DagNodeExecutionError(RuntimeError):
    """A node violated its contract or raised an exception."""

    def __init__(
        self,
        node_id: str,
        message: str,
        *,
        trace: tuple[DagTraceEntry, ...] = (),
    ) -> None:
        super().__init__(f"DAG node {node_id!r} failed: {message}")
        self.node_id = node_id
        self.trace = trace


@dataclass(frozen=True)
class DagDefinition:
    id: str
    version: str
    seed_keys: frozenset[str]
    nodes: tuple[DagNode, ...]

    def __init__(
        self,
        *,
        id: str,
        version: str,
        seed_keys: Iterable[str],
        nodes: Iterable[DagNode],
    ) -> None:
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "seed_keys", frozenset(seed_keys))
        object.__setattr__(self, "nodes", tuple(nodes))
        self._validate()

    def _validate(self) -> None:
        if not self.id or not self.version:
            raise DagDefinitionError("Graph id and version are required")

        node_ids: set[str] = set()
        provider_by_key: dict[str, str] = {}
        for node in self.nodes:
            if not node.id or not node.version:
                raise DagDefinitionError("Every node requires an id and version")
            if node.id in node_ids:
                raise DagDefinitionError(f"Duplicate node id: {node.id}")
            node_ids.add(node.id)
            for key in node.provides:
                if key in self.seed_keys:
                    raise DagDefinitionError(
                        f"Node {node.id} overwrites seed value {key!r}"
                    )
                previous = provider_by_key.get(key)
                if previous:
                    raise DagDefinitionError(
                        f"Output {key!r} has multiple providers: {previous}, {node.id}"
                    )
                provider_by_key[key] = node.id

        for node in self.nodes:
            missing = {
                key
                for key in node.requires
                if key not in self.seed_keys and key not in provider_by_key
            }
            if missing:
                raise DagDefinitionError(
                    f"Node {node.id} has requirements without providers: "
                    f"{', '.join(sorted(missing))}"
                )

        # Kahn's algorithm proves acyclicity. Dependencies through seed values
        # do not create graph edges.
        dependencies: dict[str, set[str]] = {node.id: set() for node in self.nodes}
        dependants: dict[str, set[str]] = {node.id: set() for node in self.nodes}
        for node in self.nodes:
            for key in node.requires:
                provider = provider_by_key.get(key)
                if provider:
                    dependencies[node.id].add(provider)
                    dependants[provider].add(node.id)

        ready = [node.id for node in self.nodes if not dependencies[node.id]]
        visited: list[str] = []
        while ready:
            current = ready.pop(0)
            visited.append(current)
            for dependant in tuple(dependants[current]):
                dependencies[dependant].discard(current)
                if not dependencies[dependant] and dependant not in ready:
                    ready.append(dependant)
        if len(visited) != len(self.nodes):
            cyclic = sorted(node_ids.difference(visited))
            raise DagDefinitionError(
                f"Graph contains a cycle involving: {', '.join(cyclic)}"
            )

    def ordered_nodes(self) -> tuple[DagNode, ...]:
        provider_by_key = {
            key: node.id for node in self.nodes for key in node.provides
        }
        remaining = list(self.nodes)
        available = set(self.seed_keys)
        ordered: list[DagNode] = []
        while remaining:
            progressed = False
            for node in tuple(remaining):
                if all(
                    key in available or provider_by_key.get(key) is None
                    for key in node.requires
                ):
                    ordered.append(node)
                    available.update(node.provides)
                    remaining.remove(node)
                    progressed = True
            if not progressed:  # Defensive; _validate() already detects cycles.
                raise DagDefinitionError("Graph could not be topologically ordered")
        return tuple(ordered)


class DagExecutor:
    def __init__(self, definition: DagDefinition) -> None:
        self.definition = definition
        self._ordered_nodes = definition.ordered_nodes()

    def run(
        self,
        seed: Mapping[str, Any],
        *,
        context: DagExecutionContext,
    ) -> DagRun:
        missing_seed = self.definition.seed_keys.difference(seed)
        if missing_seed:
            raise DagDefinitionError(
                "Missing graph seed values: " + ", ".join(sorted(missing_seed))
            )

        values: dict[str, Any] = dict(seed)
        trace: list[DagTraceEntry] = []
        for node in self._ordered_nodes:
            missing_inputs = sorted(node.requires.difference(values))
            if missing_inputs:
                now = _utc_now()
                trace.append(
                    DagTraceEntry(
                        node_id=node.id,
                        node_version=node.version,
                        status=NodeStatus.SKIPPED,
                        requires=tuple(sorted(node.requires)),
                        provides=tuple(sorted(node.provides)),
                        started_at=now,
                        finished_at=now,
                        duration_ms=0.0,
                        reason=(
                            "Required upstream values were not produced: "
                            + ", ".join(missing_inputs)
                        ),
                    )
                )
                continue

            started_at = _utc_now()
            started = perf_counter()
            node_inputs = MappingProxyType(
                {key: values[key] for key in node.requires}
            )
            try:
                outcome = node.evaluate(context, node_inputs)
            except Exception as exc:
                finished_at = _utc_now()
                entry = DagTraceEntry(
                    node_id=node.id,
                    node_version=node.version,
                    status=NodeStatus.FAILED,
                    requires=tuple(sorted(node.requires)),
                    provides=tuple(sorted(node.provides)),
                    started_at=started_at,
                    finished_at=finished_at,
                    duration_ms=round((perf_counter() - started) * 1000, 3),
                    reason=f"{type(exc).__name__}: {exc}",
                )
                trace.append(entry)
                raise DagNodeExecutionError(
                    node.id,
                    entry.reason,
                    trace=tuple(trace),
                ) from exc

            if not isinstance(outcome, NodeResult):
                raise DagNodeExecutionError(
                    node.id,
                    f"returned {type(outcome).__name__}, expected NodeResult",
                    trace=tuple(trace),
                )
            if outcome.status == NodeStatus.FAILED:
                raise DagNodeExecutionError(
                    node.id,
                    outcome.reason or "node returned failed status",
                    trace=tuple(trace),
                )

            output_keys = frozenset(outcome.outputs)
            if outcome.status == NodeStatus.SUCCEEDED:
                if output_keys != node.provides:
                    raise DagNodeExecutionError(
                        node.id,
                        "published outputs do not match its declaration: "
                        f"expected {sorted(node.provides)}, got {sorted(output_keys)}",
                        trace=tuple(trace),
                    )
                values.update(outcome.outputs)
            elif output_keys:
                raise DagNodeExecutionError(
                    node.id,
                    f"status {outcome.status.value!r} must not publish outputs",
                    trace=tuple(trace),
                )

            trace.append(
                DagTraceEntry(
                    node_id=node.id,
                    node_version=node.version,
                    status=outcome.status,
                    requires=tuple(sorted(node.requires)),
                    provides=tuple(sorted(node.provides)),
                    started_at=started_at,
                    finished_at=_utc_now(),
                    duration_ms=round((perf_counter() - started) * 1000, 3),
                    reason=outcome.reason,
                    provenance=dict(outcome.provenance),
                    warnings=outcome.warnings,
                )
            )

        return DagRun(
            graph_id=self.definition.id,
            graph_version=self.definition.version,
            values=MappingProxyType(values),
            trace=tuple(trace),
        )

    async def run_async(
        self,
        seed: Mapping[str, Any],
        *,
        context: DagExecutionContext,
    ) -> DagRun:
        """Execute dependency-ready nodes concurrently.

        Synchronous nodes remain valid and are evaluated in the event-loop
        thread. Provider nodes may return an awaitable and perform blocking
        lookups through ``asyncio.to_thread``. Values from one readiness layer
        are committed only after every node in that layer has completed
        successfully, so a provider failure cannot publish a partial layer.
        """
        missing_seed = self.definition.seed_keys.difference(seed)
        if missing_seed:
            raise DagDefinitionError(
                "Missing graph seed values: " + ", ".join(sorted(missing_seed))
            )

        values: dict[str, Any] = dict(seed)
        trace: list[DagTraceEntry] = []
        remaining = list(self._ordered_nodes)
        while remaining:
            available = set(values)
            ready = [node for node in remaining if node.requires.issubset(available)]
            if not ready:
                now = _utc_now()
                for node in remaining:
                    missing_inputs = sorted(node.requires.difference(available))
                    trace.append(
                        DagTraceEntry(
                            node_id=node.id,
                            node_version=node.version,
                            status=NodeStatus.SKIPPED,
                            requires=tuple(sorted(node.requires)),
                            provides=tuple(sorted(node.provides)),
                            started_at=now,
                            finished_at=now,
                            duration_ms=0.0,
                            reason=(
                                "Required upstream values were not produced: "
                                + ", ".join(missing_inputs)
                            ),
                        )
                    )
                break

            evaluations = await asyncio.gather(
                *(self._evaluate_async_node(node, context, values) for node in ready),
                return_exceptions=True,
            )
            pending_outputs: dict[str, Any] = {}
            first_error: DagNodeExecutionError | None = None
            for node, evaluation in zip(ready, evaluations):
                if isinstance(evaluation, BaseException):
                    now = _utc_now()
                    entry = DagTraceEntry(
                        node_id=node.id,
                        node_version=node.version,
                        status=NodeStatus.FAILED,
                        requires=tuple(sorted(node.requires)),
                        provides=tuple(sorted(node.provides)),
                        started_at=now,
                        finished_at=now,
                        duration_ms=0.0,
                        reason=f"{type(evaluation).__name__}: {evaluation}",
                    )
                    trace.append(entry)
                    if first_error is None:
                        first_error = DagNodeExecutionError(
                            node.id, entry.reason, trace=tuple(trace)
                        )
                        first_error.__cause__ = evaluation
                    continue

                outcome, entry = evaluation
                trace.append(entry)
                try:
                    self._validate_outcome(node, outcome)
                except Exception as exc:
                    failed = DagTraceEntry(
                        node_id=node.id,
                        node_version=node.version,
                        status=NodeStatus.FAILED,
                        requires=entry.requires,
                        provides=entry.provides,
                        started_at=entry.started_at,
                        finished_at=entry.finished_at,
                        duration_ms=entry.duration_ms,
                        reason=f"{type(exc).__name__}: {exc}",
                    )
                    trace[-1] = failed
                    if first_error is None:
                        first_error = DagNodeExecutionError(
                            node.id, failed.reason, trace=tuple(trace)
                        )
                        first_error.__cause__ = exc
                    continue
                if outcome.status == NodeStatus.SUCCEEDED:
                    pending_outputs.update(outcome.outputs)

            if first_error is not None:
                first_error.trace = tuple(trace)
                raise first_error from first_error.__cause__

            values.update(pending_outputs)
            for node in ready:
                remaining.remove(node)

        return DagRun(
            graph_id=self.definition.id,
            graph_version=self.definition.version,
            values=MappingProxyType(values),
            trace=tuple(trace),
        )

    @staticmethod
    async def _evaluate_async_node(node, context, values):
        started_at = _utc_now()
        started = perf_counter()
        node_inputs = MappingProxyType({key: values[key] for key in node.requires})
        outcome = node.evaluate(context, node_inputs)
        if inspect.isawaitable(outcome):
            outcome = await outcome
        finished_at = _utc_now()
        entry = DagTraceEntry(
            node_id=node.id,
            node_version=node.version,
            status=(outcome.status if isinstance(outcome, NodeResult) else NodeStatus.FAILED),
            requires=tuple(sorted(node.requires)),
            provides=tuple(sorted(node.provides)),
            started_at=started_at,
            finished_at=finished_at,
            duration_ms=round((perf_counter() - started) * 1000, 3),
            reason=(outcome.reason if isinstance(outcome, NodeResult) else ""),
            provenance=(dict(outcome.provenance) if isinstance(outcome, NodeResult) else {}),
            warnings=(outcome.warnings if isinstance(outcome, NodeResult) else ()),
        )
        return outcome, entry

    @staticmethod
    def _validate_outcome(node, outcome: Any) -> None:
        if not isinstance(outcome, NodeResult):
            raise TypeError(
                f"returned {type(outcome).__name__}, expected NodeResult"
            )
        if outcome.status == NodeStatus.FAILED:
            raise RuntimeError(outcome.reason or "node returned failed status")
        output_keys = frozenset(outcome.outputs)
        if outcome.status == NodeStatus.SUCCEEDED:
            if output_keys != node.provides:
                raise ValueError(
                    "published outputs do not match its declaration: "
                    f"expected {sorted(node.provides)}, got {sorted(output_keys)}"
                )
        elif output_keys:
            raise ValueError(
                f"status {outcome.status.value!r} must not publish outputs"
            )


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()
