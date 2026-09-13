"""Application workflows for single and batch variant classification.

The service coordinates classification, cache access, usage recording, batch
validation, and bounded concurrency. HTTP authentication, cookies, response
headers, quotas, and audit-log transport remain in the API controller.
"""

from __future__ import annotations

import asyncio
from collections.abc import Awaitable, Callable
from dataclasses import dataclass
import logging
import time
from typing import Any

from pydantic import ValidationError

from backend.classification_runtime import (
    ClassificationCacheRepository,
    ClassificationUsageRepository,
    classification_fingerprint,
)
from backend.contracts import (
    BatchItemResult,
    BatchRequest,
    BatchResponse,
    ClassificationResult,
    VariantRequest,
)
from backend.policy.gene import get_gene_policy
from backend.services.classification_presentation import ClassificationPresentationService
from backend.services.evidence_orchestration import (
    ClassificationCommand,
    EvidenceOrchestrationService,
)


LOGGER = logging.getLogger("ariane.variant_classification")

UncachedClassifier = Callable[
    [str, str, str, str], Awaitable[ClassificationResult]
]
CachedClassifier = Callable[
    [str, str, str, str, str],
    Awaitable[tuple[ClassificationResult, str, str]],
]
ErrorMapper = Callable[[Exception], tuple[str, str, bool]]


@dataclass(frozen=True)
class UsageContext:
    actor_id: str
    actor_type: str
    request_id: str


@dataclass(frozen=True)
class CompletedClassification:
    result: ClassificationResult
    cache_status: str
    classifier_fingerprint: str
    duration_ms: float
    policy: dict[str, str]


@dataclass(frozen=True)
class ServiceAuditEvent:
    event: str
    level: str
    fields: dict[str, Any]


@dataclass(frozen=True)
class PreparedBatch:
    valid_items: tuple[tuple[int, VariantRequest], ...]
    validation_errors: tuple[BatchItemResult, ...]
    audit_events: tuple[ServiceAuditEvent, ...]

    @property
    def classification_units(self) -> int:
        return len(self.valid_items)


@dataclass(frozen=True)
class CompletedBatch:
    response: BatchResponse
    duration_ms: float
    audit_events: tuple[ServiceAuditEvent, ...]


def _policy_metadata(gene: str) -> dict[str, str]:
    configured = get_gene_policy(gene)
    return {
        "gene": gene,
        "policy_id": configured["policy"]["runtime_policy_id"],
        "policy_version": configured["policy"]["version"],
        "reference_transcript": configured["gene_config"]["reference_transcript"],
    }


def _batch_variant_label(index: int, raw_item: object) -> str:
    if isinstance(raw_item, dict):
        gene = str(raw_item.get("gene", "")).strip().upper()
        notation = str(raw_item.get("c_notation", "")).strip()
        label = " ".join(value for value in (gene, notation) if value)
        if label:
            return label
    return f"Batch item {index + 1}"


def _batch_validation_message(exc: ValidationError) -> str:
    messages: list[str] = []
    for error in exc.errors(include_url=False):
        location = ".".join(str(value) for value in error.get("loc", ()))
        message = str(error.get("msg", "Invalid variant input"))
        if message.startswith("Value error, "):
            message = message[len("Value error, "):]
        rendered = f"{location}: {message}" if location else message
        if rendered not in messages:
            messages.append(rendered)
    return "; ".join(messages) or "Invalid variant input"


def _safe_batch_input(raw_item: object) -> dict[str, str]:
    if isinstance(raw_item, dict):
        return {
            "gene": str(raw_item.get("gene", ""))[:40],
            "c_notation": str(raw_item.get("c_notation", ""))[:200],
        }
    return {"input_type": type(raw_item).__name__}


def _apply_normalization_metadata(
    result: ClassificationResult,
    request: VariantRequest,
) -> None:
    result.reference_transcript = request.reference_transcript
    result.submitted_notation = request.submitted_notation
    result.normalization_source = request.normalization_source
    result.consequence_status = request.consequence_status
    result.normalization_provenance = request.normalization_provenance
    result.protein_consequence_explanation = request.protein_consequence_explanation


def _record_usage(
    repository: ClassificationUsageRepository | None,
    **values: Any,
) -> None:
    """Keep analytics auxiliary so storage failure cannot change a result."""
    if repository is None:
        return
    try:
        repository.record(**values)
    except Exception:
        LOGGER.exception("Classification usage event could not be stored")


async def execute_variant_classification(
    command: ClassificationCommand,
    *,
    orchestration: EvidenceOrchestrationService | None = None,
    presentation: ClassificationPresentationService | None = None,
) -> ClassificationResult:
    """Run evidence orchestration and project its result for API consumers."""
    orchestrator = orchestration or EvidenceOrchestrationService()
    presenter = presentation or ClassificationPresentationService()
    evidence = await orchestrator.orchestrate(command)
    return presenter.build(evidence)


class VariantClassificationService:
    """Execute the shared application workflow behind all API transports."""

    def __init__(self, *, batch_concurrency: int = 3):
        self.batch_concurrency = max(1, batch_concurrency)
        self._batch_semaphore = asyncio.Semaphore(self.batch_concurrency)

    async def classify_cached(
        self,
        gene: str,
        c_notation: str,
        p_notation: str,
        dup_type: str,
        reference_transcript: str,
        *,
        cache_repository: ClassificationCacheRepository | None,
        classify_uncached: UncachedClassifier,
    ) -> tuple[ClassificationResult, str, str]:
        """Return a fingerprint-bound result, using only complete cache entries."""
        fingerprint = classification_fingerprint(gene)
        cache_status = "disabled" if cache_repository is None else "miss"
        cached = None
        if cache_repository is not None:
            try:
                cached = cache_repository.get(
                    gene=gene,
                    transcript=reference_transcript,
                    c_notation=c_notation,
                    p_notation=p_notation,
                    dup_type=dup_type,
                    fingerprint=fingerprint,
                )
                if cached is not None:
                    cache_status = cached.status
            except Exception:
                LOGGER.exception(
                    "Classification cache read failed for %s:%s", gene, c_notation
                )
                cached = None
                cache_status = "read_error"
        if cached is not None and cached.result is not None:
            return cached.result, cached.status, fingerprint

        response = await classify_uncached(gene, c_notation, p_notation, dup_type)
        if cache_repository is not None:
            try:
                cache_repository.put(
                    response,
                    transcript=reference_transcript,
                    dup_type=dup_type,
                    fingerprint=fingerprint,
                )
            except Exception:
                LOGGER.exception(
                    "Classification cache write failed for %s:%s", gene, c_notation
                )
                cache_status = "write_error"
        return response, cache_status, fingerprint

    async def classify_single(
        self,
        request: VariantRequest,
        *,
        usage: UsageContext,
        usage_repository: ClassificationUsageRepository | None,
        classify_cached: CachedClassifier,
    ) -> CompletedClassification:
        """Classify one normalized request and record exactly one usage event."""
        started = time.monotonic()
        cache_status = "not_attempted"
        fingerprint = classification_fingerprint(request.gene)
        policy = _policy_metadata(request.gene)
        try:
            result, cache_status, fingerprint = await classify_cached(
                request.gene,
                request.c_notation,
                request.p_notation or "",
                request.dup_type,
                request.reference_transcript,
            )
            _apply_normalization_metadata(result, request)
        except Exception:
            _record_usage(
                usage_repository,
                actor_id=usage.actor_id,
                actor_type=usage.actor_type,
                request_id=usage.request_id,
                request_mode="single",
                gene=request.gene,
                transcript=request.reference_transcript,
                c_notation=request.c_notation,
                status="error",
                cache_status=cache_status,
                predicted_class=None,
                total_points=None,
                duration_ms=(time.monotonic() - started) * 1000,
                policy_id=policy["policy_id"],
                policy_version=policy["policy_version"],
                classifier_fingerprint=fingerprint,
            )
            raise
        duration_ms = (time.monotonic() - started) * 1000
        _record_usage(
            usage_repository,
            actor_id=usage.actor_id,
            actor_type=usage.actor_type,
            request_id=usage.request_id,
            request_mode="single",
            gene=request.gene,
            transcript=request.reference_transcript,
            c_notation=request.c_notation,
            status="completed",
            cache_status=cache_status,
            predicted_class=result.predicted_class,
            total_points=result.total_points,
            duration_ms=duration_ms,
            policy_id=policy["policy_id"],
            policy_version=policy["policy_version"],
            classifier_fingerprint=fingerprint,
        )
        return CompletedClassification(
            result=result,
            cache_status=cache_status,
            classifier_fingerprint=fingerprint,
            duration_ms=duration_ms,
            policy=policy,
        )

    @staticmethod
    def prepare_batch(request: BatchRequest) -> PreparedBatch:
        """Validate items independently so one bad variant does not reject the batch."""
        valid_items: list[tuple[int, VariantRequest]] = []
        validation_errors: list[BatchItemResult] = []
        audit_events: list[ServiceAuditEvent] = []
        for index, raw_item in enumerate(request.variants):
            try:
                item = VariantRequest.model_validate(raw_item)
            except ValidationError as exc:
                error = _batch_validation_message(exc)
                validation_errors.append(BatchItemResult(
                    index=index,
                    status="error",
                    variant=_batch_variant_label(index, raw_item),
                    error=error,
                    error_code="invalid_variant",
                    error_retryable=False,
                ))
                audit_events.append(ServiceAuditEvent(
                    event="batch_item_validation_error",
                    level="warning",
                    fields={
                        "item_index": index,
                        "input": _safe_batch_input(raw_item),
                        "error": error,
                    },
                ))
            else:
                valid_items.append((index, item))
        return PreparedBatch(
            valid_items=tuple(valid_items),
            validation_errors=tuple(validation_errors),
            audit_events=tuple(audit_events),
        )

    async def classify_batch(
        self,
        prepared: PreparedBatch,
        *,
        usage: UsageContext,
        usage_repository: ClassificationUsageRepository | None,
        classify_cached: CachedClassifier,
        map_error: ErrorMapper,
    ) -> CompletedBatch:
        """Classify valid batch items concurrently while preserving input order."""
        batch_started = time.monotonic()

        async def classify_item(
            index: int,
            item: VariantRequest,
        ) -> tuple[BatchItemResult, ServiceAuditEvent]:
            started = time.monotonic()
            cache_status = "not_attempted"
            fingerprint = classification_fingerprint(item.gene)
            async with self._batch_semaphore:
                try:
                    result, cache_status, fingerprint = await classify_cached(
                        item.gene,
                        item.c_notation,
                        item.p_notation or "",
                        item.dup_type,
                        item.reference_transcript,
                    )
                    _apply_normalization_metadata(result, item)
                    output = BatchItemResult(
                        index=index,
                        status="ok",
                        variant=f"{item.gene} {item.c_notation}",
                        result=result,
                    )
                except Exception as exc:
                    error_code, error_message, retryable = map_error(exc)
                    output = BatchItemResult(
                        index=index,
                        status="error",
                        variant=f"{item.gene} {item.c_notation}",
                        error=error_message,
                        error_code=error_code,
                        error_retryable=retryable,
                    )

            policy = _policy_metadata(item.gene)
            duration_ms = (time.monotonic() - started) * 1000
            output.cache_status = cache_status
            output.classifier_fingerprint = fingerprint
            output.policy_id = policy["policy_id"]
            output.policy_version = policy["policy_version"]
            output.duration_ms = duration_ms
            _record_usage(
                usage_repository,
                actor_id=usage.actor_id,
                actor_type=usage.actor_type,
                request_id=usage.request_id,
                request_mode="batch",
                gene=item.gene,
                transcript=item.reference_transcript,
                c_notation=item.c_notation,
                status="completed" if output.status == "ok" else "error",
                cache_status=cache_status,
                predicted_class=output.result.predicted_class if output.result else None,
                total_points=output.result.total_points if output.result else None,
                duration_ms=duration_ms,
                policy_id=policy["policy_id"],
                policy_version=policy["policy_version"],
                classifier_fingerprint=fingerprint,
            )
            event = ServiceAuditEvent(
                event=(
                    "batch_item_completed"
                    if output.status == "ok"
                    else "batch_item_error"
                ),
                level="info" if output.status == "ok" else "warning",
                fields={
                    "item_index": output.index,
                    "input": item.model_dump(mode="json"),
                    "result": {
                        "predicted_class": output.result.predicted_class,
                        "predicted_label": output.result.predicted_label,
                        "total_points": output.result.total_points,
                        "classification_cache": cache_status,
                    } if output.result else None,
                    "error": output.error,
                },
            )
            return output, event

        executions = await asyncio.gather(*(
            classify_item(index, item) for index, item in prepared.valid_items
        ))
        items = sorted(
            [*prepared.validation_errors, *(output for output, _event in executions)],
            key=lambda value: value.index,
        )
        success_count = sum(item.status == "ok" for item in items)
        response = BatchResponse(
            total=len(items),
            success_count=success_count,
            error_count=len(items) - success_count,
            results=items,
        )
        return CompletedBatch(
            response=response,
            duration_ms=(time.monotonic() - batch_started) * 1000,
            audit_events=tuple(event for _output, event in executions),
        )
