"""HTTP transport for single and batch classification workflows."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timezone
import logging

from fastapi import APIRouter, Depends, HTTPException, Request, Response

from backend.api.auth import require_public_api_key
from backend.api.routing import PairedApiRoutes
from backend.classification_dag import CLASSIFICATION_ENGINE_ID
from backend.classification_runtime import (
    ClassificationCacheRepository,
    ClassificationUsageRepository,
    resolve_usage_identity,
)
from backend.contracts import BatchRequest, BatchResponse, ClassificationResult, VariantRequest
from backend.api.public import create_public_api_router
from backend.services.classification_presentation import ClassificationPresentationService
from backend.services import (
    ClassificationCommand,
    EvidenceExecutionError,
    EvidenceOrchestrationService,
    RequiredEvidenceUnavailableError,
    UsageContext,
    VariantClassificationService,
    VariantPreparationError,
    execute_variant_classification,
)


AuditWriter = Callable[..., None]
RepositoryGetter = Callable[[], ClassificationUsageRepository | None]
CacheGetter = Callable[[], ClassificationCacheRepository | None]
OrchestrationGetter = Callable[[], EvidenceOrchestrationService]


class ClassificationApi:
    """Adapt HTTP requests to the shared classification application service."""

    def __init__(
        self,
        *,
        service: VariantClassificationService,
        orchestration: OrchestrationGetter,
        cache_repository: CacheGetter,
        usage_repository: RepositoryGetter,
        quota_repository: RepositoryGetter,
        daily_classification_limit: Callable[[], int],
        audit: AuditWriter,
    ) -> None:
        self.service = service
        self._orchestration = orchestration
        self._cache_repository = cache_repository
        self._usage_repository = usage_repository
        self._quota_repository = quota_repository
        self._daily_classification_limit = daily_classification_limit
        self.audit = audit
        self.presentation = ClassificationPresentationService()

    @staticmethod
    def map_execution_error(exc: Exception) -> tuple[str, str, bool]:
        if isinstance(exc, HTTPException):
            message = str(exc.detail)
            headers = exc.headers or {}
            explicit_code = headers.get("X-ARIANE-Error-Code")
            retryable_header = headers.get("X-ARIANE-Retryable")
            explicit_retryable = (
                retryable_header.lower() == "true"
                if retryable_header is not None
                else None
            )
            if exc.status_code == 422:
                return explicit_code or "variant_not_classifiable", message, False
            if exc.status_code == 503:
                return (
                    explicit_code or "evidence_unavailable",
                    message,
                    True if explicit_retryable is None else explicit_retryable,
                )
            return (
                explicit_code or "request_failed",
                message,
                (
                    exc.status_code in {429, 502, 504}
                    if explicit_retryable is None
                    else explicit_retryable
                ),
            )
        return "internal_error", str(exc) or type(exc).__name__, False

    async def classify_uncached(
        self,
        gene: str,
        c_notation: str,
        p_notation: str = "",
        dup_type: str = "Unknown",
    ) -> ClassificationResult:
        """Map application failures to the established HTTP error contract."""
        try:
            return await execute_variant_classification(
                ClassificationCommand(
                    gene=gene,
                    c_notation=c_notation,
                    p_notation=p_notation,
                    dup_type=dup_type,
                ),
                orchestration=self._orchestration(),
            )
        except VariantPreparationError as exc:
            raise HTTPException(
                status_code=422,
                detail=str(exc),
                headers={"X-ARIANE-Error-Code": exc.code},
            ) from exc
        except RequiredEvidenceUnavailableError as exc:
            headers = {
                "X-ARIANE-Error-Code": exc.code,
                "X-ARIANE-Retryable": str(exc.retryable).lower(),
            }
            if exc.retryable:
                headers["Retry-After"] = "5"
            raise HTTPException(
                status_code=503,
                detail=str(exc),
                headers=headers,
            ) from exc
        except EvidenceExecutionError as exc:
            raise HTTPException(status_code=503, detail=str(exc)) from exc

    async def classify_cached(
        self,
        gene: str,
        c_notation: str,
        p_notation: str = "",
        dup_type: str = "Unknown",
        reference_transcript: str = "",
    ) -> tuple[ClassificationResult, str, str]:
        async def refresh_cached(result: ClassificationResult) -> ClassificationResult:
            clinvar, clingen, _diagnostics = await self._orchestration().lookup_external(
                result.gene,
                result.c_notation,
            )
            return self.presentation.refresh_external(
                result,
                clinvar=clinvar,
                clingen=clingen,
            )

        return await self.service.classify_cached(
            gene,
            c_notation,
            p_notation,
            dup_type,
            reference_transcript,
            cache_repository=self._cache_repository(),
            classify_uncached=self.classify_uncached,
            refresh_cached=refresh_cached,
        )

    def reserve_public_api_classifications(
        self,
        request: Request,
        response: Response,
        units: int,
    ) -> None:
        """Reserve per-key work before classification; browser traffic is separate."""
        api_key_id = str(getattr(request.state, "api_key_id", "") or "").strip()
        if not api_key_id:
            return
        repository = self._quota_repository()
        if repository is None:
            raise HTTPException(
                status_code=503,
                detail={
                    "code": "api_quota_unavailable",
                    "message": "Public API quota storage is unavailable",
                    "retryable": True,
                },
                headers={"Retry-After": "60"},
            )
        try:
            reservation = repository.reserve_public_api_classifications(
                api_key_id=api_key_id,
                units=units,
                limit=self._daily_classification_limit(),
            )
        except Exception as exc:
            logging.getLogger("ariane.api_quota").exception(
                "Public API quota reservation failed for key %s", api_key_id
            )
            raise HTTPException(
                status_code=503,
                detail={
                    "code": "api_quota_unavailable",
                    "message": "Public API quota storage is unavailable",
                    "retryable": True,
                },
                headers={"Retry-After": "60"},
            ) from exc
        reset_epoch = str(int(reservation.reset_at.timestamp()))
        quota_headers = {
            "X-RateLimit-Limit": str(reservation.limit),
            "X-RateLimit-Remaining": str(reservation.remaining),
            "X-RateLimit-Reset": reset_epoch,
        }
        if not reservation.allowed:
            retry_after = max(
                1,
                int(
                    reservation.reset_at.timestamp()
                    - datetime.now(timezone.utc).timestamp()
                ),
            )
            raise HTTPException(
                status_code=429,
                detail={
                    "code": "daily_classification_quota_exceeded",
                    "message": (
                        "The API key has reached its daily classification quota. "
                        "The quota resets at 00:00 UTC."
                    ),
                    "retryable": True,
                },
                headers={**quota_headers, "Retry-After": str(retry_after)},
            )
        for name, value in quota_headers.items():
            response.headers[name] = value

    async def classify_variant(
        self,
        request_model: VariantRequest,
        request: Request,
        http_response: Response,
    ) -> ClassificationResult:
        self.reserve_public_api_classifications(request, http_response, 1)
        input_data = request_model.model_dump(mode="json")
        identity = resolve_usage_identity(request)
        identity.set_cookie(http_response)
        try:
            completed = await self.service.classify_single(
                request_model,
                usage=UsageContext(
                    actor_id=identity.actor_id,
                    actor_type=identity.actor_type,
                    request_id=request.state.request_id,
                ),
                usage_repository=self._usage_repository(),
                classify_cached=self.classify_cached,
            )
        except Exception as exc:
            self.audit(
                request,
                "classification_error",
                level="exception",
                input=input_data,
                error_type=type(exc).__name__,
                error=str(exc)[:2000],
            )
            raise
        result = completed.result
        request.state.classification_execution = {
            "policy": completed.policy,
            "classifier_fingerprint": completed.classifier_fingerprint,
            "cache_status": completed.cache_status,
            "duration_ms": completed.duration_ms,
        }
        http_response.headers["X-ARIANE-Cache-Status"] = completed.cache_status
        http_response.headers[
            "X-ARIANE-Classifier-Fingerprint"
        ] = completed.classifier_fingerprint
        self.audit(
            request,
            "classification_completed",
            input=input_data,
            result={
                "predicted_class": result.predicted_class,
                "predicted_label": result.predicted_label,
                "total_points": result.total_points,
                "evidence_direction": result.evidence_direction,
                "mixed_evidence": result.mixed_evidence,
                "pathogenic_points": result.pathogenic_points,
                "benign_points": result.benign_points,
                "classification_cache": completed.cache_status,
                "evidence_interactions": [
                    warning.model_dump(mode="json")
                    for warning in result.evidence_interactions
                ],
                "spliceai_audit": result.spliceai_audit.model_dump(mode="json")
                if result.spliceai_audit else None,
            },
        )
        return result

    async def classify_batch(
        self,
        request_model: BatchRequest,
        request: Request,
        http_response: Response,
    ) -> BatchResponse:
        identity = resolve_usage_identity(request)
        identity.set_cookie(http_response)
        prepared = self.service.prepare_batch(request_model)
        for event in prepared.audit_events:
            self.audit(request, event.event, level=event.level, **event.fields)
        if prepared.classification_units:
            self.reserve_public_api_classifications(
                request,
                http_response,
                prepared.classification_units,
            )
        completed = await self.service.classify_batch(
            prepared,
            usage=UsageContext(
                actor_id=identity.actor_id,
                actor_type=identity.actor_type,
                request_id=request.state.request_id,
            ),
            usage_repository=self._usage_repository(),
            classify_cached=self.classify_cached,
            map_error=self.map_execution_error,
        )
        for event in completed.audit_events:
            self.audit(request, event.event, level=event.level, **event.fields)
        request.state.batch_duration_ms = completed.duration_ms
        return completed.response

    def router(self) -> APIRouter:
        router = APIRouter()
        paired = PairedApiRoutes(router)
        paired.post(
            "/classify",
            response_model=ClassificationResult,
        )(
            self.classify_variant,
        )
        router.add_api_route(
            "/api/classify/batch",
            self.classify_batch,
            methods=["POST"],
            dependencies=[Depends(require_public_api_key)],
            response_model=BatchResponse,
        )
        router.include_router(create_public_api_router(
            classify_single=self.classify_variant,
            classify_batch=self.classify_batch,
            engine=CLASSIFICATION_ENGINE_ID,
            batch_concurrency=self.service.batch_concurrency,
            per_key_classifications_per_utc_day=self._daily_classification_limit(),
        ))
        return router
