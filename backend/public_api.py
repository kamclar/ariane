"""Versioned public API contract and projections.

This module contains transport models only. It never evaluates evidence or
derives a classification. Public responses are projections of the result
returned by the production classification DAG.
"""

from __future__ import annotations

from collections.abc import Awaitable, Callable
from datetime import datetime, timezone
from typing import Any, Literal

from fastapi import APIRouter, Depends, HTTPException, Request, Response
from pydantic import BaseModel, Field

from backend.gene_policy import active_genes, get_gene_policy
from backend.models import (
    BatchRequest,
    BatchResponse,
    BatchItemResult,
    ClassificationResult,
    ProteinPs1ReviewRecommendation,
    RnaReviewRecommendation,
    VariantRequest,
)
from backend.version import ARIANE_VERSION
from backend.api_auth import API_KEY_HEADER_NAME, require_public_api_key


PUBLIC_API_VERSION = "1.0"
PUBLIC_API_MAXIMUM_BATCH_ITEMS = 10
RECOMMENDED_UNCACHED_BATCH_ITEMS = 5
_REVIEW_FIELDS = (
    "rna_review",
    "splice_ps1_review",
    "protein_ps1_review",
    "initiation_review",
)


class PublicApiPolicy(BaseModel):
    gene: str
    policy_id: str
    policy_version: str
    reference_transcript: str


class PublicApiMetadata(BaseModel):
    api_version: str = PUBLIC_API_VERSION
    application_version: str = ARIANE_VERSION
    request_id: str
    generated_at: str
    classification_engine: str = "dag"
    policy: PublicApiPolicy | None = None
    classifier_fingerprint: str = ""
    cache_status: str = ""
    duration_ms: float | None = None


class PublicApiReviewItem(BaseModel):
    kind: Literal[
        "rna",
        "splice_ps1",
        "protein_ps1",
        "initiation_codon",
    ]
    recommended: bool
    priority: str
    title: str
    summary: str
    potential_branches: list[str] = Field(default_factory=list)
    details: dict[str, Any] = Field(default_factory=dict)


class PublicApiManualReview(BaseModel):
    available: bool
    recommended: bool
    affects_automatic_classification: Literal[False] = False
    notice: str = (
        "Review aids do not add criteria or points and do not change the "
        "automatic classification. An amended result requires documented "
        "manual evidence and separate backend evaluation."
    )
    items: list[PublicApiReviewItem] = Field(default_factory=list)


class PublicApiClassificationResponse(BaseModel):
    metadata: PublicApiMetadata
    classification: ClassificationResult
    manual_review: PublicApiManualReview


class PublicApiError(BaseModel):
    code: str
    message: str
    retryable: bool = False


class PublicApiErrorResponse(BaseModel):
    metadata: PublicApiMetadata
    error: PublicApiError


class PublicApiBatchItem(BaseModel):
    index: int
    status: Literal["ok", "error"]
    variant: str
    metadata: PublicApiMetadata | None = None
    classification: ClassificationResult | None = None
    manual_review: PublicApiManualReview | None = None
    error: PublicApiError | None = None


class PublicApiBatchResponse(BaseModel):
    metadata: PublicApiMetadata
    total: int
    success_count: int
    error_count: int
    results: list[PublicApiBatchItem]


class PublicApiLimits(BaseModel):
    maximum_batch_items: int
    recommended_uncached_batch_items: int
    batch_concurrency: int
    requests_per_second: int
    request_burst: int
    per_key_requests_per_minute: int
    per_key_request_burst: int


class PublicApiGene(BaseModel):
    symbol: str
    reference_transcript: str
    reference_protein: str
    policy_id: str
    policy_version: str


class PublicApiCapabilities(BaseModel):
    api_version: str = PUBLIC_API_VERSION
    application_version: str = ARIANE_VERSION
    status: Literal["beta"] = "beta"
    classification_engine: Literal["dag"] = "dag"
    supported_genes: list[PublicApiGene]
    limits: PublicApiLimits
    endpoints: list[str]
    manual_review_semantics: str
    authentication_required: Literal[True] = True
    authentication_header: str = API_KEY_HEADER_NAME


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def build_metadata(
    *,
    request_id: str,
    engine: str,
    policy: dict[str, str] | None = None,
    classifier_fingerprint: str = "",
    cache_status: str = "",
    duration_ms: float | None = None,
) -> PublicApiMetadata:
    return PublicApiMetadata(
        request_id=request_id,
        generated_at=utc_timestamp(),
        classification_engine=engine,
        policy=PublicApiPolicy(**policy) if policy else None,
        classifier_fingerprint=classifier_fingerprint,
        cache_status=cache_status,
        duration_ms=round(duration_ms, 1) if duration_ms is not None else None,
    )


def build_manual_review(result: ClassificationResult) -> PublicApiManualReview:
    definitions: tuple[
        tuple[str, RnaReviewRecommendation | ProteinPs1ReviewRecommendation | None],
        ...,
    ] = (
        ("rna", result.rna_review),
        ("splice_ps1", result.splice_ps1_review),
        ("protein_ps1", result.protein_ps1_review),
        ("initiation_codon", result.initiation_review),
    )
    items: list[PublicApiReviewItem] = []
    for kind, recommendation in definitions:
        if recommendation is None:
            continue
        visible = recommendation.recommended or bool(
            getattr(recommendation, "display", False)
        )
        if not visible:
            continue
        details = recommendation.model_dump(mode="json", exclude_none=True)
        items.append(PublicApiReviewItem(
            kind=kind,
            recommended=recommendation.recommended,
            priority=recommendation.priority,
            title=recommendation.title,
            summary=recommendation.summary,
            potential_branches=list(recommendation.potential_branches),
            details=details,
        ))
    return PublicApiManualReview(
        available=bool(items),
        recommended=any(item.recommended for item in items),
        items=items,
    )


def automatic_result(result: ClassificationResult) -> ClassificationResult:
    """Return the scored result without duplicating review aids."""
    return result.model_copy(update={field: None for field in _REVIEW_FIELDS})


def build_classification_response(
    result: ClassificationResult,
    metadata: PublicApiMetadata,
) -> PublicApiClassificationResponse:
    return PublicApiClassificationResponse(
        metadata=metadata,
        classification=automatic_result(result),
        manual_review=build_manual_review(result),
    )


def build_batch_item(
    item: BatchItemResult,
    *,
    request_id: str,
    engine: str,
) -> PublicApiBatchItem:
    if item.status == "ok" and item.result is not None:
        metadata = build_metadata(
            request_id=request_id,
            engine=engine,
            policy={
                "gene": item.result.gene,
                "policy_id": item.policy_id,
                "policy_version": item.policy_version,
                "reference_transcript": item.result.reference_transcript,
            },
            classifier_fingerprint=item.classifier_fingerprint,
            cache_status=item.cache_status,
            duration_ms=item.duration_ms,
        )
        return PublicApiBatchItem(
            index=item.index,
            status="ok",
            variant=item.variant,
            metadata=metadata,
            classification=automatic_result(item.result),
            manual_review=build_manual_review(item.result),
        )
    message = item.error or "Classification failed"
    if item.error_code == "internal_error":
        message = "Classification failed because of an internal server error"
    return PublicApiBatchItem(
        index=item.index,
        status="error",
        variant=item.variant,
        error=PublicApiError(
            code=item.error_code or "classification_failed",
            message=message,
            retryable=item.error_retryable,
        ),
    )


SingleClassificationHandler = Callable[
    [VariantRequest, Request, Response],
    Awaitable[ClassificationResult],
]
BatchClassificationHandler = Callable[
    [BatchRequest, Request, Response],
    Awaitable[BatchResponse],
]


def create_public_api_router(
    *,
    classify_single: SingleClassificationHandler,
    classify_batch: BatchClassificationHandler,
    engine: str,
    batch_concurrency: int,
    requests_per_second: int = 5,
    request_burst: int = 20,
    per_key_requests_per_minute: int = 30,
    per_key_request_burst: int = 3,
) -> APIRouter:
    """Create the versioned transport router around shared handlers."""
    router = APIRouter(prefix="/api/v1", tags=["Public API v1"])

    @router.get("/capabilities", response_model=PublicApiCapabilities)
    async def capabilities() -> PublicApiCapabilities:
        return PublicApiCapabilities(
            supported_genes=[
                PublicApiGene(
                    symbol=symbol,
                    reference_transcript=get_gene_policy(symbol)["gene_config"][
                        "reference_transcript"
                    ],
                    reference_protein=get_gene_policy(symbol)["gene_config"][
                        "reference_protein"
                    ],
                    policy_id=get_gene_policy(symbol)["policy"]["runtime_policy_id"],
                    policy_version=get_gene_policy(symbol)["policy"]["version"],
                )
                for symbol in active_genes()
            ],
            limits=PublicApiLimits(
                maximum_batch_items=PUBLIC_API_MAXIMUM_BATCH_ITEMS,
                recommended_uncached_batch_items=RECOMMENDED_UNCACHED_BATCH_ITEMS,
                batch_concurrency=batch_concurrency,
                requests_per_second=requests_per_second,
                request_burst=request_burst,
                per_key_requests_per_minute=per_key_requests_per_minute,
                per_key_request_burst=per_key_request_burst,
            ),
            endpoints=[
                "/api/v1/capabilities",
                "/api/v1/classify",
                "/api/v1/classify/batch",
            ],
            manual_review_semantics=(
                "Manual review recommendations are review aids only. They add "
                "no criteria or points and do not change the automatic "
                "classification."
            ),
        )

    @router.post(
        "/classify",
        response_model=PublicApiClassificationResponse,
        response_model_exclude_none=True,
        responses={
            401: {"model": PublicApiErrorResponse},
            422: {"model": PublicApiErrorResponse},
            503: {"model": PublicApiErrorResponse},
        },
    )
    async def classify(
        req: VariantRequest,
        request: Request,
        http_response: Response,
        _api_key_id: str = Depends(require_public_api_key),
    ) -> PublicApiClassificationResponse:
        try:
            result = await classify_single(req, request, http_response)
        except HTTPException:
            raise
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail="Classification failed because of an internal server error",
            ) from exc
        execution = request.state.classification_execution
        metadata = build_metadata(
            request_id=request.state.request_id,
            engine=engine,
            policy=execution["policy"],
            classifier_fingerprint=execution["classifier_fingerprint"],
            cache_status=execution["cache_status"],
            duration_ms=execution["duration_ms"],
        )
        return build_classification_response(result, metadata)

    @router.post(
        "/classify/batch",
        response_model=PublicApiBatchResponse,
        response_model_exclude_none=True,
        responses={
            401: {"model": PublicApiErrorResponse},
            422: {"model": PublicApiErrorResponse},
            503: {"model": PublicApiErrorResponse},
        },
    )
    async def classify_many(
        req: BatchRequest,
        request: Request,
        http_response: Response,
        _api_key_id: str = Depends(require_public_api_key),
    ) -> PublicApiBatchResponse:
        if len(req.variants) > PUBLIC_API_MAXIMUM_BATCH_ITEMS:
            raise HTTPException(
                status_code=422,
                detail={
                    "code": "invalid_request",
                    "message": (
                        "Public API v1 accepts at most "
                        f"{PUBLIC_API_MAXIMUM_BATCH_ITEMS} variants per synchronous "
                        "batch"
                    ),
                    "retryable": False,
                },
            )
        result = await classify_batch(req, request, http_response)
        metadata = build_metadata(
            request_id=request.state.request_id,
            engine=engine,
            duration_ms=request.state.batch_duration_ms,
        )
        return PublicApiBatchResponse(
            metadata=metadata,
            total=result.total,
            success_count=result.success_count,
            error_count=result.error_count,
            results=[
                build_batch_item(
                    item,
                    request_id=request.state.request_id,
                    engine=engine,
                )
                for item in result.results
            ],
        )

    return router
