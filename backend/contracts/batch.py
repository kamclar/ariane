"""Batch classification HTTP schemas."""

from __future__ import annotations

from pydantic import BaseModel, SkipValidation, field_validator

from backend.contracts.result import ClassificationResult
from backend.contracts.variant import VariantRequest


MAXIMUM_BATCH_ITEMS = 200


class BatchItemResult(BaseModel):
    index: int
    status: str
    variant: str
    error: str | None = None
    error_code: str = ""
    error_retryable: bool = False
    cache_status: str = ""
    classifier_fingerprint: str = ""
    policy_id: str = ""
    policy_version: str = ""
    duration_ms: float | None = None
    result: ClassificationResult | None = None


class BatchRequest(BaseModel):
    variants: list[SkipValidation[VariantRequest]]

    @field_validator("variants")
    @classmethod
    def validate_count(cls, value: list) -> list:
        if not value:
            raise ValueError("At least one variant required")
        if len(value) > MAXIMUM_BATCH_ITEMS:
            raise ValueError(f"Maximum {MAXIMUM_BATCH_ITEMS} variants per batch")
        return value


class BatchResponse(BaseModel):
    total: int
    success_count: int
    error_count: int
    results: list[BatchItemResult]
