"""Authenticated API for immutable manual-review records."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
from typing import Annotated

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from pydantic import BaseModel, field_validator, model_validator

from backend.api.admin import require_admin
from backend.classification_dag import DagNodeExecutionError
from backend.classification_dag.manual import execute_manual_evidence
from backend.contracts import ClassificationResult, ManualEvidenceRequest
from backend.infrastructure.review_repository import ReviewRecordRepository


router = APIRouter(prefix="/api/review-records", tags=["review records"])
logger = logging.getLogger("ariane.audit")


class ReviewRecordCreateRequest(BaseModel):
    module1_result: ClassificationResult
    manual_evidence: ManualEvidenceRequest
    reviewer_role: str

    @field_validator("reviewer_role")
    @classmethod
    def require_reviewer_role(cls, value: str) -> str:
        normalized = value.strip()
        if not normalized:
            raise ValueError("Reviewer role is required")
        return normalized

    @model_validator(mode="after")
    def require_matching_module1_context(self):
        context = self.manual_evidence.variant_context
        if context is None:
            raise ValueError("A variant context is required for a persisted review")
        if (
            context.gene != self.module1_result.gene
            or context.c_notation != self.module1_result.c_notation
            or context.p_notation != self.module1_result.p_notation
        ):
            raise ValueError(
                "The manual-evidence variant context must match the Module 1 result"
            )
        submitted = [
            (item.name, item.strength, item.points, item.applies)
            for item in self.manual_evidence.base_criteria
        ]
        module1 = [
            (item.name, item.strength, item.points, item.applies)
            for item in self.module1_result.criteria
        ]
        if submitted != module1:
            raise ValueError(
                "The manual-evidence base criteria must match the Module 1 result"
            )
        return self


class ReviewApprovalRequest(BaseModel):
    approver_name: str
    approver_role: str
    approval_comment: str = ""
    attestation_confirmed: bool = False

    @field_validator("approver_name", "approver_role")
    @classmethod
    def require_approval_identity(cls, value: str) -> str:
        normalized = value.strip()
        if not normalized:
            raise ValueError("Approver name and role are required")
        return normalized

    @model_validator(mode="after")
    def require_attestation(self):
        if not self.attestation_confirmed:
            raise ValueError(
                "Approval requires confirmation that the evidence and result were reviewed"
            )
        return self


def _repository() -> ReviewRecordRepository:
    return ReviewRecordRepository()


def _audit(request: Request, event: str, **fields) -> None:
    logger.info(json.dumps({
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "log_type": "ariane_audit",
        "request_id": getattr(request.state, "request_id", ""),
        "source_ip": request.client.host if request.client else "unknown",
        "method": request.method,
        "path": request.url.path,
        "event": event,
        **fields,
    }, ensure_ascii=True, separators=(",", ":")))


@router.post("")
async def create_review_record(
    payload: ReviewRecordCreateRequest,
    request: Request,
    authenticated_account: Annotated[str, Depends(require_admin)],
):
    manual = payload.manual_evidence
    try:
        execution = execute_manual_evidence(
            [criterion.model_dump() for criterion in manual.base_criteria],
            [criterion.model_dump() for criterion in manual.manual_criteria],
            manual.variant_context.model_dump() if manual.variant_context else None,
        )
    except DagNodeExecutionError as exc:
        cause = exc.__cause__
        detail = str(cause) if isinstance(cause, ValueError) else (
            "The manual evidence could not be recomputed for persistence"
        )
        raise HTTPException(status_code=422, detail=detail) from exc

    amended_result = {
        **execution.result,
        "assessor": manual.assessor,
        "assessed_at": manual.assessed_at,
    }
    context = manual.variant_context.model_dump()
    record = _repository().create_draft(
        variant_context=context,
        module1_result=payload.module1_result.model_dump(mode="json"),
        manual_evidence={
            "variant_context": context,
            "criteria": [
                item.model_dump(mode="json")
                for item in manual.manual_criteria
                if item.enabled
            ],
            "assessor": manual.assessor,
            "assessed_at": manual.assessed_at,
        },
        amended_result=amended_result,
        reviewer_name=manual.assessor,
        reviewer_role=payload.reviewer_role,
        review_date=manual.assessed_at,
        authenticated_account=authenticated_account,
    )
    _audit(
        request,
        "review_record_draft_created",
        record_id=record["record_id"],
        variant_key=record["variant_key"],
        version=record["version"],
        authenticated_account=authenticated_account,
    )
    return record


@router.post("/{record_id}/approve")
async def approve_review_record(
    record_id: str,
    payload: ReviewApprovalRequest,
    request: Request,
    authenticated_account: Annotated[str, Depends(require_admin)],
):
    try:
        record = _repository().approve(
            record_id,
            approver_name=payload.approver_name,
            approver_role=payload.approver_role,
            approval_comment=payload.approval_comment.strip(),
            authenticated_account=authenticated_account,
        )
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    _audit(
        request,
        "review_record_approved",
        record_id=record["record_id"],
        parent_record_id=record["parent_record_id"],
        variant_key=record["variant_key"],
        version=record["version"],
        authenticated_account=authenticated_account,
    )
    return record


@router.get("")
async def list_review_records(
    authenticated_account: Annotated[str, Depends(require_admin)],
    limit: int = Query(default=100, ge=1, le=500),
):
    return {"records": _repository().list(limit=limit)}


@router.get("/{record_id}")
async def get_review_record(
    record_id: str,
    authenticated_account: Annotated[str, Depends(require_admin)],
):
    try:
        return _repository().get(record_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
