"""Application workflow for expert-reviewed evidence."""

from __future__ import annotations

from backend.classification_dag.manual import execute_manual_evidence
from backend.contracts import (
    EvidenceInteractionWarning,
    ManualCriterionResult,
    ManualCriterionStatus,
    ManualEvidenceRequest,
    ManualEvidenceResult,
    ManualEvidenceStatusRequest,
    ManualEvidenceStatusResponse,
    Ps1ReferenceResolutionRequest,
    Ps1ReferenceResolutionResponse,
)
from backend.review.validation import manual_criterion_statuses
from backend.services.ps1_reference_resolution import resolve_ps1_reference


class ManualEvidenceService:
    """Evaluate complete manual records without exposing DAG details to HTTP routes."""

    def evaluate(self, request: ManualEvidenceRequest) -> ManualEvidenceResult:
        execution = execute_manual_evidence(
            [criterion.model_dump() for criterion in request.base_criteria],
            [criterion.model_dump() for criterion in request.manual_criteria],
            request.variant_context.model_dump() if request.variant_context else None,
        )
        result = execution.result
        return ManualEvidenceResult(
            predicted_class=result["predicted_class"],
            predicted_label=result["predicted_label"],
            total_points=result["total_points"],
            classification_note=result["classification_note"],
            manual_criteria=[
                ManualCriterionResult(**criterion) for criterion in result["manual_criteria"]
            ],
            evidence_interactions=[
                EvidenceInteractionWarning(**warning)
                for warning in result["evidence_interactions"]
            ],
            assessor=request.assessor,
            assessed_at=request.assessed_at,
        )

    def status(
        self,
        request: ManualEvidenceStatusRequest,
    ) -> ManualEvidenceStatusResponse:
        statuses = manual_criterion_statuses(
            [criterion.model_dump() for criterion in request.base_criteria],
            [criterion.model_dump() for criterion in request.manual_criteria],
            request.variant_context.model_dump() if request.variant_context else None,
        )
        return ManualEvidenceStatusResponse(
            criteria=[ManualCriterionStatus(**status) for status in statuses]
        )

    async def resolve_ps1_reference(
        self,
        request: Ps1ReferenceResolutionRequest,
    ) -> Ps1ReferenceResolutionResponse:
        result = await resolve_ps1_reference(
            request.gene,
            request.assessed_c_notation,
            request.reference_c_notation,
        )
        return Ps1ReferenceResolutionResponse(**result)
