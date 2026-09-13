"""Manual-evidence, client-diagnostic, and normalization HTTP routes."""

from __future__ import annotations

from collections.abc import Callable

from fastapi import APIRouter, HTTPException, Request, Response

from backend.api.routing import PairedApiRoutes
from backend.classification_dag import DagNodeExecutionError
from backend.contracts import (
    ClientValidationRequest,
    ManualEvidenceRequest,
    ManualEvidenceResult,
    ManualEvidenceStatusRequest,
    ManualEvidenceStatusResponse,
    Ps1ReferenceResolutionRequest,
    Ps1ReferenceResolutionResponse,
    VariantNormalizationResponse,
    VariantRequest,
)
from backend.services import ManualEvidenceService


AuditWriter = Callable[..., None]


def create_manual_router(
    *,
    service: ManualEvidenceService,
    audit: AuditWriter,
) -> APIRouter:
    router = APIRouter()
    paired = PairedApiRoutes(router)

    @paired.post("/audit/client-validation", status_code=204)
    async def client_validation_error(req: ClientValidationRequest, request: Request):
        audit(
            request,
            "client_validation_error",
            level="warning",
            form=req.form,
            input=req.input,
            error=req.error,
        )
        return Response(status_code=204)

    @paired.post("/manual-evidence/evaluate")
    async def evaluate_manual_evidence_endpoint(
        req: ManualEvidenceRequest,
        request: Request,
    ) -> ManualEvidenceResult:
        try:
            response = service.evaluate(req)
        except DagNodeExecutionError as exc:
            cause = exc.__cause__
            if isinstance(cause, ValueError):
                detail = str(cause)
                status_code = 422
            else:
                detail = (
                    "Manual evidence evaluation could not complete at internal "
                    f"step {exc.node_id}. No adjusted classification was returned."
                )
                status_code = 503
            audit(
                request,
                "manual_evidence_error",
                level="warning",
                input=req.model_dump(mode="json"),
                error=str(exc),
                trace=[entry.as_dict() for entry in exc.trace],
            )
            raise HTTPException(status_code=status_code, detail=detail) from exc

        audit(
            request,
            "manual_evidence_completed",
            input=req.model_dump(mode="json"),
            result={
                "predicted_class": response.predicted_class,
                "predicted_label": response.predicted_label,
                "total_points": response.total_points,
                "evidence_interactions": [
                    warning.model_dump(mode="json")
                    for warning in response.evidence_interactions
                ],
            },
        )
        return response

    @paired.post("/manual-evidence/status")
    async def manual_evidence_status_endpoint(
        req: ManualEvidenceStatusRequest,
    ) -> ManualEvidenceStatusResponse:
        return service.status(req)

    @paired.post("/manual-evidence/resolve-ps1-reference")
    async def resolve_ps1_reference_endpoint(
        req: Ps1ReferenceResolutionRequest,
        request: Request,
    ) -> Ps1ReferenceResolutionResponse:
        try:
            response = await service.resolve_ps1_reference(req)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        except Exception as exc:
            audit(
                request,
                "ps1_reference_resolution_error",
                level="exception",
                input=req.model_dump(mode="json"),
                error_type=type(exc).__name__,
                error=str(exc)[:2000],
            )
            raise HTTPException(
                status_code=503,
                detail="PS1 reference facts could not be resolved; no criterion was added.",
            ) from exc
        audit(
            request,
            "ps1_reference_resolved",
            input=req.model_dump(mode="json"),
            result=response.model_dump(mode="json"),
        )
        return response

    @paired.post("/normalize")
    async def normalize_variant(
        req: VariantRequest,
        request: Request,
    ) -> VariantNormalizationResponse:
        response = VariantNormalizationResponse(
            gene=req.gene,
            submitted_notation=req.submitted_notation,
            c_notation=req.c_notation,
            p_notation=req.p_notation,
            reference_transcript=req.reference_transcript,
            normalization_source=req.normalization_source,
            consequence_status=req.consequence_status,
            normalization_provenance=req.normalization_provenance,
            protein_consequence_explanation=req.protein_consequence_explanation,
            assembly=req.assembly,
        )
        audit(
            request,
            "variant_normalized",
            input={
                "gene": req.gene,
                "notation": req.submitted_notation,
                "assembly": req.assembly,
            },
            result=response.model_dump(mode="json"),
        )
        return response

    return router
