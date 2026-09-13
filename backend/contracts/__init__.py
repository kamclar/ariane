"""Stable imports for ARIANE HTTP request and response schemas."""

from backend.contracts.batch import (
    MAXIMUM_BATCH_ITEMS,
    BatchItemResult,
    BatchRequest,
    BatchResponse,
)
from backend.contracts.client import ClientValidationRequest
from backend.contracts.ps1 import (
    Ps1ReferenceResolutionRequest,
    Ps1ReferenceResolutionResponse,
    Ps1ResolvedVariant,
)
from backend.contracts.result import (
    AlphaMissenseResult,
    ClassificationResult,
    ClinicalAnnotation,
    ClinicalAnnotationPublication,
    ClinicalLrAudit,
    ClinicalLrThresholdComparison,
    CriterionResult,
    EvidenceInteractionWarning,
    ExternalComparison,
    ExternalSubmitter,
    ProteinPs1Candidate,
    ProteinPs1ReviewRecommendation,
    RnaReviewRecommendation,
    SpliceAIAudit,
    VusExplanation,
)
from backend.contracts.review import (
    ManualCriterionInput,
    ManualCriterionResult,
    ManualCriterionStatus,
    ManualEvidenceRequest,
    ManualEvidenceResult,
    ManualEvidenceStatusRequest,
    ManualEvidenceStatusResponse,
    ManualVariantContext,
)
from backend.contracts.variant import VariantNormalizationResponse, VariantRequest


CLASS_LABELS = {
    5: "Pathogenic",
    4: "Likely Pathogenic",
    3: "VUS",
    2: "Likely Benign",
    1: "Benign",
}

__all__ = (
    "AlphaMissenseResult",
    "BatchItemResult",
    "BatchRequest",
    "BatchResponse",
    "CLASS_LABELS",
    "ClientValidationRequest",
    "ClassificationResult",
    "ClinicalAnnotation",
    "ClinicalAnnotationPublication",
    "ClinicalLrAudit",
    "ClinicalLrThresholdComparison",
    "CriterionResult",
    "EvidenceInteractionWarning",
    "ExternalComparison",
    "ExternalSubmitter",
    "MAXIMUM_BATCH_ITEMS",
    "ManualCriterionInput",
    "ManualCriterionResult",
    "ManualCriterionStatus",
    "ManualEvidenceRequest",
    "ManualEvidenceResult",
    "ManualEvidenceStatusRequest",
    "ManualEvidenceStatusResponse",
    "ManualVariantContext",
    "ProteinPs1Candidate",
    "ProteinPs1ReviewRecommendation",
    "Ps1ReferenceResolutionRequest",
    "Ps1ReferenceResolutionResponse",
    "Ps1ResolvedVariant",
    "RnaReviewRecommendation",
    "SpliceAIAudit",
    "VariantNormalizationResponse",
    "VariantRequest",
    "VusExplanation",
)
