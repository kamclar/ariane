"""Manual evidence and review HTTP schemas."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from backend.contracts.result import CriterionResult, EvidenceInteractionWarning
from backend.policy.gene import active_genes


class ManualCriterionInput(BaseModel):
    code: str
    enabled: bool = False
    evidence: dict[str, Any] = Field(default_factory=dict)
    notes: str = ""
    references: list[str] = Field(default_factory=list)

    @model_validator(mode="before")
    @classmethod
    def reject_strength_override(cls, data):
        if isinstance(data, dict) and data.get("override_strength") not in {None, ""}:
            raise ValueError(
                "Manual strength overrides are not permitted; criterion strength "
                "is derived from the configured VCEP evidence thresholds"
            )
        return data

    @field_validator("code")
    @classmethod
    def validate_manual_code(cls, value):
        code = value.strip().upper()
        if code not in {
            "PS3", "PS4", "PM3", "PP1", "PP4", "BS2", "BS3", "BS4", "BP5",
            "PVS1_RNA", "BP7_RNA", "PVS1_INIT", "PS1_SPLICE", "PS1_PROTEIN",
        }:
            raise ValueError("Unsupported manually reviewed criterion")
        return code


class ManualVariantContext(BaseModel):
    gene: str
    c_notation: str
    p_notation: str

    @field_validator("gene")
    @classmethod
    def validate_context_gene(cls, value):
        gene = value.strip().upper()
        supported = set(active_genes())
        if gene not in supported:
            raise ValueError(f"Gene must be one of: {', '.join(sorted(supported))}")
        return gene

    @field_validator("c_notation", "p_notation")
    @classmethod
    def require_variant_notation(cls, value):
        if not value.strip():
            raise ValueError("Manual evidence variant context must not be empty")
        return value.strip()


class ManualEvidenceStatusRequest(BaseModel):
    base_criteria: list[CriterionResult]
    manual_criteria: list[ManualCriterionInput]
    variant_context: ManualVariantContext | None = None


class ManualCriterionStatus(BaseModel):
    code: str
    status: Literal["not_started", "incomplete", "ready"]
    message: str = ""
    suggested_strength: str | None = None


class ManualEvidenceStatusResponse(BaseModel):
    criteria: list[ManualCriterionStatus] = Field(default_factory=list)


class ManualEvidenceRequest(BaseModel):
    base_criteria: list[CriterionResult]
    manual_criteria: list[ManualCriterionInput]
    variant_context: ManualVariantContext | None = None
    assessor: str
    assessed_at: str

    @field_validator("assessor", "assessed_at")
    @classmethod
    def require_audit_value(cls, value):
        if not value.strip():
            raise ValueError("Audit fields must not be empty")
        return value.strip()

    @model_validator(mode="after")
    def require_complete_enabled_records(self):
        enabled = [item for item in self.manual_criteria if item.enabled]
        if not enabled:
            raise ValueError("Select at least one manually reviewed criterion")
        for item in enabled:
            if item.code == "PS1_PROTEIN" and self.variant_context is None:
                raise ValueError(
                    "PS1_PROTEIN requires the assessed variant context so the backend "
                    "can verify the protein consequence and nucleotide change"
                )
            if not item.notes.strip():
                raise ValueError(f"{item.code} requires evidence notes")
            if not any(reference.strip() for reference in item.references):
                raise ValueError(f"{item.code} requires at least one evidence reference")
        return self


class ManualCriterionResult(BaseModel):
    code: str
    applies: bool
    suggested_strength: str | None = None
    selected_strength: str | None = None
    points: int = 0
    reason: str = ""
    threshold_note: str = ""
    overridden: bool = False
    notes: str = ""
    references: list[str] = Field(default_factory=list)
    single_strong_likely_benign_eligible: bool = False
    single_strong_likely_benign_basis: str = ""
    independent_evidence_contribution_count: int = 0


class ManualEvidenceResult(BaseModel):
    predicted_class: int
    predicted_label: str
    total_points: int
    classification_note: str = ""
    manual_criteria: list[ManualCriterionResult]
    evidence_interactions: list[EvidenceInteractionWarning] = Field(default_factory=list)
    assessor: str
    assessed_at: str
