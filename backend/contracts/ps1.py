"""PS1 reference-resolution HTTP schemas."""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator

from backend.policy.gene import active_genes


class Ps1ReferenceResolutionRequest(BaseModel):
    gene: str
    assessed_c_notation: str
    reference_c_notation: str

    @field_validator("gene")
    @classmethod
    def validate_ps1_gene(cls, value):
        gene = value.strip().upper()
        supported = set(active_genes())
        if gene not in supported:
            raise ValueError(f"Gene must be one of: {', '.join(sorted(supported))}")
        return gene

    @field_validator("assessed_c_notation", "reference_c_notation")
    @classmethod
    def require_ps1_notation(cls, value):
        if not value.strip():
            raise ValueError("Both assessed and reference c. notations are required")
        return value.strip()


class Ps1ResolvedVariant(BaseModel):
    gene: str
    reference_transcript: str
    c_notation: str
    p_notation: str
    spliceai_score: float | None = None
    spliceai_status: str = "unavailable"
    spliceai_reason: str = ""


class Ps1ReferenceResolutionResponse(BaseModel):
    assessed: Ps1ResolvedVariant
    reference: Ps1ResolvedVariant
    same_missense_substitution: bool
    different_nucleotide_change: bool
    clinvar_status: str = "not_found"
    clinvar_error: str = ""
    clinvar_variation_id: str = ""
    clinvar_accession: str = ""
    clinvar_classification: str = ""
    clinvar_review_status: str = ""
    clinvar_stars: int = 0
    clingen_status: str = "not_found"
    clingen_error: str = ""
    clingen_caid: str = ""
    erepo_registry_status: str = "not_found"
    erepo_assertion_uuid: str = ""
    erepo_assertion_method_version: str = ""
    historical_expert_panel_warning: str = ""
    classification: str = ""
    classification_verification: str = "unresolved"
    classification_source: str = ""
    objective_ps1_checks_pass: bool = False
    review_message: str = ""
    references: list[str] = Field(default_factory=list)
