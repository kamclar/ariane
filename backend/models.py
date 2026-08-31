# ============================================================
# ARIANE data models
# ============================================================
from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Optional, List, Dict, Any, Literal
import re

from backend.modules.variant_input import normalize_variant_input
from backend.gene_policy import active_genes


class VariantRequest(BaseModel):
    gene: str
    c_notation: str
    p_notation: str = ""
    dup_type: str = "Unknown"
    assembly: Optional[Literal["GRCh37", "GRCh38"]] = None
    submitted_notation: str = ""
    reference_transcript: str = ""
    normalization_source: str = ""
    consequence_status: str = ""
    normalization_provenance: Dict[str, str] = Field(default_factory=dict)
    protein_consequence_explanation: str = ""

    @model_validator(mode="before")
    @classmethod
    def normalize_hgvs_fields(cls, data):
        if not isinstance(data, dict):
            return data
        gene = str(data.get("gene", "")).strip().upper()
        raw_notation = str(data.get("c_notation", "")).strip()
        normalized = normalize_variant_input(
            gene,
            raw_notation,
            assembly=data.get("assembly"),
            p_notation=data.get("p_notation"),
        )
        return {
            **data,
            "gene": normalized.gene,
            "submitted_notation": raw_notation,
            "c_notation": normalized.c_notation,
            "p_notation": normalized.p_notation,
            "assembly": normalized.assembly or None,
            "reference_transcript": normalized.reference_transcript,
            "normalization_source": normalized.normalization_source,
            "consequence_status": normalized.consequence_status,
            "normalization_provenance": normalized.normalization_provenance or {},
            "protein_consequence_explanation": normalized.protein_consequence_explanation,
        }

    @field_validator("gene")
    @classmethod
    def validate_gene(cls, v):
        v = v.strip().upper()
        supported = set(active_genes())
        if v not in supported:
            raise ValueError(f"Gene must be one of: {', '.join(sorted(supported))}")
        return v

    @field_validator("c_notation")
    @classmethod
    def validate_c_notation(cls, v):
        v = v.strip()
        if not v.startswith("c."):
            raise ValueError("c. notation must start with 'c.', for example c.4185G>A")
        substitution = r"^c\.[-*]?\d+(?:[+-]\d+)?[ACGT]>[ACGT]$"
        equality = r"^c\.[-*]?\d+(?:[+-]\d+)?[ACGT]?=$"
        sequence_change = r"^c\.[0-9*+_?()\-]+(?:delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+)$"
        if not any(re.fullmatch(pattern, v, re.IGNORECASE) for pattern in (
            substitution,
            equality,
            sequence_change,
        )):
            raise ValueError(
                "Unrecognised c. notation. Use HGVS format, for example "
                "c.4185G>A, c.68_69delAG, or c.212+1G>T"
            )
        return v

    @field_validator("p_notation")
    @classmethod
    def validate_p_notation(cls, v):
        v = v.strip()
        if not v:
            raise ValueError(
                "p. notation is required, for example p.(Gln1395=); "
                "use p.? when the protein consequence is unknown"
            )
        if not v.startswith("p."):
            raise ValueError("p. notation must start with 'p.', for example p.(Gln1395=)")
        # The model receives the canonical result of the HGVS engine. Keep the
        # lexical guard broad enough for extension and multi-residue delins
        # forms; biological validation and supplied-vs-derived comparison have
        # already happened in normalize_variant_input().
        protein = r"^(?:p\.\?|p\.\([A-Za-z0-9_?=*]+\))$"
        if not re.fullmatch(protein, v):
            raise ValueError(
                "Unrecognised p. notation. Use HGVS format, for example "
                "p.(Cys61Gly), p.(Gln1395=), or p.(Glu23ValfsTer17)"
            )
        return v

    @field_validator("dup_type")
    @classmethod
    def validate_dup_type(cls, v):
        value = v.strip().title()
        if value not in ("Unknown", "Tandem"):
            raise ValueError("dup_type must be Unknown or Tandem")
        return value

class CriterionResult(BaseModel):
    name: str
    applies: bool
    strength: Optional[str] = None
    points: int = 0
    reason: str = ""
    source: str = ""
    decision_path: Optional[Dict[str, Any]] = None
    single_strong_likely_benign_eligible: bool = False
    single_strong_likely_benign_basis: str = ""
    independent_evidence_contribution_count: int = 0
    likelihood_ratio_contribution_count: int = 0
    clinical_evidence_types: List[str] = Field(default_factory=list)
    distinct_clinical_evidence_type_count: int = 0


class ExternalSubmitter(BaseModel):
    scv: str
    org: str
    classification: str
    date_eval: str = ""
    is_enigma_ep: bool = False
    review_status: str = ""
    curated_status: str = ""
    comment: str = ""


class ExternalComparison(BaseModel):
    clinvar_status: str = "not_found"
    clinvar_message: str = ""
    clinvar_error: str = ""
    clinvar_classification: str = ""
    clinvar_review_status: str = ""
    clinvar_review_stars: int = 0
    clinvar_n_submitters: int = 0
    clinvar_has_conflict: bool = False
    clinvar_submitters: List[ExternalSubmitter] = []
    clingen_status: str = "not_found"
    clingen_message: str = ""
    clingen_error: str = ""
    enigma_ep_class: str = ""
    enigma_ep_source: str = ""
    erepo_evidence_codes: List[str] = []


class AlphaMissenseResult(BaseModel):
    am_score: Optional[float] = None
    am_class: str = ""        # "likely_pathogenic" | "ambiguous" | "likely_benign"


class SpliceAIAudit(BaseModel):
    status: str = ""
    score: Optional[float] = None
    scoring_profile_id: str = ""
    scoring_profile_sha256: str = ""
    genome_assembly: str = ""
    distance: Optional[int] = None
    mask: Optional[int] = None
    annotation_subset: str = ""
    aggregation: str = ""
    source: str = ""
    transcript_policy: str = ""
    selected_transcript: str = ""
    reference_transcript_score: Optional[float] = None
    max_any_transcript_score: Optional[float] = None
    max_any_transcript: str = ""
    max_delta_field: str = ""
    delta_scores: Dict[str, float] = Field(default_factory=dict)
    reference_scores: Dict[str, float] = Field(default_factory=dict)
    alternate_scores: Dict[str, float] = Field(default_factory=dict)
    grch38: str = ""
    cache_key: str = ""
    reason: str = ""


class ClinicalLrAudit(BaseModel):
    application_status: str = "not_found"
    likelihood_ratio: Optional[float] = None
    candidate_likelihood_ratio: Optional[float] = None
    code: Optional[str] = None
    strength: Optional[str] = None
    source_bundle_ids: List[str] = Field(default_factory=list)
    source_bundle_count: int = 0
    independent_source_group_count: int = 0
    clinical_evidence_types: List[str] = Field(default_factory=list)
    distinct_clinical_evidence_type_count: int = 0
    likelihood_ratio_contribution_count: int = 0
    overlap_status: str = "not_assessed"
    double_counting_risk: bool = False
    automatic_combination_allowed: bool = False
    overlap_assessment_note: str = ""
    overlap_assessment_sources: List[str] = Field(default_factory=list)
    source_components: List[Dict[str, Any]] = Field(default_factory=list)
    reason: str = ""


class EvidenceInteractionWarning(BaseModel):
    status: Literal["info", "review_required", "deduplicated", "conflict"]
    mechanism: str
    criteria: List[str] = []
    retained: List[str] = []
    suppressed: List[str] = []
    reason: str
    source: str
    source_url: str
    review_required: bool = False


class ClinicalAnnotationPublication(BaseModel):
    pmid: str
    label: str
    url: str


class ClinicalAnnotation(BaseModel):
    category: Literal["reduced_penetrance"]
    label: str
    summary: str
    evidence: str
    source: str
    source_url: str
    source_row: int
    publications: List[ClinicalAnnotationPublication] = Field(default_factory=list)
    affects_classification: bool = False


class VusExplanation(BaseModel):
    category: str = ""
    tier: str = ""
    title: str = ""
    summary: str = ""
    what_to_check: str = ""
    review_priority: str = ""


class RnaReviewRecommendation(BaseModel):
    recommended: bool = False
    priority: str = "none"
    title: str = ""
    summary: str = ""
    reasons: List[str] = Field(default_factory=list)
    what_to_test: List[str] = Field(default_factory=list)
    potential_branches: List[str] = Field(default_factory=list)
    limitations: str = ""
    reference_source: str = ""
    source_url: str = ""
    is_evidence_criterion: bool = False
    manual_review_prefill: Dict[str, Any] = Field(default_factory=dict)


class ProteinPs1Candidate(BaseModel):
    key: str = ""
    reference_id: str = ""
    gene: str = ""
    transcript: str = ""
    c_notation: str = ""
    p_notation: str = ""
    classification: str = ""
    iarc_class: Optional[int] = None
    classification_basis: str = ""
    classification_source: str = ""
    reference_status: str = "review_required"
    status_reason: str = ""
    source_dataset: str = ""
    reference_splice_evidence_status: str = "not_assessed"
    reference_splice_sources_checked: List[str] = Field(default_factory=list)


class ProteinPs1ReviewRecommendation(BaseModel):
    display: bool = False
    recommended: bool = False
    priority: str = "none"
    title: str = ""
    summary: str = ""
    reasons: List[str] = Field(default_factory=list)
    what_to_check: List[str] = Field(default_factory=list)
    potential_branches: List[str] = Field(default_factory=list)
    limitations: str = ""
    reference_source: str = ""
    source_url: str = ""
    is_evidence_criterion: bool = False
    application_status: str = "not_applied"
    candidates: List[ProteinPs1Candidate] = Field(default_factory=list)
    splice_sources_checked: List[str] = Field(default_factory=list)
    vua_splice_evidence_status: str = "not_assessed"
    vua_spliceai_score: Optional[float] = None
    reference_spliceai_scores: Dict[str, Optional[float]] = Field(default_factory=dict)
    manual_review_prefill: Dict[str, Any] = Field(default_factory=dict)


class ClassificationResult(BaseModel):
    variant: str
    gene: str
    c_notation: str
    p_notation: str = ""
    variant_type: str = "unknown"
    bp7_rna_context: Dict[str, Any] = Field(default_factory=dict)
    reference_transcript: str = ""
    submitted_notation: str = ""
    normalization_source: str = ""
    consequence_status: str = ""
    normalization_provenance: Dict[str, str] = Field(default_factory=dict)
    protein_consequence_explanation: str = ""
    predicted_class: int
    predicted_label: str = ""
    total_points: int = 0
    criteria: List[CriterionResult] = []
    excluded_criteria: List[CriterionResult] = []
    not_applicable_criteria: List[CriterionResult] = []
    warnings: List[str] = []
    external: Optional[ExternalComparison] = None
    has_functional_evidence: bool = False
    classification_note: str = ""
    evidence_direction: str = "none"
    mixed_evidence: bool = False
    pathogenic_points: int = 0
    benign_points: int = 0
    narrative: str = ""
    alphamissense: Optional[AlphaMissenseResult] = None
    spliceai_audit: Optional[SpliceAIAudit] = None
    clinical_lr_audit: Optional[ClinicalLrAudit] = None
    population_frequency_audit: Dict[str, Any] = Field(default_factory=dict)
    evidence_interactions: List[EvidenceInteractionWarning] = []
    clinical_annotations: List[ClinicalAnnotation] = Field(default_factory=list)
    vus_explanation: Optional[VusExplanation] = None
    rna_review: Optional[RnaReviewRecommendation] = None
    splice_ps1_review: Optional[RnaReviewRecommendation] = None
    protein_ps1_review: Optional[ProteinPs1ReviewRecommendation] = None
    initiation_review: Optional[RnaReviewRecommendation] = None


class VariantNormalizationResponse(BaseModel):
    gene: str
    submitted_notation: str
    c_notation: str
    p_notation: str
    reference_transcript: str
    normalization_source: str
    consequence_status: str = ""
    normalization_provenance: Dict[str, str] = Field(default_factory=dict)
    protein_consequence_explanation: str = ""
    assembly: Optional[Literal["GRCh37", "GRCh38"]] = None


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
    spliceai_score: Optional[float] = None
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
    classification: str = ""
    classification_verification: str = "unresolved"
    classification_source: str = ""
    objective_ps1_checks_pass: bool = False
    review_message: str = ""
    references: List[str] = Field(default_factory=list)


class ManualCriterionInput(BaseModel):
    code: str
    enabled: bool = False
    evidence: Dict[str, Any] = Field(default_factory=dict)
    notes: str = ""
    references: List[str] = Field(default_factory=list)

    @model_validator(mode="before")
    @classmethod
    def reject_strength_override(cls, data):
        if (
            isinstance(data, dict)
            and data.get("override_strength") not in {None, ""}
        ):
            raise ValueError(
                "Manual strength overrides are not permitted; criterion strength "
                "is derived from the configured VCEP evidence thresholds"
            )
        return data

    @field_validator("code")
    @classmethod
    def validate_manual_code(cls, v):
        code = v.strip().upper()
        if code not in {
            "PS3",
            "PS4",
            "PM3",
            "PP1",
            "PP4",
            "BS2",
            "BS3",
            "BS4",
            "BP5",
            "PVS1_RNA",
            "BP7_RNA",
            "PVS1_INIT",
            "PS1_SPLICE",
            "PS1_PROTEIN",
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


class ManualEvidenceRequest(BaseModel):
    base_criteria: List[CriterionResult]
    manual_criteria: List[ManualCriterionInput]
    variant_context: Optional[ManualVariantContext] = None
    assessor: str
    assessed_at: str

    @field_validator("assessor", "assessed_at")
    @classmethod
    def require_audit_value(cls, v):
        if not v.strip():
            raise ValueError("Audit fields must not be empty")
        return v.strip()

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
                raise ValueError(
                    f"{item.code} requires evidence notes"
                )
            if not any(reference.strip() for reference in item.references):
                raise ValueError(
                    f"{item.code} requires at least one evidence reference"
                )
        return self


class ManualCriterionResult(BaseModel):
    code: str
    applies: bool
    suggested_strength: Optional[str] = None
    selected_strength: Optional[str] = None
    points: int = 0
    reason: str = ""
    threshold_note: str = ""
    overridden: bool = False
    notes: str = ""
    references: List[str] = []
    single_strong_likely_benign_eligible: bool = False
    single_strong_likely_benign_basis: str = ""
    independent_evidence_contribution_count: int = 0


class ManualEvidenceResult(BaseModel):
    predicted_class: int
    predicted_label: str
    total_points: int
    classification_note: str = ""
    manual_criteria: List[ManualCriterionResult]
    evidence_interactions: List[EvidenceInteractionWarning] = []
    assessor: str
    assessed_at: str


CLASS_LABELS = {
    5: "Pathogenic",
    4: "Likely Pathogenic",
    3: "VUS",
    2: "Likely Benign",
    1: "Benign",
}


class BatchItemResult(BaseModel):
    index: int
    status: str                          # "ok" or "error"
    variant: str
    error: Optional[str] = None
    result: Optional[ClassificationResult] = None


class BatchRequest(BaseModel):
    variants: List[VariantRequest]

    @field_validator("variants")
    @classmethod
    def validate_count(cls, v: list) -> list:
        if not v:
            raise ValueError("At least one variant required")
        if len(v) > 200:
            raise ValueError("Maximum 200 variants per batch")
        return v


class BatchResponse(BaseModel):
    total: int
    success_count: int
    error_count: int
    results: List[BatchItemResult]


class ClientValidationRequest(BaseModel):
    form: Literal["single", "batch"]
    input: Dict[str, Any]
    error: str = Field(min_length=1, max_length=500)

    @model_validator(mode="after")
    def limit_input_size(self):
        import json

        if len(json.dumps(self.input, ensure_ascii=False)) > 5000:
            raise ValueError("Client validation input is too large")
        return self
