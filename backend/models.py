# ============================================================
# ARIANE data models
# ============================================================
from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Optional, List, Dict, Any, Literal
import re

from backend.modules.variant_input import normalize_variant_input


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
        if v not in ("BRCA1", "BRCA2"):
            raise ValueError("Gene must be BRCA1 or BRCA2")
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
    clinvar_classification: str = ""
    clinvar_review_status: str = ""
    clinvar_review_stars: int = 0
    clinvar_n_submitters: int = 0
    clinvar_has_conflict: bool = False
    clinvar_submitters: List[ExternalSubmitter] = []
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


class ClassificationResult(BaseModel):
    variant: str
    gene: str
    c_notation: str
    p_notation: str = ""
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
    population_frequency_audit: Dict[str, Any] = Field(default_factory=dict)
    evidence_interactions: List[EvidenceInteractionWarning] = []
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


class ManualCriterionInput(BaseModel):
    code: str
    enabled: bool = False
    evidence: Dict[str, Any] = Field(default_factory=dict)
    override_strength: Optional[str] = None
    notes: str = ""
    references: List[str] = Field(default_factory=list)

    @field_validator("code")
    @classmethod
    def validate_manual_code(cls, v):
        code = v.strip().upper()
        if code not in {
            "PS4",
            "PM3",
            "PP1",
            "PP4",
            "BS2",
            "BS4",
            "PVS1_RNA",
            "BP7_RNA",
            "PVS1_INIT",
            "PS1_SPLICE",
            "PS1_PROTEIN",
        }:
            raise ValueError("Unsupported manually reviewed criterion")
        return code


class ManualEvidenceRequest(BaseModel):
    base_criteria: List[CriterionResult]
    manual_criteria: List[ManualCriterionInput]
    assessor: str
    assessed_at: str

    @field_validator("assessor", "assessed_at")
    @classmethod
    def require_audit_value(cls, v):
        if not v.strip():
            raise ValueError("Audit fields must not be empty")
        return v.strip()


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
