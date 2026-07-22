# ============================================================
# ARIANE data models
# ============================================================
from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Optional, List, Dict, Any, Literal
import re

from backend.modules.hgvs import split_combined_hgvs


class VariantRequest(BaseModel):
    gene: str
    c_notation: str
    p_notation: Optional[str] = None
    dup_type: str = "Unknown"

    @model_validator(mode="before")
    @classmethod
    def normalize_hgvs_fields(cls, data):
        if not isinstance(data, dict):
            return data
        gene = str(data.get("gene", "")).strip().upper()
        raw_c = str(data.get("c_notation", "")).strip()
        transcript_match = re.match(r"^(NM_\d+\.\d+):(c\..+)$", raw_c, re.IGNORECASE)
        if transcript_match:
            transcript = transcript_match.group(1).upper()
            expected = {"BRCA1": "NM_007294.4", "BRCA2": "NM_000059.4"}.get(gene)
            if expected and transcript != expected:
                raise ValueError(
                    f"Transcript {transcript} does not match {gene}; expected {expected}"
                )
            raw_c = transcript_match.group(2)
        c_notation, p_notation = split_combined_hgvs(
            raw_c,
            data.get("p_notation"),
        )
        return {**data, "c_notation": c_notation, "p_notation": p_notation or None}

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
        if v is None:
            return None
        v = v.strip()
        if not v:
            return None
        if not v.startswith("p."):
            raise ValueError("p. notation must start with 'p.', for example p.(Gln1395=)")
        protein = (
            r"^p\.\((?:\?|[A-Z][a-z]{2}\d+"
            r"(?:_[A-Z][a-z]{2}\d+)?"
            r"(?:=|\?|Ter|[A-Z][a-z]{2}|del|dup|fs(?:Ter)?\d*|"
            r"[A-Z][a-z]{2}fs(?:Ter)?\d*|delins[A-Z][a-z]{2}|ins[A-Z][a-z]{2})"
            r")\)$"
        )
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


class ClassificationResult(BaseModel):
    variant: str
    gene: str
    c_notation: str
    p_notation: str = ""
    predicted_class: int
    predicted_label: str = ""
    total_points: int = 0
    criteria: List[CriterionResult] = []
    warnings: List[str] = []
    external: Optional[ExternalComparison] = None
    has_functional_evidence: bool = False
    classification_note: str = ""
    narrative: str = ""
    alphamissense: Optional[AlphaMissenseResult] = None
    vus_explanation: Optional[VusExplanation] = None
    rna_review: Optional[RnaReviewRecommendation] = None
    splice_ps1_review: Optional[RnaReviewRecommendation] = None
    initiation_review: Optional[RnaReviewRecommendation] = None


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
            "BS2",
            "BS4",
            "PVS1_RNA",
            "BP7_RNA",
            "PVS1_INIT",
            "PS1_SPLICE",
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
