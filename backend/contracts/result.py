"""Classification result HTTP schemas."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field


class CriterionResult(BaseModel):
    name: str
    applies: bool
    strength: str | None = None
    points: int = 0
    reason: str = ""
    source: str = ""
    decision_path: dict[str, Any] | None = None
    table9_audit: dict[str, Any] | None = None
    single_strong_likely_benign_eligible: bool = False
    single_strong_likely_benign_basis: str = ""
    independent_evidence_contribution_count: int = 0
    likelihood_ratio_contribution_count: int = 0
    clinical_evidence_types: list[str] = Field(default_factory=list)
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
    clinvar_candidate_ids: list[str] = Field(default_factory=list)
    clinvar_classification: str = ""
    clinvar_review_status: str = ""
    clinvar_review_stars: int = 0
    clinvar_n_submitters: int = 0
    clinvar_has_conflict: bool = False
    clinvar_submitters: list[ExternalSubmitter] = Field(default_factory=list)
    clingen_status: str = "not_found"
    clingen_message: str = ""
    clingen_error: str = ""
    enigma_ep_class: str = ""
    enigma_ep_source: str = ""
    erepo_evidence_codes: list[str] = Field(default_factory=list)
    erepo_guideline_versions: list[str] = Field(default_factory=list)
    erepo_cspec_ids: list[str] = Field(default_factory=list)
    erepo_assertion_id: str = ""
    erepo_registry_status: str = "not_found"
    erepo_registry_assertion_uuid: str = ""
    erepo_registry_method_version: str = ""
    historical_expert_panel_warning: str = ""


class AlphaMissenseResult(BaseModel):
    am_score: float | None = None
    am_class: str = ""


class SpliceAIAudit(BaseModel):
    status: str = ""
    score: float | None = None
    required_for_classification: bool = False
    retryable: bool = False
    reference_lookup_required: bool = False
    reference_lookup_complete: bool = True
    reference_variant_statuses: dict[str, Any] = Field(default_factory=dict)
    scoring_profile_id: str = ""
    scoring_profile_sha256: str = ""
    genome_assembly: str = ""
    distance: int | None = None
    mask: int | None = None
    annotation_subset: str = ""
    aggregation: str = ""
    source: str = ""
    transcript_policy: str = ""
    selected_transcript: str = ""
    reference_transcript_score: float | None = None
    max_any_transcript_score: float | None = None
    max_any_transcript: str = ""
    max_delta_field: str = ""
    delta_scores: dict[str, float] = Field(default_factory=dict)
    reference_scores: dict[str, float] = Field(default_factory=dict)
    alternate_scores: dict[str, float] = Field(default_factory=dict)
    grch38: str = ""
    cache_key: str = ""
    reason: str = ""


class ClinicalLrThresholdComparison(BaseModel):
    status: str = "not_available"
    source_label: str = ""
    source_code: str | None = None
    source_strength: str | None = None
    vcep_policy: str = ""
    vcep_label: str = ""
    vcep_code: str | None = None
    vcep_strength: str | None = None
    combined_lr: float | None = None
    threshold_operator: str = ""
    threshold_value: float | None = None
    threshold_rule: str = ""
    reason: str = ""


class ClinicalLrAudit(BaseModel):
    application_status: str = "not_found"
    likelihood_ratio: float | None = None
    likelihood_ratio_status: Literal[
        "not_found", "available", "source_reported_zero", "unavailable_conflict"
    ] = "not_found"
    candidate_likelihood_ratio: float | None = None
    code: str | None = None
    strength: str | None = None
    source_bundle_ids: list[str] = Field(default_factory=list)
    source_bundle_count: int = 0
    independent_source_group_count: int = 0
    clinical_evidence_types: list[str] = Field(default_factory=list)
    distinct_clinical_evidence_type_count: int = 0
    likelihood_ratio_contribution_count: int = 0
    overlap_status: str = "not_assessed"
    double_counting_risk: bool = False
    source_reported_overlap_caveat: bool = False
    automatic_combination_allowed: bool = False
    data_release: str = ""
    overlap_assessment_note: str = ""
    overlap_assessment_sources: list[str] = Field(default_factory=list)
    source_components: list[dict[str, Any]] = Field(default_factory=list)
    threshold_comparison: ClinicalLrThresholdComparison = Field(
        default_factory=ClinicalLrThresholdComparison
    )
    reason: str = ""


class EvidenceInteractionWarning(BaseModel):
    status: Literal["info", "review_required", "deduplicated", "conflict"]
    mechanism: str
    criteria: list[str] = Field(default_factory=list)
    retained: list[str] = Field(default_factory=list)
    suppressed: list[str] = Field(default_factory=list)
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
    publications: list[ClinicalAnnotationPublication] = Field(default_factory=list)
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
    reasons: list[str] = Field(default_factory=list)
    what_to_test: list[str] = Field(default_factory=list)
    potential_branches: list[str] = Field(default_factory=list)
    limitations: str = ""
    reference_source: str = ""
    source_url: str = ""
    is_evidence_criterion: bool = False
    manual_review_prefill: dict[str, Any] = Field(default_factory=dict)


class ProteinPs1Candidate(BaseModel):
    key: str = ""
    reference_id: str = ""
    gene: str = ""
    transcript: str = ""
    c_notation: str = ""
    p_notation: str = ""
    classification: str = ""
    iarc_class: int | None = None
    classification_basis: str = ""
    classification_source: str = ""
    reference_status: str = "review_required"
    status_reason: str = ""
    source_dataset: str = ""
    reference_splice_evidence_status: str = "not_assessed"
    reference_splice_sources_checked: list[str] = Field(default_factory=list)


class ProteinPs1ReviewRecommendation(BaseModel):
    display: bool = False
    recommended: bool = False
    priority: str = "none"
    title: str = ""
    summary: str = ""
    reasons: list[str] = Field(default_factory=list)
    what_to_check: list[str] = Field(default_factory=list)
    potential_branches: list[str] = Field(default_factory=list)
    limitations: str = ""
    reference_source: str = ""
    source_url: str = ""
    is_evidence_criterion: bool = False
    application_status: str = "not_applied"
    candidates: list[ProteinPs1Candidate] = Field(default_factory=list)
    splice_sources_checked: list[str] = Field(default_factory=list)
    vua_splice_evidence_status: str = "not_assessed"
    vua_spliceai_score: float | None = None
    reference_spliceai_scores: dict[str, float | None] = Field(default_factory=dict)
    manual_review_prefill: dict[str, Any] = Field(default_factory=dict)


class ClassificationResult(BaseModel):
    variant: str
    gene: str
    c_notation: str
    p_notation: str = ""
    variant_type: str = "unknown"
    bp7_rna_context: dict[str, Any] = Field(default_factory=dict)
    reference_transcript: str = ""
    submitted_notation: str = ""
    normalization_source: str = ""
    consequence_status: str = ""
    normalization_provenance: dict[str, str] = Field(default_factory=dict)
    protein_consequence_explanation: str = ""
    predicted_class: int
    predicted_label: str = ""
    total_points: int = 0
    criteria: list[CriterionResult] = Field(default_factory=list)
    excluded_criteria: list[CriterionResult] = Field(default_factory=list)
    not_applicable_criteria: list[CriterionResult] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    external: ExternalComparison | None = None
    has_functional_evidence: bool = False
    classification_note: str = ""
    evidence_direction: str = "none"
    mixed_evidence: bool = False
    pathogenic_points: int = 0
    benign_points: int = 0
    narrative: str = ""
    alphamissense: AlphaMissenseResult | None = None
    spliceai_audit: SpliceAIAudit | None = None
    clinical_lr_audit: ClinicalLrAudit | None = None
    population_frequency_audit: dict[str, Any] = Field(default_factory=dict)
    evidence_interactions: list[EvidenceInteractionWarning] = Field(default_factory=list)
    clinical_annotations: list[ClinicalAnnotation] = Field(default_factory=list)
    vus_explanation: VusExplanation | None = None
    rna_review: RnaReviewRecommendation | None = None
    splice_ps1_review: RnaReviewRecommendation | None = None
    protein_ps1_review: ProteinPs1ReviewRecommendation | None = None
    initiation_review: RnaReviewRecommendation | None = None
