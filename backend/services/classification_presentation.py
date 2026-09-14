"""Project structured DAG output into the stable public classification model."""

from __future__ import annotations

from backend.contracts import (
    AlphaMissenseResult,
    CLASS_LABELS,
    ClassificationResult,
    ClinicalAnnotation,
    ClinicalLrAudit,
    CriterionResult,
    EvidenceInteractionWarning,
    ExternalComparison,
    ExternalSubmitter,
    ProteinPs1ReviewRecommendation,
    RnaReviewRecommendation,
    SpliceAIAudit,
    VusExplanation,
)
from backend.lookups.clinvar import clinvar_review_stars
from backend.criteria.bp7_rna import evaluate_bp7_rna_variant_context
from backend.presentation.criterion_order import sorted_criterion_items
from backend.reference_data.enigma_rules import clinical_annotations_for_variant
from backend.reference_data.erepo_vcep import lookup_erepo_vcep_assertion
from backend.presentation.external import external_comparison
from backend.presentation.narrative import generate_narrative
from backend.presentation.vus_explanation import explain_vus
from backend.services.evidence_orchestration import (
    OrchestratedEvidence,
    refreshed_external_warnings,
)


def _criterion_models(values: dict, *, applies: bool) -> list[CriterionResult]:
    return [
        CriterionResult(
            name=name,
            applies=applies,
            strength=criterion.get("strength"),
            points=criterion.get("points", 0) if applies else 0,
            reason=criterion.get("reason", ""),
            source=criterion.get("source", ""),
            decision_path=criterion.get("decision_path"),
            table9_audit=criterion.get("table9_audit"),
            single_strong_likely_benign_eligible=criterion.get(
                "single_strong_likely_benign_eligible", False
            ),
            single_strong_likely_benign_basis=criterion.get(
                "single_strong_likely_benign_basis", ""
            ),
            independent_evidence_contribution_count=criterion.get(
                "independent_evidence_contribution_count", 0
            ),
            likelihood_ratio_contribution_count=criterion.get(
                "likelihood_ratio_contribution_count", 0
            ),
            clinical_evidence_types=criterion.get("clinical_evidence_types", []),
            distinct_clinical_evidence_type_count=criterion.get(
                "distinct_clinical_evidence_type_count", 0
            ),
        )
        for name, criterion in sorted_criterion_items(values)
    ]


def _external_status_message(source: str, value: dict) -> str:
    status = str(value.get("status") or "unavailable")
    if status == "ok":
        return f"{source} returned an exact record for the assessed variant."
    if status == "not_found":
        return f"No exact record was found in {source} for the assessed variant."
    if status == "ambiguous":
        return (
            f"{source} returned more than one possible record for the assessed "
            "variant. None was selected. Candidate IDs are retained in the audit data."
        )
    return f"{source} is unavailable for this request (status: {status})."


def _historical_expert_panel_warning(
    *,
    local_erepo: dict,
    local_erepo_record: dict,
    clinvar_stars: int,
    enigma_submission: dict,
) -> str:
    if local_erepo.get("status") in {
        "historical_vcep_assertion",
        "unversioned_vcep_assertion",
    }:
        version = str(local_erepo_record.get("assertion_method_version") or "not recorded")
        return (
            "An ENIGMA expert-panel assertion is present, but its recorded "
            f"specification version is {version}, not the active v1.2. It is "
            "shown for context and is not used as a current VCEP assertion."
        )

    if (
        clinvar_stars != 3
        or not enigma_submission
        or local_erepo.get("status") == "current_vcep_assertion"
    ):
        return ""

    comment = str(enigma_submission.get("comment") or "")
    historical_multifactorial = any(
        phrase in comment.lower()
        for phrase in (
            "multifactorial likelihood",
            "posterior probability",
            "iarc class",
        )
    )
    date_evaluated = str(enigma_submission.get("date_eval") or "").strip()
    if historical_multifactorial:
        date_text = f" evaluated on {date_evaluated}" if date_evaluated else ""
        return (
            f"The ClinVar ENIGMA expert-panel assertion{date_text} records an "
            "IARC class based on multifactorial posterior probability. It is not "
            "a criterion-by-criterion assertion under the active VCEP v1.2 and "
            "is shown as historical context."
        )

    return (
        "ClinVar contains a three-star ENIGMA expert-panel assertion, but "
        "the active local ClinGen ERepo snapshot has no matching current "
        "v1.2 assertion. It is shown for context only."
    )


def _evidence_display_flags(criteria: list[CriterionResult]) -> tuple[bool, bool]:
    """Keep protein-function and RNA evidence labels distinct in the UI."""
    has_table9_functional_evidence = any(
        criterion.applies
        and criterion.name in {"PS3", "BS3"}
        and criterion.table9_audit is not None
        for criterion in criteria
    )
    has_curated_rna_evidence = any(
        criterion.applies and criterion.name == "PVS1_RNA"
        for criterion in criteria
    )
    return has_table9_functional_evidence, has_curated_rna_evidence


def _external_model_from_values(
    *,
    gene: str,
    c_notation: str,
    predicted_class: int,
    clinvar: dict,
    clingen: dict,
) -> ExternalComparison:
    local_erepo = lookup_erepo_vcep_assertion(gene, c_notation)
    local_erepo_record = dict(local_erepo.get("record") or {})
    comparison = external_comparison(
        gene,
        c_notation,
        predicted_class,
        clinvar,
        clingen,
        local_erepo,
    )
    submitters = [
        ExternalSubmitter(
            scv=item.get("scv", ""),
            org=item.get("org", ""),
            classification=item.get("class") or "",
            date_eval=item.get("date_eval", ""),
            is_enigma_ep=item.get("is_enigma_ep", False),
            review_status=item.get("review", ""),
            curated_status=(
                "ENIGMA expert-panel ClinVar assertion"
                if item.get("is_enigma_ep", False)
                and "expert panel" in str(item.get("review") or "").lower()
                else "ENIGMA-labelled ClinVar submission"
                if item.get("is_enigma_ep", False)
                else ""
            ),
            comment=item.get("comment", "")[:200],
        )
        for item in clinvar.get("submissions", [])
        if clinvar.get("status") == "ok"
    ]
    aggregate = clinvar.get("aggregate", {})
    review_status = aggregate.get("review_status", "")
    stars = clinvar_review_stars(review_status)
    enigma_submission = dict(clinvar.get("enigma_submission") or {})
    warning = _historical_expert_panel_warning(
        local_erepo=local_erepo,
        local_erepo_record=local_erepo_record,
        clinvar_stars=stars,
        enigma_submission=enigma_submission,
    )
    difference_message = ""
    if comparison.get("match") is False:
        prefix = (
            "Historical ENIGMA expert-panel classification"
            if warning
            else "ENIGMA expert-panel classification"
        )
        difference_message = (
            f"{prefix}: {comparison.get('enigma_class', '')}; current ARIANE "
            "v1.2 automated result: "
            f"{CLASS_LABELS.get(predicted_class, '')}."
        )
    return ExternalComparison(
        clinvar_status=str(clinvar.get("status") or "unavailable"),
        clinvar_message=_external_status_message("ClinVar", clinvar),
        clinvar_error=str(clinvar.get("error") or "")[:500],
        clinvar_candidate_ids=[
            str(value) for value in clinvar.get("candidate_ids", []) if value
        ],
        clinvar_classification=aggregate.get("classification", ""),
        clinvar_review_status=review_status,
        clinvar_review_stars=stars,
        clinvar_n_submitters=aggregate.get("n_submitters", 0),
        clinvar_has_conflict=clinvar.get("has_conflict", False),
        clinvar_submitters=submitters,
        clingen_status=str(clingen.get("status") or "unavailable"),
        clingen_message=_external_status_message("ClinGen ERepo", clingen),
        clingen_error=str(clingen.get("error") or "")[:500],
        enigma_ep_class=comparison.get("enigma_class", ""),
        enigma_ep_source=comparison.get("enigma_source", ""),
        erepo_evidence_codes=[
            str(item.get("code") or "")
            for item in clingen.get("evidence_codes", [])
            if item.get("code")
        ],
        erepo_guideline_versions=[
            str(value) for value in clingen.get("guideline_versions", []) if value
        ],
        erepo_cspec_ids=[
            str(value) for value in clingen.get("cspec_ids", []) if value
        ],
        erepo_assertion_id=str(clingen.get("assertion_id") or ""),
        erepo_registry_status=str(local_erepo.get("status") or "not_found"),
        erepo_registry_assertion_uuid=str(local_erepo_record.get("uuid") or ""),
        erepo_registry_method_version=str(
            local_erepo_record.get("assertion_method_version") or ""
        ),
        historical_expert_panel_warning=warning,
        expert_panel_difference_message=difference_message,
    )


def _external_model(evidence: OrchestratedEvidence) -> ExternalComparison:
    """Compatibility wrapper for a fully orchestrated classification."""
    return _external_model_from_values(
        gene=evidence.variant.gene,
        c_notation=evidence.variant.c_notation,
        predicted_class=evidence.result["predicted_class"],
        clinvar=dict(evidence.clinvar),
        clingen=dict(evidence.clingen),
    )


class ClassificationPresentationService:
    """Build the API model without acquiring evidence or changing classification."""

    def refresh_external(
        self,
        result: ClassificationResult,
        *,
        clinvar: dict,
        clingen: dict,
    ) -> ClassificationResult:
        """Attach current read-only comparisons to a cached Module 1 result."""
        external = _external_model_from_values(
            gene=result.gene,
            c_notation=result.c_notation,
            predicted_class=result.predicted_class,
            clinvar=clinvar,
            clingen=clingen,
        )
        warnings = refreshed_external_warnings(
            result.warnings,
            clinvar=clinvar,
            clingen=clingen,
        )
        return result.model_copy(update={"external": external, "warnings": warnings})

    def build(self, evidence: OrchestratedEvidence) -> ClassificationResult:
        result = evidence.result
        artifacts = evidence.artifacts
        normalized = evidence.normalized_input
        variant = evidence.variant
        criteria = _criterion_models(result["criteria"], applies=True)
        excluded = _criterion_models(
            result.get("excluded_criteria", {}),
            applies=False,
        )
        not_applicable = _criterion_models(
            result.get("not_applicable_criteria", {}),
            applies=False,
        )
        spliceai_score = artifacts.get("spliceai_score")
        bayesdel_score = artifacts.get("bayesdel_score")
        alphamissense = artifacts.get("alphamissense")
        splice_status = artifacts.get("spliceai_status") or {}
        gnomad_data = artifacts.get("gnomad_data")
        vus_explanation = explain_vus(result)
        narrative = generate_narrative(
            gene=variant.gene,
            c_notation=variant.c_notation,
            p_notation=variant.p_notation,
            variant_type=evidence.variant_type,
            result=result,
            spliceai_score=spliceai_score,
            bayesdel_score=bayesdel_score,
            alphamissense=alphamissense,
        )
        bp7_rna_context = evaluate_bp7_rna_variant_context(
            {
                "gene": variant.gene,
                "c_notation": variant.c_notation,
                "p_notation": variant.p_notation,
            },
            [criterion.model_dump() for criterion in criteria],
        )
        has_table9_functional_evidence, has_curated_rna_evidence = (
            _evidence_display_flags(criteria)
        )
        return ClassificationResult(
            variant=result["variant"],
            gene=variant.gene,
            c_notation=variant.c_notation,
            p_notation=variant.p_notation,
            variant_type=evidence.variant_type,
            bp7_rna_context=bp7_rna_context,
            reference_transcript=normalized.reference_transcript,
            normalization_source=normalized.normalization_source,
            consequence_status=normalized.consequence_status,
            normalization_provenance=normalized.normalization_provenance or {},
            protein_consequence_explanation=normalized.protein_consequence_explanation,
            predicted_class=result["predicted_class"],
            predicted_label=CLASS_LABELS.get(result["predicted_class"], ""),
            total_points=result["total_points"],
            criteria=criteria,
            excluded_criteria=excluded,
            not_applicable_criteria=not_applicable,
            warnings=result["warnings"],
            external=_external_model(evidence),
            has_functional_evidence=result.get("has_functional_evidence", False),
            has_table9_functional_evidence=has_table9_functional_evidence,
            has_curated_rna_evidence=has_curated_rna_evidence,
            classification_note=result.get("classification_note", ""),
            evidence_direction=result.get("evidence_direction", "none"),
            mixed_evidence=result.get("mixed_evidence", False),
            pathogenic_points=result.get("pathogenic_points", 0),
            benign_points=result.get("benign_points", 0),
            narrative=narrative,
            alphamissense=(
                AlphaMissenseResult(
                    am_score=alphamissense.get("am_score"),
                    am_class=alphamissense.get("am_class", ""),
                )
                if alphamissense
                else None
            ),
            vus_explanation=(
                VusExplanation(**vus_explanation) if vus_explanation else None
            ),
            rna_review=(
                RnaReviewRecommendation(**result["rna_review"])
                if result.get("rna_review")
                else None
            ),
            splice_ps1_review=(
                RnaReviewRecommendation(**result["splice_ps1_review"])
                if result.get("splice_ps1_review")
                else None
            ),
            protein_ps1_review=(
                ProteinPs1ReviewRecommendation(**result["protein_ps1_review"])
                if result.get("protein_ps1_review")
                else None
            ),
            initiation_review=(
                RnaReviewRecommendation(**result["initiation_review"])
                if result.get("initiation_review")
                else None
            ),
            spliceai_audit=(
                SpliceAIAudit(
                    **{
                        field: splice_status.get(field)
                        for field in SpliceAIAudit.model_fields
                        if splice_status.get(field) is not None
                    }
                )
                if splice_status
                else None
            ),
            clinical_lr_audit=(
                ClinicalLrAudit(**artifacts["clinical_lr_result"])
                if artifacts.get("clinical_lr_result", {}).get("application_status")
                != "not_found"
                else None
            ),
            population_frequency_audit=(
                gnomad_data.get("population_frequency_audit", {})
                if gnomad_data
                else {}
            ),
            evidence_interactions=[
                EvidenceInteractionWarning(**warning)
                for warning in result.get("evidence_interactions", [])
            ],
            clinical_annotations=[
                ClinicalAnnotation(**annotation)
                for annotation in clinical_annotations_for_variant(
                    variant.gene,
                    variant.c_notation,
                )
            ],
        )
