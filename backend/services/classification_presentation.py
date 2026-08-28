"""Project structured DAG output into the stable public classification model."""

from __future__ import annotations

from backend.models import (
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
from backend.modules.criterion_order import sorted_criterion_items
from backend.services.evidence_orchestration import OrchestratedEvidence


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
        return f"{source} returned an exact record."
    if status == "not_found":
        return f"No exact record was found in {source}."
    if status == "ambiguous":
        candidates = ", ".join(str(item) for item in value.get("candidate_ids", []))
        suffix = f" Candidate IDs: {candidates}." if candidates else ""
        return f"{source} returned an ambiguous match; no record was selected.{suffix}"
    return f"{source} is unavailable for this request (status: {status})."


def _external_model(evidence: OrchestratedEvidence) -> ExternalComparison:
    from backend.lookups.clinvar import clinvar_review_stars
    from backend.modules.external import external_comparison

    clinvar = evidence.clinvar
    clingen = evidence.clingen
    result = evidence.result
    variant = evidence.variant
    comparison = external_comparison(
        variant.gene,
        variant.c_notation,
        result["predicted_class"],
        clinvar,
        clingen,
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
                "ClinGen/ENIGMA curated submitter"
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
    return ExternalComparison(
        clinvar_status=str(clinvar.get("status") or "unavailable"),
        clinvar_message=_external_status_message("ClinVar", clinvar),
        clinvar_error=str(clinvar.get("error") or "")[:500],
        clinvar_classification=aggregate.get("classification", ""),
        clinvar_review_status=review_status,
        clinvar_review_stars=clinvar_review_stars(review_status),
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
    )


class ClassificationPresentationService:
    """Build the API model without acquiring evidence or changing classification."""

    def build(self, evidence: OrchestratedEvidence) -> ClassificationResult:
        from backend.modules.bp7_rna import evaluate_bp7_rna_variant_context
        from backend.modules.enigma_rules import clinical_annotations_for_variant
        from backend.modules.narrative import generate_narrative
        from backend.modules.vus_explanation import explain_vus

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
