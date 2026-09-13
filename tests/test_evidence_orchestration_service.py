"""Boundaries and contracts of the classification application service."""

import asyncio
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest

from backend.services.variant_classification_service import execute_variant_classification
from backend.services.evidence_orchestration import (
    ClassificationCommand,
    EvidenceOrchestrationService,
    RequiredEvidenceUnavailableError,
    VariantPreparationError,
)
from backend.variant_processing.variant_input import NormalizedVariantInput


def _complete_population_result():
    policy = {
        "policy_id": "enigma_brca_v1_2",
        "frequency_criteria": {
            "pm2": {
                "required_absence_dataset_runtime_keys": [
                    "v2_1_non_cancer",
                    "v3_1_non_cancer",
                ]
            }
        },
    }
    return {
        "status": "pm2_coverage_method_unresolved",
        "policy_id": "enigma_brca_v1_2",
        "classification_policy": policy,
        "datasets": {
            "v2_1_non_cancer": {"status": "absent"},
            "v3_1_non_cancer": {"status": "absent"},
        },
    }


def _complete_artifacts(**overrides):
    result = {
        "gnomad_data": _complete_population_result(),
        "exon_cnv_result": None,
        "bayesdel_score": 0.438,
        "bayesdel_status": {"status": "ok", "reason": "test score"},
    }
    result.update(overrides)
    return result


def test_orchestrator_prepares_policy_bound_normalized_variant():
    service = EvidenceOrchestrationService()
    normalized, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
    )

    assert normalized.reference_transcript == "NM_007294.4"
    assert variant.gene == "BRCA1"
    assert variant.c_notation == "c.303T>G"
    assert variant.p_notation == "p.(Tyr101Ter)"
    assert variant.variant_type == "nonsense"


def test_orchestrator_routes_dna_delins_with_missense_consequence_to_figure1a():
    service = EvidenceOrchestrationService()
    normalized, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.131_132delinsCT")
    )

    assert normalized.p_notation == "p.(Cys44Ser)"
    assert variant.variant_type == "missense"


def test_orchestrator_rejects_invalid_input_before_provider_planning():
    service = EvidenceOrchestrationService()
    with pytest.raises(VariantPreparationError):
        service.prepare(ClassificationCommand("BRCA1", "c.181A>C"))


def test_orchestrator_rejects_unresolved_coding_delins_with_explicit_reason():
    service = EvidenceOrchestrationService()
    unresolved = NormalizedVariantInput(
        gene="BRCA1",
        submitted_notation="c.922_923delinsGC",
        c_notation="c.922_923delinsGC",
        p_notation="p.?",
        reference_transcript="NM_007294.4",
        normalization_source="test source",
        consequence_status="unresolved",
    )

    with patch(
        "backend.services.evidence_orchestration.normalize_variant_input",
        return_value=unresolved,
    ), pytest.raises(VariantPreparationError) as error:
        service.prepare(ClassificationCommand("BRCA1", "c.922_923delinsGC"))

    assert "protein consequence" in str(error.value)
    assert "ENIGMA PTC or Figure 1A branch" in str(error.value)
    assert "No classification was returned" in str(error.value)


def test_required_spliceai_timeout_stops_figure1a_classification():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.5366C>T", "p.(Ala1789Val)")
    )

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            _complete_artifacts(**{
                "spliceai_score": None,
                "spliceai_status": {
                    "status": "api_error",
                    "reason": "TimeoutError: The read operation timed out",
                    "retryable": True,
                },
            }),
        )

    assert error.value.code == "spliceai_temporarily_unavailable"
    assert error.value.retryable is True


def test_required_spliceai_timeout_for_dna_delins_reports_source_and_no_result():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.131_132delinsCT")
    )

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            _complete_artifacts(**{
                "spliceai_score": None,
                "spliceai_status": {
                    "status": "api_error",
                    "reason": "TimeoutError: The read operation timed out",
                    "retryable": True,
                },
            }),
        )

    assert error.value.source == "SpliceAI"
    assert error.value.status == "api_error"
    assert error.value.code == "spliceai_temporarily_unavailable"
    assert error.value.retryable is True
    assert "TimeoutError: The read operation timed out" in str(error.value)
    assert "No classification was returned" in str(error.value)


def test_missing_spliceai_does_not_block_non_figure1a_ptc_classification():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
    )

    service._require_complete_classification_evidence(
        variant,
        _complete_artifacts(**{
            "spliceai_score": None,
            "spliceai_status": {
                "status": "not_applicable",
                "required_for_classification": False,
            },
        }),
    )


def test_incomplete_ps1_reference_spliceai_stops_classification():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.5217T>A", "p.(Asp1739Glu)")
    )

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            _complete_artifacts(**{
                "spliceai_score": 0.01,
                "spliceai_status": {
                    "status": "ok",
                    "score": 0.01,
                    "reference_lookup_required": True,
                    "reference_lookup_complete": False,
                    "reference_variant_statuses": {
                        "c.5217T>G": {
                            "status": "api_error",
                            "score": None,
                            "reason": "TimeoutError: The read operation timed out",
                            "retryable": True,
                        },
                    },
                },
            }),
        )

    assert error.value.code == "spliceai_temporarily_unavailable"
    assert error.value.retryable is True


def test_required_bayesdel_timeout_stops_in_domain_figure1a_classification():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.5366C>T", "p.(Ala1789Val)")
    )

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            _complete_artifacts(
                spliceai_score=0.03,
                spliceai_status={"status": "ok", "score": 0.03},
                bayesdel_score=None,
                bayesdel_status={
                    "status": "api_error",
                    "reason": "TimeoutError: The read operation timed out",
                },
            ),
        )

    assert error.value.source == "BayesDel_noAF"
    assert error.value.code == "bayesdel_temporarily_unavailable"
    assert error.value.retryable is True


def test_bayesdel_is_not_required_after_spliceai_pp3_branch_is_resolved():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.5366C>T", "p.(Ala1789Val)")
    )

    service._require_complete_classification_evidence(
        variant,
        _complete_artifacts(
            spliceai_score=0.20,
            spliceai_status={"status": "ok", "score": 0.20},
            bayesdel_score=None,
            bayesdel_status={"status": "api_error", "reason": "timeout"},
        ),
    )


def test_bayesdel_is_not_required_outside_functional_domain():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.3247A>C", "p.(Met1083Leu)")
    )

    service._require_complete_classification_evidence(
        variant,
        _complete_artifacts(
            spliceai_score=0.03,
            spliceai_status={"status": "ok", "score": 0.03},
            bayesdel_score=None,
            bayesdel_status={"status": "api_error", "reason": "timeout"},
        ),
    )


def test_incomplete_required_gnomad_dataset_stops_ptc_classification():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
    )
    population = _complete_population_result()
    population["datasets"]["v3_1_non_cancer"] = {
        "status": "no_coordinates",
        "errors": ["GRCh38 coordinates were unavailable"],
    }

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            _complete_artifacts(
                gnomad_data=population,
                spliceai_score=None,
                spliceai_status={"status": "not_applicable"},
            ),
        )

    assert error.value.source == "gnomAD"
    assert error.value.status == "no_coordinates"
    assert error.value.code == "population_evidence_unavailable"
    assert "v3_1_non_cancer" in str(error.value)


def test_one_found_gnomad_dataset_does_not_hide_failure_of_the_other():
    service = EvidenceOrchestrationService()
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
    )
    population = _complete_population_result()
    population["status"] = "found"
    population["datasets"]["v2_1_non_cancer"] = {"status": "found"}
    population["datasets"]["v3_1_non_cancer"] = {"status": "cache_missing"}

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            _complete_artifacts(
                gnomad_data=population,
                spliceai_score=None,
                spliceai_status={"status": "not_applicable"},
            ),
        )

    assert error.value.code == "population_evidence_unavailable"
    assert error.value.status == "cache_missing"


def test_unresolved_appendix_g_path_stops_exon_cnv_classification():
    service = EvidenceOrchestrationService()
    _, base_variant = service.prepare(
        ClassificationCommand("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
    )
    variant = replace(
        base_variant,
        c_notation="deletion of exon 13",
        p_notation="p.?",
        variant_type="exon_deletion",
    )

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            {
                "spliceai_score": None,
                "spliceai_status": {"status": "not_applicable"},
                "exon_cnv_result": {
                    "is_indel": True,
                    "pm2_applicability": "unavailable",
                    "reason": "The exon interval could not be resolved",
                },
                "bayesdel_score": None,
                "bayesdel_status": {"status": "not_applicable"},
            },
        )

    assert error.value.code == "structural_population_evidence_unavailable"


def test_classification_service_composes_orchestration_and_presentation():
    marker = object()
    expected = object()

    class StubOrchestration:
        async def orchestrate(self, command):
            assert command.gene == "BRCA1"
            return marker

    class StubPresentation:
        def build(self, evidence):
            assert evidence is marker
            return expected

    result = asyncio.run(execute_variant_classification(
        ClassificationCommand("BRCA1", "c.303T>G"),
        orchestration=StubOrchestration(),
        presentation=StubPresentation(),
    ))
    assert result is expected


def test_main_module_contains_no_evidence_or_presentation_implementation():
    main_path = Path(__file__).resolve().parents[1] / "backend" / "main.py"
    source = main_path.read_text(encoding="utf-8")
    forbidden = (
        "execute_classification_request",
        "normalize_variant_input",
        "clinvar_lookup",
        "clingen_erepo_lookup",
        "generate_narrative",
        "sorted_criterion_items",
    )
    assert not any(name in source for name in forbidden)
