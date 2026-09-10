"""Boundaries and contracts of the classification application service."""

import asyncio
from pathlib import Path

import pytest

from backend.classification_dag import ClassifierEngineMode
from backend.services.variant_classification_service import execute_variant_classification
from backend.services.evidence_orchestration import (
    ClassificationCommand,
    EvidenceOrchestrationService,
    RequiredEvidenceUnavailableError,
    VariantPreparationError,
)


def test_orchestrator_prepares_policy_bound_normalized_variant():
    service = EvidenceOrchestrationService(engine_mode=ClassifierEngineMode.DAG)
    normalized, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
    )

    assert normalized.reference_transcript == "NM_007294.4"
    assert variant.gene == "BRCA1"
    assert variant.c_notation == "c.303T>G"
    assert variant.p_notation == "p.(Tyr101Ter)"
    assert variant.variant_type == "nonsense"


def test_orchestrator_rejects_invalid_input_before_provider_planning():
    service = EvidenceOrchestrationService(engine_mode=ClassifierEngineMode.DAG)
    with pytest.raises(VariantPreparationError):
        service.prepare(ClassificationCommand("BRCA1", "c.181A>C"))


def test_required_spliceai_timeout_stops_figure1a_classification():
    service = EvidenceOrchestrationService(engine_mode=ClassifierEngineMode.DAG)
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.5366C>T", "p.(Ala1789Val)")
    )

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            {
                "spliceai_score": None,
                "spliceai_status": {
                    "status": "api_error",
                    "reason": "TimeoutError: The read operation timed out",
                    "retryable": True,
                },
            },
        )

    assert error.value.code == "spliceai_temporarily_unavailable"
    assert error.value.retryable is True


def test_missing_spliceai_does_not_block_non_figure1a_ptc_classification():
    service = EvidenceOrchestrationService(engine_mode=ClassifierEngineMode.DAG)
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
    )

    service._require_complete_classification_evidence(
        variant,
        {
            "spliceai_score": None,
            "spliceai_status": {
                "status": "not_applicable",
                "required_for_classification": False,
            },
        },
    )


def test_incomplete_ps1_reference_spliceai_stops_classification():
    service = EvidenceOrchestrationService(engine_mode=ClassifierEngineMode.DAG)
    _, variant = service.prepare(
        ClassificationCommand("BRCA1", "c.5217T>A", "p.(Asp1739Glu)")
    )

    with pytest.raises(RequiredEvidenceUnavailableError) as error:
        service._require_complete_classification_evidence(
            variant,
            {
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
            },
        )

    assert error.value.code == "spliceai_temporarily_unavailable"
    assert error.value.retryable is True


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
        engine_mode=ClassifierEngineMode.DAG,
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
