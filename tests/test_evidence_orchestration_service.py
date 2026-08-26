"""Boundaries and contracts of the classification application service."""

import asyncio
from pathlib import Path

import pytest

from backend.classification_dag import ClassifierEngineMode
from backend.services.variant_classification_service import execute_variant_classification
from backend.services.evidence_orchestration import (
    ClassificationCommand,
    EvidenceOrchestrationService,
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
