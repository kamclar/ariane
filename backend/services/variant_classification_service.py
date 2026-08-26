"""Application use case for classifying one sequence variant."""

from __future__ import annotations

from backend.classification_dag import ClassifierEngineMode
from backend.models import ClassificationResult
from backend.services.classification_presentation import ClassificationPresentationService
from backend.services.evidence_orchestration import (
    ClassificationCommand,
    EvidenceOrchestrationService,
)


async def execute_variant_classification(
    command: ClassificationCommand,
    *,
    engine_mode: ClassifierEngineMode,
    orchestration: EvidenceOrchestrationService | None = None,
    presentation: ClassificationPresentationService | None = None,
) -> ClassificationResult:
    """Run evidence orchestration and project its result for API consumers."""
    orchestrator = orchestration or EvidenceOrchestrationService(
        engine_mode=engine_mode
    )
    presenter = presentation or ClassificationPresentationService()
    evidence = await orchestrator.orchestrate(command)
    return presenter.build(evidence)
