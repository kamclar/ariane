from __future__ import annotations

import asyncio

from backend.contracts import BatchRequest, ClassificationResult, VariantRequest
from backend.services.variant_classification_service import (
    UsageContext,
    VariantClassificationService,
)


def _result() -> ClassificationResult:
    return ClassificationResult(
        variant="BRCA1 c.4185G>A p.(Gln1395=)",
        gene="BRCA1",
        c_notation="c.4185G>A",
        p_notation="p.(Gln1395=)",
        predicted_class=3,
        predicted_label="VUS",
        total_points=0,
    )


class _UsageRepository:
    def __init__(self) -> None:
        self.events: list[dict] = []

    def record(self, **values) -> None:
        self.events.append(values)


def test_single_workflow_adds_request_metadata_and_records_usage() -> None:
    request = VariantRequest(gene="BRCA1", c_notation="c.4185G>A")
    usage_repository = _UsageRepository()
    service = VariantClassificationService()

    async def classify_cached(*_args):
        return _result(), "hit", "fingerprint-1"

    completed = asyncio.run(service.classify_single(
        request,
        usage=UsageContext(
            actor_id="visitor-1",
            actor_type="visitor",
            request_id="request-1",
        ),
        usage_repository=usage_repository,
        classify_cached=classify_cached,
    ))

    assert completed.cache_status == "hit"
    assert completed.classifier_fingerprint == "fingerprint-1"
    assert completed.result.reference_transcript == "NM_007294.4"
    assert completed.result.submitted_notation == "c.4185G>A"
    assert completed.policy["policy_id"] == "ENIGMA_BRCA_VCEP_1.2"
    assert len(usage_repository.events) == 1
    assert usage_repository.events[0]["status"] == "completed"
    assert usage_repository.events[0]["request_mode"] == "single"


def test_batch_preparation_keeps_valid_items_and_describes_invalid_items() -> None:
    request = BatchRequest(variants=[
        {"gene": "BRCA1", "c_notation": "c.4185G>A"},
        {"gene": "BRCA1", "c_notation": "not-hgvs"},
    ])

    prepared = VariantClassificationService.prepare_batch(request)

    assert prepared.classification_units == 1
    assert prepared.valid_items[0][0] == 0
    assert prepared.validation_errors[0].index == 1
    assert prepared.validation_errors[0].error_code == "invalid_variant"
    assert prepared.audit_events[0].event == "batch_item_validation_error"
