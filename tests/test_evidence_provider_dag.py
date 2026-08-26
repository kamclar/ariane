"""Production evidence-provider DAG regression tests."""

import asyncio
from dataclasses import replace
from threading import Barrier
from types import SimpleNamespace

import pytest

from backend.classification_dag import (
    ClassificationRequest,
    DagNodeExecutionError,
    EvidenceStatus,
    NormalizedVariant,
    ProviderDependencies,
    execute_classification_request,
)


def _request() -> ClassificationRequest:
    return ClassificationRequest(
        variant=NormalizedVariant(
            gene="BRCA1",
            reference_transcript="NM_007294.4",
            c_notation="c.5366C>T",
            p_notation="p.(Ala1789Val)",
            variant_type="missense",
            submitted_notation="BRCA1 c.5366C>T",
            normalization_source="test",
        )
    )


def _dependencies(*, spliceai_score=0.03, bayesdel_score=0.438):
    resolved = SimpleNamespace(
        status="ok",
        source="test-coordinate-provider",
        warnings=[],
    )

    def ps1_result(*args, **kwargs):
        return {
            "applies": False,
            "strength": None,
            "points": 0,
            "reason": "No matching reference in the test provider.",
            "reference_status": "not_found",
            "application_status": "not_applied",
            "review_required": False,
            "blocking_reasons": [],
            "candidates": [],
            "vua_splice_sources_checked": [],
            "vua_splice_evidence_status": "none_identified",
        }

    return ProviderDependencies(
        resolve_variant=lambda gene, c: resolved,
        get_grch37=lambda resolved_map, gene, c: {
            "chrom": "17", "pos": 1, "ref": "C", "alt": "T",
        },
        get_grch38=lambda resolved_map, gene, c: {
            "chrom": "17", "pos": 2, "ref": "C", "alt": "T",
        },
        spliceai_lookup=lambda gene, c: spliceai_score,
        spliceai_status=lambda gene, c: {
            "status": "ok" if spliceai_score is not None else "api_error",
            "score": spliceai_score,
            "reason": "test SpliceAI provider",
            "source": "test-spliceai",
            "scoring_profile_id": "test-profile",
        },
        bayesdel_lookup=lambda gene, c: (bayesdel_score, None),
        bayesdel_status=lambda gene, c: {
            "status": "ok" if bayesdel_score is not None else "no_score",
            "reason": "test BayesDel provider",
        },
        gnomad_lookup=lambda **kwargs: {},
        clinical_lr_lookup=lambda gene, c: {"applies": False},
        exon_cnv_lookup=lambda gene, c: None,
        residue_lookup=lambda gene, p: {},
        ps1_candidate_lookup=lambda gene, c, p, variant_type: [],
        splice_source_lookup=lambda gene, c, table9: {
            "status": "none_identified", "sources_checked": [],
        },
        select_ps1_spliceai=lambda score, table9: (score, "test-spliceai"),
        ps1_lookup=ps1_result,
    )


def test_provider_dag_acquires_evidence_and_classifies_without_preloaded_values():
    execution = asyncio.run(execute_classification_request(
        _request(), dependencies=_dependencies()
    ))

    assert execution.graph_version == "4.0.0-gene-policy-provider-dag"
    assert execution.result["criteria"]["PP3"]["points"] == 1
    assert execution.provider_artifacts["spliceai_score"] == 0.03
    assert execution.provider_artifacts["bayesdel_score"] == 0.438
    assert execution.audit_record()["provider_evidence"]["spliceai"]["status"] == "available"
    provider_ids = {
        entry.node_id for entry in execution.trace
        if entry.node_id.startswith("provider.")
    }
    assert {
        "provider.coordinates",
        "provider.spliceai",
        "provider.bayesdel",
        "provider.gnomad",
        "provider.enigma.table9",
        "provider.clinical_lr",
        "provider.protein_ps1",
    }.issubset(provider_ids)


def test_unavailable_spliceai_remains_unavailable_and_does_not_enable_bayesdel_pp3():
    execution = asyncio.run(execute_classification_request(
        _request(), dependencies=_dependencies(spliceai_score=None)
    ))

    assert "PP3" not in execution.result["criteria"]
    assert execution.provider_artifacts["spliceai_score"] is None
    splice_trace = next(
        entry for entry in execution.trace if entry.node_id == "provider.spliceai"
    )
    assert splice_trace.provenance["evidence_status"] == EvidenceStatus.UNAVAILABLE.value
    assert any("SpliceAI is unavailable" in warning for warning in execution.result["warnings"])


def test_independent_remote_providers_run_concurrently():
    barrier = Barrier(2, timeout=2)
    dependencies = _dependencies()

    def coordinate_lookup(gene, c):
        barrier.wait()
        return SimpleNamespace(
            status="ok", source="test-coordinate-provider", warnings=[]
        )

    def bayesdel_lookup(gene, c):
        barrier.wait()
        return 0.438, None

    execution = asyncio.run(execute_classification_request(
        _request(),
        dependencies=replace(
            dependencies,
            resolve_variant=coordinate_lookup,
            bayesdel_lookup=bayesdel_lookup,
        ),
    ))
    assert execution.result["criteria"]["PP3"]["points"] == 1


def test_local_provider_contract_failure_stops_classification():
    def broken_gnomad(**kwargs):
        raise RuntimeError("validated gnomAD dataset could not be read")

    with pytest.raises(DagNodeExecutionError) as error:
        asyncio.run(execute_classification_request(
            _request(),
            dependencies=replace(_dependencies(), gnomad_lookup=broken_gnomad),
        ))
    assert error.value.node_id == "provider.gnomad"
    assert "validated gnomAD dataset could not be read" in str(error.value)
