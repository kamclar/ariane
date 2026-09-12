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
        erepo_pvs1_rna_lookup=lambda gene, c: {
            "status": "not_in_registry",
            "record": None,
            "reason": "No exact approved record in the test registry.",
        },
    )


def test_provider_dag_acquires_evidence_and_classifies_without_preloaded_values():
    execution = asyncio.run(execute_classification_request(
        _request(), dependencies=_dependencies()
    ))

    assert execution.graph_version == "4.1.0-gene-policy-provider-dag"
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
        "provider.enigma.erepo_pvs1_rna",
        "provider.clinical_lr",
        "provider.protein_ps1",
    }.issubset(provider_ids)


def test_approved_erepo_pvs1_rna_is_applied_and_replaces_figure1a_prediction():
    from backend.modules.erepo_pvs1_rna import lookup_erepo_pvs1_rna
    from backend.services.classification_completeness import first_required_evidence_gap

    request = ClassificationRequest(
        variant=NormalizedVariant(
            gene="BRCA1",
            reference_transcript="NM_007294.4",
            c_notation="c.5332G>A",
            p_notation="p.(Asp1778Asn)",
            variant_type="missense",
            submitted_notation="BRCA1 c.5332G>A",
            normalization_source="test",
        )
    )
    dependencies = replace(
        _dependencies(spliceai_score=0.25, bayesdel_score=0.50),
        erepo_pvs1_rna_lookup=lookup_erepo_pvs1_rna,
    )

    execution = asyncio.run(execute_classification_request(
        request,
        dependencies=dependencies,
    ))

    assert execution.result["criteria"]["PVS1_RNA"]["strength"] == "Very Strong"
    assert execution.result["criteria"]["PVS1_RNA"]["points"] == 8
    assert "PP3" not in execution.result["criteria"]
    audit = execution.audit_record()["provider_evidence"]["erepo_pvs1_rna"]
    assert audit["status"] == "available"
    assert audit["source_checksum"]
    assert execution.provider_artifacts["spliceai_status"]["status"] == "not_applicable"
    assert execution.provider_artifacts["spliceai_status"][
        "replaced_by_curated_pvs1_rna"
    ] is True
    completeness_gap = first_required_evidence_gap(
        request.variant,
        execution.provider_artifacts,
    )
    assert completeness_gap is None or completeness_gap.source != "SpliceAI"

    invalid_artifacts = dict(execution.provider_artifacts)
    invalid_curated = dict(invalid_artifacts["erepo_pvs1_rna_result"])
    invalid_record = dict(invalid_curated["record"])
    invalid_record["guideline_version"] = "1.1.0"
    invalid_curated["record"] = invalid_record
    invalid_artifacts["erepo_pvs1_rna_result"] = invalid_curated
    invalid_gap = first_required_evidence_gap(request.variant, invalid_artifacts)
    assert invalid_gap is not None
    assert invalid_gap.source == "SpliceAI"


def test_variant_absent_from_partial_erepo_registry_is_not_treated_as_no_evidence():
    from backend.modules.erepo_pvs1_rna import lookup_erepo_pvs1_rna

    result = lookup_erepo_pvs1_rna("BRCA1", "c.4185G>A")

    assert result["status"] == "not_in_registry"
    assert result["coverage_status"] == "partial_verified_seed"
    assert "does not mean" in result["reason"]


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


def test_provider_does_not_query_spliceai_when_automatic_path_does_not_require_it():
    calls = []
    dependencies = replace(
        _dependencies(),
        spliceai_lookup=lambda gene, c: calls.append((gene, c)),
    )
    request = ClassificationRequest(
        variant=replace(
            _request().variant,
            p_notation="p.(Ala1789Ter)",
            variant_type="nonsense",
        )
    )

    execution = asyncio.run(execute_classification_request(
        request,
        dependencies=dependencies,
    ))

    assert calls == []
    assert execution.provider_artifacts["spliceai_score"] is None
    assert execution.provider_artifacts["spliceai_status"]["status"] == "not_applicable"
    assert execution.provider_artifacts["spliceai_status"][
        "required_for_classification"
    ] is False


def test_provider_audits_incomplete_ps1_reference_spliceai_lookup():
    assessed = "c.5366C>T"
    reference = "c.5366C>G"

    def score(gene, c_notation):
        return 0.03 if c_notation == assessed else None

    def status(gene, c_notation):
        if c_notation == assessed:
            return {"status": "ok", "score": 0.03, "source": "test-spliceai"}
        return {
            "status": "api_error",
            "score": None,
            "reason": "TimeoutError: The read operation timed out",
        }

    dependencies = replace(
        _dependencies(),
        spliceai_lookup=score,
        spliceai_status=status,
        ps1_candidate_lookup=lambda *args: [reference],
    )
    execution = asyncio.run(execute_classification_request(
        _request(), dependencies=dependencies
    ))
    audit = execution.provider_artifacts["spliceai_status"]

    assert audit["status"] == "ok"
    assert audit["reference_lookup_required"] is True
    assert audit["reference_lookup_complete"] is False
    assert audit["reference_variant_statuses"][reference]["retryable"] is True


def test_unavailable_bayesdel_is_audited_and_does_not_create_protein_prediction():
    execution = asyncio.run(execute_classification_request(
        _request(), dependencies=_dependencies(bayesdel_score=None)
    ))

    assert "PP3" not in execution.result["criteria"]
    assert "BP4" not in execution.result["criteria"]
    assert execution.provider_artifacts["bayesdel_score"] is None
    assert (
        execution.audit_record()["provider_evidence"]["bayesdel"]["status"]
        == EvidenceStatus.UNAVAILABLE.value
    )
    assert any(
        "BayesDel_noAF not available" in warning
        for warning in execution.result["warnings"]
    )


def test_unavailable_coordinates_do_not_query_population_as_if_variant_were_absent():
    dependencies = replace(
        _dependencies(),
        resolve_variant=lambda gene, c_notation: None,
        get_grch37=lambda resolved, gene, c_notation: None,
        get_grch38=lambda resolved, gene, c_notation: None,
    )
    execution = asyncio.run(execute_classification_request(
        _request(), dependencies=dependencies
    ))

    assert execution.provider_artifacts["coordinate_status"] == "unavailable"
    assert execution.provider_artifacts["grch37"] is None
    assert execution.provider_artifacts["grch38"] is None
    population = execution.audit_record()["provider_evidence"]["gnomad"]
    assert population["status"] == EvidenceStatus.NOT_APPLICABLE.value
    assert "coordinates are unavailable" in population["reason"]
    assert not any(
        code.startswith(("BA1", "BS1", "PM2"))
        for code in execution.result["criteria"]
    )


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
