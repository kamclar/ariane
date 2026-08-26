"""Regression tests for the fail-closed classification DAG foundation."""

from dataclasses import dataclass
from itertools import product

import pytest

from backend.classification_dag import (
    ClassificationInputs,
    DagDefinition,
    DagDefinitionError,
    DagExecutionContext,
    DagExecutor,
    DagNodeExecutionError,
    EvidenceBundle,
    EvidenceItem,
    EvidenceStatus,
    NodeResult,
    NodeStatus,
    ClassifierEngineMode,
    execute_classification,
    get_configured_engine_mode,
)
from backend.classification_dag.runtime import compare_classification_results
from backend.classification_dag.policy import classify_by_enigma_combination
from backend.modules.classifier import evaluate_variant
from backend.modules.table9 import table9_lookup_ps3_bs3


CONTEXT = DagExecutionContext(
    run_id="test-run",
    variant_key="BRCA1:c.5366C>T",
    policy_id="ENIGMA_BRCA_VCEP_1.2",
)


@dataclass(frozen=True)
class AddOneNode:
    id: str = "add_one"
    version: str = "1"
    requires: frozenset[str] = frozenset({"input"})
    provides: frozenset[str] = frozenset({"intermediate"})

    def evaluate(self, context, inputs):
        return NodeResult.succeeded({"intermediate": inputs["input"] + 1})


@dataclass(frozen=True)
class DoubleNode:
    id: str = "double"
    version: str = "1"
    requires: frozenset[str] = frozenset({"intermediate"})
    provides: frozenset[str] = frozenset({"output"})

    def evaluate(self, context, inputs):
        return NodeResult.succeeded({"output": inputs["intermediate"] * 2})


def test_graph_topologically_orders_nodes_and_records_trace():
    graph = DagDefinition(
        id="test",
        version="1",
        seed_keys={"input"},
        nodes=(DoubleNode(), AddOneNode()),
    )
    run = DagExecutor(graph).run({"input": 4}, context=CONTEXT)

    assert run.values["output"] == 10
    assert [entry.node_id for entry in run.trace] == ["add_one", "double"]
    assert all(entry.status.value == "succeeded" for entry in run.trace)


def test_graph_rejects_missing_provider():
    with pytest.raises(DagDefinitionError, match="without providers"):
        DagDefinition(
            id="test",
            version="1",
            seed_keys=set(),
            nodes=(DoubleNode(),),
        )


def test_graph_rejects_duplicate_output_provider():
    with pytest.raises(DagDefinitionError, match="multiple providers"):
        DagDefinition(
            id="test",
            version="1",
            seed_keys={"input"},
            nodes=(AddOneNode(id="first"), AddOneNode(id="second")),
        )


def test_graph_rejects_cycle():
    @dataclass(frozen=True)
    class Left:
        id: str = "left"
        version: str = "1"
        requires: frozenset[str] = frozenset({"right_value"})
        provides: frozenset[str] = frozenset({"left_value"})

        def evaluate(self, context, inputs):
            return NodeResult.succeeded({"left_value": 1})

    @dataclass(frozen=True)
    class Right:
        id: str = "right"
        version: str = "1"
        requires: frozenset[str] = frozenset({"left_value"})
        provides: frozenset[str] = frozenset({"right_value"})

        def evaluate(self, context, inputs):
            return NodeResult.succeeded({"right_value": 1})

    with pytest.raises(DagDefinitionError, match="cycle"):
        DagDefinition(
            id="test",
            version="1",
            seed_keys=set(),
            nodes=(Left(), Right()),
        )


def test_node_contract_violation_fails_closed():
    @dataclass(frozen=True)
    class WrongOutput:
        id: str = "wrong"
        version: str = "1"
        requires: frozenset[str] = frozenset({"input"})
        provides: frozenset[str] = frozenset({"declared"})

        def evaluate(self, context, inputs):
            return NodeResult.succeeded({"different": 1})

    graph = DagDefinition(
        id="test", version="1", seed_keys={"input"}, nodes=(WrongOutput(),)
    )
    with pytest.raises(DagNodeExecutionError, match="published outputs"):
        DagExecutor(graph).run({"input": 1}, context=CONTEXT)


def test_unavailable_upstream_evidence_is_not_treated_as_not_met():
    @dataclass(frozen=True)
    class UnavailableEvidence:
        id: str = "evidence"
        version: str = "1"
        requires: frozenset[str] = frozenset({"input"})
        provides: frozenset[str] = frozenset({"evidence"})

        def evaluate(self, context, inputs):
            return NodeResult.unavailable("dataset checksum validation failed")

    @dataclass(frozen=True)
    class Criterion:
        id: str = "criterion"
        version: str = "1"
        requires: frozenset[str] = frozenset({"evidence"})
        provides: frozenset[str] = frozenset({"decision"})

        def evaluate(self, context, inputs):
            raise AssertionError("must not run without evidence")

    graph = DagDefinition(
        id="test",
        version="1",
        seed_keys={"input"},
        nodes=(UnavailableEvidence(), Criterion()),
    )
    run = DagExecutor(graph).run({"input": 1}, context=CONTEXT)
    assert "decision" not in run.values
    assert [entry.status for entry in run.trace] == [
        NodeStatus.UNAVAILABLE,
        NodeStatus.SKIPPED,
    ]
    assert "checksum validation failed" in run.trace[0].reason


def _representative_inputs():
    return ClassificationInputs(
        gene="BRCA1",
        variant_type="missense",
        c_notation="c.5366C>T",
        p_notation="p.(Ala1789Val)",
        spliceai_score=0.05,
        bayesdel_score=0.5,
        table9_result=table9_lookup_ps3_bs3("BRCA1", "c.5366C>T"),
        reference_transcript="NM_007294.4",
        submitted_notation="BRCA1 c.5366C>T",
        normalization_source="test reference engine",
    )


def _oracle_kwargs(inputs):
    return {
        "gene": inputs.gene,
        "variant_type": inputs.variant_type,
        "p_notation": inputs.p_notation,
        "c_notation": inputs.c_notation,
        "spliceai_score": inputs.spliceai_score,
        "bayesdel_score": inputs.bayesdel_score,
        "gnomad_data": inputs.gnomad_data,
        "table9_result": inputs.table9_result,
        "pp4_bp5_result": inputs.pp4_bp5_result,
        "ps1_result": inputs.ps1_result,
        "exon_cnv_result": inputs.exon_cnv_result,
        "residue_info": inputs.residue_info,
        "dup_type": inputs.dup_type,
    }


def test_classification_inputs_create_typed_variant_and_evidence_bundle():
    inputs = _representative_inputs()
    variant = inputs.normalized_variant()
    evidence = inputs.evidence_bundle()

    assert variant.variant_key == "BRCA1:c.5366C>T"
    assert variant.reference_transcript == "NM_007294.4"
    assert evidence.get("spliceai").status == EvidenceStatus.AVAILABLE
    assert evidence.get("spliceai").value == 0.05
    assert evidence.get("gnomad").status == EvidenceStatus.NOT_PROVIDED


def test_evidence_bundle_rejects_duplicate_ids():
    duplicate = EvidenceItem(
        id="spliceai",
        kind="splice_prediction",
        status=EvidenceStatus.AVAILABLE,
        value=0.05,
    )
    with pytest.raises(ValueError, match="duplicate evidence ids"):
        EvidenceBundle((duplicate, duplicate))


def test_single_strong_benign_requires_explicit_table3_eligibility():
    bp5_single_lr = {
        "BP5": {
            "applies": True,
            "strength": "Strong",
            "points": -4,
            "single_strong_likely_benign_eligible": False,
        }
    }
    single_lr_result = classify_by_enigma_combination(bp5_single_lr, -4)
    assert single_lr_result[:2] == (3, "VUS")
    assert "multiple evidence contributions" in single_lr_result[2]

    bp5_multiple_evidence = {
        "BP5": {
            **bp5_single_lr["BP5"],
            "single_strong_likely_benign_eligible": True,
        }
    }
    assert classify_by_enigma_combination(bp5_multiple_evidence, -4)[:2] == (
        2,
        "Likely Benign",
    )


def test_single_very_strong_benign_is_likely_benign():
    criterion = {
        "BS4": {
            "applies": True,
            "strength": "Very Strong",
            "points": -8,
        }
    }
    assert classify_by_enigma_combination(criterion, -8)[:2] == (
        2,
        "Likely Benign",
    )


def test_two_very_strong_pathogenic_criteria_are_pathogenic():
    criteria = {
        "PVS1": {
            "applies": True,
            "strength": "Very Strong",
            "points": 8,
        },
        "PP4": {
            "applies": True,
            "strength": "Very Strong",
            "points": 8,
        },
        "PM2": {
            "applies": True,
            "strength": "Supporting",
            "points": 1,
        },
    }

    assert classify_by_enigma_combination(criteria, 17)[:2] == (
        5,
        "Pathogenic",
    )


def test_increasing_strong_pathogenic_evidence_to_very_strong_cannot_lower_class():
    with_strong_pp4 = {
        "PVS1": {"strength": "Very Strong", "points": 8},
        "PP4": {"strength": "Strong", "points": 4},
    }
    with_very_strong_pp4 = {
        "PVS1": {"strength": "Very Strong", "points": 8},
        "PP4": {"strength": "Very Strong", "points": 8},
    }

    assert classify_by_enigma_combination(with_strong_pp4, 12)[:2] == (
        5,
        "Pathogenic",
    )
    assert classify_by_enigma_combination(with_very_strong_pp4, 16)[:2] == (
        5,
        "Pathogenic",
    )


def test_pathogenic_table3_combinations_are_monotonic_for_stronger_evidence():
    levels = (
        ("Supporting", 1),
        ("Moderate", 2),
        ("Strong", 4),
        ("Very Strong", 8),
    )
    for item_count in range(1, 6):
        for level_indexes in product(range(len(levels)), repeat=item_count):
            criteria = {
                f"P{index}": {
                    "strength": levels[level_index][0],
                    "points": levels[level_index][1],
                }
                for index, level_index in enumerate(level_indexes)
            }
            original_class = classify_by_enigma_combination(
                criteria, sum(item["points"] for item in criteria.values())
            )[0]
            for index, level_index in enumerate(level_indexes):
                if level_index == len(levels) - 1:
                    continue
                upgraded = {name: dict(item) for name, item in criteria.items()}
                upgraded[f"P{index}"] = {
                    "strength": levels[level_index + 1][0],
                    "points": levels[level_index + 1][1],
                }
                upgraded_class = classify_by_enigma_combination(
                    upgraded, sum(item["points"] for item in upgraded.values())
                )[0]
                assert upgraded_class >= original_class
            for strength, points in levels:
                extended = {
                    **criteria,
                    f"P{item_count}": {"strength": strength, "points": points},
                }
                extended_class = classify_by_enigma_combination(
                    extended, sum(item["points"] for item in extended.values())
                )[0]
                assert extended_class >= original_class


def test_benign_table3_combinations_are_monotonic_for_stronger_evidence():
    levels = (
        ("Supporting", -1),
        ("Moderate", -2),
        ("Strong", -4),
        ("Very Strong", -8),
    )
    for item_count in range(1, 6):
        for level_indexes in product(range(len(levels)), repeat=item_count):
            criteria = {
                f"B{index}": {
                    "strength": levels[level_index][0],
                    "points": levels[level_index][1],
                    "single_strong_likely_benign_eligible": False,
                }
                for index, level_index in enumerate(level_indexes)
            }
            original_class = classify_by_enigma_combination(
                criteria, sum(item["points"] for item in criteria.values())
            )[0]
            for index, level_index in enumerate(level_indexes):
                if level_index == len(levels) - 1:
                    continue
                upgraded = {name: dict(item) for name, item in criteria.items()}
                upgraded[f"B{index}"] = {
                    "strength": levels[level_index + 1][0],
                    "points": levels[level_index + 1][1],
                    "single_strong_likely_benign_eligible": False,
                }
                upgraded_class = classify_by_enigma_combination(
                    upgraded, sum(item["points"] for item in upgraded.values())
                )[0]
                assert upgraded_class <= original_class
            for strength, points in levels:
                extended = {
                    **criteria,
                    f"B{item_count}": {
                        "strength": strength,
                        "points": points,
                        "single_strong_likely_benign_eligible": False,
                    },
                }
                extended_class = classify_by_enigma_combination(
                    extended, sum(item["points"] for item in extended.values())
                )[0]
                assert extended_class <= original_class


def test_dag_preserves_bp5_provenance_and_rejects_single_lr_class2_shortcut():
    from backend.modules.pp4_bp5 import evaluate_pp4_bp5

    clinical_lr = evaluate_pp4_bp5("BRCA1", "c.1005C>A")
    execution = execute_classification(
        ClassificationInputs(
            gene="BRCA1",
            variant_type="unknown",
            c_notation="c.1005C>A",
            p_notation="p.(?)",
            pp4_bp5_result=clinical_lr,
        ),
        mode="dag",
    )

    bp5 = execution.result["criteria"]["BP5"]
    assert bp5["likelihood_ratio_contribution_count"] == 1
    assert bp5["single_strong_likely_benign_eligible"] is False
    assert execution.result["predicted_class"] == 3
    assert "multiple evidence contributions" in execution.result["classification_note"]


def test_dag_mode_preserves_current_clinical_result_exactly():
    inputs = _representative_inputs()
    legacy = evaluate_variant(**_oracle_kwargs(inputs))
    execution = execute_classification(inputs, mode="dag")

    assert compare_classification_results(legacy, execution.result) == ()
    assert execution.graph_id == "ariane.vcep.classification"
    assert [entry.node_id for entry in execution.trace] == [
        "contract.classification_inputs",
        "provider.enigma.table9",
        "context.spliceai.provenance",
        "rule.population_frequency",
        "rule.exon_cnv.population",
        "rule.functional.table9",
        "rule.pvs1_pm5",
        "rule.clinical_lr",
        "rule.protein_ps1",
        "rule.bioinformatic.figure1a",
        "policy.evidence_interactions",
        "policy.enigma_combination",
        "review.manual_triage",
        "contract.variant_assertion",
        "projection.public_result",
    ]


def test_table9_provider_records_version_and_runtime_checksum():
    inputs = ClassificationInputs(
        gene="BRCA1",
        variant_type="missense",
        c_notation="c.509G>A",
        p_notation="p.(Arg170Gln)",
        reference_transcript="NM_007294.4",
    )
    execution = execute_classification(inputs, mode="dag")
    provider_trace = next(
        entry for entry in execution.trace
        if entry.node_id == "provider.enigma.table9"
    )

    assert provider_trace.status == NodeStatus.SUCCEEDED
    assert provider_trace.provenance["source_id"] == "enigma-v1.2-table9"
    assert provider_trace.provenance["source_version"] == "1.2.0"
    assert len(provider_trace.provenance["source_checksum"]) == 64
    assert provider_trace.provenance["reviewed"] is True


def test_table9_provider_rejects_conflicting_pre_dag_value():
    inputs = ClassificationInputs(
        gene="BRCA1",
        variant_type="missense",
        c_notation="c.509G>A",
        p_notation="p.(Arg170Gln)",
        table9_result={"reviewed": False, "applies": False},
    )

    with pytest.raises(DagNodeExecutionError, match="Table 9 evidence conflicts"):
        execute_classification(inputs, mode="dag")


@pytest.mark.parametrize(
    ("gene", "variant_type", "c_notation", "p_notation"),
    (
        ("BRCA1", "nonsense", "c.303T>G", "p.(Tyr101Ter)"),
        ("BRCA1", "frameshift", "c.5266dup", "p.(Gln1756ProfsTer74)"),
        ("BRCA1", "synonymous", "c.4185G>A", "p.(Gln1395=)"),
        ("BRCA1", "intronic", "c.548-9A>G", "p.(?)"),
        (
            "BRCA2",
            "exon_deletion",
            "c.(793+1_794-1)_(1909+1_1910-1)del",
            "p.(?)",
        ),
    ),
)
def test_dag_contract_accepts_supported_variant_families(
    gene, variant_type, c_notation, p_notation
):
    inputs = ClassificationInputs(
        gene=gene,
        variant_type=variant_type,
        c_notation=c_notation,
        p_notation=p_notation,
        table9_result=table9_lookup_ps3_bs3(gene, c_notation),
    )
    direct = evaluate_variant(**_oracle_kwargs(inputs))
    dag = execute_classification(inputs, mode="dag")
    assert compare_classification_results(direct, dag.result) == ()


def test_unknown_runtime_mode_is_rejected(monkeypatch):
    monkeypatch.setenv("ARIANE_CLASSIFIER_ENGINE", "fallback")
    with pytest.raises(ValueError, match="Invalid ARIANE_CLASSIFIER_ENGINE"):
        execute_classification(_representative_inputs())


def test_native_dag_is_default_and_does_not_require_legacy_evaluator(monkeypatch):
    monkeypatch.delenv("ARIANE_CLASSIFIER_ENGINE", raising=False)
    assert get_configured_engine_mode() == ClassifierEngineMode.DAG
    execution = execute_classification(_representative_inputs())
    assert execution.engine_mode == ClassifierEngineMode.DAG
    assert all("legacy." not in entry.node_id for entry in execution.trace)


def test_legacy_runtime_mode_is_no_longer_available():
    with pytest.raises(ValueError, match="Invalid ARIANE_CLASSIFIER_ENGINE"):
        execute_classification(_representative_inputs(), mode="legacy")
