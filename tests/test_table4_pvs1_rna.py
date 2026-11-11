"""Regression coverage for exact allele-specific Table 4 RNA decisions."""

from __future__ import annotations

from collections import Counter
import re

from backend.criteria.pvs1 import evaluate_pvs1
from backend.criteria.pvs1_rna import resolve_pvs1_rna
from backend.reference_data.st2_pvs1_rna import lookup_st2_pvs1_rna_evidence
from backend.reference_data.table4 import (
    load_table4_data,
    lookup_table4_pvs1_rna,
    parse_pvs1_code_strength,
)


def _input_hgvs(table4_key: str) -> str:
    """Restore the substitution separator omitted by the Table 4 workbook."""
    return re.sub(r"([ACGT])([ACGT])$", r"\1>\2", table4_key)


def test_every_exact_table4_rna_rule_has_a_deterministic_runtime_decision():
    data = load_table4_data()
    observed = Counter()

    for gene, rules in data["splice_rules"].items():
        for source_key, rule in rules.items():
            source_code = str(rule.get("pvs1_code") or "")
            if "(RNA)" not in source_code:
                continue

            result = lookup_table4_pvs1_rna(gene, _input_hgvs(source_key))
            strength, points, requires_rna = parse_pvs1_code_strength(source_code)

            assert requires_rna is True
            assert result["status"] == (
                "eligible" if strength is not None else "not_applicable_rna"
            )
            assert result["record"]["source_code"] == source_code
            assert result["record"]["strength"] == strength
            assert result["record"]["points"] == points
            assert result["record"]["decision_scope"] == "exact_allele"
            assert result["source_version"] == "1.2.0"
            assert result["source_checksum"]

            decision = resolve_pvs1_rna(
                gene,
                _input_hgvs(source_key),
                erepo_registry_result=None,
                table4_evidence_result=result,
                st2_evidence_result=None,
            ).decision
            assert decision["applies"] is (strength is not None)
            assert decision["strength"] == strength
            assert decision["points"] == points
            assert decision["prediction_replaced"] is True
            assert decision["application_status"] == (
                "applied_from_table4"
                if strength is not None
                else "not_applicable_from_table4"
            )
            observed[source_code] += 1

    assert observed == {
        "PVS1 (RNA)": 126,
        "PVS1_Strong (RNA)": 3,
        "PVS1_Moderate (RNA)": 2,
        "PVS1_Supporting (RNA)": 11,
        "PVS1_N/A (RNA)": 2,
    }


def test_complex_cohort_examples_use_the_final_table4_rna_weight():
    for c_notation, expected_note in (
        ("c.476-2A>G", ""),
        ("c.9118-2A>G", "PTC_NMD observed"),
    ):
        result = lookup_table4_pvs1_rna("BRCA2", c_notation)

        assert result["status"] == "eligible"
        assert result["record"]["source_code"] == "PVS1 (RNA)"
        assert result["record"]["strength"] == "Very Strong"
        assert result["record"]["points"] == 8
        assert result["record"]["notes"] == expected_note


def test_table4_rna_not_applicable_rows_do_not_acquire_points():
    for gene, c_notation in (
        ("BRCA1", "c.594-2A>C"),
        ("BRCA2", "c.794-2A>G"),
    ):
        result = lookup_table4_pvs1_rna(gene, c_notation)

        assert result["status"] == "not_applicable_rna"
        assert result["record"]["source_code"] == "PVS1_N/A (RNA)"
        assert result["record"]["strength"] is None
        assert result["record"]["points"] == 0


def test_rna_source_precedence_is_defined_once():
    selected = resolve_pvs1_rna(
        "BRCA1",
        "c.1A>G",
        erepo_registry_result={"status": "ambiguous", "record": {}},
        table4_evidence_result={"status": "eligible", "record": {}},
        st2_evidence_result={"status": "eligible", "record": {}},
    )
    assert selected.source == "erepo"
    assert selected.evidence_item_id == "erepo_pvs1_rna"

    selected = resolve_pvs1_rna(
        "BRCA2",
        "c.476-2A>G",
        erepo_registry_result={"status": "not_in_registry"},
        table4_evidence_result=lookup_table4_pvs1_rna(
            "BRCA2",
            "c.476-2A>G",
        ),
        st2_evidence_result=None,
    )
    assert selected.source == "table4"
    assert selected.evidence_item_id == "table4_pvs1_rna"

    selected = resolve_pvs1_rna(
        "BRCA1",
        "c.4185G>A",
        erepo_registry_result={"status": "not_in_registry"},
        table4_evidence_result={"status": "not_in_source"},
        st2_evidence_result=lookup_st2_pvs1_rna_evidence(
            "BRCA1",
            "c.4185G>A",
        ),
    )
    assert selected.source == "st2"
    assert selected.evidence_item_id == "st2_pvs1_rna"


def test_table4_rna_not_applicable_replaces_the_prediction_branch():
    selected = resolve_pvs1_rna(
        "BRCA1",
        "c.594-2A>C",
        erepo_registry_result=None,
        table4_evidence_result=lookup_table4_pvs1_rna(
            "BRCA1",
            "c.594-2A>C",
        ),
        st2_evidence_result=None,
    )

    assert selected.source == "table4"
    assert selected.replaces_prediction is True


def test_predictive_pvs1_branch_points_to_the_exact_rna_decision():
    result = evaluate_pvs1(
        "BRCA2",
        "splice_site",
        "p.?",
        c_notation="c.476-2A>G",
    )

    assert result["applies"] is False
    assert result["requires_rna"] is True
    assert "evaluated by the PVS1_RNA branch" in result["reason"]
    assert "outside the automated" not in result["reason"]
