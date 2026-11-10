import json
from pathlib import Path

import pytest

from backend.criteria.pvs1 import evaluate_pvs1
from backend.domain.classification import NormalizedVariant
from backend.reference_data.table4 import (
    table4_lookup_deletion,
    table4_lookup_duplication,
    table4_lookup_initiation,
    table4_lookup_pvs1_ptc,
    table4_lookup_splice,
)
from backend.services.classification_completeness import (
    first_required_evidence_gap,
)
from scripts.build_enigma_table4_snapshot import (
    EXPECTED_HEADERS,
    SNAPSHOT_DATE,
    build_snapshot,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "backend" / "data" / "enigma_table4.json"
SOURCE = (
    ROOT
    / "docs"
    / "enigma"
    / "v1.2"
    / "source"
    / "Specifications_Table4_V1.2.xlsx"
)


def _stored_table4() -> dict:
    return json.loads(DATA.read_text(encoding="utf-8"))


def test_generated_runtime_sections_match_the_workbook() -> None:
    stored = _stored_table4()
    rebuilt = build_snapshot(SOURCE)

    for section in (
        "exon_ranges",
        "initiation_rules",
        "ptc_rules",
        "splice_rules",
        "deletion_rules",
        "duplication_rules",
        "critical_boundaries",
    ):
        assert stored[section] == rebuilt[section], section
    assert rebuilt["generated"] == SNAPSHOT_DATE


def test_builder_rejects_an_unexpected_workbook_layout(tmp_path: Path) -> None:
    import openpyxl

    workbook = openpyxl.load_workbook(SOURCE)
    sheet = workbook["Table 4 - Annotated Exons"]
    first_column = next(iter(EXPECTED_HEADERS))
    sheet.cell(row=1, column=int(first_column) + 1, value="Unexpected header")
    changed_source = tmp_path / "changed-table4.xlsx"
    workbook.save(changed_source)

    with pytest.raises(
        RuntimeError,
        match="Unexpected Table 4 workbook columns",
    ):
        build_snapshot(changed_source)


def test_table4_runtime_rule_counts() -> None:
    table4 = _stored_table4()

    assert sum(len(rules) for rules in table4["splice_rules"].values()) == 575
    assert sum(len(rules) for rules in table4["ptc_rules"].values()) == 48
    assert sum(len(rules) for rules in table4["deletion_rules"].values()) == 50
    assert sum(
        len(arrangements)
        for gene_rules in table4["duplication_rules"].values()
        for arrangements in gene_rules.values()
    ) == 100
    assert len(table4["initiation_rules"]) == 2


def test_all_generated_rules_are_available_to_runtime_lookups() -> None:
    table4 = _stored_table4()

    for gene, rules in table4["splice_rules"].items():
        for notation, expected in rules.items():
            actual = table4_lookup_splice(gene, notation)
            assert actual["found"], f"{gene} {notation}"
            assert actual["pvs1_code"] == expected["pvs1_code"]
            assert actual["exon"] == expected["exon"]
            assert actual["notes"] == expected["notes"]

    for gene, rules in table4["deletion_rules"].items():
        for exon, expected in rules.items():
            actual = table4_lookup_deletion(gene, exon)
            assert actual["found"], f"{gene} {exon} deletion"
            assert actual["pvs1_code"] == expected["pvs1_code"]
            assert actual["notes"] == expected["notes"]

    for gene, rules in table4["duplication_rules"].items():
        for exon, arrangements in rules.items():
            for arrangement, expected in arrangements.items():
                actual = table4_lookup_duplication(gene, exon, arrangement)
                assert actual["found"], f"{gene} {exon} {arrangement}"
                assert actual["pvs1_code"] == expected["pvs1_code"]
                assert actual["notes"] == expected["notes"]

    for gene, expected in table4["initiation_rules"].items():
        actual = table4_lookup_initiation(gene)
        assert actual["found"], gene
        assert actual["pvs1_code"] == expected["pvs1_code"]
        assert actual["exon"] == expected["exon"]
        assert actual["notes"] == expected["notes"]


def test_all_ptc_exons_are_available_to_runtime_lookup() -> None:
    table4 = _stored_table4()

    for gene, rules in table4["ptc_rules"].items():
        critical = table4["critical_boundaries"].get(gene, {})
        for exon, expected in rules.items():
            cds_pos = table4["exon_ranges"][gene][exon][0]
            if exon == critical.get("exon"):
                first_aa = critical["boundary_aa"] + 1
                cds_pos = first_aa * 3 - 2
                termination_aa = None
            else:
                first_aa = (cds_pos + 2) // 3
                termination_aa = first_aa

            actual = table4_lookup_pvs1_ptc(
                gene,
                cds_pos,
                first_altered_aa=first_aa,
                termination_aa=termination_aa,
            )
            assert actual["exon"] == exon
            assert actual["pvs1_code"] == expected["pvs1_code"]
            assert actual["notes"] == expected.get("notes", "")


def test_both_critical_boundary_branches_match_the_snapshot() -> None:
    table4 = _stored_table4()

    for gene, boundary in table4["critical_boundaries"].items():
        boundary_aa = boundary["boundary_aa"]
        at_boundary = table4_lookup_pvs1_ptc(
            gene,
            boundary_aa * 3 - 2,
            first_altered_aa=boundary_aa,
            termination_aa=boundary_aa,
        )
        after_boundary = table4_lookup_pvs1_ptc(
            gene,
            (boundary_aa + 1) * 3 - 2,
            first_altered_aa=boundary_aa + 1,
            termination_aa=None,
        )

        assert at_boundary["pvs1_code"] == boundary["at_or_before"]["pvs1_code"]
        assert at_boundary["pm5_code"] == boundary["at_or_before"]["pm5_code"]
        assert after_boundary["pvs1_code"] == boundary["after"]["pvs1_code"]
        assert after_boundary["pm5_code"] is None


def test_brca1_last_exon_uses_the_complete_cds_range() -> None:
    table4 = _stored_table4()
    assert table4["exon_ranges"]["BRCA1"]["E23(24)"] == [5468, 5592]

    result = table4_lookup_pvs1_ptc(
        "BRCA1",
        5589,
        first_altered_aa=1863,
        termination_aa=1863,
    )
    assert result["exon"] == "E23(24)"
    assert result["pvs1_code"] == "PVS1_N/A"
    assert "CDS position" not in result["reason"]


def test_unbound_workbook_notes_reach_the_relevant_rules() -> None:
    table4 = _stored_table4()
    donor = table4["splice_rules"]["BRCA1"]["c.301+1GA"]["notes"]
    assert "Conflicting functional/clinical data" in donor

    for gene, exon in (
        ("BRCA1", "E10(11)"),
        ("BRCA2", "E4"),
        ("BRCA2", "E8"),
        ("BRCA2", "E12"),
    ):
        note = table4["deletion_rules"][gene][exon]["notes"]
        assert "intermediate risk variants" in note


def test_unknown_stop_before_last_exon_boundary_requires_review() -> None:
    result = table4_lookup_pvs1_ptc(
        "BRCA1",
        5556,
        first_altered_aa=1853,
        termination_aa=None,
    )
    assert result["pvs1_code"] is None
    assert result["review_required"] is True
    assert result["application_status"] == "review_required"

    criterion = evaluate_pvs1(
        "BRCA1",
        "frameshift",
        "p.(Tyr1853Aspfs)",
        "c.5556_5560del",
    )
    assert criterion["applies"] is False
    assert criterion["review_required"] is True
    assert criterion["application_status"] == "review_required"

    variant = NormalizedVariant(
        gene="BRCA1",
        reference_transcript="NM_007294.4",
        c_notation="c.5556_5560del",
        p_notation="p.(Tyr1853Aspfs)",
        variant_type="frameshift",
    )
    gap = first_required_evidence_gap(variant, {})
    assert gap is not None
    assert gap.source == "protein consequence"
    assert gap.status == "manual_review_required"


def test_first_change_after_last_exon_boundary_is_not_ambiguous() -> None:
    result = table4_lookup_pvs1_ptc(
        "BRCA1",
        5572,
        first_altered_aa=1858,
        termination_aa=None,
    )
    assert result["review_required"] is False
    assert result["pvs1_code"] == "PVS1_N/A"
    assert result["application_status"] == "not_applicable"


def test_layered_rna_rules_and_source_gaps_remain_visible() -> None:
    brca1_rna = table4_lookup_splice("BRCA1", "c.81-1G>C")
    brca2_rna = table4_lookup_splice("BRCA2", "c.8487+2T>C")
    assert brca1_rna["pvs1_code"] == "PVS1 (RNA)"
    assert brca2_rna["pvs1_code"] == "PVS1 (RNA)"

    assert table4_lookup_splice("BRCA2", "c.6937+2T>A")["found"] is False
    assert table4_lookup_splice("BRCA2", "c.8487+2T>A")["found"] is False
