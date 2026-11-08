"""General UTR-SNV projection through the pinned BRCA transcript alignments."""

import json
from pathlib import Path

from backend.domain.transcript_projection import (
    parse_transcript_alignment,
    project_utr_substitution,
)
from backend.lookups.coordinates import resolve_variant, validate_coordinate_source_manifest
from backend.population_frequency import PopulationFrequencyService
from backend.population_frequency.criteria import evaluate_frequency_criteria
from backend.reference_data.table9 import table9_lookup_ps3_bs3
from tests.dag_test_support import classify_with_dag


ROOT = Path(__file__).resolve().parents[1]
ALIGNMENTS = json.loads(
    (ROOT / "data" / "coordinates" / "brca_transcript_alignments.json").read_text(
        encoding="utf-8"
    )
)


def _transcript_sequences() -> dict[str, str]:
    records: dict[str, list[str]] = {}
    accession = ""
    for line in (
        ROOT / "data" / "reference" / "panel" / "fasta" / "transcripts.fa"
    ).read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            accession = line[1:].split()[0]
            records[accession] = []
        else:
            records[accession].append(line.strip())
    return {key: "".join(value) for key, value in records.items()}


SEQUENCES = _transcript_sequences()


def _alternate(base: str) -> str:
    return "A" if base != "A" else "C"


def test_known_utr_examples_project_to_both_assemblies():
    validate_coordinate_source_manifest()

    brca2 = resolve_variant("BRCA2", "c.*841A>G")
    assert brca2.status == "ok"
    assert brca2.source == "versioned_transcript_alignment"
    assert brca2.grch37.variant_id() == "13-32973748-A-G"
    assert brca2.grch38.variant_id() == "13-32399611-A-G"

    brca1 = resolve_variant("BRCA1", "c.-86C>T")
    assert brca1.status == "ok"
    assert brca1.grch37.variant_id() == "17-41277354-G-A"
    assert brca1.grch38.variant_id() == "17-43125337-G-A"


def test_every_reference_utr_position_is_covered_on_both_assemblies():
    for gene, raw_record in ALIGNMENTS["genes"].items():
        alignment = parse_transcript_alignment(gene, raw_record)
        sequence = SEQUENCES[alignment.transcript]
        positions = [
            (index, f"-{alignment.cds_start_index - index}")
            for index in range(alignment.cds_start_index)
        ] + [
            (index, f"*{index - alignment.cds_end_index + 1}")
            for index in range(alignment.cds_end_index, len(sequence))
        ]
        for index, position in positions:
            ref = sequence[index]
            c_notation = f"c.{position}{ref}>{_alternate(ref)}"
            for assembly in ("GRCh37", "GRCh38"):
                assert project_utr_substitution(
                    c_notation,
                    transcript_sequence=sequence,
                    transcript_alignment=alignment,
                    assembly=assembly,
                ) is not None, (gene, c_notation, assembly)


def test_projection_rejects_a_reference_allele_mismatch():
    alignment = parse_transcript_alignment("BRCA2", ALIGNMENTS["genes"]["BRCA2"])

    assert project_utr_substitution(
        "c.*841C>G",
        transcript_sequence=SEQUENCES[alignment.transcript],
        transcript_alignment=alignment,
        assembly="GRCh38",
    ) is None


def test_exact_table9_utr_record_is_used_without_figure1a_inference():
    result = classify_with_dag(
        "BRCA1",
        "5utr",
        "p.?",
        "c.-10A>G",
        spliceai_score=0.99,
        bayesdel_score=0.99,
        table9_result=table9_lookup_ps3_bs3("BRCA1", "c.-10A>G"),
    )

    assert result["criteria"]["BS3"]["applies"] is True
    assert result["criteria"]["BS3"]["strength"] == "Strong"
    assert result["criteria"]["BS3"]["decision_path"]["branch_id"] == "utr-exact"
    assert not any(code in result["criteria"] for code in ("PP3", "BP4", "BP7"))


def test_real_common_utr_snv_can_reach_appendix_g_ba1():
    resolved = resolve_variant("BRCA2", "c.*841A>G")
    data = PopulationFrequencyService.load_default().get_frequencies(
        gene="BRCA2",
        c_notation="c.*841A>G",
        grch37=resolved.grch37,
        grch38=resolved.grch38,
    )

    criteria = evaluate_frequency_criteria(
        data,
        "3utr",
        gene="BRCA2",
        c_notation="c.*841A>G",
    )
    assert criteria["BA1"]["strength"] == "Stand-alone"
    assert criteria["BA1"]["points"] == -99

