"""Cohort generation uses tiny fake FASTA objects, never model inference."""

from collections import defaultdict
from copy import deepcopy
from types import SimpleNamespace

import pytest

from scripts.spliceai_distance_cases import build_cohort, variant_kind


class FakeSequence:
    def __init__(self, sequence):
        self.seq = sequence

    def __getitem__(self, item):
        return FakeSequence(self.seq[item])


def inputs():
    # Each strand has long introns to exercise all +/- radius boundaries.
    sequence = "ACGT" * 12500
    annotator = SimpleNamespace(
        genes=["BRCA1_tx", "BRCA2_tx"], chroms=["chr17", "chr13"], strands=["-", "+"],
        tx_starts=[101, 101], tx_ends=[40100, 40100],
        exon_starts=[[101, 15101, 30101], [101, 15101, 30101]],
        exon_ends=[[200, 15200, 40100], [200, 15200, 40100]],
        ref_fasta={"chr17": FakeSequence(sequence), "chr13": FakeSequence(sequence)},
    )
    profile = {"profile_id": "test", "reference_transcripts": {
        "BRCA1": {"ensembl": "BRCA1_tx"}, "BRCA2": {"ensembl": "BRCA2_tx"}}}
    fixture = {"profile_id": "test", "cases": []}
    return annotator, profile, fixture


def cohort(**kwargs):
    return build_cohort(*inputs(), seed=23, positions_per_region=3, **kwargs)


def test_reproducible_sampling_and_ids_all_three_alts_per_locus():
    result = cohort()
    assert result == cohort()
    alternate_seed = build_cohort(*inputs(), seed=24, positions_per_region=3)
    first = {case["variant"]: case["id"] for case in result["cases"]}
    second = {case["variant"]: case["id"] for case in alternate_seed["cases"]}
    assert all(first[variant] == second[variant] for variant in first.keys() & second.keys())
    first_loci = [row["random_locus_ids"] for row in result["coverage"]["regions"]]
    assert first_loci != [row["random_locus_ids"] for row in alternate_seed["coverage"]["regions"]]
    clusters = defaultdict(list)
    for case in result["cases"]:
        if "random" in case["selections"]:
            clusters[case["locus_id"]].append(case)
    assert len(clusters) == 30
    for cases in clusters.values():
        assert len(cases) == 3
        ref = cases[0]["variant"].split("-")[2]
        assert {case["variant"].split("-")[3] for case in cases} == set("ACGT") - {ref}
    assert result["coverage"]["random_snv_cases"] == 90
    assert "not IID" in result["design"]["statistical_unit"]


def test_every_exon_intron_full_transcript_and_transcription_order():
    result = cohort(include_indels=False)
    rows = result["coverage"]["regions"]
    assert len(rows) == 10
    for row in rows:
        assert row["random_loci"] == 3
        assert row["random_snv_cases"] == 9
        assert row["random_locus_shortfall"] == 0
    brca1_exon1 = next(row for row in rows if row["id"] == "BRCA1:exon:1")
    brca2_exon1 = next(row for row in rows if row["id"] == "BRCA2:exon:1")
    assert brca1_exon1["start_1based"] == 30101
    assert brca2_exon1["start_1based"] == 101
    positions = {int(case["variant"].split("-")[1]) for case in result["cases"]}
    assert {101, 40100} <= positions
    assert min(positions) == 101 and max(positions) == 40100
    assert result["coverage"]["variant_kinds"].keys() == {"snv"}


def test_all_refs_genomic_forward_and_simple_indel_anchors():
    annotator, _, _ = inputs()
    result = cohort()
    kinds = defaultdict(set)
    for case in result["cases"]:
        chrom, position, ref, alt = case["variant"].split("-")
        position = int(position)
        assert annotator.ref_fasta[chrom][position - 1:position - 1 + len(ref)].seq == ref
        assert 101 <= position <= position + len(ref) - 1 <= 40100
        assert variant_kind(ref, alt) == case["variant_kind"]
        if case["variant_kind"] != "snv":
            assert ref[0] == alt[0]
            assert case["selection"] == "targeted"
            assert case["selections"] == ["targeted"]
            kinds[case["variant_kind"]].add(abs(len(ref) - len(alt)))
    assert kinds == {"insertion": {1, 2, 3, 6}, "deletion": {1, 2, 3, 6}}
    assert sum(result["coverage"]["variant_kinds"][kind] for kind in kinds) <= 128


def test_radius_stress_has_signed_genomic_offsets_and_not_random_labels():
    result = cohort(include_indels=False)
    offsets_by_gene = defaultdict(set)
    for case in result["cases"]:
        for label in case["strata"]:
            if label.startswith("splice_site:"):
                _, kind, coordinate, _, offset = label.split(":")
                assert int(coordinate) - int(case["variant"].split("-")[1]) == int(offset)
                offsets_by_gene[case["gene"]].add(int(offset))
                assert "targeted" in case["selections"]
                if int(coordinate) == 15200:
                    assert kind == ("acceptor" if case["gene"] == "BRCA1" else "donor")
    assert offsets_by_gene == {gene: {-5001, -5000, -4999, 4999, 5000, 5001}
                               for gene in ("BRCA1", "BRCA2")}


def test_duplicate_controls_union_labels_keep_expected_values_and_stable_id():
    annotator, profile, fixture = inputs()
    original = cohort(include_indels=False)
    selected = next(case for case in original["cases"] if "random" in case["selections"])
    fixture["cases"] = [{key: selected[key] for key in ("variant", "gene", "transcript")}]
    fixture["cases"][0].update(id="fixture_case", delta_scores={"DS_AG": 0.123})
    fixture["cases"].append(dict(fixture["cases"][0], id="duplicate_fixture_case"))
    unchanged_fixture = deepcopy(fixture)
    result = build_cohort(annotator, profile, fixture, seed=23, positions_per_region=3, include_indels=False)
    assert len(result["cases"]) == len(original["cases"])
    control = next(case for case in result["cases"] if case["id"] == selected["id"])
    assert control["control_ids"] == ["duplicate_fixture_case", "fixture_case"]
    assert control["selection"] == "control"
    assert {"control", "random"} <= set(control["selections"])
    assert {"validation_control", "stratified_random"} <= set(control["strata"])
    assert control["delta_scores"] == {"DS_AG": 0.123}
    assert result["coverage"]["random_loci"] == original["coverage"]["random_loci"]
    assert fixture == unchanged_fixture


def test_short_regions_ambiguous_bases_are_reported_not_replaced():
    annotator, profile, fixture = inputs()
    annotator.tx_starts = [1, 1]
    annotator.tx_ends = [5, 5]
    annotator.exon_starts = [[1, 5], [1, 5]]
    annotator.exon_ends = [[2, 5], [2, 5]]
    annotator.ref_fasta = {"chr17": FakeSequence("ANCGT"), "chr13": FakeSequence("ANCGT")}
    result = build_cohort(annotator, profile, fixture, seed=1, positions_per_region=3)
    assert result["coverage"]["random_loci"] == 8
    assert result["coverage"]["random_locus_shortfall"] == 10
    assert sum(row["ambiguous_reference_loci"] for row in result["coverage"]["regions"]) == 2
    assert all("N" not in case["variant"].split("-")[2] for case in result["cases"])


@pytest.mark.parametrize("ref,alt", [("AC", "GT"), ("AC", "AG"), ("AC", "AAG"),
                                    ("CA", "A"), ("A", "CA"), ("N", "A"), ("A", "A")])
def test_reject_complex_or_invalid_alleles(ref, alt):
    with pytest.raises(ValueError):
        variant_kind(ref, alt)


def test_reject_bad_reference_bounds_transcript_and_fixture_profile():
    annotator, profile, fixture = inputs()
    fixture["cases"] = [{"id": "bad_ref", "variant": "chr17-101-C-A",
                         "gene": "BRCA1", "transcript": "BRCA1_tx"}]
    with pytest.raises(ValueError, match="Reference mismatch"):
        build_cohort(annotator, profile, fixture, seed=1, positions_per_region=1)
    fixture["cases"][0]["variant"] = "chr17-40100-TA-T"
    with pytest.raises(ValueError, match="outside transcript bounds"):
        build_cohort(annotator, profile, fixture, seed=1, positions_per_region=1)
    fixture["cases"] = []
    fixture["profile_id"] = "incorrect"
    with pytest.raises(ValueError, match="different profile"):
        build_cohort(annotator, profile, fixture, seed=1, positions_per_region=1)
    fixture["profile_id"] = "test"
    annotator.exon_starts[0][1] = 100
    with pytest.raises(ValueError, match="outside transcript bounds"):
        build_cohort(annotator, profile, fixture, seed=1, positions_per_region=1)


def test_reject_truncated_fasta_missing_pinned_transcript_and_invalid_sample_size():
    annotator, profile, fixture = inputs()
    annotator.ref_fasta["chr17"] = FakeSequence("ACGT")
    with pytest.raises(ValueError, match="Truncated"):
        build_cohort(annotator, profile, fixture, seed=1, positions_per_region=1)
    annotator, profile, fixture = inputs()
    annotator.genes[0] = "other"
    with pytest.raises(ValueError, match="pinned transcript"):
        build_cohort(annotator, profile, fixture, seed=1, positions_per_region=1)
    for count in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            build_cohort(*inputs(), seed=1, positions_per_region=count)
