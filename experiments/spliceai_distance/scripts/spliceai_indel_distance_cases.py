"""Deterministic targeted simple-indel cohort for SpliceAI distance validation.

The pinned upstream ``spliceai.utils.get_delta_scores`` implementation scores
SNVs and simple VCF-style insertions/deletions only.  It emits missing values
when both REF and ALT contain more than one base, so MNVs and complex delins are
deliberately outside this cohort rather than being decomposed or reinterpreted.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
import hashlib
import random

from experiments.spliceai_distance.scripts.spliceai_distance_cases import (
    BASES,
    _load_transcripts,
    _regions,
    variant_kind,
)


INDEL_LENGTHS = (1, 2, 3, 6)
REGION_ANCHORS = 3


def _even_positions(start: int, end: int, count: int) -> list[int]:
    """Return stable distributed integer positions, without duplicates."""
    if start > end:
        return []
    candidates = list(range(start, end + 1))
    if len(candidates) <= count:
        return candidates
    indexes = [index * (len(candidates) - 1) // (count - 1) for index in range(count)]
    return [candidates[index] for index in indexes]


def _inserted_sequence(anchor: str, length: int, pattern_index: int) -> tuple[str, str]:
    """Generate explicit low-complexity and mixed insertion stress patterns."""
    mode = pattern_index % 3
    if mode == 0:
        return anchor * length, "homopolymer"
    if mode == 1:
        partner = BASES[(BASES.index(anchor) + 1) % len(BASES)]
        return (anchor + partner) * ((length + 1) // 2), "alternating"
    start = (BASES.index(anchor) + pattern_index) % len(BASES)
    sequence = "".join(BASES[(start + offset) % len(BASES)] for offset in range(length))
    return sequence, "mixed"


def build_targeted_indel_cohort(annotator, profile, fixture, *, seed: int):
    """Build about 780 simple indels across every BRCA1/2 exon and intron.

    Three distributed anchors in each canonical exon/intron contribute one
    insertion and one deletion.  Both sides of every internal exon boundary add
    one insertion and one deletion.  Lengths 1, 2, 3 and 6 are rotated across
    loci, providing both frame-shifting and frame-preserving changes.  Fixture
    SNVs are retained solely as pinned scorer controls and are reported as such.
    """
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    if fixture.get("profile_id") != profile.get("profile_id"):
        raise ValueError("Validation fixture belongs to a different profile")

    transcripts = _load_transcripts(annotator, profile)
    by_name = {transcript.name: transcript for transcript in transcripts}
    regions = {transcript.name: _regions(transcript) for transcript in transcripts}
    case_map = {}
    skipped = Counter()

    def add(transcript, position, ref, alt, selection, strata, control=None):
        kind = variant_kind(ref, alt)
        if transcript.reference(position, len(ref)) != ref:
            raise ValueError(f"Reference mismatch for {transcript.name}:{position}:{ref}")
        variant = f"{transcript.chrom}-{position}-{ref}-{alt}"
        key = (transcript.name, variant)
        anchor_region = next(
            region for region in regions[transcript.name]
            if region["start_1based"] <= position <= region["end_1based"]
        )
        labels = set(strata) | {anchor_region["id"], anchor_region["kind"], kind}
        if key not in case_map:
            digest = hashlib.sha256(f"{transcript.name}|{variant}".encode()).hexdigest()[:24]
            case_map[key] = {
                "id": f"{transcript.gene}_{kind}_{digest}",
                "variant": variant,
                "gene": transcript.gene,
                "transcript": transcript.name,
                "variant_kind": kind,
                "locus_id": f"{transcript.name}:{position}",
                "strata": [],
                "selections": [],
                "region_id": anchor_region["id"],
            }
        case = case_map[key]
        case["strata"] = sorted(set(case["strata"]) | labels)
        case["selections"] = sorted(set(case["selections"]) | {selection})
        case["selection"] = "control" if "control" in case["selections"] else "targeted"
        if control is not None:
            case["control_ids"] = sorted(set(case.get("control_ids", [])) | {control["id"]})
            for field, value in control.items():
                if field not in ("id", "variant", "gene", "transcript"):
                    if field in case and case[field] != value:
                        raise ValueError(f"Conflicting control field: {field}")
                    case[field] = deepcopy(value)

    def simple_indel_pair(transcript, position, length, selection, labels, pattern_index):
        if position < transcript.start or position + length > transcript.end:
            skipped["allele_outside_transcript"] += 1
            return
        anchor = transcript.reference(position)
        deleted_ref = transcript.reference(position, length + 1)
        if anchor not in BASES or set(deleted_ref) - set(BASES):
            skipped["ambiguous_reference"] += 1
            return
        # Length modulo three is factual for every allele.  Calling an intronic
        # allele "frameshifting" would incorrectly assert a protein consequence.
        frame_label = ("length_multiple_of_three" if length % 3 == 0
                       else "length_not_multiple_of_three")
        common = [*labels, f"indel_length:{length}", frame_label]
        inserted, pattern = _inserted_sequence(anchor, length, pattern_index)
        inserted = inserted[:length]
        add(transcript, position, anchor, anchor + inserted, selection,
            [*common, f"insertion_pattern:{pattern}"])
        deletion_complexity = len(set(deleted_ref[1:]))
        add(transcript, position, deleted_ref, anchor, selection,
            [*common, f"deleted_sequence_unique_bases:{deletion_complexity}"])

    pattern_index = 0
    for transcript in transcripts:
        transcript_regions = regions[transcript.name]
        for region_index, region in enumerate(transcript_regions):
            # Keep the largest tested deletion entirely inside the region.
            anchors = _even_positions(
                region["start_1based"], region["end_1based"] - max(INDEL_LENGTHS),
                REGION_ANCHORS,
            )
            if len(anchors) < REGION_ANCHORS:
                skipped["region_anchor_shortfall"] += REGION_ANCHORS - len(anchors)
            region["targeted_anchor_positions"] = anchors
            for anchor_index, position in enumerate(anchors):
                length = INDEL_LENGTHS[(region_index * REGION_ANCHORS + anchor_index) % len(INDEL_LENGTHS)]
                simple_indel_pair(
                    transcript, position, length, "targeted",
                    ["region_systematic", f"region_anchor:{anchor_index + 1}"], pattern_index,
                )
                pattern_index += 1

        # Test both genomic sides of every canonical intron.  The left anchor
        # deletes into the intron; the right anchor deletes into the next exon.
        for intron_index, ((_, left_end), (right_start, _)) in enumerate(
            zip(transcript.exons, transcript.exons[1:]), start=1
        ):
            for side_index, (position, side) in enumerate((
                (left_end, "left_exon_to_intron"),
                (right_start - 1, "right_intron_to_exon"),
            )):
                length = INDEL_LENGTHS[(intron_index * 2 + side_index) % len(INDEL_LENGTHS)]
                simple_indel_pair(
                    transcript, position, length, "targeted",
                    ["splice_boundary_stress", f"boundary_side:{side}",
                     f"canonical_intron:{intron_index}"], pattern_index,
                )
                pattern_index += 1

    for control in fixture["cases"]:
        transcript = by_name.get(control["transcript"])
        if transcript is None or transcript.gene != control["gene"]:
            raise ValueError("Control does not use the pinned gene/transcript")
        chrom, position, ref, alt = control["variant"].split("-")
        if chrom != transcript.chrom:
            raise ValueError("Control chromosome does not match the pinned transcript")
        add(transcript, int(position), ref, alt, "control", ["validation_control"], control)

    cases = sorted(case_map.values(), key=lambda case: (
        case["selection"] != "control", case["gene"],
        int(case["variant"].split("-")[1]), case["variant"], case["transcript"],
    ))
    controls = [case for case in cases if case["selection"] == "control"]
    targeted = [case for case in cases if case["selection"] != "control"]
    random.Random(f"{seed}|targeted-indel-execution-order").shuffle(targeted)
    cases = controls + targeted

    region_rows = [region for transcript in transcripts for region in regions[transcript.name]]
    for region in region_rows:
        anchored = [case for case in cases if case["region_id"] == region["id"]]
        region["total_cases"] = len(anchored)
        region["total_loci"] = len({case["locus_id"] for case in anchored})

    return {
        "schema_version": 1,
        "design": {
            "seed": seed,
            "profile_id": profile.get("profile_id"),
            "reference_transcripts": {item.gene: item.name for item in transcripts},
            "coordinates": "Loaded Annotator coordinates are 1-based inclusive; alleles are genomic-forward",
            "scope": "Simple VCF-style insertions/deletions across every canonical BRCA1/2 exon and intron",
            "region_anchors": REGION_ANCHORS,
            "indel_lengths": list(INDEL_LENGTHS),
            "sampling": "Deterministic engineering coverage, not a prevalence sample",
            "statistical_unit": "Transcript/locus cluster; alleles at the same locus are correlated",
            "complex_allele_policy": "MNV and complex delins excluded because pinned upstream get_delta_scores emits missing values when both REF and ALT lengths exceed one",
            "representation": "Genomic-forward, simple left-anchored VCF alleles; deliberately includes low-complexity insertion stress cases",
        },
        "cases": cases,
        "coverage": {
            "total_cases": len(cases),
            "targeted_indel_cases": sum(case["selection"] == "targeted" for case in cases),
            "control_cases": sum(case["selection"] == "control" for case in cases),
            "total_loci": len({case["locus_id"] for case in cases}),
            "variant_kinds": dict(sorted(Counter(case["variant_kind"] for case in cases).items())),
            "length_modulo_three": {
                label: sum(label in case["strata"] for case in cases)
                for label in ("length_multiple_of_three", "length_not_multiple_of_three")
            },
            "regions": region_rows,
            "skipped_target_requests": dict(sorted(skipped.items())),
        },
    }
