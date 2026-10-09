"""Deterministic, transcript-wide cohorts without importing SpliceAI or models.

The input is a *loaded* SpliceAI Annotator: its transcript/exon coordinates are
1-based inclusive (the annotation TSV starts are converted by Annotator).
All REF/ALT alleles and offsets use forward genomic orientation, on both strands.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import random


BASES = "ACGT"
RADIUS_OFFSETS = (-5001, -5000, -4999, 4999, 5000, 5001)
BOUNDARY_OFFSETS = (-2, -1, 0, 1, 2)
INDEL_LENGTHS = (1, 2, 3, 6)


@dataclass(frozen=True)
class Transcript:
    gene: str
    name: str
    chrom: str
    strand: str
    start: int
    end: int
    exons: tuple
    sequence: str

    def reference(self, position, length=1):
        if position < self.start or position + length - 1 > self.end or length < 1:
            raise ValueError(f"Allele outside transcript bounds: {self.name}:{position}")
        return self.sequence[position - self.start:position - self.start + length]


def _load_transcripts(annotator, profile):
    names = [str(name) for name in annotator.genes]
    transcripts = []
    for gene in ("BRCA1", "BRCA2"):
        name = profile["reference_transcripts"][gene]["ensembl"]
        if names.count(name) != 1:
            raise ValueError(f"Expected exactly one pinned transcript: {name}")
        index = names.index(name)
        chrom, strand = str(annotator.chroms[index]), str(annotator.strands[index])
        start, end = int(annotator.tx_starts[index]), int(annotator.tx_ends[index])
        starts = tuple(map(int, annotator.exon_starts[index]))
        ends = tuple(map(int, annotator.exon_ends[index]))
        if strand not in ("+", "-") or start < 1 or end < start:
            raise ValueError(f"Invalid transcript coordinates/strand: {name}")
        if not starts or len(starts) != len(ends):
            raise ValueError(f"Invalid exon arrays: {name}")
        exons = tuple(zip(starts, ends))
        if any(left < start or right > end or left > right for left, right in exons):
            raise ValueError(f"Exon outside transcript bounds: {name}")
        if any(right >= next_left for (_, right), (next_left, _) in zip(exons, exons[1:])):
            raise ValueError(f"Exons must be sorted and non-overlapping: {name}")
        if exons[0][0] != start or exons[-1][1] != end:
            raise ValueError(f"Full transcript must start/end at terminal exons: {name}")
        record = annotator.ref_fasta[chrom][start - 1:end]
        sequence = str(getattr(record, "seq", record)).upper()
        if len(sequence) != end - start + 1:
            raise ValueError(f"Truncated reference FASTA for {name}")
        transcripts.append(Transcript(gene, name, chrom, strand, start, end, exons, sequence))
    return transcripts


def _regions(transcript):
    intervals = [("exon", left, right) for left, right in transcript.exons]
    intervals += [("intron", right + 1, next_left - 1)
                  for (_, right), (next_left, _) in zip(transcript.exons, transcript.exons[1:])
                  if next_left > right + 1]
    intervals.sort(key=lambda item: item[1], reverse=transcript.strand == "-")
    numbers = Counter()
    regions = []
    for kind, start, end in intervals:
        numbers[kind] += 1
        regions.append({
            "id": f"{transcript.gene}:{kind}:{numbers[kind]}",
            "gene": transcript.gene, "transcript": transcript.name,
            "kind": kind, "number": numbers[kind],
            "start_1based": start, "end_1based": end, "length": end - start + 1,
        })
    return regions


def variant_kind(ref, alt):
    """Accept only SNVs and simple left-anchored insertions/deletions."""
    if not ref or not alt or set(ref + alt) - set(BASES) or ref == alt:
        raise ValueError("Alleles must be distinct, nonempty canonical DNA")
    if len(ref) == len(alt) == 1:
        return "snv"
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt:
        return "deletion"
    if len(ref) == 1 and len(alt) > 1 and alt[0] == ref:
        return "insertion"
    raise ValueError("Unsupported complex/MNV or non-left-anchored allele")


def _evenly_spaced(items, limit):
    if len(items) <= limit:
        return list(items)
    return [items[index * (len(items) - 1) // (limit - 1)] for index in range(limit)]


def build_cohort(annotator, profile, fixture, *, seed: int,
                 positions_per_region: int, include_indels: bool = True):
    """Build random locus clusters, targeted stress alleles, and fixture controls.

    Random sampling is uniform without replacement within each exon/intron among
    canonical reference bases; all three SNV alternatives share a locus_id. Small
    or ambiguous regions report a shortfall. Counts are not an IID allele sample
    and target/control variants must not enter random-sample denominators.
    """
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    if (isinstance(positions_per_region, bool) or not isinstance(positions_per_region, int)
            or positions_per_region < 1):
        raise ValueError("positions_per_region must be a positive integer")
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
        anchor_region = next(region for region in regions[transcript.name]
                             if region["start_1based"] <= position <= region["end_1based"])
        labels = set(strata) | {anchor_region["id"], anchor_region["kind"], kind}
        if key not in case_map:
            digest = hashlib.sha256(f"{transcript.name}|{variant}".encode()).hexdigest()[:24]
            case_map[key] = {
                "id": f"{transcript.gene}_{kind}_{digest}", "variant": variant,
                "gene": transcript.gene, "transcript": transcript.name,
                "variant_kind": kind, "locus_id": f"{transcript.name}:{position}",
                "strata": [], "selections": [], "region_id": anchor_region["id"],
            }
        case = case_map[key]
        case["strata"] = sorted(set(case["strata"]) | labels)
        case["selections"] = sorted(set(case["selections"]) | {selection})
        case["selection"] = next(item for item in ("control", "random", "targeted")
                                 if item in case["selections"])
        if control is not None:
            case["control_ids"] = sorted(set(case.get("control_ids", [])) | {control["id"]})
            # Preserve original expected score fields for the runner's control checks.
            for field, value in control.items():
                if field not in ("id", "variant", "gene", "transcript"):
                    if field in case and case[field] != value:
                        raise ValueError(f"Conflicting control field: {field}")
                    case[field] = deepcopy(value)

    def snvs(transcript, position, selection, strata):
        if not transcript.start <= position <= transcript.end:
            skipped["targeted_loci_outside_transcript"] += 1
            return
        ref = transcript.reference(position)
        if ref not in BASES:
            skipped["targeted_loci_ambiguous_reference"] += 1
            return
        for alt in BASES:
            if ref != alt:
                add(transcript, position, ref, alt, selection, strata)

    for transcript in transcripts:
        for region in regions[transcript.name]:
            eligible = [position for position in range(region["start_1based"], region["end_1based"] + 1)
                        if transcript.reference(position) in BASES]
            # Independent streams keep one region's draw stable when other regions change.
            rng = random.Random(f"{seed}|{transcript.name}|{region['id']}")
            selected = sorted(rng.sample(eligible, min(positions_per_region, len(eligible))))
            region.update({"eligible_loci": len(eligible),
                           "ambiguous_reference_loci": region["length"] - len(eligible),
                           "requested_random_loci": positions_per_region,
                           "random_locus_ids": [f"{transcript.name}:{position}" for position in selected],
                           "random_loci": len(selected),
                           "random_locus_shortfall": positions_per_region - len(selected)})
            for position in selected:
                snvs(transcript, position, "random", ["stratified_random"])
        for exon_start, exon_end in transcript.exons:
            for boundary in (exon_start, exon_end):
                for offset in BOUNDARY_OFFSETS:
                    snvs(transcript, boundary + offset, "targeted", [
                        "boundary_stress", f"boundary_offset:{offset:+d}",
                        f"boundary_coordinate:{boundary}",
                    ])
        for (_, left_end), (right_start, _) in zip(transcript.exons, transcript.exons[1:]):
            for site, site_kind in ((left_end, "donor" if transcript.strand == "+" else "acceptor"),
                                    (right_start, "acceptor" if transcript.strand == "+" else "donor")):
                for offset in RADIUS_OFFSETS:
                    snvs(transcript, site - offset, "targeted", [
                        "radius_stress", f"splice_site_kind:{site_kind}",
                        f"splice_site_offset:{offset:+d}", f"splice_site_coordinate:{site}",
                        f"splice_site:{site_kind}:{site}:offset:{offset:+d}",
                    ])
        if include_indels:
            # Four distributed exons and four introns, alternating boundary and
            # midpoint anchors; 8 * 4 * 2 = at most 64 indels per gene.
            for kind in ("exon", "intron"):
                chosen = _evenly_spaced([region for region in regions[transcript.name]
                                        if region["kind"] == kind], 4)
                for index, region in enumerate(chosen):
                    position = (region["start_1based"] if index % 2 == 0 else
                                (region["start_1based"] + region["end_1based"]) // 2)
                    anchor = transcript.reference(position)
                    if anchor not in BASES:
                        skipped["indel_ambiguous_anchor"] += 1
                        continue
                    for length in INDEL_LENGTHS:
                        labels = ["indel_stress", f"indel_length:{length}",
                                  "indel_boundary_anchor" if index % 2 == 0 else "indel_midpoint_anchor"]
                        # A defined, reference-independent synthetic insert; never
                        # reverse-complemented on the negative-strand transcript.
                        inserted = ("ACGTAC" * 2)[:length]
                        add(transcript, position, anchor, anchor + inserted, "targeted", labels)
                        if position + length > transcript.end:
                            skipped["indel_deletion_outside_transcript"] += 1
                            continue
                        ref = transcript.reference(position, length + 1)
                        if set(ref) - set(BASES):
                            skipped["indel_deletion_ambiguous_reference"] += 1
                            continue
                        add(transcript, position, ref, anchor, "targeted", labels)

    for control in fixture["cases"]:
        transcript = by_name.get(control["transcript"])
        if transcript is None or transcript.gene != control["gene"]:
            raise ValueError("Control does not use the pinned gene/transcript")
        chrom, pos, ref, alt = control["variant"].split("-")
        if chrom != transcript.chrom:
            raise ValueError("Control chromosome does not match the pinned transcript")
        add(transcript, int(pos), ref, alt, "control", ["validation_control"], control)

    cases = sorted(case_map.values(), key=lambda case: (
        case["gene"], int(case["variant"].split("-")[1]), case["variant"], case["transcript"]))
    # Validate known controls first, then exercise both strands and indel
    # alignment early. Shuffle the remainder reproducibly, so bounded runs
    # do not visit only the beginning of BRCA1. This changes no membership.
    controls = [case for case in cases if "control" in case["selections"]]
    remaining = [case for case in cases if "control" not in case["selections"]]
    random.Random(f"{seed}|execution-order").shuffle(remaining)
    prefix = list(controls)
    for gene in sorted({case["gene"] for case in remaining}):
        for kind in ("deletion", "insertion"):
            candidate = next((case for case in remaining
                              if case["gene"] == gene and case["variant_kind"] == kind), None)
            if candidate is not None:
                prefix.append(candidate)
                remaining.remove(candidate)
    cases = prefix + remaining
    region_rows = [region for transcript in transcripts for region in regions[transcript.name]]
    for region in region_rows:
        anchored = [case for case in cases if case["region_id"] == region["id"]]
        region["total_cases"] = len(anchored)
        region["total_loci"] = len({case["locus_id"] for case in anchored})
        region["random_snv_cases"] = sum("random" in case["selections"] and case["variant_kind"] == "snv"
                                         for case in anchored)
    return {
        "schema_version": 1,
        "design": {
            "seed": seed, "positions_per_region": positions_per_region,
            "include_indels": include_indels, "profile_id": profile.get("profile_id"),
            "reference_transcripts": {item.gene: item.name for item in transcripts},
            "coordinates": "Loaded Annotator coordinates are 1-based inclusive; alleles are genomic-forward",
            "region_numbering": "5-prime to 3-prime transcript order; coverage uses allele anchor position",
            "scope": "Full pinned transcripts: every exon (including UTR sequence) and intron; no CDS restriction",
            "sampling": "Uniform canonical loci without replacement within each region; all 3 SNV ALT alleles per locus",
            "statistical_unit": "Transcript/locus cluster, not ALT allele; stratified finite-population sample, not IID pooled variants",
            "targeted_policy": "Boundary/radius/indel stress and controls are not random observations; selections may overlap",
            "radius_offsets": list(RADIUS_OFFSETS),
            "radius_offset_definition": "splice-site exonic boundary coordinate minus variant anchor, genomic orientation",
            "boundary_offsets": list(BOUNDARY_OFFSETS),
            "indel_lengths": list(INDEL_LENGTHS) if include_indels else [],
            "indel_limit_per_gene": 64 if include_indels else 0,
            "indel_representation": "Simple genomic left-anchored alleles; not repeat-left-normalized; no MNV/delins",
        },
        "cases": cases,
        "coverage": {
            "total_cases": len(cases), "total_loci": len({case["locus_id"] for case in cases}),
            "random_loci": sum(region["random_loci"] for region in region_rows),
            "random_snv_cases": sum(region["random_snv_cases"] for region in region_rows),
            "random_locus_shortfall": sum(region["random_locus_shortfall"] for region in region_rows),
            "variant_kinds": dict(sorted(Counter(case["variant_kind"] for case in cases).items())),
            "selection_memberships": {selection: sum(selection in case["selections"] for case in cases)
                                      for selection in ("random", "targeted", "control")},
            "regions": region_rows, "skipped_target_requests": dict(sorted(skipped.items())),
        },
    }
