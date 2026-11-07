"""Pure projection of exact transcript indels onto validated base maps.

The projector never guesses across an uncovered reference interval.  It only
returns an allele when every affected genomic base, the VCF anchor, and all
left-normalization context are present in the supplied versioned maps.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import re
from typing import Any, Mapping

from backend.domain.indel import indel_operation


_DNA = re.compile(r"[ACGT]+", re.IGNORECASE)
_SNV = re.compile(r"^(BRCA[12]):c\.(\d+)([ACGT])>[ACGT]$", re.IGNORECASE)
_INTRONIC_SNV = re.compile(
    r"^(BRCA[12]):c\.(\d+)([+-]\d+)([ACGT])>[ACGT]$", re.IGNORECASE
)
_SNAPSHOT_COORD = re.compile(
    r"(?:chr)?([^:]+):(\d+):([ACGT]+)>([ACGT]+)", re.IGNORECASE
)
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@lru_cache(maxsize=1)
def _parser():
    from hgvs.parser import Parser

    return Parser()


def _parse_variant(value: str):
    from hgvs.exceptions import HGVSParseError

    try:
        return _parser().parse_hgvs_variant(value)
    except (HGVSParseError, ValueError):
        return None


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class ReferenceBase:
    chrom: str
    pos: int
    transcript_ref: str
    genomic_ref: str


@dataclass(frozen=True)
class AssemblyBaseMap:
    strand: int
    by_transcript_position: Mapping[tuple[int, int], ReferenceBase]
    by_genomic_position: Mapping[int, str]


@dataclass(frozen=True)
class ProjectedAllele:
    chrom: str
    pos: int
    ref: str
    alt: str


def _parse_snapshot_coordinate(value: object) -> tuple[str, int, str] | None:
    match = _SNAPSHOT_COORD.fullmatch(str(value or ""))
    if not match:
        return None
    return match.group(1), int(match.group(2)), match.group(3).upper()


def build_assembly_base_maps(
    coding_snapshot: Mapping[str, Mapping[str, Any]],
    coordinate_maps: Mapping[str, Mapping[str, Mapping[str, Any]]],
) -> dict[tuple[str, str], AssemblyBaseMap]:
    """Build per-base maps from already validated local coordinate sources."""
    raw: dict[tuple[str, str], dict[tuple[int, int], ReferenceBase]] = {}

    def add(
        gene: str,
        assembly: str,
        base: int,
        offset: int,
        transcript_ref: str,
        chrom: str,
        pos: int,
        genomic_ref: str,
    ) -> None:
        key = (gene.upper(), assembly)
        item = ReferenceBase(
            chrom=str(chrom).removeprefix("chr"),
            pos=int(pos),
            transcript_ref=transcript_ref.upper(),
            genomic_ref=genomic_ref.upper(),
        )
        previous = raw.setdefault(key, {}).setdefault((base, offset), item)
        if previous != item:
            raise RuntimeError(
                f"Conflicting local coordinates for {gene} c.{base}{offset:+d}/{assembly}"
            )

    for key, record in coding_snapshot.items():
        match = _SNV.fullmatch(key)
        if not match:
            continue
        gene, base_text, transcript_ref = match.groups()
        for assembly, field in (("GRCh37", "grch37"), ("GRCh38", "grch38")):
            parsed = _parse_snapshot_coordinate(record.get(field))
            if parsed is None:
                continue
            chrom, pos, genomic_ref = parsed
            add(gene, assembly, int(base_text), 0, transcript_ref, chrom, pos, genomic_ref)

    for records in coordinate_maps.values():
        for key, record in records.items():
            match = _INTRONIC_SNV.fullmatch(key)
            if not match:
                continue
            gene, base_text, offset_text, transcript_ref = match.groups()
            for assembly, field in (("GRCh37", "grch37"), ("GRCh38", "grch38")):
                value = record.get(field)
                if not isinstance(value, Mapping):
                    continue
                add(
                    gene,
                    assembly,
                    int(base_text),
                    int(offset_text),
                    transcript_ref,
                    str(value.get("chrom") or ""),
                    int(value.get("pos") or 0),
                    str(value.get("ref") or ""),
                )

    result: dict[tuple[str, str], AssemblyBaseMap] = {}
    for map_key, positions in raw.items():
        zero_offset = sorted(
            (base, item.pos) for (base, offset), item in positions.items() if offset == 0
        )
        steps = [
            right_pos - left_pos
            for (left_base, left_pos), (right_base, right_pos)
            in zip(zero_offset, zero_offset[1:])
            if right_base == left_base + 1 and abs(right_pos - left_pos) == 1
        ]
        if not steps:
            raise RuntimeError(
                f"Cannot infer transcript strand for {map_key[0]}/{map_key[1]}"
            )
        strand = 1 if steps.count(1) > steps.count(-1) else -1
        by_genomic: dict[int, str] = {}
        for item in positions.values():
            expected = (
                item.transcript_ref
                if strand == 1
                else _reverse_complement(item.transcript_ref)
            )
            if item.genomic_ref != expected:
                raise RuntimeError(
                    f"Reference-base mismatch for {map_key[0]}/{map_key[1]} at {item.pos}"
                )
            previous = by_genomic.setdefault(item.pos, item.genomic_ref)
            if previous != item.genomic_ref:
                raise RuntimeError(
                    f"Conflicting genomic reference at "
                    f"{map_key[0]}/{map_key[1]}:{item.pos}"
                )
        result[map_key] = AssemblyBaseMap(strand, positions, by_genomic)
    return result


def _normalized_vcf_allele(
    reference: Mapping[int, str],
    chrom: str,
    edit_pos: int,
    deleted: str,
    inserted: str,
) -> ProjectedAllele | None:
    deleted = deleted.upper()
    inserted = inserted.upper()
    while deleted and inserted and deleted[-1] == inserted[-1]:
        deleted = deleted[:-1]
        inserted = inserted[:-1]
    while deleted and inserted and deleted[0] == inserted[0]:
        deleted = deleted[1:]
        inserted = inserted[1:]
        edit_pos += 1
    if not deleted and not inserted:
        return None

    if not deleted or not inserted:
        sequence = inserted or deleted
        # Rotate pure insertions/deletions to the leftmost representation.  A
        # missing preceding base means normalization cannot be proven complete.
        while True:
            previous = reference.get(edit_pos - 1)
            if previous is None:
                return None
            if previous != sequence[-1]:
                break
            sequence = previous + sequence[:-1]
            edit_pos -= 1
        if inserted:
            inserted, deleted = sequence, ""
        else:
            deleted, inserted = sequence, ""
        anchor = reference.get(edit_pos - 1)
        if anchor is None:
            return None
        return ProjectedAllele(
            chrom=chrom,
            pos=edit_pos - 1,
            ref=anchor + deleted,
            alt=anchor + inserted,
        )

    return ProjectedAllele(chrom=chrom, pos=edit_pos, ref=deleted, alt=inserted)


def project_exact_indel(
    gene: str,
    transcript: str,
    c_notation: str,
    assembly_map: AssemblyBaseMap,
) -> ProjectedAllele | None:
    """Project one exact small indel, returning ``None`` when proof is incomplete."""
    operation = indel_operation(c_notation)
    if operation not in {"del", "dup", "ins", "delins"}:
        return None
    variant = _parse_variant(f"{transcript}:{c_notation}")
    if variant is None:
        return None
    start = variant.posedit.pos.start
    end = variant.posedit.pos.end
    if (
        start.base is None
        or end.base is None
        or getattr(start.datum, "name", "") != "CDS_START"
        or getattr(end.datum, "name", "") != "CDS_START"
        or int(start.base) < 1
        or int(end.base) < 1
    ):
        return None
    start_key = (int(start.base), int(start.offset or 0))
    end_key = (int(end.base), int(end.offset or 0))
    first = assembly_map.by_transcript_position.get(start_key)
    last = assembly_map.by_transcript_position.get(end_key)
    if first is None or last is None or first.chrom != last.chrom:
        return None

    edit = variant.posedit.edit
    explicit_ref = getattr(edit, "ref", None)
    explicit_ref = (
        explicit_ref.upper()
        if isinstance(explicit_ref, str) and _DNA.fullmatch(explicit_ref)
        else None
    )
    explicit_alt = getattr(edit, "alt", None)
    explicit_alt = (
        explicit_alt.upper()
        if isinstance(explicit_alt, str) and _DNA.fullmatch(explicit_alt)
        else None
    )
    strand = assembly_map.strand
    reference = assembly_map.by_genomic_position

    if operation == "ins":
        if explicit_alt is None or last.pos - first.pos != strand:
            return None
        genomic_alt = explicit_alt if strand == 1 else _reverse_complement(explicit_alt)
        edit_pos = first.pos + 1 if strand == 1 else last.pos + 1
        return _normalized_vcf_allele(
            reference, first.chrom, edit_pos, "", genomic_alt
        )

    low, high = sorted((first.pos, last.pos))
    genomic_bases = [reference.get(pos) for pos in range(low, high + 1)]
    if any(base is None for base in genomic_bases):
        return None
    genomic_ref = "".join(str(base) for base in genomic_bases)
    transcript_ref = genomic_ref if strand == 1 else _reverse_complement(genomic_ref)
    if explicit_ref is not None and explicit_ref != transcript_ref:
        return None

    if operation == "dup":
        edit_pos = high + 1 if strand == 1 else low
        return _normalized_vcf_allele(
            reference, first.chrom, edit_pos, "", genomic_ref
        )

    genomic_alt = ""
    if operation == "delins":
        if explicit_alt is None:
            return None
        genomic_alt = explicit_alt if strand == 1 else _reverse_complement(explicit_alt)
    return _normalized_vcf_allele(
        reference, first.chrom, low, genomic_ref, genomic_alt
    )


def expected_transcript_reference(
    transcript: str,
    c_notation: str,
    assembly_map: AssemblyBaseMap,
) -> str | None:
    """Return the exact local reference span for an indel when fully covered."""
    if indel_operation(c_notation) not in {"del", "dup", "delins"}:
        return None
    variant = _parse_variant(f"{transcript}:{c_notation}")
    if variant is None:
        return None
    start = variant.posedit.pos.start
    end = variant.posedit.pos.end
    if (
        start.base is None
        or end.base is None
        or getattr(start.datum, "name", "") != "CDS_START"
        or getattr(end.datum, "name", "") != "CDS_START"
        or int(start.base) < 1
        or int(end.base) < 1
    ):
        return None
    first = assembly_map.by_transcript_position.get(
        (int(start.base), int(start.offset or 0))
    )
    last = assembly_map.by_transcript_position.get(
        (int(end.base), int(end.offset or 0))
    )
    if first is None or last is None or first.chrom != last.chrom:
        return None
    low, high = sorted((first.pos, last.pos))
    bases = [assembly_map.by_genomic_position.get(pos) for pos in range(low, high + 1)]
    if any(base is None for base in bases):
        return None
    genomic_reference = "".join(str(base) for base in bases)
    return (
        genomic_reference
        if assembly_map.strand == 1
        else _reverse_complement(genomic_reference)
    )
