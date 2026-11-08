"""Pure projection of UTR substitutions through versioned transcript alignments.

The projection is deliberately limited to exact single-nucleotide substitutions.
It validates the submitted reference allele against the pinned transcript sequence
and returns ``None`` whenever the requested position is outside a covered exon.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Mapping

from backend.domain.small_variant_projection import ProjectedAllele


_UTR_SUBSTITUTION = re.compile(
    r"^c\.(?P<position>-\d+|\*\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$",
    re.IGNORECASE,
)
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def is_exact_utr_substitution(c_notation: str) -> bool:
    """Return whether the HGVS edit is an exact 5' or 3' UTR SNV."""
    return _UTR_SUBSTITUTION.fullmatch((c_notation or "").strip()) is not None


@dataclass(frozen=True)
class TranscriptExonAlignment:
    transcript_start: int
    transcript_end: int
    genomic_start: int
    genomic_end: int


@dataclass(frozen=True)
class TranscriptAssemblyAlignment:
    accession: str
    chrom: str
    strand: int
    exons: tuple[TranscriptExonAlignment, ...]


@dataclass(frozen=True)
class TranscriptAlignment:
    gene: str
    transcript: str
    transcript_length: int
    cds_start_index: int
    cds_end_index: int
    sequence_sha256: str
    assemblies: Mapping[str, TranscriptAssemblyAlignment]


def _required_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    return value


def parse_transcript_alignment(
    gene: str,
    record: Mapping[str, Any],
) -> TranscriptAlignment:
    """Validate one immutable alignment snapshot record."""
    transcript = str(record.get("transcript") or "")
    transcript_length = _required_int(
        record.get("transcript_length"), f"{gene} transcript_length"
    )
    cds_start_index = _required_int(
        record.get("cds_start_index"), f"{gene} cds_start_index"
    )
    cds_end_index = _required_int(
        record.get("cds_end_index"), f"{gene} cds_end_index"
    )
    if not (
        transcript
        and transcript_length > 0
        and 0 <= cds_start_index < cds_end_index <= transcript_length
    ):
        raise ValueError(f"{gene} transcript alignment metadata is invalid")

    raw_assemblies = record.get("assemblies")
    if not isinstance(raw_assemblies, Mapping) or set(raw_assemblies) != {
        "GRCh37",
        "GRCh38",
    }:
        raise ValueError(f"{gene} must define GRCh37 and GRCh38 alignments")

    assemblies: dict[str, TranscriptAssemblyAlignment] = {}
    for assembly, raw_alignment in raw_assemblies.items():
        if not isinstance(raw_alignment, Mapping):
            raise ValueError(f"{gene}/{assembly} alignment is invalid")
        strand_symbol = str(raw_alignment.get("strand") or "")
        strand = 1 if strand_symbol == "+" else -1 if strand_symbol == "-" else 0
        accession = str(raw_alignment.get("accession") or "")
        chrom = str(raw_alignment.get("chrom") or "").removeprefix("chr")
        raw_exons = raw_alignment.get("exons")
        if not accession or not chrom or strand == 0 or not isinstance(raw_exons, list):
            raise ValueError(f"{gene}/{assembly} alignment metadata is invalid")

        exons: list[TranscriptExonAlignment] = []
        expected_transcript_start = 0
        for index, raw_exon in enumerate(raw_exons, 1):
            if not isinstance(raw_exon, Mapping):
                raise ValueError(f"{gene}/{assembly} exon {index} is invalid")
            exon = TranscriptExonAlignment(
                transcript_start=_required_int(
                    raw_exon.get("transcript_start"),
                    f"{gene}/{assembly} exon {index} transcript_start",
                ),
                transcript_end=_required_int(
                    raw_exon.get("transcript_end"),
                    f"{gene}/{assembly} exon {index} transcript_end",
                ),
                genomic_start=_required_int(
                    raw_exon.get("genomic_start"),
                    f"{gene}/{assembly} exon {index} genomic_start",
                ),
                genomic_end=_required_int(
                    raw_exon.get("genomic_end"),
                    f"{gene}/{assembly} exon {index} genomic_end",
                ),
            )
            if (
                exon.transcript_start != expected_transcript_start
                or exon.transcript_end <= exon.transcript_start
                or exon.genomic_start < 1
                or exon.genomic_end < exon.genomic_start
                or exon.transcript_end - exon.transcript_start
                != exon.genomic_end - exon.genomic_start + 1
            ):
                raise ValueError(f"{gene}/{assembly} exon {index} is inconsistent")
            expected_transcript_start = exon.transcript_end
            exons.append(exon)
        if expected_transcript_start != transcript_length:
            raise ValueError(
                f"{gene}/{assembly} alignment length does not match the transcript"
            )
        assemblies[str(assembly)] = TranscriptAssemblyAlignment(
            accession=accession,
            chrom=chrom,
            strand=strand,
            exons=tuple(exons),
        )

    return TranscriptAlignment(
        gene=gene,
        transcript=transcript,
        transcript_length=transcript_length,
        cds_start_index=cds_start_index,
        cds_end_index=cds_end_index,
        sequence_sha256=str(record.get("sequence_sha256") or "").lower(),
        assemblies=assemblies,
    )


def utr_transcript_index(
    c_notation: str,
    *,
    cds_start_index: int,
    cds_end_index: int,
) -> tuple[int, str, str] | None:
    """Return the zero-based transcript index and submitted alleles."""
    match = _UTR_SUBSTITUTION.fullmatch((c_notation or "").strip())
    if not match:
        return None
    position = match.group("position")
    distance = int(position[1:])
    if distance < 1:
        return None
    transcript_index = (
        cds_end_index + distance - 1
        if position.startswith("*")
        else cds_start_index - distance
    )
    return transcript_index, match.group("ref").upper(), match.group("alt").upper()


def project_utr_substitution(
    c_notation: str,
    *,
    transcript_sequence: str,
    transcript_alignment: TranscriptAlignment,
    assembly: str,
) -> ProjectedAllele | None:
    """Project a 5' or 3' UTR SNV without guessing across an alignment gap."""
    parsed = utr_transcript_index(
        c_notation,
        cds_start_index=transcript_alignment.cds_start_index,
        cds_end_index=transcript_alignment.cds_end_index,
    )
    if parsed is None:
        return None
    transcript_index, ref, alt = parsed
    sequence = transcript_sequence.upper()
    if (
        len(sequence) != transcript_alignment.transcript_length
        or transcript_index < 0
        or transcript_index >= len(sequence)
        or sequence[transcript_index] != ref
        or ref == alt
    ):
        return None
    alignment = transcript_alignment.assemblies.get(assembly)
    if alignment is None:
        return None
    exon = next(
        (
            item
            for item in alignment.exons
            if item.transcript_start <= transcript_index < item.transcript_end
        ),
        None,
    )
    if exon is None:
        return None
    offset = transcript_index - exon.transcript_start
    genomic_position = (
        exon.genomic_start + offset
        if alignment.strand == 1
        else exon.genomic_end - offset
    )
    genomic_ref = ref if alignment.strand == 1 else ref.translate(_COMPLEMENT)
    genomic_alt = alt if alignment.strand == 1 else alt.translate(_COMPLEMENT)
    return ProjectedAllele(
        chrom=alignment.chrom,
        pos=genomic_position,
        ref=genomic_ref,
        alt=genomic_alt,
    )
