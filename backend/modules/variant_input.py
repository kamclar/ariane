"""Fail-closed normalization of user variant descriptions."""
from __future__ import annotations

import re
from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Optional

from backend.config import TRANSCRIPTS
from backend.lookups.indels import load_indel_snapshot, lookup_indel_snapshot
from backend.lookups.precomputed import (
    load_classification_snapshot_index,
    load_classification_snapshot_metadata,
)
from backend.modules.hgvs import (
    normalize_protein_notation,
    protein_notations_compatible,
    split_combined_hgvs,
)


_ASSEMBLIES = {"GRCH37": "GRCh37", "GRCH38": "GRCh38"}
_TRANSCRIPT_C = re.compile(
    r"^(NM_\d+(?:\.\d+)?)\s*:\s*(c\..+)$", re.IGNORECASE
)
_LEADING_C_SEPARATOR = re.compile(r"^:\s*(c\..+)$", re.IGNORECASE)
_GENE_QUALIFIED_INPUT = re.compile(
    r"^((?!CHR[0-9XYM]+(?:\b|:))[A-Z][A-Z0-9-]*)"
    r"\s*(?::\s*|\s+)(.+)$",
    re.IGNORECASE,
)
_TRANSCRIPT_ACCESSION_GENE = {
    transcript.split(".", 1)[0]: gene for gene, transcript in TRANSCRIPTS.items()
}
_GENOMIC_PATTERNS = (
    re.compile(
        r"^(?:chr)?(?P<chrom>[0-9XYM]+)[:\s-]+(?P<pos>\d+)"
        r"[:\s-]+(?P<ref>[ACGT]+)[>:\/\s-]+(?P<alt>[ACGT]+)$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?:chr)?(?P<chrom>[0-9XYM]+):(?P<pos>\d+)"
        r"[:\s]+(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)$",
        re.IGNORECASE,
    ),
)


@dataclass(frozen=True)
class NormalizedVariantInput:
    gene: str
    submitted_notation: str
    c_notation: str
    p_notation: str
    reference_transcript: str
    normalization_source: str
    consequence_status: str = ""
    normalization_provenance: dict[str, str] | None = None
    protein_consequence_explanation: str = ""
    assembly: str = ""
    genomic_notation: str = ""


def normalize_assembly(value: Optional[str]) -> str:
    raw = (value or "").strip().upper()
    raw = raw.replace("HG19", "GRCH37").replace("HG38", "GRCH38")
    if not raw:
        return ""
    normalized = _ASSEMBLIES.get(raw)
    if not normalized:
        raise ValueError("Genome assembly must be GRCh37 or GRCh38")
    return normalized


def explicit_gene_from_notation(notation: str) -> Optional[str]:
    """Return an explicitly supplied active gene, rejecting unsupported identifiers."""
    submitted = (notation or "").strip()
    gene_from_prefix = None
    gene_match = _GENE_QUALIFIED_INPUT.fullmatch(submitted)
    if gene_match:
        gene_from_prefix = gene_match.group(1).upper()
        if gene_from_prefix not in TRANSCRIPTS:
            raise ValueError(
                f"No reference transcript is configured for gene {gene_from_prefix}"
            )
        submitted = gene_match.group(2).strip()

    gene_from_transcript = None
    transcript_match = _TRANSCRIPT_C.fullmatch(submitted)
    if transcript_match:
        transcript = transcript_match.group(1).upper()
        accession = transcript.split(".", 1)[0]
        gene_from_transcript = _TRANSCRIPT_ACCESSION_GENE.get(accession)

    if (
        gene_from_prefix
        and gene_from_transcript
        and gene_from_prefix != gene_from_transcript
    ):
        raise ValueError(
            f"Conflicting gene identifiers: {gene_from_prefix} does not match "
            f"transcript {transcript_match.group(1).upper()} ({gene_from_transcript})"
        )
    return gene_from_prefix or gene_from_transcript


def _parse_genomic(value: str) -> Optional[tuple[str, int, str, str]]:
    compact = re.sub(r"\s+", " ", value.strip())
    for pattern in _GENOMIC_PATTERNS:
        match = pattern.fullmatch(compact)
        if match:
            return (
                match.group("chrom").upper().removeprefix("CHR"),
                int(match.group("pos")),
                match.group("ref").upper(),
                match.group("alt").upper(),
            )
    return None


def _snapshot_coordinate(record: dict, assembly: str) -> Optional[tuple[str, int, str, str]]:
    value = record.get(assembly.lower())
    if isinstance(value, dict):
        return (
            str(value.get("chrom", "")).removeprefix("chr"),
            int(value["pos"]),
            str(value.get("ref", "")).upper(),
            str(value.get("alt", "")).upper(),
        )
    if isinstance(value, str):
        return _parse_genomic(value)
    return None


@lru_cache(maxsize=1)
def _genomic_reverse_index() -> dict[tuple[str, str, str, int, str, str], tuple[dict, ...]]:
    candidates: dict[tuple[str, str, str, int, str, str], list[dict]] = {}

    def add(record: dict, source: str) -> None:
        gene = str(record.get("gene", "")).upper()
        c_notation = str(
            record.get("canonical_c_notation") or record.get("c_notation") or ""
        )
        p_notation = normalize_protein_notation(str(record.get("p_notation") or ""))
        transcript = str(record.get("reference_transcript") or TRANSCRIPTS.get(gene, ""))
        if not (gene and c_notation and p_notation and transcript):
            return
        for assembly in ("GRCh37", "GRCh38"):
            coords = _snapshot_coordinate(record, assembly)
            if not coords:
                continue
            chrom, pos, ref, alt = coords
            key = (gene, assembly, chrom, pos, ref, alt)
            item = {
                "c_notation": c_notation,
                "p_notation": p_notation,
                "reference_transcript": transcript,
                "source": source,
            }
            if item not in candidates.setdefault(key, []):
                candidates[key].append(item)

    metadata = load_classification_snapshot_metadata()
    snv_source = "BRCA coding SNV snapshot"
    if metadata.get("created"):
        snv_source += f" ({metadata['created']})"
    for record in load_classification_snapshot_index().values():
        add(record, snv_source)
    indels, _aliases = load_indel_snapshot()
    for record in indels.values():
        add(record, "Normalized BRCA indel snapshot")
    return {key: tuple(values) for key, values in candidates.items()}


def _from_c_notation(
    gene: str, submitted: str, c_notation: str, supplied_p: str
) -> NormalizedVariantInput:
    from backend.modules.hgvs_engine import derive_protein_consequence
    from backend.modules.reference_validation import validate_reference_allele

    # Validate the stated transcript reference before attempting consequence
    # lookup. This also prevents a wrong-reference SNV from being accepted with
    # a user-supplied protein description.
    validate_reference_allele(gene, c_notation)

    consequence = derive_protein_consequence(gene, c_notation)
    canonical_c = consequence.canonical_c_notation
    engine_p = consequence.p_notation

    # Existing snapshots remain independent regression/cross-check sources.
    # They are not used as a hidden runtime fallback for sequence derivation.
    indel = lookup_indel_snapshot(gene, c_notation)
    if indel:
        snapshot_p = normalize_protein_notation(str(indel.get("p_notation") or ""))
        snapshot_source = "Normalized BRCA indel snapshot"
    else:
        record = load_classification_snapshot_index().get(f"{gene}:{canonical_c}")
        snapshot_p = normalize_protein_notation(
            str(record.get("p_notation") or "") if record else ""
        )
        snapshot_source = "BRCA coding SNV snapshot" if record else ""

    if (
        snapshot_p
        and snapshot_p != "p.?"
        and engine_p != "p.?"
        and not protein_notations_compatible(snapshot_p, engine_p)
    ):
        raise ValueError(
            f"Validated source conflict for {gene} {canonical_c}: local HGVS engine "
            f"gives {engine_p}, while {snapshot_source} gives {snapshot_p}. "
            "The variant was not classified."
        )

    normalized_supplied_p = normalize_protein_notation(supplied_p)
    if (
        engine_p != "p.?"
        and normalized_supplied_p
        and not protein_notations_compatible(normalized_supplied_p, engine_p)
    ):
        raise ValueError(
            f"Protein consequence mismatch for {gene} {canonical_c}: "
            f"{TRANSCRIPTS[gene]} gives {engine_p}, not {normalized_supplied_p}"
        )
    if engine_p == "p.?" and normalized_supplied_p and normalized_supplied_p != "p.?":
        raise ValueError(
            f"Protein consequence cannot be verified from {TRANSCRIPTS[gene]} for "
            f"{canonical_c}; supplied {normalized_supplied_p} was not accepted."
        )
    p_notation = engine_p
    protein_explanation = ""
    if p_notation == "p.?":
        protein_explanation = (
            "The protein consequence cannot be determined from DNA notation alone. "
            "The variant is outside the translated coding sequence, may affect splicing, "
            "or has uncertain breakpoints; transcript, RNA evidence, or breakpoint evidence "
            "is required."
        )

    return NormalizedVariantInput(
        gene=gene,
        submitted_notation=submitted,
        c_notation=canonical_c,
        p_notation=p_notation,
        reference_transcript=TRANSCRIPTS[gene],
        normalization_source="Local biocommons.hgvs engine with checksum-verified panel reference",
        consequence_status=consequence.consequence_status,
        normalization_provenance=consequence.provenance,
        protein_consequence_explanation=protein_explanation,
    )


def normalize_variant_input(
    gene: str,
    notation: str,
    *,
    assembly: Optional[str] = None,
    p_notation: Optional[str] = None,
) -> NormalizedVariantInput:
    gene = (gene or "").strip().upper()
    if gene not in TRANSCRIPTS:
        raise ValueError(f"No reference transcript is configured for gene {gene or '(missing)'}")
    submitted = (notation or "").strip()
    if not submitted:
        raise ValueError("Variant description is required")
    explicit_gene = explicit_gene_from_notation(submitted)
    if explicit_gene:
        gene = explicit_gene

    # A gene or canonical transcript explicitly included in the notation is
    # authoritative. A bare c. description continues to use the separately
    # selected gene. Never guess a gene from a coordinate number or sequence.
    leading_separator = _LEADING_C_SEPARATOR.fullmatch(submitted)
    if leading_separator:
        submitted = leading_separator.group(1)

    gene_match = _GENE_QUALIFIED_INPUT.fullmatch(submitted)
    if gene_match:
        submitted_gene = gene_match.group(1).upper()
        gene = submitted_gene
        submitted = gene_match.group(2).strip()

    transcript_match = _TRANSCRIPT_C.fullmatch(submitted)
    if transcript_match:
        transcript = transcript_match.group(1).upper()
        accession = transcript.split(".", 1)[0]
        transcript_gene = _TRANSCRIPT_ACCESSION_GENE.get(accession)
        if transcript_gene is None:
            raise ValueError(f"Unsupported reference transcript {transcript}")
        if gene_match and transcript_gene != gene:
            raise ValueError(
                f"Conflicting gene identifiers: {gene} does not match "
                f"transcript {transcript} ({transcript_gene})"
            )
        gene = transcript_gene
        if "." not in transcript:
            raise ValueError(
                f"Transcript accession must include an explicit version; "
                f"expected {TRANSCRIPTS[gene]} for {gene}"
            )
        if transcript != TRANSCRIPTS[gene]:
            raise ValueError(
                f"Transcript {transcript} does not match {gene}; expected {TRANSCRIPTS[gene]}"
            )
        submitted = transcript_match.group(2)

    c_notation, combined_p = split_combined_hgvs(submitted, p_notation)
    if c_notation.lower().startswith("c."):
        return _from_c_notation(gene, notation.strip(), c_notation, combined_p)

    genomic = _parse_genomic(submitted)
    if not genomic:
        raise ValueError(
            "Unrecognised variant description. Enter reference-transcript c. HGVS "
            "or a genomic variant such as chr17:43124096:T>G."
        )
    normalized_assembly = normalize_assembly(assembly)
    if not normalized_assembly:
        raise ValueError(
            "Genome assembly is required for genomic coordinates; select GRCh37 or GRCh38"
        )
    chrom, pos, ref, alt = genomic
    key = (gene, normalized_assembly, chrom, pos, ref, alt)
    matches = _genomic_reverse_index().get(key, ())
    if not matches:
        raise ValueError(
            f"No validated local {gene} mapping was found for "
            f"{normalized_assembly} chr{chrom}:{pos}:{ref}>{alt}. "
            "The variant was not classified."
        )
    canonical = {
        (item["c_notation"], item["p_notation"], item["reference_transcript"])
        for item in matches
    }
    if len(canonical) != 1:
        descriptions = ", ".join(
            f"{item['reference_transcript']}:{item['c_notation']} {item['p_notation']}"
            for item in matches
        )
        raise ValueError(
            f"Ambiguous genomic variant: multiple reference-transcript descriptions "
            f"match {normalized_assembly} chr{chrom}:{pos}:{ref}>{alt}: {descriptions}"
        )
    item = matches[0]
    normalized_c = _from_c_notation(gene, notation.strip(), item["c_notation"], "")
    provenance = dict(normalized_c.normalization_provenance or {})
    provenance["genomic_mapping_source"] = item["source"]
    return replace(
        normalized_c,
        normalization_source=(
            f"{normalized_c.normalization_source}; genomic mapping: {item['source']}"
        ),
        normalization_provenance=provenance,
        assembly=normalized_assembly,
        genomic_notation=f"chr{chrom}:{pos}:{ref}>{alt}",
    )
