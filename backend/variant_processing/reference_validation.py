"""Reference-transcript allele validation backed by versioned local datasets."""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from backend.domain.indel import indel_operation
from backend.domain.small_variant_projection import (
    AssemblyBaseMap,
    build_assembly_base_maps,
    expected_transcript_reference,
)
from backend.policy.catalog import TRANSCRIPTS
from backend.reference_data.classification_snapshot import load_classification_snapshot_index


INTRONIC_COORDINATES_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "coordinates"
    / "brca_intronic_snv_coordinates.json"
)
_SNV_RE = re.compile(r"^c\.([0-9*+-]+)([ACGT])>([ACGT])$", re.IGNORECASE)
_DNA_RE = re.compile(r"[ACGT]+", re.IGNORECASE)


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


@lru_cache(maxsize=1)
def _reference_bases() -> dict[tuple[str, str], str]:
    """Build transcript-position reference bases from generated local maps."""
    bases: dict[tuple[str, str], str] = {}

    def add_keys(keys) -> None:
        for key in keys:
            try:
                gene, notation = key.split(":", 1)
            except ValueError:
                continue
            match = _SNV_RE.fullmatch(notation)
            if not match:
                continue
            position, reference = match.group(1), match.group(2).upper()
            map_key = (gene.upper(), position)
            previous = bases.get(map_key)
            if previous is not None and previous != reference:
                raise RuntimeError(
                    f"Conflicting reference bases in local datasets for {gene} c.{position}: "
                    f"{previous} and {reference}"
                )
            bases[map_key] = reference

    add_keys(load_classification_snapshot_index().keys())
    if not INTRONIC_COORDINATES_PATH.exists():
        raise RuntimeError(
            f"Reference validation dataset is missing: {INTRONIC_COORDINATES_PATH}"
        )
    try:
        with INTRONIC_COORDINATES_PATH.open(encoding="utf-8") as handle:
            add_keys(json.load(handle).keys())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Reference validation dataset could not be loaded: {type(exc).__name__}: {exc}"
        ) from exc
    return bases


@lru_cache(maxsize=1)
def _grch38_indel_base_maps() -> dict[tuple[str, str], AssemblyBaseMap]:
    try:
        with INTRONIC_COORDINATES_PATH.open(encoding="utf-8") as handle:
            intronic = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Reference validation dataset could not be loaded: {type(exc).__name__}: {exc}"
        ) from exc
    if not isinstance(intronic, dict):
        raise RuntimeError("Reference validation dataset must contain an object")
    return build_assembly_base_maps(
        load_classification_snapshot_index(),
        {"brca_intronic_snv_map": intronic},
    )


def validate_reference_allele(gene: str, c_notation: str) -> None:
    """Reject an explicit reference allele that disagrees with local references."""
    match = _SNV_RE.fullmatch((c_notation or "").strip())
    if match:
        position, supplied = match.group(1), match.group(2).upper()
        expected = _reference_bases().get((gene.strip().upper(), position))
        if expected is None:
            raise ValueError(
                f"Reference allele could not be verified for {gene} c.{position} "
                "using the installed reference-transcript datasets; classification was not run."
            )
        if supplied != expected:
            raise ValueError(
                f"Reference allele mismatch for {gene} {c_notation}: position c.{position} "
                f"is {expected} in the configured reference transcript, not {supplied}."
            )
        return

    operation = indel_operation(c_notation)
    if operation not in {"del", "dup", "delins"}:
        return
    normalized_gene = gene.strip().upper()
    transcript = TRANSCRIPTS.get(normalized_gene)
    if not transcript:
        return
    variant = _parse_variant(f"{transcript}:{c_notation}")
    if variant is None:
        return
    supplied_reference = getattr(variant.posedit.edit, "ref", None)
    if not (
        isinstance(supplied_reference, str)
        and _DNA_RE.fullmatch(supplied_reference)
    ):
        return
    assembly_map = _grch38_indel_base_maps().get((normalized_gene, "GRCh38"))
    if assembly_map is None:
        return
    expected_reference = expected_transcript_reference(
        transcript, c_notation, assembly_map
    )
    if expected_reference is not None and supplied_reference.upper() != expected_reference:
        raise ValueError(
            f"Reference allele mismatch for {gene} {c_notation}: the explicit "
            f"{operation} sequence has {len(supplied_reference)} bp, but the "
            f"configured transcript interval contains {len(expected_reference)} bp "
            f"({expected_reference}). Classification was not run."
        )
