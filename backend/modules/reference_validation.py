"""Reference-transcript allele validation backed by versioned local datasets."""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from backend.lookups.precomputed import load_classification_snapshot_index


INTRONIC_COORDINATES_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "coordinates"
    / "brca_intronic_snv_coordinates.json"
)
_SNV_RE = re.compile(r"^c\.([0-9*+-]+)([ACGT])>([ACGT])$", re.IGNORECASE)


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


def validate_reference_allele(gene: str, c_notation: str) -> None:
    """Reject an SNV whose stated reference disagrees with the BRCA transcript."""
    match = _SNV_RE.fullmatch((c_notation or "").strip())
    if not match:
        return
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
