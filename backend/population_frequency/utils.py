"""Small value and coordinate helpers shared by population-frequency modules."""

from __future__ import annotations

from typing import Any


def as_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def as_int(value: Any) -> int | None:
    try:
        if value is None or value == "":
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def strip_chr(chrom: Any) -> str:
    return str(chrom).replace("chr", "", 1)


def add_chr(chrom: Any) -> str:
    value = str(chrom)
    return value if value.startswith("chr") else "chr" + value


def coordinate_value(coords: Any, key: str) -> Any:
    return coords.get(key) if isinstance(coords, dict) else getattr(coords, key)


def variant_id_from_coords(
    coords: Any | None,
    with_chr: bool | None = None,
) -> str | None:
    if not coords:
        return None
    try:
        chrom = coordinate_value(coords, "chrom")
        pos = coordinate_value(coords, "pos")
        ref = coordinate_value(coords, "ref")
        alt = coordinate_value(coords, "alt")
        if with_chr is True:
            chrom = add_chr(chrom)
        elif with_chr is False:
            chrom = strip_chr(chrom)
        return f"{chrom}-{pos}-{ref}-{alt}"
    except (AttributeError, KeyError, TypeError):
        return None


def reference_span_from_coords(
    coords: Any | None,
) -> tuple[str | None, int | None, int | None]:
    """Return the genomic span represented by the VCF REF allele."""
    if not coords:
        return None, None, None
    try:
        chrom = coordinate_value(coords, "chrom")
        pos = int(coordinate_value(coords, "pos"))
        ref = coordinate_value(coords, "ref")
        ref_length = max(len(str(ref or "")), 1)
        return str(chrom), pos, pos + ref_length - 1
    except (AttributeError, KeyError, TypeError, ValueError):
        return None, None, None
