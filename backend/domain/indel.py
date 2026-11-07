"""Pure HGVS indel-shape and size helpers shared across backend layers."""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Any

_OPERATION = re.compile(r"(delins|del|dup|ins)", re.IGNORECASE)
_LITERAL_SEQUENCE = re.compile(r"[ACGT]+", re.IGNORECASE)


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


def indel_operation(c_notation: str | None) -> str | None:
    """Return the HGVS indel operation, without inferring event size."""
    match = _OPERATION.search(c_notation or "")
    return match.group(1).lower() if match else None


def is_indel_allele(c_notation: str | None) -> bool:
    return indel_operation(c_notation) is not None


def assess_indel_size_for_transcript(transcript: str, c_notation: str) -> dict[str, Any]:
    """Measure an exactly described c. HGVS indel without sequence lookups."""
    operation = indel_operation(c_notation)
    if operation is None:
        return {
            "is_indel": False,
            "operation": None,
            "status": "not_indel",
            "affected_bp": None,
        }

    if not transcript:
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reason": "no approved reference transcript is configured for the gene",
        }
    variant = _parse_variant(f"{transcript}:{c_notation}")
    if variant is None:
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reason": "the c. HGVS does not define an exactly measurable sequence change",
        }

    start = variant.posedit.pos.start
    end = variant.posedit.pos.end
    if start.base is None or end.base is None:
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reason": "uncertain breakpoints do not define an exact event size",
        }

    edit = variant.posedit.edit
    explicit_reference = getattr(edit, "ref", None)
    explicit_reference_bp = (
        len(explicit_reference)
        if isinstance(explicit_reference, str)
        and _LITERAL_SEQUENCE.fullmatch(explicit_reference)
        else None
    )
    inserted = getattr(edit, "alt", None)
    inserted_bp = (
        len(inserted)
        if isinstance(inserted, str) and _LITERAL_SEQUENCE.fullmatch(inserted)
        else None
    )

    start_base = int(start.base)
    end_base = int(end.base)
    start_offset = int(start.offset or 0)
    end_offset = int(end.offset or 0)
    reference_span: int | None = None
    if start_base == end_base and end_offset >= start_offset:
        # This covers exact intronic ranges such as c.317-19_317-12del as
        # well as a single intronic position. Both endpoints share one HGVS
        # anchor, so their inclusive distance is unambiguous.
        reference_span = end_offset - start_offset + 1
    elif start_offset == 0 and end_offset == 0 and end_base >= start_base:
        reference_span = end_base - start_base + 1
    elif explicit_reference_bp is not None:
        # An explicitly stated deleted/duplicated sequence establishes the
        # event size even when transcript structure is needed to measure the
        # endpoints (for example an intron-to-exon event).
        reference_span = explicit_reference_bp

    if (
        reference_span is not None
        and explicit_reference_bp is not None
        and reference_span != explicit_reference_bp
    ):
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reference_span_bp": reference_span,
            "explicit_reference_bp": explicit_reference_bp,
            "inserted_bp": inserted_bp,
            "reason": (
                "the explicit reference sequence length conflicts with the "
                "HGVS breakpoint span"
            ),
        }

    if reference_span is None and operation != "ins":
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "inserted_bp": inserted_bp,
            "reason": (
                "the breakpoints require transcript structure and no explicit "
                "reference sequence establishes the event size"
            ),
        }

    affected_bp: int | None
    if operation in {"del", "dup"}:
        affected_bp = reference_span
    elif operation == "ins":
        affected_bp = inserted_bp
    else:
        if inserted_bp is not None:
            return {
                "is_indel": True,
                "operation": operation,
                "status": "exact_components",
                "affected_bp": None,
                "reference_span_bp": reference_span,
                "inserted_bp": inserted_bp,
            }
        affected_bp = None

    if affected_bp is None:
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reference_span_bp": reference_span,
            "inserted_bp": inserted_bp,
            "reason": "the inserted sequence length is not explicit in the c. HGVS",
        }
    return {
        "is_indel": True,
        "operation": operation,
        "status": "exact",
        "affected_bp": affected_bp,
        "reference_span_bp": reference_span,
        "inserted_bp": inserted_bp,
    }
