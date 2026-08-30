"""Structured indel-size assessment for the ENIGMA Appendix G PM2 path."""

from __future__ import annotations

import re
from typing import Any

from hgvs.exceptions import HGVSParseError
from hgvs.parser import Parser

from backend.config import TRANSCRIPTS


_PARSER = Parser()
_OPERATION = re.compile(r"(delins|del|dup|ins)", re.IGNORECASE)
_LITERAL_SEQUENCE = re.compile(r"[ACGT]+", re.IGNORECASE)


def indel_operation(c_notation: str | None) -> str | None:
    """Return the HGVS indel operation, without inferring event size."""
    match = _OPERATION.search(c_notation or "")
    return match.group(1).lower() if match else None


def is_indel_allele(c_notation: str | None) -> bool:
    return indel_operation(c_notation) is not None


def assess_indel_size(gene: str, c_notation: str) -> dict[str, Any]:
    """Measure an exactly described c. HGVS indel without sequence lookups.

    For delins alleles, the deleted and inserted lengths are retained
    separately. If they lie on opposite sides of the Appendix G boundary, this
    parser does not invent a single event-size convention. Symbolic, intronic,
    or uncertain descriptions remain unknown and require structural review.
    """
    operation = indel_operation(c_notation)
    if operation is None:
        return {
            "is_indel": False,
            "operation": None,
            "status": "not_indel",
            "affected_bp": None,
        }

    transcript = TRANSCRIPTS.get((gene or "").upper())
    if not transcript:
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reason": "no approved reference transcript is configured for the gene",
        }
    try:
        variant = _PARSER.parse_hgvs_variant(f"{transcript}:{c_notation}")
    except (HGVSParseError, ValueError):
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reason": "the c. HGVS does not define an exactly measurable sequence change",
        }

    start = variant.posedit.pos.start
    end = variant.posedit.pos.end
    if (
        start.base is None
        or end.base is None
        or start.offset != 0
        or end.offset != 0
        or int(end.base) < int(start.base)
    ):
        return {
            "is_indel": True,
            "operation": operation,
            "status": "unknown",
            "affected_bp": None,
            "reason": "intronic or uncertain breakpoints do not define an exact event size",
        }

    reference_span = int(end.base) - int(start.base) + 1
    edit = variant.posedit.edit
    inserted = getattr(edit, "alt", None)
    inserted_bp = (
        len(inserted)
        if isinstance(inserted, str) and _LITERAL_SEQUENCE.fullmatch(inserted)
        else None
    )

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
            "reason": (
                "the inserted sequence length is not explicit in the c. HGVS"
            ),
        }
    return {
        "is_indel": True,
        "operation": operation,
        "status": "exact",
        "affected_bp": affected_bp,
        "reference_span_bp": reference_span,
        "inserted_bp": inserted_bp,
    }
