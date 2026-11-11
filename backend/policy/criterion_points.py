"""Point values for ACMG/AMP evidence strengths."""

from __future__ import annotations


_POINTS_BY_STRENGTH = {
    "Supporting": 1,
    "Moderate": 2,
    "Strong": 4,
    "Very Strong": 8,
}


def points_for_strength(strength: str | None, *, benign: bool = False) -> int:
    if strength is None:
        return 0
    points = _POINTS_BY_STRENGTH[strength]
    return -points if benign else points
