"""Canonical display order for ACMG/AMP and ENIGMA evidence codes."""

import re
from typing import Any, Dict, Iterable, Tuple


_GROUP_ORDER = {
    "PVS": 0,
    "PS": 1,
    "PM": 2,
    "PP": 3,
    "BA": 4,
    "BS": 5,
    "BP": 6,
}


def criterion_sort_key(code: str) -> Tuple[int, int, int, str]:
    """Sort PVS, PS, PM, PP, BA, BS and BP codes as listed in guidelines.

    Internal qualifiers such as ``PM5_PTC`` and ``BS1_Supporting`` remain next
    to their base criterion. Unknown codes are retained at the end.
    """
    normalized = str(code or "").strip().upper()
    match = re.match(r"^(PVS|PS|PM|PP|BA|BS|BP)(\d+)", normalized)
    if not match:
        return (99, 99, 1, normalized)
    prefix, number = match.groups()
    base = f"{prefix}{number}"
    return (
        _GROUP_ORDER[prefix],
        int(number),
        0 if normalized == base else 1,
        normalized,
    )


def sorted_criterion_items(
    criteria: Dict[str, Any],
) -> Iterable[Tuple[str, Any]]:
    """Return mapping items in canonical criterion order."""
    return sorted(criteria.items(), key=lambda item: criterion_sort_key(item[0]))
