import re
from typing import Optional


def normalize_protein_notation(value: Optional[str]) -> str:
    notation = (value or "").strip()
    if not notation:
        return ""
    if notation.startswith("(p.") and notation.endswith(")"):
        notation = notation[1:-1]
    if notation.startswith("p.") and not notation.startswith("p.("):
        notation = f"p.({notation[2:]})"
    # HGVS permits both Ter and * for a termination codon. Keep one canonical
    # representation so equivalent user and source notations compare equal.
    notation = re.sub(r"\*(?=\d|\))", "Ter", notation)
    return notation


def split_combined_hgvs(c_notation: str, p_notation: Optional[str] = None) -> tuple[str, str]:
    c_value = (c_notation or "").strip()
    p_value = normalize_protein_notation(p_notation)
    if p_value:
        return c_value, p_value

    match = re.match(r"^(c\.\S+)\s+(\(?p\..+\)?)$", c_value, re.IGNORECASE)
    if not match:
        return c_value, ""
    return match.group(1), normalize_protein_notation(match.group(2))
