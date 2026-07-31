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


def protein_notations_compatible(supplied: Optional[str], canonical: Optional[str]) -> bool:
    """Return whether a supplied p. description agrees with the canonical one.

    HGVS frameshift descriptions are often reported in the legacy abbreviated
    form ``p.(Cys1225fs)``. Accept that form only when its original amino acid
    and position match a validated full consequence such as
    ``p.(Cys1225SerfsTer10)``. The canonical notation remains the output.
    """
    supplied_normalized = normalize_protein_notation(supplied)
    canonical_normalized = normalize_protein_notation(canonical)
    if supplied_normalized == canonical_normalized:
        return True

    abbreviated = re.fullmatch(
        r"p\.\(([A-Z][a-z]{2})(\d+)fs\)", supplied_normalized
    )
    full = re.fullmatch(
        r"p\.\(([A-Z][a-z]{2})(\d+)[A-Z][a-z]{2}fs(?:Ter\d+|\?)?\)",
        canonical_normalized,
    )
    return bool(abbreviated and full and abbreviated.groups() == full.groups())


def split_combined_hgvs(c_notation: str, p_notation: Optional[str] = None) -> tuple[str, str]:
    c_value = (c_notation or "").strip()
    p_value = normalize_protein_notation(p_notation)
    if p_value:
        return c_value, p_value

    match = re.match(r"^(c\.\S+)\s+(\(?p\..+\)?)$", c_value, re.IGNORECASE)
    if not match:
        return c_value, ""
    return match.group(1), normalize_protein_notation(match.group(2))
