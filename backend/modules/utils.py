# functional domains as defined by ENIGMA
# these are amino acid positions

from typing import Optional
import re
from backend.gene_policy import functional_domains


_PROTEIN_INTERVAL_RE = re.compile(
    r"p\.\(?[A-Za-z*]{1,3}(\d+)(?:_[A-Za-z*]{1,3}(\d+))?"
)


def get_amino_acid_interval(p_notation: str) -> Optional[tuple[int, int]]:
    """Return the complete explicitly affected amino-acid interval.

    The second number is read only from an HGVS residue range. Numbers in a
    downstream frameshift termination, for example ``fsTer10``, are not
    mistaken for absolute protein coordinates.
    """
    match = _PROTEIN_INTERVAL_RE.search(p_notation or "")
    if not match:
        return None
    start = int(match.group(1))
    end = int(match.group(2) or start)
    return (min(start, end), max(start, end))


def get_amino_acid_position(p_notation: str) -> Optional[int]:
    """
    Extract the amino acid position from protein notation.
    Examples:
        p.(Arg170Gln) -> 170
        p.(Cys1225fs) -> 1225
        p.(Val2050del) -> 2050

    For an explicit multi-residue interval this returns its first position.
    Domain-overlap decisions must use ``get_amino_acid_interval`` instead.
    """
    interval = get_amino_acid_interval(p_notation)
    return interval[0] if interval else None


def overlapping_functional_domains(
    gene: str,
    amino_acid_interval: tuple[int, int],
) -> tuple[str, ...]:
    """Return every ENIGMA domain overlapping any part of the interval."""
    start, end = amino_acid_interval
    return tuple(
        domain_name
        for domain_name, (domain_start, domain_end) in functional_domains(gene).items()
        if start <= domain_end and end >= domain_start
    )

def is_in_functional_domain(gene: str, aa_position: int) -> tuple:
    """
    Check if an amino acid position is in a functional domain.
    Returns (is_in_domain, domain_name)
    """
    domains = overlapping_functional_domains(gene, (aa_position, aa_position))
    return (bool(domains), domains[0] if domains else None)


def get_cds_position_from_c_notation(c_notation: str) -> Optional[int]:
    """
    Extract CDS position from coding variant notation.
    c.509G>A -> 509, c.628C>T -> 628
    Does not handle intronic (c.8953+2T>C) - returns None for those.
    """
    import re
    match = re.match(r'c\.(-?\d+)', c_notation)
    if match:
        return int(match.group(1))
    return None


def get_intron_offset_from_c_notation(c_notation: str) -> Optional[tuple]:
    """
    Extract intron offset from intronic notation.
    c.8953+2T>C -> (8953, +2)
    c.794-1G>A -> (794, -1)
    """
    import re
    match = re.match(r'c\.(-?\d+)([+-])(\d+)', c_notation)
    if match:
        pos = int(match.group(1))
        sign = 1 if match.group(2) == '+' else -1
        offset = sign * int(match.group(3))
        return (pos, offset)
    return None
