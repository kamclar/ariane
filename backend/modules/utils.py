# functional domains as defined by ENIGMA
# these are amino acid positions

from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse
from backend.gene_policy import functional_domains

def get_amino_acid_position(p_notation: str) -> Optional[int]:
    """
    Extract the amino acid position from protein notation.
    Examples:
        p.(Arg170Gln) -> 170
        p.(Cys1225fs) -> 1225
        p.(Val2050del) -> 2050

    This is a bit crude but works for most cases..
    """
    import re
    # look for pattern like Arg170 or Cys1225
    match = re.search(r'[A-Z][a-z]{2}(\d+)', p_notation)
    if match:
        return int(match.group(1))
    return None

def is_in_functional_domain(gene: str, aa_position: int) -> tuple:
    """
    Check if an amino acid position is in a functional domain.
    Returns (is_in_domain, domain_name)
    """
    for domain_name, (start, end) in functional_domains(gene).items():
        if start <= aa_position <= end:
            return (True, domain_name)

    return (False, None)


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
