from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse

from backend.modules.utils import (
    get_amino_acid_position,
    is_in_functional_domain,
)

def evaluate_bp1(
    gene: str,
    variant_type: str,
    p_notation: str,
    spliceai_score: Optional[float] = None
) -> Dict:
    """
    Evaluate BP1 criterion (variant outside functional domain).

    BP1_Strong: silent/missense/in-frame variants outside functional domain
                AND no splicing prediction (SpliceAI <= 0.1)

    NOTE: spliceai_score MUST be passed explicitly. A default of 0 would
    silently treat "score not available" as "no splicing predicted",
    which is wrong - BP1 needs a confirmed low score.
    """
    result = {
        "applies": False,
        "strength": None,
        "points": 0,
        "reason": ""
    }

    # BP1 only applies to certain variant types
    # BP1_Strong applies to missense, synonymous, AND inframe insertion/deletion/delins
    # outside functional domain with SpliceAI <= 0.1
    # Source: ENIGMA VCEP v1.2 - BP1 criteria specification
    applicable_types = [
        "missense",
        "synonymous", "silent",  # silent is an alias for synonymous
        "inframe_deletion", "inframe_insertion", "inframe_delins", "delins"
    ]
    if variant_type not in applicable_types:
        result["reason"] = f"BP1 not applicable for {variant_type} variants"
        return result

    # check splicing prediction
    # None means score unavailable - cannot confirm no splice effect -> BP1 does not apply
    if spliceai_score is None:
        result["reason"] = "SpliceAI score not available - cannot confirm no splice effect, BP1 not applied"
        return result
    if spliceai_score > 0.1:
        result["reason"] = f"SpliceAI score {spliceai_score:.3f} > 0.1 - possible splicing effect"
        return result

    # check if in functional domain
    aa_pos = get_amino_acid_position(p_notation)
    if aa_pos is None:
        result["reason"] = "Could not determine amino acid position"
        return result

    in_domain, domain_name = is_in_functional_domain(gene, aa_pos)

    if in_domain:
        result["reason"] = f"Variant at aa {aa_pos} is inside {domain_name} domain"
        return result

    # variant is outside functional domain and no splicing effect
    result["applies"] = True
    result["strength"] = "Strong"
    result["points"] = -4  # benign evidence
    result["reason"] = f"Variant at aa {aa_pos} is outside functional domains, no splicing predicted"

    return result

if __name__ == "__main__":
    # test it - now we must pass spliceai_score explicitly
    print("Testing BP1 evaluation:")
    print(f"  Missense at aa 170 (outside domain), SpliceAI=0.03: "
          f"{evaluate_bp1('BRCA1', 'missense', 'p.(Arg170Gln)', spliceai_score=0.03)}")
    print(f"  Missense at aa 170 (outside domain), SpliceAI=None:  "
          f"{evaluate_bp1('BRCA1', 'missense', 'p.(Arg170Gln)', spliceai_score=None)}")
    print(f"  Missense at aa 50 (in RING domain), SpliceAI=0.02:   "
          f"{evaluate_bp1('BRCA1', 'missense', 'p.(Arg50Gln)', spliceai_score=0.02)}")
    print(f"  Missense at aa 1700 (in BRCT domain), SpliceAI=0.02: "
          f"{evaluate_bp1('BRCA1', 'missense', 'p.(Arg1700Gln)', spliceai_score=0.02)}")
