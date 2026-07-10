# BayesDel_noAF thresholds are gene-specific per ENIGMA VCEP v1.2
# BRCA1: PP3 >= 0.28, BP4 <= 0.15
# BRCA2: PP3 >= 0.30, BP4 <= 0.18
# Using generic thresholds for both genes would be wrong.

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
from backend.lookups.spliceai import (
    normalize_variant_type,
    spliceai_predicts_splice_effect,
    variant_type_allows_spliceai_pp3,
)

BAYESDEL_THRESHOLDS = {
    "BRCA1": {"pp3": 0.28, "bp4": 0.15},
    "BRCA2": {"pp3": 0.30, "bp4": 0.18},
}


def evaluate_pp3_bp4(
    gene: str,
    variant_type: str,
    p_notation: str,
    bayesdel_score: Optional[float] = None,
    spliceai_score: Optional[float] = None
) -> Dict:
    # evaluate PP3 and BP4 using BayesDel_noAF and SpliceAI
    #
    # v1.5.2 SpliceAI guardrail:
    # PP3 from SpliceAI is restricted to silent/synonymous, missense/in-frame,
    # and intronic variants. It is NOT applied to nonsense/PTC, frameshift,
    # exon-level deletion, or canonical splice-site variants, where PVS1/RNA
    # logic evaluates the loss-of-function/splicing mechanism.
    #
    # BP4 requires confirmed low SpliceAI (<= 0.1), never missing SpliceAI.
    # For missense/in-frame variants in a functional domain it also requires
    # BayesDel_noAF <= the gene-specific BP4 threshold.
    results = {}

    vtype = normalize_variant_type(variant_type)

    thresholds = BAYESDEL_THRESHOLDS.get(gene, {"pp3": 0.28, "bp4": 0.15})
    pp3_threshold = thresholds["pp3"]
    bp4_threshold = thresholds["bp4"]

    aa_pos = get_amino_acid_position(p_notation)
    in_domain = False
    domain_name = None
    if aa_pos:
        in_domain, domain_name = is_in_functional_domain(gene, aa_pos)

    protein_prediction_types = {
        "missense",
        "inframe_deletion", "inframe_insertion", "inframe_delins", "delins",
    }

    # ---- PP3 ----
    # SpliceAI branch. This must not be "any variant type".
    if (
        variant_type_allows_spliceai_pp3(vtype)
        and spliceai_predicts_splice_effect(spliceai_score)
    ):
        results["PP3"] = {
            "applies": True,
            "strength": "Supporting",
            "points": 1,
            "reason": f"SpliceAI {spliceai_score:.3f} >= 0.2 - predicted splice effect"
        }

    # BayesDel branch: missense/in-frame in functional domain only.
    # Only if SpliceAI did not already trigger PP3 (same criterion cannot stack).
    elif vtype in protein_prediction_types:
        if in_domain and bayesdel_score is not None and bayesdel_score >= pp3_threshold:
            results["PP3"] = {
                "applies": True,
                "strength": "Supporting",
                "points": 1,
                "reason": (
                    f"BayesDel_noAF {bayesdel_score:.3f} >= {pp3_threshold} ({gene}-specific threshold) "
                    f"in functional domain ({domain_name})"
                )
            }
        elif in_domain and bayesdel_score is None:
            results["PP3"] = {
                "applies": False,
                "reason": f"BayesDel_noAF not available for this variant (in domain: {domain_name})"
            }
        elif not in_domain and bayesdel_score is not None and bayesdel_score >= pp3_threshold:
            results["PP3"] = {
                "applies": False,
                "reason": (
                    f"BayesDel_noAF {bayesdel_score:.3f} >= {pp3_threshold} ({gene}-specific threshold) "
                    f"but variant is outside functional domains: PP3 from BayesDel applies only inside "
                    f"a clinically important functional domain per ENIGMA VCEP v1.2. "
                    f"For variants outside domain, BP1 applies if SpliceAI <= 0.1."
                )
            }

    # ---- BP4 ----
    if vtype in protein_prediction_types:
        if in_domain:
            if bayesdel_score is not None and spliceai_score is not None:
                if bayesdel_score <= bp4_threshold and spliceai_score <= 0.1:
                    results["BP4"] = {
                        "applies": True,
                        "strength": "Supporting",
                        "points": -1,
                        "reason": (
                            f"BayesDel_noAF {bayesdel_score:.3f} <= {bp4_threshold} ({gene}-specific threshold) AND "
                            f"SpliceAI {spliceai_score:.3f} <= 0.1 "
                            f"(in domain: {domain_name})"
                        )
                    }
            elif bayesdel_score is None:
                results["BP4"] = {
                    "applies": False,
                    "reason": "BayesDel_noAF not available - cannot evaluate BP4 for missense/in-frame variant in domain"
                }
            elif bayesdel_score is not None and spliceai_score is None and bayesdel_score <= bp4_threshold:
                results["BP4"] = {
                    "applies": False,
                    "reason": (
                        f"BayesDel_noAF {bayesdel_score:.3f} <= {bp4_threshold} ({gene}-specific threshold) "
                        f"but SpliceAI not available: BP4 requires confirmed SpliceAI <= 0.1 per ENIGMA VCEP v1.2"
                    )
                }

    elif vtype in ["synonymous", "silent"]:
        if in_domain and spliceai_score is not None and spliceai_score <= 0.1:
            results["BP4"] = {
                "applies": True,
                "strength": "Supporting",
                "points": -1,
                "reason": f"Silent variant in domain ({domain_name}), SpliceAI {spliceai_score:.3f} <= 0.1"
            }

    elif vtype == "intronic":
        if spliceai_score is not None and spliceai_score <= 0.1:
            results["BP4"] = {
                "applies": True,
                "strength": "Supporting",
                "points": -1,
                "reason": f"Intronic variant, SpliceAI {spliceai_score:.3f} <= 0.1"
            }

    return results

if __name__ == "__main__":
    print("Testing PP3/BP4 with BayesDel and SpliceAI:")
    print("=" * 60)
    print(f"SpliceAI 0.5 missense:             {evaluate_pp3_bp4('BRCA1', 'missense', 'p.(Arg170Gln)', spliceai_score=0.5)}")
    print(f"SpliceAI 0.5 nonsense:             {evaluate_pp3_bp4('BRCA1', 'nonsense', 'p.(Gln210Ter)', spliceai_score=0.5)}")
    print(f"BayesDel 0.35, in BRCT domain:     {evaluate_pp3_bp4('BRCA1', 'missense', 'p.(Arg1700Gln)', bayesdel_score=0.35, spliceai_score=0.02)}")
    print(f"BayesDel 0.10, in BRCT domain:     {evaluate_pp3_bp4('BRCA1', 'missense', 'p.(Arg1700Gln)', bayesdel_score=0.10, spliceai_score=0.03)}")
    print(f"BayesDel 0.35, outside domain:     {evaluate_pp3_bp4('BRCA1', 'missense', 'p.(Arg500Gln)', bayesdel_score=0.35, spliceai_score=0.02)}")
    print(f"BayesDel 0.20 (uninformative):     {evaluate_pp3_bp4('BRCA1', 'missense', 'p.(Arg1700Gln)', bayesdel_score=0.20, spliceai_score=0.03)}")
    print(f"No scores available:               {evaluate_pp3_bp4('BRCA1', 'missense', 'p.(Arg170Gln)')}")
