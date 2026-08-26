from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse
from backend.modules.utils import get_intron_offset_from_c_notation
from backend.modules.decision_trace import figure1a_path, step
from backend.gene_policy import spliceai_thresholds

def evaluate_bp7(
    variant_type: str,
    spliceai_score: Optional[float] = None,
    in_domain: bool = False,
    bp4_met: bool = False,
    c_notation: str = "",
    *,
    gene: str,
) -> Dict:
    # evaluate BP7 criterion per ENIGMA VCEP v1.2
    #
    # ENIGMA: "Following convention, BP7 is applied in addition to BP4"
    # The logic depends on whether the variant is inside or outside a functional domain:
    #
    # Silent variant INSIDE functional domain:
    #   -> BP4_Supporting requires BayesDel <= threshold AND SpliceAI <= 0.1
    #   -> if BP4 met, also add BP7_Supporting
    #   -> if BP4 not met (e.g. BayesDel unavailable), BP7 does not apply alone
    #
    # Silent variant OUTSIDE functional domain:
    #   -> BP1_Strong applies (handled in evaluate_bp1, not here)
    #   -> BP7 does not additionally apply for out-of-domain silent variants
    #
    # So BP7 is only relevant for synonymous variants INSIDE a domain where BP4 is also met.
    result = {
        "applies": False,
        "strength": None,
        "points": 0,
        "reason": ""
    }
    splice_low = spliceai_thresholds(gene)["bp4"]

    if variant_type.lower() == "intronic":
        if not bp4_met:
            result["reason"] = "Intronic variant but BP4 not met - BP7 not applied"
            return result
        intron_info = get_intron_offset_from_c_notation(c_notation)
        if intron_info is None:
            result["reason"] = "Could not determine intronic offset - BP7 not applied"
            return result
        _, offset = intron_info
        if offset < 0 and offset > -21:
            result["reason"] = f"Intronic offset {offset:+d} is inside conserved acceptor motif - BP7 not applied"
            return result
        if offset > 0 and offset < 7:
            result["reason"] = f"Intronic offset {offset:+d} is inside conserved donor motif - BP7 not applied"
            return result
        result["applies"] = True
        result["strength"] = "Supporting"
        result["points"] = -1
        result["reason"] = (
            f"Intronic variant at offset {offset:+d}, outside conserved donor/acceptor motif, "
            "BP4 met (BP7 applied in addition to BP4 per ENIGMA)"
        )
        result["decision_path"] = figure1a_path(
            branch_id="intronic",
            criterion="BP7",
            outcome_node="int-bp4-bp7",
            steps=[
                step("int-canonical-site", "Donor/acceptor ±1,2 position?", "no", f"Intronic offset {offset:+d}"),
                step("int-splice-impact", "Predicted impact on splicing?", "no_impact", f"SpliceAI {spliceai_score:.3f} ≤ {splice_low}"),
                step("int-deep", "Deep intronic, at or beyond +7 or -21?", "yes", f"Intronic offset {offset:+d}"),
            ],
        )
        return result

    # BP7 also applies to synonymous variants inside a functional domain.
    if variant_type.lower() not in ["synonymous", "silent"]:
        result["reason"] = f"BP7 not applicable for {variant_type} variants"
        return result

    # BP7 is only relevant inside functional domains
    # silent variant outside domain -> BP1_Strong handles it, not BP7
    if not in_domain:
        result["reason"] = "Silent variant outside functional domain - BP1_Strong applies instead, not BP7"
        return result

    # inside domain: BP7 requires BP4 to be met first
    if not bp4_met:
        result["reason"] = "Silent variant in domain but BP4 not met - BP7 not applied"
        return result

    # SpliceAI must be confirmed <= 0.1 (BP4 already requires this, but verify)
    if spliceai_score is None:
        result["reason"] = "SpliceAI not available - cannot confirm no splice effect, BP7 not applied"
        return result
    if spliceai_score > splice_low:
        result["reason"] = f"SpliceAI {spliceai_score:.3f} > {splice_low} - possible splice effect, BP7 not applied"
        return result

    # silent variant inside domain, BP4 met, SpliceAI <= 0.1 -> BP7 in addition to BP4
    result["applies"] = True
    result["strength"] = "Supporting"
    result["points"] = -1
    result["reason"] = (
        f"Silent variant in functional domain, BP4 met, SpliceAI {spliceai_score:.3f} <= {splice_low} "
        "(BP7 applied in addition to BP4 per ENIGMA convention)"
    )
    result["decision_path"] = figure1a_path(
        branch_id="synonymous",
        criterion="BP7",
        outcome_node="syn-bp4-bp7",
        steps=[
            step("syn-splice-impact", "Predicted impact on splicing?", "no_impact", f"SpliceAI {spliceai_score:.3f} ≤ {splice_low}"),
            step("syn-domain", "Inside functional domain?", "yes", "Variant is inside an ENIGMA functional domain"),
        ],
    )

    return result


if __name__ == "__main__":
    # Test BP7
    print("Testing BP7 evaluation:")
    print("=" * 60)
    print(f"Synonymous outside domain, SpliceAI=0.05:            {evaluate_bp7('synonymous', 0.05, in_domain=False, bp4_met=False)}")
    print(f"Synonymous in domain, BP4 met, SpliceAI=0.05:     {evaluate_bp7('synonymous', 0.05, in_domain=True, bp4_met=True)}")
    print(f"Synonymous in domain, BP4 not met, SpliceAI=0.05: {evaluate_bp7('synonymous', 0.05, in_domain=True, bp4_met=False)}")
    print(f"Synonymous in domain, BP4 met, SpliceAI=0.3:      {evaluate_bp7('synonymous', 0.3, in_domain=True, bp4_met=True)}")
    print(f"Missense, SpliceAI=0.05:                          {evaluate_bp7('missense', 0.05)}")
