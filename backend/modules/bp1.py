from typing import Optional, Dict

from backend.modules.utils import (
    get_amino_acid_interval,
    overlapping_functional_domains,
)
from backend.modules.decision_trace import figure1a_path, step
from backend.gene_policy import spliceai_thresholds

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
    splice_low = spliceai_thresholds(gene)["bp4"]

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
    if spliceai_score > splice_low:
        result["reason"] = f"SpliceAI score {spliceai_score:.3f} > {splice_low} - possible splicing effect"
        return result

    # Check the complete explicitly affected protein interval. A boundary-
    # spanning in-frame change is not outside a domain merely because its first
    # residue lies outside that domain.
    protein_interval = get_amino_acid_interval(p_notation)
    if protein_interval is None:
        result["reason"] = "Could not determine the complete amino acid interval"
        return result

    domains = overlapping_functional_domains(gene, protein_interval)
    interval_label = (
        str(protein_interval[0])
        if protein_interval[0] == protein_interval[1]
        else f"{protein_interval[0]}-{protein_interval[1]}"
    )

    if domains:
        result["reason"] = (
            f"Variant amino acid interval {interval_label} overlaps "
            f"{', '.join(domains)} domain"
        )
        return result

    # variant is outside functional domain and no splicing effect
    result["applies"] = True
    result["strength"] = "Strong"
    result["points"] = -4  # benign evidence
    result["single_strong_likely_benign_eligible"] = True
    result["single_strong_likely_benign_basis"] = (
        "ENIGMA Table 3 footnote: BP1 Strong combines variant type, "
        "position outside a functional domain, and a no-impact splicing prediction"
    )
    result["independent_evidence_contribution_count"] = 3
    result["reason"] = (
        f"Variant amino acid interval {interval_label} is outside functional "
        "domains, no splicing predicted"
    )
    branch_id = "synonymous" if variant_type in {"synonymous", "silent"} else "missense-inframe"
    splice_node = "syn-splice-impact" if branch_id == "synonymous" else "mi-splice-impact"
    domain_node = "syn-domain" if branch_id == "synonymous" else "mi-domain-after-low"
    outcome_node = "syn-bp1" if branch_id == "synonymous" else "mi-bp1"
    result["decision_path"] = figure1a_path(
        branch_id=branch_id,
        criterion="BP1",
        outcome_node=outcome_node,
        steps=[
            step(splice_node, "Predicted impact on splicing?", "no_impact", f"SpliceAI {spliceai_score:.3f} ≤ {splice_low}"),
            step(
                domain_node,
                "Inside functional domain?",
                "no",
                f"Amino-acid interval {interval_label} is outside ENIGMA functional domains",
            ),
        ],
    )

    return result
