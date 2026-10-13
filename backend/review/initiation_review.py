"""Informational recommendation for initiation-codon PVS1 review."""

from typing import Dict
from backend.policy.gene import vcep_specification
from backend.reference_data.table4 import table4_lookup_initiation


INITIATION_REFERENCE_SOURCE = (
    "Rule/weighting: ENIGMA BRCA1/2 VCEP v1.2 Specifications Table 4 and "
    "Appendix D initiation-codon PVS1 flowchart."
)


def evaluate_initiation_review(
    gene: str,
    variant_type: str,
    *,
    c_notation: str = "",
    p_notation: str = "",
    reference_transcript: str = "",
) -> Dict:
    """Flag Met1/start-loss variants for curated PVS1_INIT review."""
    if (variant_type or "").lower() != "initiation_codon":
        return _empty(gene)

    table4 = table4_lookup_initiation(gene)
    prefill = {
        "met1_loss_confirmed": True,
        "assessed_gene": gene,
        "assessed_c_notation": c_notation,
        "assessed_p_notation": p_notation,
        "reference_transcript": reference_transcript,
        "table4_code": table4.get("pvs1_code") or "",
        "table4_exon": table4.get("exon") or "",
        "source_citation": "ENIGMA BRCA1/2 VCEP v1.2 Specifications Table 4",
        "source_references": [
            INITIATION_REFERENCE_SOURCE,
            vcep_specification(gene)["url"],
        ],
    }
    if table4.get("found"):
        appendix_figure = table4["appendix_figure"]
        nearest_start = table4["nearest_alternative_start"]
        prefill.update({
            "curated_strength": table4["pvs1_strength"],
            "alternative_start_assessed": table4["alternative_start_assessed"],
            "nearest_alternative_start": nearest_start,
            "upstream_pathogenic_evidence": table4["upstream_pathogenic_evidence"],
            "functional_domain_impact": table4["functional_domain_impact"],
            "alternative_start_conclusion": table4["alternative_start_conclusion"],
            "appendix_figure": appendix_figure,
            "initiation_flowchart_rationale": (
                f"ENIGMA Table 4 assigns {table4['pvs1_code']} "
                f"({table4['pvs1_strength']}) to the {gene} start-codon branch "
                f"in {table4['exon']}. {appendix_figure} records no known "
                f"alternative start codon, identifies {nearest_start} as the "
                "closest potential in-frame start, and reports at least one "
                "pathogenic variant upstream. N-terminal functional impact is "
                "not a decision node for this initiation-codon branch."
            ),
        })
        prefill["source_references"].insert(
            1,
            f"ENIGMA BRCA1/2 VCEP v1.2 {appendix_figure}",
        )

    return {
        "recommended": True,
        "priority": "high",
        "title": "Initiation-codon PVS1 review recommended",
        "summary": (
            "This appears to be a Met1/start-loss variant. ARIANE does not score "
            "the ENIGMA initiation-codon PVS1 flowchart automatically; use the "
            "structured PVS1_INIT manual branch if the curated review supports it."
        ),
        "reasons": [
            "The variant type was inferred as initiation_codon from the HGVS/protein notation.",
            "ENIGMA initiation-codon PVS1 strength depends on alternative start codons, upstream P/LP evidence, and expected N-terminal functional impact.",
        ],
        "what_to_test": [
            "Confirm that the variant affects the canonical Met1/start codon on the relevant BRCA transcript.",
            "Assess whether an in-frame downstream alternative start codon is available.",
            "Check whether pathogenic or likely pathogenic variants occur upstream of the nearest alternative start.",
            "Assess whether the predicted N-terminal loss affects an established functional region.",
            "Select the PVS1_INIT strength from the ENIGMA initiation-codon flowchart and document the rationale.",
        ],
        "potential_branches": ["PVS1_INIT"],
        "limitations": (
            "This is a review aid only. It adds no points and does not replace the "
            "curated ENIGMA initiation-codon flowchart assessment."
        ),
        "reference_source": INITIATION_REFERENCE_SOURCE,
        "source_url": vcep_specification(gene)["url"],
        "is_evidence_criterion": False,
        "manual_review_prefill": prefill,
    }


def _empty(gene: str) -> Dict:
    return {
        "recommended": False,
        "priority": "none",
        "title": "",
        "summary": "",
        "reasons": [],
        "what_to_test": [],
        "potential_branches": [],
        "limitations": "",
        "reference_source": "",
        "source_url": vcep_specification(gene)["url"],
        "is_evidence_criterion": False,
        "manual_review_prefill": {},
    }
