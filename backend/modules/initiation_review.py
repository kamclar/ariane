"""Informational recommendation for initiation-codon PVS1 review."""

from typing import Dict


ENIGMA_CSPEC_URL = "https://cspec.genome.network/cspec/ui/svi/doc/GN097"
INITIATION_REFERENCE_SOURCE = (
    "Rule/weighting: ENIGMA BRCA1/2 VCEP v1.2 Specifications Table 4 and "
    "Appendix D initiation-codon PVS1 flowchart."
)


def evaluate_initiation_review(variant_type: str) -> Dict:
    """Flag Met1/start-loss variants for curated PVS1_INIT review."""
    if (variant_type or "").lower() != "initiation_codon":
        return _empty()

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
        "source_url": ENIGMA_CSPEC_URL,
        "is_evidence_criterion": False,
    }


def _empty() -> Dict:
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
        "source_url": ENIGMA_CSPEC_URL,
        "is_evidence_criterion": False,
    }
