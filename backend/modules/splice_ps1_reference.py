"""Splice PS1 pilot reference candidates for manual review support."""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional


REFERENCE_PATH = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "splice_ps1_reference_set.json"
)


def load_splice_ps1_reference_candidates(
    gene: Optional[str] = None,
) -> Dict[str, Any]:
    """Return compact pilot candidates for UI prefill, not automatic scoring."""
    if not REFERENCE_PATH.exists():
        return {
            "curation_status": "missing",
            "description": "Splice PS1 pilot reference set is not available.",
            "usage_requirements": [],
            "candidates": [],
        }

    data = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    requested_gene = gene.upper() if gene else None
    variants: List[Dict[str, Any]] = []

    for record in data.get("variants", []):
        record_gene = record.get("gene", "")
        if requested_gene and record_gene != requested_gene:
            continue
        classification = record.get("classification", "")
        suggested_strength = {
            "Pathogenic": "Strong",
            "Likely Pathogenic": "Moderate",
        }.get(classification, "")
        source_file = record.get("source_file", "")
        source_sheet = record.get("source_sheet", "")
        source_row = record.get("source_row", "")
        variants.append(
            {
                "key": f"{record_gene}|{record.get('reference_variant', '')}|{source_row}",
                "gene": record_gene,
                "reference_variant": record.get("reference_variant", ""),
                "p_notation": record.get("p_notation", ""),
                "classification": classification,
                "classification_source": record.get("classification_source", ""),
                "prefill_strength_suggestion": suggested_strength,
                "prefill_strength_basis": "Conservative prefill from reference classification only; reviewer must confirm ENIGMA Appendix J/Table 17.",
                "splice_event_label": record.get("splice_event_label", ""),
                "event_type": record.get("event_type", ""),
                "assay_result_category": record.get("assay_result_category", ""),
                "variant_context": record.get("variant_context", ""),
                "prior_probability": record.get("prior_probability", ""),
                "source_file": source_file,
                "source_sheet": source_sheet,
                "source_row": source_row,
                "source_label": f"{source_file} / {source_sheet} row {source_row}",
                "source_url": record.get("source_url", ""),
                "curation_status": record.get("curation_status", ""),
                "curation_note": record.get("curation_note", ""),
            }
        )

    return {
        "curation_status": data.get("curation_status", ""),
        "description": data.get("description", ""),
        "usage_requirements": data.get("usage_requirements", []),
        "candidates": variants,
    }
