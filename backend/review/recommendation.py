"""Shared result shape for manual-review recommendations."""

from typing import Any

from backend.policy.gene import vcep_specification


def empty_review_recommendation(gene: str) -> dict[str, Any]:
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
