"""User-facing recommendation for non-automatic protein-level PS1 candidates."""

from __future__ import annotations

from typing import Any, Dict
from backend.gene_policy import (
    resolve_policy_gene,
    spliceai_thresholds,
    vcep_specification,
)


def evaluate_protein_ps1_review(
    ps1_result: Dict[str, Any] | None, *, gene: str | None = None
) -> Dict[str, Any]:
    policy_gene = resolve_policy_gene(gene)
    splice_low = spliceai_thresholds(policy_gene)["bp4"]
    ps1_result = ps1_result or {}
    if ps1_result.get("application_status") in {"reference_ineligible", "not_applicable"}:
        candidates = ps1_result.get("candidates", [])
        candidate_names = ", ".join(
            f"{item.get('gene')} {item.get('c_notation')} {item.get('p_notation')}"
            for item in candidates
        )
        reasons = [str(ps1_result.get("reason") or "Protein PS1 is not applicable")]
        if candidate_names:
            reasons.insert(
                0,
                "ST7 classifies the different reference variant(s) "
                f"{candidate_names}; it does not classify the variant currently under assessment.",
            )
        reasons.extend(
            item.get("status_reason", "")
            for item in candidates
            if item.get("status_reason")
        )
        return {
            **_empty(policy_gene),
            "display": True,
            "title": "Protein PS1 not applicable",
            "summary": (
                "ST7 contains a different P/LP reference variant with the same "
                "protein consequence. The recorded ENIGMA splice conditions exclude "
                "using that reference for protein-level PS1. No points were added."
            ),
            "reasons": list(dict.fromkeys(reasons)),
            "limitations": "An excluded reference cannot be confirmed through manual PS1 review.",
            "application_status": ps1_result.get("application_status"),
            "candidates": candidates,
            "splice_sources_checked": ps1_result.get("vua_splice_sources_checked", []),
            "vua_splice_evidence_status": ps1_result.get(
                "vua_splice_evidence_status", "not_assessed"
            ),
            "vua_spliceai_score": ps1_result.get("vua_spliceai_score"),
            "reference_spliceai_scores": ps1_result.get("reference_spliceai_scores", {}),
        }
    if not ps1_result.get("review_required"):
        return {
            **_empty(policy_gene),
            "application_status": ps1_result.get(
                "application_status", "not_applied"
            ),
        }
    candidates = ps1_result.get("candidates", [])
    candidate_names = ", ".join(
        f"{item.get('gene')} {item.get('c_notation')} {item.get('p_notation')}"
        for item in candidates
    )
    reasons = list(ps1_result.get("blocking_reasons", []))
    if candidate_names:
        reasons.insert(0, f"Matching P/LP reference candidate(s): {candidate_names}")
    return {
        "display": True,
        "recommended": True,
        "priority": "medium",
        "title": "Protein PS1 review candidate",
        "summary": (
            "ST7 contains a different P/LP reference variant with the same missense "
            "substitution; it does not classify the variant under assessment. The "
            "complete ENIGMA PS1 eligibility record is unavailable, so PS1 was not scored."
        ),
        "reasons": reasons,
        "what_to_check": [
            "Confirm that the reference P/LP classification was assigned using ENIGMA/ClinGen VCEP specifications.",
            "Confirm the same normalized missense substitution and a different nucleotide change.",
            f"Confirm SpliceAI <= {splice_low} for both variants.",
            "Check the defined RNA/splice evidence sources and confirm that no damaging splice effect is recorded for either variant.",
            "If the reference classification is known to use PS1, exclude a direct reciprocal dependency.",
        ],
        "potential_branches": ["PS1 (protein)"],
        "limitations": (
            "ST7 is an official ENIGMA reference dataset and a trusted candidate "
            "source, but it does not itself record every PS1-specific VCEP and "
            "splice eligibility field. This recommendation adds no points."
        ),
        "reference_source": "ENIGMA Supplementary Table 7 v1.2 candidate reference set",
        "source_url": vcep_specification(policy_gene)["url"],
        "is_evidence_criterion": False,
        "application_status": ps1_result.get(
            "application_status", "manual_review_required"
        ),
        "candidates": candidates,
        "splice_sources_checked": ps1_result.get(
            "vua_splice_sources_checked", []
        ),
        "vua_splice_evidence_status": ps1_result.get(
            "vua_splice_evidence_status", "not_assessed"
        ),
        "vua_spliceai_score": ps1_result.get("vua_spliceai_score"),
        "reference_spliceai_scores": ps1_result.get("reference_spliceai_scores", {}),
    }


def _empty(gene: str) -> Dict[str, Any]:
    return {
        "display": False,
        "recommended": False,
        "priority": "none",
        "title": "",
        "summary": "",
        "reasons": [],
        "what_to_check": [],
        "potential_branches": [],
        "limitations": "",
        "reference_source": "",
        "source_url": vcep_specification(gene)["url"],
        "is_evidence_criterion": False,
        "application_status": "not_applied",
        "candidates": [],
        "splice_sources_checked": [],
        "vua_splice_evidence_status": "not_assessed",
        "vua_spliceai_score": None,
        "reference_spliceai_scores": {},
    }
