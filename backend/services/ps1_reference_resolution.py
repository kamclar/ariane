"""Resolve objective facts for a manually entered protein-level PS1 reference."""

from __future__ import annotations

import asyncio
from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional

from backend.modules.variant_input import normalize_variant_input
from backend.modules.variant_type import infer_variant_type


P_LP = {"Pathogenic", "Likely Pathogenic"}


def _clinvar_review_stars(review_status: str) -> int:
    status = (review_status or "").strip().lower()
    if "no assertion criteria" in status or "no classification" in status:
        return 0
    if "practice guideline" in status:
        return 4
    if "expert panel" in status:
        return 3
    if "multiple submitters" in status and "no conflicts" in status:
        return 2
    if "criteria provided" in status:
        return 1
    return 0


def _classification_label(value: Any) -> str:
    text = str(value or "").strip().lower()
    return {
        "pathogenic": "Pathogenic",
        "likely pathogenic": "Likely Pathogenic",
    }.get(text, "")


@dataclass(frozen=True)
class Ps1ReferenceDependencies:
    spliceai_lookup: Callable[[str, str], Optional[float]]
    spliceai_status: Callable[[str, str], Dict[str, Any]]
    clinvar_lookup: Callable[[str, str], Dict[str, Any]]
    clingen_lookup: Callable[[str, str], Dict[str, Any]]

    @classmethod
    def production(cls) -> "Ps1ReferenceDependencies":
        from backend.lookups.clingen import clingen_erepo_lookup
        from backend.lookups.clinvar import clinvar_lookup
        from backend.lookups.spliceai import get_spliceai_score, get_spliceai_status

        return cls(
            spliceai_lookup=get_spliceai_score,
            spliceai_status=get_spliceai_status,
            clinvar_lookup=clinvar_lookup,
            clingen_lookup=clingen_erepo_lookup,
        )


def _resolved_variant(normalized, score, status):
    status = dict(status or {})
    return {
        "gene": normalized.gene,
        "reference_transcript": normalized.reference_transcript,
        "c_notation": normalized.c_notation,
        "p_notation": normalized.p_notation,
        "spliceai_score": score,
        "spliceai_status": str(status.get("status") or ("ok" if score is not None else "unavailable")),
        "spliceai_reason": str(status.get("reason") or ""),
    }


def _safe_external_lookup(lookup, gene: str, c_notation: str) -> Dict[str, Any]:
    try:
        return dict(lookup(gene, c_notation) or {})
    except Exception as exc:
        return {
            "status": "api_error",
            "error": f"{type(exc).__name__}: {exc}",
        }


async def resolve_ps1_reference(
    gene: str,
    assessed_c_notation: str,
    reference_c_notation: str,
    *,
    dependencies: Ps1ReferenceDependencies | None = None,
) -> Dict[str, Any]:
    """Normalize both variants and collect PS1 facts without assigning points."""
    deps = dependencies or Ps1ReferenceDependencies.production()
    assessed = normalize_variant_input(gene, assessed_c_notation)
    reference = normalize_variant_input(gene, reference_c_notation)
    if assessed.gene != reference.gene:
        raise ValueError(
            f"Reference variant gene {reference.gene} does not match assessed variant gene {assessed.gene}"
        )
    assessed_type = infer_variant_type(assessed.c_notation, assessed.p_notation)
    reference_type = infer_variant_type(reference.c_notation, reference.p_notation)
    if assessed_type != "missense" or reference_type != "missense":
        raise ValueError("Protein-level PS1 requires two normalized missense variants")

    assessed_score, reference_score, clinvar, clingen = await asyncio.gather(
        asyncio.to_thread(deps.spliceai_lookup, assessed.gene, assessed.c_notation),
        asyncio.to_thread(deps.spliceai_lookup, reference.gene, reference.c_notation),
        asyncio.to_thread(
            _safe_external_lookup,
            deps.clinvar_lookup,
            reference.gene,
            reference.c_notation,
        ),
        asyncio.to_thread(
            _safe_external_lookup,
            deps.clingen_lookup,
            reference.gene,
            reference.c_notation,
        ),
    )
    assessed_status = deps.spliceai_status(assessed.gene, assessed.c_notation)
    reference_status = deps.spliceai_status(reference.gene, reference.c_notation)

    clinvar = dict(clinvar or {})
    clingen = dict(clingen or {})
    aggregate = dict(clinvar.get("aggregate") or {})
    stars = _clinvar_review_stars(str(aggregate.get("review_status") or ""))
    aggregate_class = _classification_label(aggregate.get("classification"))
    enigma_submission = dict(clinvar.get("enigma_submission") or {})
    enigma_class = _classification_label(enigma_submission.get("class"))
    erepo_class = _classification_label(clingen.get("classification"))

    classification = ""
    verification = "unresolved"
    source = ""
    if enigma_class in P_LP:
        classification = enigma_class
        verification = "external_vcep_assertion"
        source = "ClinVar ENIGMA expert-panel assertion"
        if enigma_submission.get("scv"):
            source += f" {enigma_submission['scv']}"
    elif clingen.get("status") == "ok" and erepo_class in P_LP:
        classification = erepo_class
        verification = "external_vcep_assertion"
        source = "ClinGen Evidence Repository ENIGMA assertion"
        if clingen.get("caid"):
            source += f" {clingen['caid']}"
    elif aggregate_class in P_LP:
        classification = aggregate_class
        verification = "historical_classification_only"
        source = f"ClinVar aggregate, {stars} star" + ("s" if stars != 1 else "")

    same_missense = assessed.p_notation == reference.p_notation
    different_nucleotide = assessed.c_notation != reference.c_notation
    objective_checks_pass = (
        verification == "external_vcep_assertion"
        and same_missense
        and different_nucleotide
        and assessed_score is not None
        and reference_score is not None
        and assessed_score <= 0.1
        and reference_score <= 0.1
    )
    if verification == "external_vcep_assertion":
        review_message = (
            "An ENIGMA VCEP P/LP assertion was found. Complete the defined RNA/splice "
            "source check and reciprocal PS1 dependency review before submitting PS1."
        )
    elif aggregate_class in P_LP:
        review_message = (
            f"ClinVar reports {aggregate_class} with {stars} star"
            f"{'s' if stars != 1 else ''}, but this is not a verified ENIGMA VCEP "
            "assertion. It can support candidate review only and adds no PS1 points."
        )
    else:
        review_message = (
            "No Pathogenic or Likely Pathogenic ENIGMA VCEP assertion was found for "
            "the entered reference variant."
        )
    if not same_missense:
        review_message += " The variants do not have the same normalized missense substitution."
    if not different_nucleotide:
        review_message += " PS1 requires a different nucleotide change."
    if assessed_score is None or reference_score is None:
        review_message += " At least one required SpliceAI result is unavailable."
    if clinvar.get("status") == "api_error":
        review_message += " ClinVar was unavailable, so its classification could not be checked."
    if clingen.get("status") == "api_error":
        review_message += " ClinGen ERepo was unavailable, so its VCEP assertion could not be checked."

    variation_id = str(clinvar.get("variation_id") or "")
    accession = str(clinvar.get("accession") or "")
    references = []
    if variation_id:
        references.append(f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/")
    if clingen.get("status") == "ok":
        references.append("https://erepo.clinicalgenome.org/evrepo/")

    return {
        "assessed": _resolved_variant(assessed, assessed_score, assessed_status),
        "reference": _resolved_variant(reference, reference_score, reference_status),
        "same_missense_substitution": same_missense,
        "different_nucleotide_change": different_nucleotide,
        "clinvar_status": str(clinvar.get("status") or "not_found"),
        "clinvar_error": str(clinvar.get("error") or ""),
        "clinvar_variation_id": variation_id,
        "clinvar_accession": accession,
        "clinvar_classification": aggregate_class,
        "clinvar_review_status": str(aggregate.get("review_status") or ""),
        "clinvar_stars": stars,
        "clingen_status": str(clingen.get("status") or "not_found"),
        "clingen_error": str(clingen.get("error") or ""),
        "clingen_caid": str(clingen.get("caid") or ""),
        "classification": classification,
        "classification_verification": verification,
        "classification_source": source,
        "objective_ps1_checks_pass": objective_checks_pass,
        "review_message": review_message,
        "references": references,
    }
