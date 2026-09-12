"""Resolve objective facts for a manually entered protein-level PS1 reference."""

from __future__ import annotations

import asyncio
from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional

from backend.lookups import clingen, clinvar, spliceai
from backend.modules.erepo_vcep import lookup_erepo_vcep_assertion
from backend.modules.ps1 import lookup_ps1_reference_variant
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
    erepo_registry_lookup: Callable[[str, str], Dict[str, Any]]
    registry_lookup: Callable[[str, str], Optional[Dict[str, Any]]]

    @classmethod
    def production(cls) -> "Ps1ReferenceDependencies":
        return cls(
            spliceai_lookup=spliceai.get_spliceai_score,
            spliceai_status=spliceai.get_spliceai_status,
            clinvar_lookup=clinvar.clinvar_lookup,
            clingen_lookup=clingen.clingen_erepo_lookup,
            erepo_registry_lookup=lookup_erepo_vcep_assertion,
            registry_lookup=lookup_ps1_reference_variant,
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
    registry_reference = dict(
        deps.registry_lookup(reference.gene, reference.c_notation) or {}
    )
    erepo_registry_result = dict(
        deps.erepo_registry_lookup(reference.gene, reference.c_notation) or {}
    )
    erepo_registry_record = dict(erepo_registry_result.get("record") or {})
    aggregate = dict(clinvar.get("aggregate") or {})
    stars = _clinvar_review_stars(str(aggregate.get("review_status") or ""))
    aggregate_class = _classification_label(aggregate.get("classification"))
    enigma_submission = dict(clinvar.get("enigma_submission") or {})
    enigma_class = _classification_label(enigma_submission.get("class"))
    local_erepo_class = _classification_label(erepo_registry_record.get("classification"))
    registry_class = _classification_label(registry_reference.get("classification"))

    classification = ""
    verification = "unresolved"
    source = ""
    if (
        erepo_registry_result.get("status") == "current_vcep_assertion"
        and local_erepo_class in P_LP
    ):
        classification = local_erepo_class
        verification = "external_vcep_assertion"
        source = (
            "ClinGen Evidence Repository ENIGMA BRCA1/2 VCEP v1.2 assertion "
            f"{erepo_registry_record.get('uuid', '')}"
        ).strip()
    elif (
        registry_reference.get("classification_basis")
        == "enigma_st7_v1_2_reference_set"
        and registry_class in P_LP
    ):
        classification = registry_class
        verification = "enigma_st7_v1_2_reference_set"
        source = "ENIGMA Supplementary Table 7 v1.2"
        original_source = str(registry_reference.get("classification_source") or "").strip()
        if original_source:
            source += f"; source recorded in ST7: {original_source}"
    elif aggregate_class in P_LP:
        # An aggregate ClinVar conclusion is candidate-discovery context only.
        # It is not the underlying VCEP evidence record and must not prefill
        # classification fields used to award PS1.
        classification = ""
        verification = "unresolved"
        source = ""

    non_current_erepo = erepo_registry_result.get("status") in {
        "historical_vcep_assertion",
        "unversioned_vcep_assertion",
    }
    historical_three_star = bool(
        stars == 3 and enigma_class in P_LP and verification == "unresolved"
    )
    historical_warning = ""
    if non_current_erepo:
        version = str(erepo_registry_record.get("assertion_method_version") or "unknown")
        historical_warning = (
            "ClinGen ERepo contains an ENIGMA expert-panel assertion for this variant, "
            f"but its recorded specification version is {version}, not the active v1.2. "
            "It is shown for context "
            "but is not used as an automatic v1.2 PS1 classification basis."
        )
    elif historical_three_star:
        historical_warning = (
            "ClinVar contains a three-star ENIGMA expert-panel assertion, but the "
            "active local ClinGen ERepo snapshot has no matching current v1.2 "
            "assertion. It is shown for context and does not qualify the PS1 "
            "reference automatically."
        )

    same_missense = assessed.p_notation == reference.p_notation
    different_nucleotide = assessed.c_notation != reference.c_notation
    registry_conditions_pass = (
        verification != "enigma_st7_v1_2_reference_set"
        or (
            registry_reference.get("reference_status") == "approved"
            and registry_reference.get("reference_splice_evidence_status")
            in {"none_identified", "normal"}
            and registry_reference.get("classification_ps1_dependency_used") is False
        )
    )
    objective_checks_pass = (
        verification == "external_vcep_assertion"
        or verification == "enigma_st7_v1_2_reference_set"
    ) and (
        same_missense
        and different_nucleotide
        and assessed_score is not None
        and reference_score is not None
        and assessed_score <= 0.1
        and reference_score <= 0.1
        and registry_conditions_pass
    )
    if verification == "external_vcep_assertion":
        review_message = (
            "A current ENIGMA BRCA1/2 VCEP v1.2 P/LP assertion was found in the "
            "checksum-validated local ClinGen ERepo snapshot. Complete the defined RNA/splice "
            "source check and reciprocal PS1 dependency review before submitting PS1."
        )
    elif verification == "enigma_st7_v1_2_reference_set":
        review_message = (
            "An ENIGMA ST7 v1.2 P/LP reference classification was found and is "
            "accepted as the classification basis. ARIANE also checked its recorded "
            "protein branch and defined RNA/splice sources."
        )
    elif aggregate_class in P_LP:
        review_message = (
            f"ClinVar reports an aggregate {aggregate_class} classification with "
            f"{stars} review star{'s' if stars != 1 else ''}. This is discovery "
            "context only. It did not prefill the reference classification, "
            "verification, or classification source and adds no PS1 points."
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
        review_message += (
            " The live ClinGen ERepo comparison was unavailable. The validated local "
            "snapshot remained the only source used for VCEP eligibility."
        )

    variation_id = str(clinvar.get("variation_id") or "")
    accession = str(clinvar.get("accession") or "")
    references = []
    if variation_id:
        references.append(f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/")
    if erepo_registry_record.get("erepo_url"):
        references.append(str(erepo_registry_record["erepo_url"]))

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
        "erepo_registry_status": str(
            erepo_registry_result.get("status") or "not_found"
        ),
        "erepo_assertion_uuid": str(erepo_registry_record.get("uuid") or ""),
        "erepo_assertion_method_version": str(
            erepo_registry_record.get("assertion_method_version") or ""
        ),
        "historical_expert_panel_warning": historical_warning,
        "classification": classification,
        "classification_verification": verification,
        "classification_source": source,
        "reference_registry_status": str(
            registry_reference.get("reference_status") or "not_found"
        ),
        "reference_confirmed_splice_status": str(
            registry_reference.get("reference_splice_evidence_status") or "not_assessed"
        ),
        "reference_splice_sources_checked": list(
            registry_reference.get("reference_splice_sources_checked") or []
        ),
        "reference_classification_used_ps1": (
            "no"
            if registry_reference.get("classification_ps1_dependency_used") is False
            else "yes"
            if registry_reference.get("classification_ps1_dependency_used") is True
            else "unknown"
        ),
        "objective_ps1_checks_pass": objective_checks_pass,
        "review_message": review_message,
        "references": references,
    }
