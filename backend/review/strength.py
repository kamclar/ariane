"""Derive manual-criterion strength from structured evidence and VCEP policy."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable, Dict, Mapping, Optional, Sequence

from backend.criteria.bp7_rna import evaluate_bp7_rna_variant_context
from backend.policy.gene import clinical_lr_thresholds, resolve_policy_gene, spliceai_thresholds
from backend.reference_data.ps1_protein import lookup_ps1_reference_variant
from backend.reference_data.ps1_splice_evidence import DEFINED_SOURCES as PS1_SPLICE_SOURCES
from backend.review.definitions import ENIGMA_RECOGNISED_PP4_SOURCES, MANUAL_CRITERIA
from backend.variant_processing.variant_input import normalize_variant_input
from backend.variant_processing.variant_type import infer_variant_type

def _number(evidence: Mapping[str, Any], key: str) -> Optional[float]:
    value = evidence.get(key)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def evaluate_bs4_likelihood_ratio(
    evidence: Mapping[str, Any],
) -> tuple[Optional[float], bool, int]:
    """Return the BS4 LR and Table 3 single-Strong eligibility provenance."""
    aggregate = _number(dict(evidence), "likelihood_ratio")
    components = evidence.get("likelihood_ratio_components")
    if not components:
        return aggregate, False, 0
    if not isinstance(components, list):
        raise ValueError("BS4 likelihood-ratio components must be a list")

    component_lrs: list[float] = []
    independence_groups: list[str] = []
    for index, component in enumerate(components, start=1):
        if not isinstance(component, Mapping):
            raise ValueError(f"BS4 LR component {index} must be a structured record")
        raw_lr = component.get("likelihood_ratio")
        try:
            lr = float(raw_lr) if raw_lr is not None else None
        except (TypeError, ValueError):
            raise ValueError(f"BS4 LR component {index} requires a numeric LR") from None
        if lr is None:
            raise ValueError(f"BS4 LR component {index} requires a numeric LR")
        if lr < 0:
            raise ValueError(f"BS4 LR component {index} cannot be negative")
        source = str(component.get("source") or "").strip()
        group = str(component.get("independence_group") or "").strip()
        if not source or not group:
            raise ValueError(
                f"BS4 LR component {index} requires a source and independence group"
            )
        component_lrs.append(lr)
        independence_groups.append(group)

    if len(independence_groups) != len(set(independence_groups)):
        raise ValueError("BS4 LR components must use distinct independence groups")

    combined = math.prod(component_lrs)
    if aggregate is not None and not math.isclose(
        aggregate, combined, rel_tol=1e-6, abs_tol=1e-12
    ):
        raise ValueError(
            "BS4 reported combined LR does not equal the product of its LR components"
        )
    return combined, len(component_lrs) >= 2, len(component_lrs)


def _pp4_value_and_scale(evidence: Mapping[str, Any]) -> tuple[Optional[float], str]:
    """Read the current PP4 input while retaining legacy LR audit records."""
    if evidence.get("clinical_lr_value") not in (None, ""):
        return _number(evidence, "clinical_lr_value"), str(
            evidence.get("clinical_lr_scale") or "lr"
        ).strip().lower()
    return _number(evidence, "combined_clinical_lr"), "lr"


def _pp4_strength(evidence: Mapping[str, Any], gene: str) -> Optional[str]:
    value, scale = _pp4_value_and_scale(evidence)
    if value is None or value < 0 or scale not in {"lr", "log10_lr", "acmg_points"}:
        return None
    pp4 = clinical_lr_thresholds(gene)["pp4"]
    lr_thresholds = (
        pp4["supporting_min_inclusive"],
        pp4["moderate_min_inclusive"],
        pp4["strong_min_inclusive"],
        pp4["very_strong_min_inclusive"],
    )
    thresholds = {
        "lr": lr_thresholds,
        "log10_lr": tuple(math.log10(value) for value in lr_thresholds),
        "acmg_points": (1.0, 2.0, 4.0, 8.0),
    }[scale]
    supporting, moderate, strong, very_strong = thresholds
    if value >= very_strong:
        return "Very Strong"
    if value >= strong:
        return "Strong"
    if value >= moderate:
        return "Moderate"
    if value >= supporting:
        return "Supporting"
    return None


def _bp5_strength(evidence: Mapping[str, Any], gene: str) -> Optional[str]:
    value, scale = _pp4_value_and_scale(evidence)
    if value is None or scale not in {"lr", "log10_lr", "acmg_points"}:
        return None
    bp5 = clinical_lr_thresholds(gene)["bp5"]
    lr_thresholds = (
        bp5["supporting_max_inclusive"],
        bp5["moderate_max_inclusive"],
        bp5["strong_max_inclusive"],
        bp5["very_strong_max_inclusive"],
    )
    thresholds = {
        "lr": lr_thresholds,
        "log10_lr": tuple(math.log10(item) for item in lr_thresholds),
        "acmg_points": (-1.0, -2.0, -4.0, -8.0),
    }[scale]
    supporting, moderate, strong, very_strong = thresholds
    if value <= very_strong:
        return "Very Strong"
    if value <= strong:
        return "Strong"
    if value <= moderate:
        return "Moderate"
    if value <= supporting:
        return "Supporting"
    return None


def _pp4_source_is_reviewed(evidence: Mapping[str, Any]) -> bool:
    status = str(evidence.get("source_review_status") or "unreviewed").strip().lower()
    if status == "enigma_recognised":
        return (
            str(evidence.get("source_pmid") or "").strip()
            in ENIGMA_RECOGNISED_PP4_SOURCES
        )
    if status == "other_reviewed":
        return all(
            bool((evidence.get(field) or "").strip())
            for field in ("source_citation", "source_reviewed_by", "source_review_rationale")
        )
    return False


def _pp4_source_is_recorded(evidence: Mapping[str, Any]) -> bool:
    status = str(evidence.get("source_review_status") or "unreviewed").strip().lower()
    if status == "enigma_recognised":
        return (
            str(evidence.get("source_pmid") or "").strip()
            in ENIGMA_RECOGNISED_PP4_SOURCES
        )
    if status == "other_reviewed":
        return _pp4_source_is_reviewed(evidence)
    if status == "unreviewed":
        return bool((evidence.get("source_citation") or "").strip())
    return False


@dataclass(frozen=True)
class _SuggestionContext:
    gene: str
    variant: Mapping[str, Any] | None
    base_criteria: Sequence[Mapping[str, Any]]


def _curated_strength(code: str, evidence: Mapping[str, Any]) -> Optional[str]:
    strength = evidence.get("curated_strength")
    return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None


def _suggest_clinical_lr(
    code: str, evidence: Dict[str, Any], context: _SuggestionContext
) -> Optional[str]:
    if not (evidence.get("clinical_data_summary") or "").strip():
        return None
    if not _pp4_source_is_reviewed(evidence):
        return None
    return (
        _pp4_strength(evidence, context.gene)
        if code == "PP4"
        else _bp5_strength(evidence, context.gene)
    )


def _suggest_functional(
    code: str, evidence: Dict[str, Any], _context: _SuggestionContext
) -> Optional[str]:
    expected_conclusion = "abnormal" if code == "PS3" else "normal"
    required_text_fields = (
        "assay_name",
        "source_citation",
        "calibration_summary",
        "variant_result_summary",
        "functional_reviewed_by",
    )
    if (
        evidence.get("assay_scope") not in {"protein_only", "combined_mrna_protein"}
        or evidence.get("functional_conclusion") != expected_conclusion
        or evidence.get("calibration_status") != "reviewed_under_enigma_vcep"
        or evidence.get("pathogenic_and_benign_controls_confirmed") is not True
        or any(
            not str(evidence.get(field) or "").strip() for field in required_text_fields
        )
    ):
        return None
    return _curated_strength(code, evidence)


def _suggest_rna(
    code: str, evidence: Dict[str, Any], context: _SuggestionContext
) -> Optional[str]:
    if (
        evidence.get("assay_scope") != "mrna_only"
        or not (evidence.get("transcript_accession") or "").strip()
        or not (evidence.get("tissue_or_cell_type") or "").strip()
        or evidence.get("nmd_assessed") not in {"yes", "no", "not_applicable"}
    ):
        return None
    if code == "PVS1_RNA":
        if evidence.get("rna_conclusion") != "damaging":
            return None
        if evidence.get("functional_transcript_remaining") not in {
            "absent_or_minimal",
            "reduced",
        }:
            return None
        return _curated_strength(code, evidence)
    if evidence.get("rna_conclusion") != "no_damaging_effect":
        return None
    eligibility = evaluate_bp7_rna_variant_context(
        context.variant, context.base_criteria
    )
    return "Strong" if eligibility["eligible"] else None


def _suggest_initiation(
    code: str, evidence: Dict[str, Any], _context: _SuggestionContext
) -> Optional[str]:
    if evidence.get("met1_loss_confirmed") is not True:
        return None
    if evidence.get("alternative_start_assessed") not in {"yes", "no"}:
        return None
    if evidence.get("upstream_pathogenic_evidence") not in {
        "yes",
        "no",
        "not_applicable",
    }:
        return None
    if evidence.get("functional_domain_impact") not in {
        "yes",
        "no",
        "uncertain",
        "not_applicable",
    }:
        return None
    if not (evidence.get("nearest_alternative_start") or "").strip():
        return None
    if not (evidence.get("initiation_flowchart_rationale") or "").strip():
        return None
    return _curated_strength(code, evidence)


def _suggest_ps1_splice(
    code: str, evidence: Dict[str, Any], _context: _SuggestionContext
) -> Optional[str]:
    required_text_fields = (
        "reference_variant",
        "reference_classification_source",
        "vua_splice_event",
        "reference_splice_event",
        "ps1_splice_rationale",
    )
    if any(not (evidence.get(field) or "").strip() for field in required_text_fields):
        return None
    if evidence.get("reference_classification") not in {
        "Pathogenic",
        "Likely Pathogenic",
    }:
        return None
    if evidence.get("same_splice_event_confirmed") is not True:
        return None
    if evidence.get("prediction_strength_comparison") not in {"similar", "stronger"}:
        return None
    return _curated_strength(code, evidence)


def _same_missense_change(
    evidence: Mapping[str, Any], context: _SuggestionContext
) -> bool:
    if context.variant is None:
        return True
    try:
        assessed = normalize_variant_input(
            context.gene,
            str(context.variant.get("c_notation") or ""),
            p_notation=str(context.variant.get("p_notation") or ""),
        )
        reference = normalize_variant_input(
            context.gene,
            str(evidence.get("reference_variant") or ""),
            p_notation=str(evidence.get("reference_p_notation") or ""),
        )
    except ValueError:
        return False
    return bool(
        infer_variant_type(assessed.c_notation, assessed.p_notation) == "missense"
        and infer_variant_type(reference.c_notation, reference.p_notation) == "missense"
        and assessed.p_notation == reference.p_notation
        and assessed.c_notation != reference.c_notation
    )


def _approved_st7_reference(
    evidence: Mapping[str, Any], context: _SuggestionContext
) -> bool:
    try:
        normalized = normalize_variant_input(
            context.gene,
            str(evidence.get("reference_variant") or ""),
            p_notation=str(evidence.get("reference_p_notation") or ""),
        )
        reference = lookup_ps1_reference_variant(context.gene, normalized.c_notation)
    except ValueError:
        return False
    return bool(
        reference
        and reference.get("classification_basis") == "enigma_st7_v1_2_reference_set"
        and reference.get("reference_status") == "approved"
        and reference.get("classification") == evidence.get("reference_classification")
        and reference.get("p_notation") == normalized.p_notation
        and reference.get("reference_splice_evidence_status")
        in {"none_identified", "normal"}
        and reference.get("classification_ps1_dependency_used") is False
    )


def _ps1_protein_splicing_is_clear(
    evidence: Dict[str, Any], context: _SuggestionContext
) -> bool:
    if (
        evidence.get("same_missense_confirmed") is not True
        or evidence.get("different_nucleotide_change_confirmed") is not True
        or evidence.get("splice_source_check_completed") is not True
    ):
        return False
    sources = evidence.get("splice_sources_checked")
    if not isinstance(sources, list) or not set(PS1_SPLICE_SOURCES).issubset(
        {str(source).strip() for source in sources}
    ):
        return False
    accepted_statuses = {"none_identified", "normal"}
    if (
        evidence.get("vua_confirmed_splice_status") not in accepted_statuses
        or evidence.get("reference_confirmed_splice_status") not in accepted_statuses
    ):
        return False
    vua_score = _number(evidence, "vua_spliceai_score")
    reference_score = _number(evidence, "reference_spliceai_score")
    splice_low = spliceai_thresholds(context.gene)["bp4"]
    if (
        vua_score is None
        or reference_score is None
        or vua_score > splice_low
        or reference_score > splice_low
    ):
        return False
    dependency = evidence.get("reference_classification_used_ps1")
    if dependency not in {"no", "yes"}:
        return False
    return dependency == "no" or (
        bool(str(evidence.get("reference_ps1_dependency_reference") or "").strip())
        and evidence.get("direct_reciprocal_dependency_excluded") is True
    )


def _suggest_ps1_protein(
    _code: str, evidence: Dict[str, Any], context: _SuggestionContext
) -> Optional[str]:
    required_text_fields = (
        "reference_variant",
        "reference_p_notation",
        "classification_source",
        "ps1_protein_rationale",
    )
    if any(
        not str(evidence.get(field) or "").strip() for field in required_text_fields
    ):
        return None
    if not _same_missense_change(evidence, context):
        return None
    classification = evidence.get("reference_classification")
    if classification not in {"Pathogenic", "Likely Pathogenic"}:
        return None
    verification = evidence.get("classification_verification")
    if verification not in {
        "external_vcep_assertion",
        "locally_recurated_under_enigma_vcep",
        "enigma_st7_v1_2_reference_set",
    }:
        return None
    if verification == "enigma_st7_v1_2_reference_set" and not _approved_st7_reference(
        evidence, context
    ):
        return None
    if not _ps1_protein_splicing_is_clear(evidence, context):
        return None
    return "Strong" if classification == "Pathogenic" else "Moderate"


def _suggest_ps4(
    _code: str, evidence: Dict[str, Any], _context: _SuggestionContext
) -> Optional[str]:
    p_value = _number(evidence, "p_value")
    odds_ratio = _number(evidence, "odds_ratio")
    lower_ci = _number(evidence, "lower_ci")
    if (
        p_value is not None
        and odds_ratio is not None
        and lower_ci is not None
        and p_value <= 0.05
        and odds_ratio >= 4
        and lower_ci > 2
        and evidence.get("case_control_country_matched") is True
        and evidence.get("case_control_ethnicity_matched") is True
    ):
        return "Strong"
    return None


def _suggest_cooccurrence(
    code: str, evidence: Dict[str, Any], _context: _SuggestionContext
) -> Optional[str]:
    if evidence.get("cooccurring_variant_classification_basis") != "vcep_specifications":
        return None
    if code == "PM3" and evidence.get("vua_benign_population_review") != "does_not_meet":
        return None
    points = _number(evidence, "evidence_points")
    if points is None or points < 1:
        return None
    if points >= 4:
        return "Strong"
    if points >= 2:
        return "Moderate"
    return "Supporting"


def _suggest_segregation(
    code: str, evidence: Dict[str, Any], context: _SuggestionContext
) -> Optional[str]:
    likelihood_ratio = _number(evidence, "likelihood_ratio")
    if code == "BS4":
        likelihood_ratio, _, _ = evaluate_bs4_likelihood_ratio(evidence)
    if likelihood_ratio is None or likelihood_ratio < 0:
        return None
    thresholds = clinical_lr_thresholds(context.gene)
    if code == "PP1":
        pp4 = thresholds["pp4"]
        if likelihood_ratio >= pp4["very_strong_min_inclusive"] and evidence.get(
            "very_strong_effect_basis"
        ) in {
            "predicted_protein",
            "predicted_splicing",
            "experimental_protein",
            "experimental_splicing",
        }:
            return "Very Strong"
        if likelihood_ratio >= pp4["strong_min_inclusive"]:
            return "Strong"
        if likelihood_ratio >= pp4["moderate_min_inclusive"]:
            return "Moderate"
        if likelihood_ratio >= pp4["supporting_min_inclusive"]:
            return "Supporting"
        return None
    bp5 = thresholds["bp5"]
    if likelihood_ratio <= bp5["very_strong_max_inclusive"]:
        return "Very Strong"
    if likelihood_ratio <= bp5["strong_max_inclusive"]:
        return "Strong"
    if likelihood_ratio <= bp5["moderate_max_inclusive"]:
        return "Moderate"
    if likelihood_ratio <= bp5["supporting_max_inclusive"]:
        return "Supporting"
    return None


_SuggestionHandler = Callable[
    [str, Dict[str, Any], _SuggestionContext], Optional[str]
]

_SUGGESTION_HANDLERS: dict[str, _SuggestionHandler] = {
    "PP4": _suggest_clinical_lr,
    "BP5": _suggest_clinical_lr,
    "PS3": _suggest_functional,
    "BS3": _suggest_functional,
    "PVS1_RNA": _suggest_rna,
    "BP7_RNA": _suggest_rna,
    "PVS1_INIT": _suggest_initiation,
    "PS1_SPLICE": _suggest_ps1_splice,
    "PS1_PROTEIN": _suggest_ps1_protein,
    "PS4": _suggest_ps4,
    "PM3": _suggest_cooccurrence,
    "BS2": _suggest_cooccurrence,
    "PP1": _suggest_segregation,
    "BS4": _suggest_segregation,
}


def suggest_strength(
    code: str,
    evidence: Dict[str, Any],
    *,
    variant_context: Mapping[str, Any] | None = None,
    base_criteria: Sequence[Mapping[str, Any]] = (),
) -> Optional[str]:
    handler = _SUGGESTION_HANDLERS.get(code)
    if handler is None:
        return None
    context = _SuggestionContext(
        gene=resolve_policy_gene(
            str((variant_context or {}).get("gene") or "").strip().upper() or None
        ),
        variant=variant_context,
        base_criteria=base_criteria,
    )
    return handler(code, evidence, context)


def _clinical_lr_record_is_complete(evidence: Mapping[str, Any]) -> bool:
    value, scale = _pp4_value_and_scale(evidence)
    source_status = str(
        evidence.get("source_review_status") or "unreviewed"
    ).strip().lower()
    return bool(
        value is not None
        and (value >= 0 if scale == "lr" else True)
        and scale in {"lr", "log10_lr", "acmg_points"}
        and source_status
        in {"enigma_recognised", "other_reviewed", "unreviewed"}
        and _pp4_source_is_recorded(evidence)
        and str(evidence.get("clinical_data_summary") or "").strip()
    )

