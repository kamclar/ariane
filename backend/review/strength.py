"""Derive manual-criterion strength from structured evidence and VCEP policy."""

from __future__ import annotations

import math
from typing import Any, Dict, Mapping, Optional, Sequence

from backend.criteria.bp7_rna import evaluate_bp7_rna_variant_context
from backend.criteria.ps1 import lookup_ps1_reference_variant
from backend.policy.gene import clinical_lr_thresholds, resolve_policy_gene, spliceai_thresholds
from backend.reference_data.ps1_splice_evidence import DEFINED_SOURCES as PS1_SPLICE_SOURCES
from backend.review.definitions import ENIGMA_RECOGNISED_PP4_SOURCES, MANUAL_CRITERIA
from backend.variant_processing.variant_input import normalize_variant_input
from backend.variant_processing.variant_type import infer_variant_type

def _number(evidence: Dict[str, Any], key: str) -> Optional[float]:
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
        try:
            lr = float(component.get("likelihood_ratio"))
        except (TypeError, ValueError):
            raise ValueError(f"BS4 LR component {index} requires a numeric LR") from None
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


def _pp4_value_and_scale(evidence: Dict[str, Any]) -> tuple[Optional[float], str]:
    """Read the current PP4 input while retaining legacy LR audit records."""
    if evidence.get("clinical_lr_value") not in (None, ""):
        return _number(evidence, "clinical_lr_value"), str(
            evidence.get("clinical_lr_scale") or "lr"
        ).strip().lower()
    return _number(evidence, "combined_clinical_lr"), "lr"


def _pp4_strength(evidence: Dict[str, Any], gene: str) -> Optional[str]:
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


def _bp5_strength(evidence: Dict[str, Any], gene: str) -> Optional[str]:
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


def _pp4_source_is_reviewed(evidence: Dict[str, Any]) -> bool:
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


def _pp4_source_is_recorded(evidence: Dict[str, Any]) -> bool:
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


def suggest_strength(
    code: str,
    evidence: Dict[str, Any],
    *,
    variant_context: Mapping[str, Any] | None = None,
    base_criteria: Sequence[Mapping[str, Any]] = (),
) -> Optional[str]:
    gene = resolve_policy_gene(
        str((variant_context or {}).get("gene") or "").strip().upper() or None
    )
    clinical_thresholds = clinical_lr_thresholds(gene)
    pp4_thresholds = clinical_thresholds["pp4"]
    bp5_thresholds = clinical_thresholds["bp5"]
    splice_low = spliceai_thresholds(gene)["bp4"]
    if code in {"PP4", "BP5"}:
        data_summary = (evidence.get("clinical_data_summary") or "").strip()
        if not data_summary or not _pp4_source_is_reviewed(evidence):
            return None
        return (
            _pp4_strength(evidence, gene)
            if code == "PP4"
            else _bp5_strength(evidence, gene)
        )

    if code in {"PS3", "BS3"}:
        expected_conclusion = "abnormal" if code == "PS3" else "normal"
        required_text_fields = (
            "assay_name",
            "source_citation",
            "calibration_summary",
            "variant_result_summary",
            "functional_reviewed_by",
        )
        if (
            evidence.get("assay_scope") not in {
                "protein_only",
                "combined_mrna_protein",
            }
            or evidence.get("functional_conclusion") != expected_conclusion
            or evidence.get("calibration_status") != "reviewed_under_enigma_vcep"
            or evidence.get("pathogenic_and_benign_controls_confirmed") is not True
            or any(
                not str(evidence.get(field) or "").strip()
                for field in required_text_fields
            )
        ):
            return None
        strength = evidence.get("curated_strength")
        return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

    if code in {"PVS1_RNA", "BP7_RNA"}:
        assay_scope = evidence.get("assay_scope")
        transcript_accession = (evidence.get("transcript_accession") or "").strip()
        tissue = (evidence.get("tissue_or_cell_type") or "").strip()
        nmd = evidence.get("nmd_assessed")
        if (
            assay_scope != "mrna_only"
            or not transcript_accession
            or not tissue
            or nmd not in {"yes", "no", "not_applicable"}
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
            strength = evidence.get("curated_strength")
            return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

        if evidence.get("rna_conclusion") != "no_damaging_effect":
            return None
        context_result = evaluate_bp7_rna_variant_context(
            variant_context, base_criteria
        )
        if not context_result["eligible"]:
            return None
        return "Strong"

    if code == "PVS1_INIT":
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
        }:
            return None
        nearest_start = (evidence.get("nearest_alternative_start") or "").strip()
        rationale = (evidence.get("initiation_flowchart_rationale") or "").strip()
        if not nearest_start or not rationale:
            return None
        strength = evidence.get("curated_strength")
        return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

    if code == "PS1_SPLICE":
        required_text_fields = [
            "reference_variant",
            "reference_classification_source",
            "vua_splice_event",
            "reference_splice_event",
            "ps1_splice_rationale",
        ]
        if any(not (evidence.get(field) or "").strip() for field in required_text_fields):
            return None
        if evidence.get("reference_classification") not in {
            "Pathogenic",
            "Likely Pathogenic",
        }:
            return None
        if evidence.get("same_splice_event_confirmed") is not True:
            return None
        if evidence.get("prediction_strength_comparison") not in {
            "similar",
            "stronger",
        }:
            return None
        strength = evidence.get("curated_strength")
        return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

    if code == "PS1_PROTEIN":
        required_text_fields = [
            "reference_variant",
            "reference_p_notation",
            "classification_source",
            "ps1_protein_rationale",
        ]
        if any(not str(evidence.get(field) or "").strip() for field in required_text_fields):
            return None
        if variant_context is not None:
            try:
                assessed = normalize_variant_input(
                    gene,
                    str(variant_context.get("c_notation") or ""),
                    p_notation=str(variant_context.get("p_notation") or ""),
                )
                reference = normalize_variant_input(
                    gene,
                    str(evidence.get("reference_variant") or ""),
                    p_notation=str(evidence.get("reference_p_notation") or ""),
                )
            except ValueError:
                return None
            if (
                infer_variant_type(assessed.c_notation, assessed.p_notation) != "missense"
                or infer_variant_type(reference.c_notation, reference.p_notation) != "missense"
                or assessed.p_notation != reference.p_notation
                or assessed.c_notation == reference.c_notation
            ):
                return None
        if evidence.get("reference_classification") not in {
            "Pathogenic", "Likely Pathogenic"
        }:
            return None
        classification_verification = evidence.get("classification_verification")
        if classification_verification not in {
            "external_vcep_assertion",
            "locally_recurated_under_enigma_vcep",
            "enigma_st7_v1_2_reference_set",
        }:
            return None
        if classification_verification == "enigma_st7_v1_2_reference_set":
            try:
                normalized_reference = normalize_variant_input(
                    gene,
                    str(evidence.get("reference_variant") or ""),
                    p_notation=str(evidence.get("reference_p_notation") or ""),
                )
                registry_reference = lookup_ps1_reference_variant(
                    gene,
                    normalized_reference.c_notation,
                )
            except ValueError:
                return None
            if not registry_reference:
                return None
            registry_conditions = (
                registry_reference.get("classification_basis")
                == "enigma_st7_v1_2_reference_set"
                and registry_reference.get("reference_status") == "approved"
                and registry_reference.get("classification")
                == evidence.get("reference_classification")
                and registry_reference.get("p_notation")
                == normalized_reference.p_notation
                and registry_reference.get("reference_splice_evidence_status")
                in {"none_identified", "normal"}
                and registry_reference.get("classification_ps1_dependency_used") is False
            )
            if not registry_conditions:
                return None
        if (
            evidence.get("same_missense_confirmed") is not True
            or evidence.get("different_nucleotide_change_confirmed") is not True
            or evidence.get("splice_source_check_completed") is not True
        ):
            return None
        sources_checked = evidence.get("splice_sources_checked")
        if (
            not isinstance(sources_checked, list)
            or not set(PS1_SPLICE_SOURCES).issubset(
                {str(source).strip() for source in sources_checked}
            )
        ):
            return None
        if evidence.get("vua_confirmed_splice_status") not in {
            "none_identified", "normal"
        } or evidence.get("reference_confirmed_splice_status") not in {
            "none_identified", "normal"
        }:
            return None
        vua_score = _number(evidence, "vua_spliceai_score")
        reference_score = _number(evidence, "reference_spliceai_score")
        if (
            vua_score is None or reference_score is None
            or vua_score > splice_low or reference_score > splice_low
        ):
            return None
        dependency_status = evidence.get("reference_classification_used_ps1")
        if dependency_status not in {"no", "yes"}:
            return None
        if dependency_status == "yes":
            if (
                not str(evidence.get("reference_ps1_dependency_reference") or "").strip()
                or evidence.get("direct_reciprocal_dependency_excluded") is not True
            ):
                return None
        return (
            "Strong"
            if evidence["reference_classification"] == "Pathogenic"
            else "Moderate"
        )

    if code == "PS4":
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

    if code in {"PM3", "BS2"}:
        if evidence.get("cooccurring_variant_classification_basis") != "vcep_specifications":
            return None
        if (
            code == "PM3"
            and evidence.get("vua_benign_population_review") != "does_not_meet"
        ):
            return None
        points = _number(evidence, "evidence_points")
        if points is None or points < 1:
            return None
        if points >= 4:
            return "Strong"
        if points >= 2:
            return "Moderate"
        return "Supporting"

    likelihood_ratio = _number(evidence, "likelihood_ratio")
    if code == "BS4":
        likelihood_ratio, _, _ = evaluate_bs4_likelihood_ratio(evidence)
    if likelihood_ratio is None or likelihood_ratio < 0:
        return None
    if code == "PP1":
        if likelihood_ratio >= pp4_thresholds["very_strong_min_inclusive"] and evidence.get("very_strong_effect_basis") in {
            "predicted_protein",
            "predicted_splicing",
            "experimental_protein",
            "experimental_splicing",
        }:
            return "Very Strong"
        if likelihood_ratio >= pp4_thresholds["strong_min_inclusive"]:
            return "Strong"
        if likelihood_ratio >= pp4_thresholds["moderate_min_inclusive"]:
            return "Moderate"
        if likelihood_ratio >= pp4_thresholds["supporting_min_inclusive"]:
            return "Supporting"
    elif code == "BS4":
        if likelihood_ratio <= bp5_thresholds["very_strong_max_inclusive"]:
            return "Very Strong"
        if likelihood_ratio <= bp5_thresholds["strong_max_inclusive"]:
            return "Strong"
        if likelihood_ratio <= bp5_thresholds["moderate_max_inclusive"]:
            return "Moderate"
        if likelihood_ratio <= bp5_thresholds["supporting_max_inclusive"]:
            return "Supporting"
    return None


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

