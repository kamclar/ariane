"""Decide whether provider evidence is complete enough to publish a result.

Provider nodes preserve unavailable and ambiguous values so their provenance is
auditable.  This module is the application-level publication gate.  It keeps a
missing required value from being interpreted as a completed negative result.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from backend.classification_dag.domain import NormalizedVariant
from backend.gene_policy import spliceai_thresholds
from backend.modules.spliceai_policy import (
    spliceai_failure_is_retryable,
    spliceai_required_for_classification,
)
from backend.modules.utils import (
    get_amino_acid_interval,
    overlapping_functional_domains,
)


_PROTEIN_PREDICTION_TYPES = frozenset({
    "missense",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_delins",
})

# These are completed observations from one required gnomAD dataset.  They do
# not all support a criterion, but none represents a failed lookup.
_COMPLETE_GNOMAD_DATASET_STATUSES = frozenset({
    "found",
    "absent",
    "absent_in_non_founder_populations",
    "filtered_record",
})

_COMPLETE_STRUCTURAL_PM2_STATUSES = frozenset({
    "applied",
    "not_applicable",
    "not_met",
})


@dataclass(frozen=True)
class RequiredEvidenceGap:
    source: str
    status: str
    reason: str
    retryable: bool = False


def _spliceai_gap(
    variant: NormalizedVariant,
    artifacts: Mapping[str, Any],
) -> RequiredEvidenceGap | None:
    if not spliceai_required_for_classification(variant.variant_type):
        return None

    splice_status = dict(artifacts.get("spliceai_status") or {})
    score = artifacts.get("spliceai_score")
    status = str(splice_status.get("status") or "unavailable")
    assessed_complete = status == "ok" and score is not None
    reference_required = bool(splice_status.get("reference_lookup_required"))
    reference_complete = bool(
        splice_status.get("reference_lookup_complete", True)
    )
    if assessed_complete and (not reference_required or reference_complete):
        return None

    failed_status: Mapping[str, Any] = splice_status
    if assessed_complete:
        failed_references = [
            item_status
            for item_status in dict(
                splice_status.get("reference_variant_statuses") or {}
            ).values()
            if str(item_status.get("status") or "") != "ok"
            or item_status.get("score") is None
        ]
        failed_status = failed_references[0] if failed_references else {}
        status = str(failed_status.get("status") or "unavailable")
        reason = str(
            failed_status.get("reason")
            or "A candidate protein PS1 reference has no SpliceAI result"
        )
    else:
        reason = str(
            splice_status.get("reason")
            or "The configured SpliceAI source returned no score"
        )
    retryable = bool(failed_status.get("retryable")) or (
        spliceai_failure_is_retryable(status, reason)
    )
    return RequiredEvidenceGap("SpliceAI", status, reason, retryable)


def _gnomad_gap(
    variant: NormalizedVariant,
    artifacts: Mapping[str, Any],
) -> RequiredEvidenceGap | None:
    if variant.variant_type.lower() in {"exon_deletion", "exon_duplication"}:
        return None

    gnomad = artifacts.get("gnomad_data")
    if not isinstance(gnomad, Mapping):
        return RequiredEvidenceGap(
            "gnomAD",
            "not_queried",
            "The population-frequency provider returned no structured result",
        )

    policy = gnomad.get("classification_policy")
    if not isinstance(policy, Mapping) or (
        policy.get("policy_id") != gnomad.get("policy_id")
    ):
        return RequiredEvidenceGap(
            "gnomAD",
            str(gnomad.get("status") or "policy_unavailable"),
            "No matching active gene-specific population policy was available",
        )

    frequency_policy = policy.get("frequency_criteria")
    pm2_policy = (
        frequency_policy.get("pm2")
        if isinstance(frequency_policy, Mapping)
        else None
    )
    required_datasets = (
        pm2_policy.get("required_absence_dataset_runtime_keys")
        if isinstance(pm2_policy, Mapping)
        else None
    )
    if not isinstance(required_datasets, list) or not required_datasets:
        return RequiredEvidenceGap(
            "gnomAD",
            "policy_unavailable",
            "The active population policy does not define its required datasets",
        )

    datasets = gnomad.get("datasets")
    if not isinstance(datasets, Mapping):
        datasets = {}
    for dataset_key in required_datasets:
        item = datasets.get(dataset_key)
        item_status = str(
            item.get("status") if isinstance(item, Mapping) else "not_queried"
        )
        if item_status not in _COMPLETE_GNOMAD_DATASET_STATUSES:
            errors = item.get("errors") if isinstance(item, Mapping) else None
            detail = (
                str(errors[0])
                if isinstance(errors, list) and errors
                else "the required dataset lookup did not complete"
            )
            return RequiredEvidenceGap(
                "gnomAD",
                item_status or "unavailable",
                f"Population dataset {dataset_key} is incomplete: {detail}",
            )

    return None


def _structural_population_gap(
    artifacts: Mapping[str, Any],
) -> RequiredEvidenceGap | None:
    value = artifacts.get("exon_cnv_result")
    if not isinstance(value, Mapping):
        return RequiredEvidenceGap(
            "ENIGMA Appendix G population evidence",
            "not_queried",
            "The structural population provider returned no structured result",
        )
    status = str(value.get("pm2_applicability") or "unavailable")
    if status in _COMPLETE_STRUCTURAL_PM2_STATUSES:
        return None
    return RequiredEvidenceGap(
        "ENIGMA Appendix G population evidence",
        status,
        str(
            value.get("reason")
            or "The structural population decision path did not complete"
        ),
    )


def _bayesdel_gap(
    variant: NormalizedVariant,
    artifacts: Mapping[str, Any],
) -> RequiredEvidenceGap | None:
    if variant.variant_type.lower() not in _PROTEIN_PREDICTION_TYPES:
        return None

    spliceai_score = artifacts.get("spliceai_score")
    if spliceai_score is None:
        # The SpliceAI gap is reported first and prevents this branch from being
        # interpreted without its preceding Figure 1A decision.
        return None
    if float(spliceai_score) >= spliceai_thresholds(variant.gene)["pp3"]:
        return None

    protein_interval = get_amino_acid_interval(variant.p_notation)
    if protein_interval is None:
        return RequiredEvidenceGap(
            "protein consequence",
            "unresolved_interval",
            "The complete amino-acid interval required by Figure 1A could not be determined",
        )
    if not overlapping_functional_domains(variant.gene, protein_interval):
        return None

    status_details = dict(artifacts.get("bayesdel_status") or {})
    status = str(status_details.get("status") or "unavailable")
    score = artifacts.get("bayesdel_score")
    if status == "ok" and score is not None:
        return None
    reason = str(
        status_details.get("reason")
        or "The configured BayesDel_noAF source returned no unambiguous score"
    )
    return RequiredEvidenceGap(
        "BayesDel_noAF",
        status,
        reason,
        retryable=status == "api_error",
    )


def first_required_evidence_gap(
    variant: NormalizedVariant,
    artifacts: Mapping[str, Any],
) -> RequiredEvidenceGap | None:
    """Return the first source gap that prevents publishing a classification."""
    gap = _spliceai_gap(variant, artifacts)
    if gap is not None:
        return gap

    is_exon_cnv = variant.variant_type.lower() in {
        "exon_deletion", "exon_duplication",
    }
    if not is_exon_cnv:
        gap = _gnomad_gap(variant, artifacts)
        if gap is not None:
            return gap

    structural = artifacts.get("exon_cnv_result")
    structural_path_required = is_exon_cnv or (
        isinstance(structural, Mapping) and structural.get("is_indel") is True
    )
    if structural_path_required:
        gap = _structural_population_gap(artifacts)
        if gap is not None:
            return gap

    return _bayesdel_gap(variant, artifacts)
