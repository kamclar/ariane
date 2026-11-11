"""Pure ENIGMA population-frequency criteria.

This module assigns BA1, BS1 and PM2 only from evidence supplied by the
population-frequency provider. It performs no file, network or registry I/O.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from backend.population_frequency.coverage import frequency_depth_ok, frequency_qc_ok
from backend.population_frequency.lookup import scored_frequency_label
from backend.domain.indel import is_indel_allele
from backend.population_frequency.utils import as_float, as_int


def pm2_not_applicable_decision(
    variant_type: str,
    *,
    gene: str | None = None,
    c_notation: str | None = None,
    policy: Mapping[str, Any] | None = None,
    appendix_g_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Return the policy-defined PM2 N/A decision without querying gnomAD."""
    if policy is None:
        return None
    frequency_policy = policy.get("frequency_criteria") or {}
    pm2_policy = frequency_policy.get("pm2") or {}
    excluded_types = {
        str(value).lower() for value in pm2_policy.get("excluded_variant_types", [])
    }
    c_allele_is_indel = is_indel_allele(c_notation)
    if variant_type.lower() not in excluded_types and not c_allele_is_indel:
        return None
    if (
        not appendix_g_evidence
        or appendix_g_evidence.get("pm2_applicability") != "not_applicable"
    ):
        return None
    return {
        "applies": False,
        "strength": None,
        "points": 0,
        "reason": str(appendix_g_evidence.get("reason") or "PM2 is not applicable."),
        "source": str(
            appendix_g_evidence.get("source")
            or frequency_policy.get("source_url", "")
        ),
    }


def _founder_exception(gnomad_data: Mapping[str, Any]) -> Mapping[str, Any]:
    founder = gnomad_data.get("founder_exception")
    if isinstance(founder, Mapping):
        return founder
    return {
        "status": "unavailable",
        "is_pathogenic_founder": None,
        "reason": "the population evidence does not contain a founder-exception result",
    }


def _ba1_bs1_not_applied(reason: str, **details: Any) -> dict[str, Any]:
    return {
        "_gnomad_info": {
            "applies": False,
            "reason": f"BA1/BS1 not applied: {reason}",
            **details,
        }
    }


@dataclass(frozen=True)
class _FrequencyPolicy:
    source: Mapping[str, Any]
    ba1: Mapping[str, Any]
    bs1_strong: Mapping[str, Any]
    bs1_supporting: Mapping[str, Any]
    pm2: Mapping[str, Any]
    ba1_threshold: float
    bs1_strong_threshold: float
    bs1_supporting_threshold: float
    frequency_depth_threshold: float
    pm2_depth_threshold: float
    minimum_observations: int


def _parse_frequency_policy(policy: Mapping[str, Any]) -> _FrequencyPolicy | None:
    frequency = policy.get("frequency_criteria") or {}
    ba1 = frequency.get("ba1") or {}
    bs1 = frequency.get("bs1") or {}
    bs1_strong = bs1.get("strong") or {}
    bs1_supporting = bs1.get("supporting") or {}
    pm2 = frequency.get("pm2") or {}
    observations = frequency.get("outbred_observation_count") or {}
    values = (
        as_float(ba1.get("threshold")),
        as_float(bs1_strong.get("threshold")),
        as_float(bs1_supporting.get("lower_threshold")),
        as_float(ba1.get("minimum_mean_depth")),
        as_float(pm2.get("minimum_mean_depth")),
        as_int(observations.get("minimum_observations_for_ba1_bs1")),
    )
    if any(value is None for value in values):
        return None
    return _FrequencyPolicy(
        source=frequency,
        ba1=ba1,
        bs1_strong=bs1_strong,
        bs1_supporting=bs1_supporting,
        pm2=pm2,
        ba1_threshold=float(values[0]),
        bs1_strong_threshold=float(values[1]),
        bs1_supporting_threshold=float(values[2]),
        frequency_depth_threshold=float(values[3]),
        pm2_depth_threshold=float(values[4]),
        minimum_observations=int(values[5]),
    )


def _founder_screen(
    gnomad_data: Mapping[str, Any],
    policy: _FrequencyPolicy,
    max_af: float,
    metric_note: str,
    af_pct: str,
) -> tuple[dict[str, Any] | None, str]:
    founder = _founder_exception(gnomad_data)
    if founder.get("is_pathogenic_founder") is True:
        if max_af > policy.ba1_threshold:
            code, excluded_policy, threshold = "BA1", policy.ba1, policy.ba1_threshold
        elif max_af > policy.bs1_strong_threshold:
            code = "BS1_Strong"
            excluded_policy = policy.bs1_strong
            threshold = policy.bs1_strong_threshold
        else:
            code = "BS1_Supporting"
            excluded_policy = policy.bs1_supporting
            threshold = policy.bs1_supporting_threshold
        reason = (
            f"gnomAD {metric_note} {af_pct} exceeds the "
            f"{code.replace('_', ' ')} threshold {threshold * 100:g}%, but "
            f"{code.split('_')[0]} was not applied and added no points: "
            "ENIGMA v1.2 excludes well-established pathogenic founder variants; "
            f"{founder.get('reason')}"
        )
        return {
            "_excluded_criteria": {
                code: {
                    "applies": False,
                    "strength": excluded_policy.get("strength"),
                    "points": 0,
                    "reason": reason,
                    "source": policy.source.get("source_url", ""),
                }
            },
            "_gnomad_info": {
                "applies": False,
                "reason": reason,
                "founder_exception": founder,
            },
        }, ""
    if founder.get("status") not in {"not_listed", "reviewed_not_found"}:
        return _ba1_bs1_not_applied(
            "the pathogenic founder exception screen is unavailable or invalid; "
            f"{founder.get('reason')}",
            founder_exception=founder,
        ), ""
    note = (
        "; no pathogenic-founder exception match in snapshot "
        f"{founder.get('snapshot_version') or 'unknown'}"
    )
    return None, note


def _applied_frequency_criterion(
    policy: _FrequencyPolicy,
    max_af: float,
    metric_note: str,
    af_pct: str,
    founder_note: str,
) -> dict[str, Any]:
    policy_note = (
        "; founder and other non-scoring population groups excluded per "
        f"{policy.source.get('source') or 'the active gene-specific policy'}"
    )
    if max_af > policy.ba1_threshold:
        return {
            "BA1": {
                "applies": True,
                "strength": policy.ba1["strength"],
                "points": policy.ba1["points"],
                "reason": (
                    f"gnomAD {metric_note} {af_pct} > {policy.ba1_threshold * 100:g}% "
                    f"- Stand-alone Benign{policy_note}{founder_note}"
                ),
            }
        }
    if max_af > policy.bs1_strong_threshold:
        return {
            "BS1_Strong": {
                "applies": True,
                "strength": policy.bs1_strong["strength"],
                "points": policy.bs1_strong["points"],
                "reason": (
                    f"gnomAD {metric_note} {af_pct} > "
                    f"{policy.bs1_strong_threshold * 100:g}%{policy_note}{founder_note}"
                ),
            }
        }
    return {
        "BS1_Supporting": {
            "applies": True,
            "strength": policy.bs1_supporting["strength"],
            "points": policy.bs1_supporting["points"],
            "reason": (
                f"gnomAD {metric_note} {af_pct} > "
                f"{policy.bs1_supporting_threshold * 100:g}%{policy_note}{founder_note}"
            ),
        }
    }


def _evaluate_ba1_bs1(
    gnomad_data: Mapping[str, Any], policy: _FrequencyPolicy
) -> tuple[dict[str, Any], bool]:
    criteria: dict[str, Any] = {}
    max_af = as_float(gnomad_data.get("max_af"))
    metric = gnomad_data.get("frequency_metric") or "frequency"
    if max_af is not None and metric != "faf95":
        criteria.update(
            _ba1_bs1_not_applied(
                "ENIGMA-compatible non-cancer FAF95 is unavailable; "
                f"{metric} cannot be used as a fallback"
            )
        )
        max_af = None

    if max_af is not None:
        counts = [
            as_int(dataset.get("non_founder_allele_count"))
            for dataset in gnomad_data.get("datasets", {}).values()
            if dataset.get("status") == "found"
            and as_float(dataset.get("max_af")) == max_af
        ]
        if max_af > policy.bs1_supporting_threshold and (
            not counts or any(count is None for count in counts)
        ):
            return _ba1_bs1_not_applied(
                "the outbred observation count is missing from the approved gnomAD result"
            ), True
        if max_af > policy.bs1_supporting_threshold and all(
            count < policy.minimum_observations
            for count in counts
            if count is not None
        ):
            return _ba1_bs1_not_applied(
                "a single observation in an ENIGMA outbred population is not informative"
            ), True
        if max_af > policy.bs1_supporting_threshold:
            if not frequency_qc_ok(gnomad_data):
                return _ba1_bs1_not_applied(
                    "the contributing gnomAD record did not pass dataset QC filters"
                ), True
            if not frequency_depth_ok(gnomad_data, policy.frequency_depth_threshold):
                return _ba1_bs1_not_applied(
                    "active policy requires mean read depth "
                    f">= {policy.frequency_depth_threshold:g}"
                ), True
            metric_note = scored_frequency_label(gnomad_data)
            af_pct = f"{max_af * 100:.6g}%"
            founder_result, founder_note = _founder_screen(
                gnomad_data, policy, max_af, metric_note, af_pct
            )
            if founder_result is not None:
                return founder_result, True
            return _applied_frequency_criterion(
                policy, max_af, metric_note, af_pct, founder_note
            ), True
        if gnomad_data.get("found"):
            return criteria, True

    if gnomad_data.get("found") and max_af is None:
        failed_qc = any(
            dataset.get("status") == "found"
            and dataset.get("quality_filter_passed") is False
            for dataset in gnomad_data.get("datasets", {}).values()
        )
        if "_gnomad_info" not in criteria:
            reason = (
                "the record failed dataset QC filters"
                if failed_qc
                else "ENIGMA-compatible non-cancer FAF95 is unavailable"
            )
            criteria.update(
                _ba1_bs1_not_applied(
                    f"the variant is present in gnomAD, but {reason}"
                )
            )
    return criteria, False


def _pm2_status_reason(status: str, depth_threshold: float) -> str | None:
    return {
        "cache_missing": "local gnomAD cache missing or incomplete - PM2 not applied",
        "cache_untrusted": (
            "local gnomAD cache is not an approved snapshot - frequency criteria not applied"
        ),
        "partial": "local gnomAD lookup partial - PM2 not applied",
        "no_coordinates": "No genomic coordinates for required gnomAD lookup - PM2 not applied",
        "outside_cached_region": "Variant outside cached panel gnomAD regions - PM2 not applied",
        "absent_without_sufficient_coverage": (
            "Absent from local gnomAD cache but coverage mean depth is below "
            f"{depth_threshold:g} or missing - PM2 not applied"
        ),
        "not_queried": "gnomAD not queried - PM2 not applied",
        "absent_v2_only": (
            "gnomAD v2.1.1 absence confirmed but v3.1.2 coverage insufficient - PM2 not applied"
        ),
        "policy_unavailable": (
            "Gene-specific gnomAD policy is unavailable - frequency criteria not applied"
        ),
        "coverage_release_incompatible": (
            "PM2 not applied: at least one coverage source is not proven "
            "compatible with its frequency callset"
        ),
        "pm2_coverage_method_unresolved": (
            "PM2 not applied: ENIGMA does not define the width of the coverage "
            "region around the variant"
        ),
    }.get(status)


def _evaluate_pm2(
    gnomad_data: Mapping[str, Any],
    variant_type: str,
    c_notation: str | None,
    policy: _FrequencyPolicy,
    appendix_g_evidence: Mapping[str, Any] | None,
) -> dict[str, Any]:
    is_indel = (
        variant_type.lower()
        in {str(value).lower() for value in policy.pm2.get("excluded_variant_types", [])}
        or is_indel_allele(c_notation)
    )
    if is_indel:
        structural = appendix_g_evidence or {}
        if str(structural.get("pm2_applicability") or "unavailable") == "applied":
            return {}
        return {
            "PM2": {
                "applies": False,
                "strength": None,
                "points": 0,
                "reason": str(
                    structural.get("reason")
                    or "PM2 is unavailable: the ENIGMA Appendix G indel-size and structural population path was not evaluated."
                ),
                "source": str(
                    structural.get("source") or policy.source.get("source_url", "")
                ),
            }
        }
    coverage_method = gnomad_data.get("pm2_coverage_method")
    if (
        not isinstance(coverage_method, Mapping)
        or coverage_method.get("automatic_assignment_allowed") is not True
    ):
        return {
            "PM2": {
                "applies": False,
                "strength": None,
                "points": 0,
                "reason": (
                    "PM2 not applied: the coverage-region method is not approved "
                    "for automatic ENIGMA classification; "
                    + str(
                        (coverage_method or {}).get("reason")
                        or "no explicit coverage-method review is present"
                    )
                ),
                "source": str(policy.source.get("source_url", "")),
            }
        }
    if gnomad_data.get("pm2_absence_established"):
        datasets_note = gnomad_data.get("pm2_datasets_note", "v2.1.1 + v3.1.2")
        founder_note = (
            "; observations confined to excluded founder/non-scoring populations "
            "were not treated as presence in an ENIGMA outbred population"
            if gnomad_data.get("founder_context_only_observed")
            else ""
        )
        return {
            "PM2_Supporting": {
                "applies": True,
                "strength": policy.pm2["strength"],
                "points": policy.pm2["points"],
                "reason": (
                    f"Absent from gnomAD {datasets_note} non-cancer callset(s) "
                    "in all policy-defined non-founder populations with mean depth "
                    f">= {policy.pm2_depth_threshold:g} across the variant reference span"
                    f"{founder_note}"
                ),
            }
        }
    reason = _pm2_status_reason(
        str(gnomad_data.get("status") or "not_queried"),
        policy.pm2_depth_threshold,
    )
    return (
        {
            "PM2": {
                "applies": False,
                "strength": None,
                "points": 0,
                "reason": reason,
            }
        }
        if reason
        else {}
    )


def evaluate_frequency_criteria(
    gnomad_data: Mapping[str, Any],
    variant_type: str,
    gene: str | None = None,
    c_notation: str | None = None,
    policy: Mapping[str, Any] | None = None,
    appendix_g_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate BA1, BS1 and PM2 from validated provider evidence."""
    if policy is None:
        candidate = gnomad_data.get("classification_policy")
        policy = candidate if isinstance(candidate, Mapping) else None
    if policy is None or policy.get("policy_id") != gnomad_data.get("policy_id"):
        return {
            "_gnomad_info": {
                "applies": False,
                "reason": (
                    "Frequency criteria not applied: no matching active "
                    f"gene-specific gnomAD policy for {gene or 'unspecified gene'}"
                ),
            }
        }
    parsed_policy = _parse_frequency_policy(policy)
    if parsed_policy is None:
        return {
            "_gnomad_info": {
                "applies": False,
                "reason": "Frequency criteria not applied: active policy is incomplete",
            }
        }

    criteria: dict[str, Any] = {}
    pm2_not_applicable = pm2_not_applicable_decision(
        variant_type,
        gene=gene,
        c_notation=c_notation,
        policy=policy,
        appendix_g_evidence=appendix_g_evidence,
    )
    if pm2_not_applicable:
        criteria["PM2"] = pm2_not_applicable

    frequency_criteria, terminal = _evaluate_ba1_bs1(gnomad_data, parsed_policy)
    criteria.update(frequency_criteria)
    if terminal or pm2_not_applicable:
        return criteria
    criteria.update(
        _evaluate_pm2(
            gnomad_data,
            variant_type,
            c_notation,
            parsed_policy,
            appendix_g_evidence,
        )
    )
    return criteria
