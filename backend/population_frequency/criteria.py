"""Pure ENIGMA population-frequency criteria.

This module assigns BA1, BS1 and PM2 only from evidence supplied by the
population-frequency provider. It performs no file, network or registry I/O.
"""

from __future__ import annotations

from typing import Any, Mapping

from backend.population_frequency.coverage import frequency_depth_ok, frequency_qc_ok
from backend.population_frequency.lookup import scored_frequency_label
from backend.population_frequency.indel_size import is_indel_allele
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


def evaluate_frequency_criteria(
    gnomad_data: Mapping[str, Any],
    variant_type: str,
    gene: str | None = None,
    c_notation: str | None = None,
    policy: Mapping[str, Any] | None = None,
    appendix_g_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate BA1, BS1 and PM2 from validated provider evidence."""
    criteria: dict[str, Any] = {}
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

    frequency_policy = policy.get("frequency_criteria") or {}
    ba1_policy = frequency_policy.get("ba1") or {}
    bs1_policy = frequency_policy.get("bs1") or {}
    bs1_strong_policy = bs1_policy.get("strong") or {}
    bs1_supporting_policy = bs1_policy.get("supporting") or {}
    pm2_policy = frequency_policy.get("pm2") or {}
    observation_policy = frequency_policy.get("outbred_observation_count") or {}
    ba1_threshold = as_float(ba1_policy.get("threshold"))
    bs1_strong_threshold = as_float(bs1_strong_policy.get("threshold"))
    bs1_supporting_threshold = as_float(bs1_supporting_policy.get("lower_threshold"))
    frequency_depth_threshold = as_float(ba1_policy.get("minimum_mean_depth"))
    pm2_depth_threshold = as_float(pm2_policy.get("minimum_mean_depth"))
    minimum_ba1_bs1_observations = as_int(
        observation_policy.get("minimum_observations_for_ba1_bs1")
    )
    required_values = (
        ba1_threshold,
        bs1_strong_threshold,
        bs1_supporting_threshold,
        frequency_depth_threshold,
        pm2_depth_threshold,
        minimum_ba1_bs1_observations,
    )
    if any(value is None for value in required_values):
        return {
            "_gnomad_info": {
                "applies": False,
                "reason": "Frequency criteria not applied: active policy is incomplete",
            }
        }

    # The checks above make the numeric policy values non-null.
    ba1_threshold = float(ba1_threshold)
    bs1_strong_threshold = float(bs1_strong_threshold)
    bs1_supporting_threshold = float(bs1_supporting_threshold)
    frequency_depth_threshold = float(frequency_depth_threshold)
    pm2_depth_threshold = float(pm2_depth_threshold)
    minimum_ba1_bs1_observations = int(minimum_ba1_bs1_observations)

    pm2_not_applicable = pm2_not_applicable_decision(
        variant_type,
        gene=gene,
        c_notation=c_notation,
        policy=policy,
        appendix_g_evidence=appendix_g_evidence,
    )
    if pm2_not_applicable:
        criteria["PM2"] = pm2_not_applicable

    status = str(gnomad_data.get("status") or "not_queried")
    max_af = as_float(gnomad_data.get("max_af"))
    metric = gnomad_data.get("frequency_metric") or "frequency"
    if max_af is not None and metric != "faf95":
        criteria["_gnomad_info"] = {
            "applies": False,
            "reason": (
                "BA1/BS1 not applied: ENIGMA-compatible non-cancer FAF95 "
                f"is unavailable; {metric} cannot be used as a fallback"
            ),
        }
        max_af = None

    if max_af is not None:
        af_pct = f"{max_af * 100:.6g}%"
        metric_note = scored_frequency_label(gnomad_data)
        population_policy_note = (
            "; founder and other non-scoring population groups excluded per "
            f"{frequency_policy.get('source') or 'the active gene-specific policy'}"
        )
        contributing_observation_counts = [
            as_int(dataset.get("non_founder_allele_count"))
            for dataset in gnomad_data.get("datasets", {}).values()
            if dataset.get("status") == "found"
            and as_float(dataset.get("max_af")) == max_af
        ]
        if max_af > bs1_supporting_threshold and (
            not contributing_observation_counts
            or any(count is None for count in contributing_observation_counts)
        ):
            criteria["_gnomad_info"] = {
                "applies": False,
                "reason": (
                    "BA1/BS1 not applied: the outbred observation count is "
                    "missing from the approved gnomAD result"
                ),
            }
            return criteria
        if max_af > bs1_supporting_threshold and all(
            count < minimum_ba1_bs1_observations
            for count in contributing_observation_counts
            if count is not None
        ):
            criteria["_gnomad_info"] = {
                "applies": False,
                "reason": (
                    "BA1/BS1 not applied: a single observation in an ENIGMA "
                    "outbred population is not informative"
                ),
            }
            return criteria

        if max_af > bs1_supporting_threshold:
            if not frequency_qc_ok(gnomad_data):
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": (
                        "BA1/BS1 not applied: the contributing gnomAD record "
                        "did not pass dataset QC filters"
                    ),
                }
                return criteria
            if not frequency_depth_ok(gnomad_data, frequency_depth_threshold):
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": (
                        "BA1/BS1 not applied: active policy requires mean read "
                        f"depth >= {frequency_depth_threshold:g}"
                    ),
                }
                return criteria

            founder = _founder_exception(gnomad_data)
            if founder.get("is_pathogenic_founder") is True:
                if max_af > ba1_threshold:
                    excluded_code = "BA1"
                    excluded_policy = ba1_policy
                    threshold = ba1_threshold
                elif max_af > bs1_strong_threshold:
                    excluded_code = "BS1_Strong"
                    excluded_policy = bs1_strong_policy
                    threshold = bs1_strong_threshold
                else:
                    excluded_code = "BS1_Supporting"
                    excluded_policy = bs1_supporting_policy
                    threshold = bs1_supporting_threshold
                exclusion_reason = (
                    f"gnomAD {metric_note} {af_pct} exceeds the "
                    f"{excluded_code.replace('_', ' ')} threshold "
                    f"{threshold * 100:g}%, but {excluded_code.split('_')[0]} was "
                    "not applied and added no points: ENIGMA v1.2 excludes "
                    "well-established pathogenic founder variants; "
                    f"{founder.get('reason')}"
                )
                criteria["_excluded_criteria"] = {
                    excluded_code: {
                        "applies": False,
                        "strength": excluded_policy.get("strength"),
                        "points": 0,
                        "reason": exclusion_reason,
                        "source": frequency_policy.get("source_url", ""),
                    }
                }
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": exclusion_reason,
                    "founder_exception": founder,
                }
                return criteria
            if founder.get("status") == "unavailable":
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": (
                        "BA1/BS1 not applied: the pathogenic founder exception "
                        f"could not be checked; {founder.get('reason')}"
                    ),
                    "founder_exception": founder,
                }
                return criteria
            founder_note = (
                "; pathogenic-founder exception checked against snapshot "
                f"{founder.get('snapshot_version') or 'unknown'}"
            )
            if max_af > ba1_threshold:
                criteria["BA1"] = {
                    "applies": True,
                    "strength": ba1_policy["strength"],
                    "points": ba1_policy["points"],
                    "reason": (
                        f"gnomAD {metric_note} {af_pct} > {ba1_threshold * 100:g}% "
                        "- Stand-alone Benign"
                        f"{population_policy_note}{founder_note}"
                    ),
                }
                return criteria
            if max_af > bs1_strong_threshold:
                criteria["BS1_Strong"] = {
                    "applies": True,
                    "strength": bs1_strong_policy["strength"],
                    "points": bs1_strong_policy["points"],
                    "reason": (
                        f"gnomAD {metric_note} {af_pct} > "
                        f"{bs1_strong_threshold * 100:g}%"
                        f"{population_policy_note}{founder_note}"
                    ),
                }
                return criteria
            criteria["BS1_Supporting"] = {
                "applies": True,
                "strength": bs1_supporting_policy["strength"],
                "points": bs1_supporting_policy["points"],
                "reason": (
                    f"gnomAD {metric_note} {af_pct} > "
                    f"{bs1_supporting_threshold * 100:g}%"
                    f"{population_policy_note}{founder_note}"
                ),
            }
            return criteria
        if gnomad_data.get("found"):
            return criteria

    if gnomad_data.get("found") and max_af is None:
        failed_qc = any(
            dataset.get("status") == "found"
            and dataset.get("quality_filter_passed") is False
            for dataset in gnomad_data.get("datasets", {}).values()
        )
        criteria.setdefault(
            "_gnomad_info",
            {
                "applies": False,
                "reason": (
                    "BA1/BS1 not applied: the variant is present in gnomAD, but "
                    + (
                        "the record failed dataset QC filters"
                        if failed_qc
                        else "ENIGMA-compatible non-cancer FAF95 is unavailable"
                    )
                ),
            },
        )

    if pm2_not_applicable:
        return criteria
    is_pm2_indel = (
        variant_type.lower()
        in {str(value).lower() for value in pm2_policy.get("excluded_variant_types", [])}
        or is_indel_allele(c_notation)
    )
    if is_pm2_indel:
        structural_status = str(
            (appendix_g_evidence or {}).get("pm2_applicability") or "unavailable"
        )
        if structural_status == "applied":
            return criteria
        criteria["PM2"] = {
            "applies": False,
            "strength": None,
            "points": 0,
            "reason": str(
                (appendix_g_evidence or {}).get("reason")
                or "PM2 is unavailable: the ENIGMA Appendix G indel-size and structural population path was not evaluated."
            ),
            "source": str(
                (appendix_g_evidence or {}).get("source")
                or frequency_policy.get("source_url", "")
            ),
        }
        return criteria
    if gnomad_data.get("pm2_absence_established"):
        datasets_note = gnomad_data.get("pm2_datasets_note", "v2.1.1 + v3.1.2")
        founder_only_note = (
            "; observations confined to excluded founder/non-scoring populations "
            "were not treated as presence in an ENIGMA outbred population"
            if gnomad_data.get("founder_context_only_observed")
            else ""
        )
        criteria["PM2_Supporting"] = {
            "applies": True,
            "strength": pm2_policy["strength"],
            "points": pm2_policy["points"],
            "reason": (
                f"Absent from gnomAD {datasets_note} non-cancer callset(s) "
                "in all policy-defined non-founder populations with mean depth "
                f">= {pm2_depth_threshold:g} across the variant reference span"
                f"{founder_only_note}"
            ),
        }
        return criteria

    reason_by_status = {
        "cache_missing": "local gnomAD cache missing or incomplete - PM2 not applied",
        "cache_untrusted": (
            "local gnomAD cache is not an approved snapshot - frequency criteria not applied"
        ),
        "partial": "local gnomAD lookup partial - PM2 not applied",
        "no_coordinates": "No genomic coordinates for required gnomAD lookup - PM2 not applied",
        "outside_cached_region": "Variant outside cached panel gnomAD regions - PM2 not applied",
        "absent_without_sufficient_coverage": (
            "Absent from local gnomAD cache but coverage mean depth is below "
            f"{pm2_depth_threshold:g} or missing - PM2 not applied"
        ),
        "not_queried": "gnomAD not queried - PM2 not applied",
        "absent_v2_only": (
            "gnomAD v2.1.1 absence confirmed but v3.1.2 coverage insufficient - PM2 not applied"
        ),
        "policy_unavailable": (
            "Gene-specific gnomAD policy is unavailable - frequency criteria not applied"
        ),
    }
    if status in reason_by_status:
        criteria["PM2"] = {
            "applies": False,
            "strength": None,
            "points": 0,
            "reason": reason_by_status[status],
        }
    return criteria
