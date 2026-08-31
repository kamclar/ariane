"""Read-only lookup of population-frequency evidence from a bound repository."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from backend.population_frequency.coverage import aggregate_coverage, lookup_coverage_by_position
from backend.population_frequency.models import GnomadRepository
from backend.population_frequency.policy import (
    GNOMAD_LOCAL_DATASET_CONFIG,
    PM2_COVERAGE_METHOD_REVIEW,
    classification_policy_for_gene,
)
from backend.population_frequency.utils import as_float, as_int, coordinate_value, strip_chr, variant_id_from_coords


def coords_in_cached_region(repository: GnomadRepository, coords: Any | None, build: str) -> bool:
    if not coords:
        return False
    try:
        chrom = coordinate_value(coords, "chrom")
        pos = int(coordinate_value(coords, "pos"))
    except (AttributeError, KeyError, TypeError, ValueError):
        return False
    build_regions = (repository.metadata.get("regions", {}) or {}).get(build, {}) or {}
    padding = int(repository.metadata.get("region_padding_bp") or 0)
    if not build_regions:
        return False
    chrom_no = strip_chr(chrom)
    return any(
        chrom_no == strip_chr(region.get("chrom"))
        and int(region.get("start")) - padding <= pos <= int(region.get("end")) + padding
        for region in build_regions.values()
    )


def dataset_extraction_ok(
    repository: GnomadRepository,
    dataset_names: Sequence[str],
    coords: Any | None,
    build: str,
) -> bool:
    if not coords:
        return False
    try:
        chrom_no = strip_chr(coordinate_value(coords, "chrom"))
    except (AttributeError, KeyError, TypeError):
        return False
    for item in repository.metadata.get("extraction_log", []) or []:
        if item.get("dataset") in dataset_names and item.get("status") == "ok":
            item_chrom = item.get("chrom")
            if item_chrom is None or strip_chr(item_chrom) == chrom_no:
                return True
    return False


def _empty_callset() -> dict[str, Any]:
    return {
        "available": False, "ac": None, "an": None, "af": None,
        "ac_hom": None, "filters": [], "popmax_pop": None,
        "popmax_af": None, "faf95_max": None, "faf_any_max": None,
    }


def _record_frequency_value(
    record: Mapping[str, Any], scored_ancestries: Sequence[str]
) -> tuple[float | None, str | None, str | None]:
    values = {
        code: as_float((record.get("faf95_by_ancestry") or {}).get(code))
        for code in scored_ancestries
    }
    if not values or any(value is None for value in values.values()):
        return None, None, None
    population, value = max(values.items(), key=lambda item: item[1])
    return value, "faf95", population


def _record_to_callset_summary(
    record: Mapping[str, Any], scored_ancestries: Sequence[str]
) -> dict[str, Any]:
    frequency_value, _, population = _record_frequency_value(record, scored_ancestries)
    non_founder_ac = record.get("non_founder_ac_by_ancestry") or {}
    return {
        "available": True,
        "ac": as_int(record.get("ac")),
        "an": as_int(record.get("an")),
        "af": as_float(record.get("af")),
        "ac_hom": as_int(record.get("nhomalt")),
        "filters": [] if record.get("filter") in (None, ".", "PASS") else [record.get("filter")],
        "popmax_pop": record.get("popmax_pop"),
        "popmax_af": as_float(record.get("popmax_af")),
        "faf95_max": frequency_value,
        "faf95_pop": population,
        "faf95_scope": record.get("faf95_scope"),
        "faf95_method": record.get("faf95_method"),
        "faf95_by_ancestry": record.get("faf95_by_ancestry") or {},
        "non_founder_ac_by_ancestry": non_founder_ac,
        "non_founder_allele_count": sum(as_int(value) or 0 for value in non_founder_ac.values()),
        "non_founder_an_by_ancestry": record.get("non_founder_an_by_ancestry") or {},
        "non_founder_observed": record.get("non_founder_observed"),
        "excluded_population_context": record.get("excluded_population_context") or {},
        "faf_any_max": as_float(record.get("faf_any_max")),
    }


def population_frequency_audit(gnomad_data: Mapping[str, Any]) -> dict[str, Any]:
    policy = gnomad_data.get("frequency_policy") or {}
    scored_codes = list(policy.get("scored_non_founder_ancestries") or [])
    datasets = []
    for dataset_key, result in gnomad_data.get("datasets", {}).items():
        config = GNOMAD_LOCAL_DATASET_CONFIG.get(dataset_key, {})
        callset = result.get(config.get("callset", "")) or {}
        coverage = result.get("coverage") or {}
        scored_faf = callset.get("faf95_by_ancestry") or {}
        scored_ac = callset.get("non_founder_ac_by_ancestry") or {}
        scored_an = callset.get("non_founder_an_by_ancestry") or {}
        excluded = callset.get("excluded_population_context") or {}
        datasets.append({
            "dataset_key": dataset_key,
            "label": result.get("label") or config.get("label", ""),
            "status": result.get("status"),
            "scored_max_faf95": result.get("max_af"),
            "scored_max_population": callset.get("faf95_pop"),
            "scored_non_founder_populations": [
                {"code": code, "ac": scored_ac.get(code), "an": scored_an.get(code),
                 "faf95": scored_faf.get(code), "used_for_scoring": True}
                for code in scored_codes
            ],
            "excluded_population_context": [
                {"code": code, **value} for code, value in excluded.items()
            ],
            "coverage": {
                "mean_depth": coverage.get("mean_depth"),
                "measurement_passes_threshold": coverage.get(
                    "measurement_passes_threshold"
                ),
                "classification_compatible": coverage.get(
                    "classification_compatible"
                ),
                "compatibility_status": coverage.get("compatibility_status"),
                "compatibility_reason": coverage.get("compatibility_reason"),
                "coverage_scope": coverage.get("coverage_scope"),
            },
        })
    return {
        "status": gnomad_data.get("status"),
        "policy_id": gnomad_data.get("policy_id"),
        "policy_source": policy.get("source", ""),
        "policy_source_url": policy.get("source_url", ""),
        "scoring_rule": policy.get("scoring_rule", ""),
        "scored_non_founder_population_codes": scored_codes,
        "excluded_populations_are_context_only": True,
        "founder_context_only_observed": any(
            item.get("status") == "absent_in_non_founder_populations"
            for item in gnomad_data.get("datasets", {}).values()
        ),
        "pm2_coverage_method": gnomad_data.get("pm2_coverage_method") or {},
        "datasets": datasets,
    }


def scored_frequency_label(gnomad_data: Mapping[str, Any]) -> str:
    maximum = as_float(gnomad_data.get("max_af"))
    if maximum is None:
        return "non-cancer FAF95"
    matches = []
    for dataset_key, result in gnomad_data.get("datasets", {}).items():
        if as_float(result.get("max_af")) != maximum:
            continue
        config = GNOMAD_LOCAL_DATASET_CONFIG.get(dataset_key, {})
        callset = result.get(config.get("callset", "")) or {}
        population = str(callset.get("faf95_pop") or "unknown").upper()
        label = result.get("label") or config.get("label", dataset_key)
        matches.append(f"{population}, {label}")
    return "non-cancer FAF95 in " + ("; ".join(matches) or "ENIGMA non-founder populations")


def query_gnomad_dataset(
    repository: GnomadRepository,
    variant_id: str | None,
    coords: Any | None,
    config: Mapping[str, Any],
    coverage_threshold: float,
    scored_ancestries: Sequence[str],
) -> dict[str, Any]:
    result = {
        "status": "no_coordinates" if not variant_id or not coords else "not_queried",
        "variant_id": variant_id, "dataset": None, "label": config["label"],
        "assembly": config["assembly"], "found": None,
        "exomes": _empty_callset(), "genomes": _empty_callset(),
        "max_af": None, "frequency_metric": None, "quality_filter_passed": None,
        "filtered_records": [], "coverage": None, "errors": [],
    }
    if not variant_id or not coords:
        return result
    if not repository.variants:
        result["status"] = "cache_missing"
        result["errors"].append("local gnomAD cache not loaded")
        return result
    if repository.frequency_status != "approved_snapshot":
        result["status"] = "cache_untrusted"
        result["errors"].append(
            f"gnomAD cache mode {repository.frequency_status!r} is not approved for classification"
        )
        return result
    if not coords_in_cached_region(repository, coords, config["assembly"]):
        result["status"] = "outside_cached_region"
        result["errors"].append("variant coordinates outside cached panel regions")
        return result
    if not dataset_extraction_ok(repository, config["dataset_names"], coords, config["assembly"]):
        result["status"] = "dataset_not_available"
        result["errors"].append(
            "dataset extraction was not successful or is not documented in cache metadata"
        )
        return result

    candidate_keys = (
        variant_id, variant_id_from_coords(coords, with_chr=False),
        variant_id_from_coords(coords, with_chr=True),
    )
    all_records = []
    seen: set[int] = set()
    for key in candidate_keys:
        if key and key in repository.variants:
            for record in repository.variants[key]:
                record_id = id(record)
                if record_id not in seen:
                    all_records.append(record)
                    seen.add(record_id)
    dataset_records = [
        record for record in all_records
        if record.get("dataset") in config["dataset_names"]
        and record.get("build") == config["assembly"]
    ]
    if dataset_records:
        result["database_record_found"] = True
        result["dataset"] = dataset_records[0].get("dataset")
        result["filtered_records"] = [
            {"variant_id": record.get("variant_id"), "filter": record.get("filter")}
            for record in dataset_records if record.get("filter") not in (None, ".", "PASS")
        ]
        passing = [record for record in dataset_records if record.get("filter") in (None, ".", "PASS")]
        result["quality_filter_passed"] = bool(passing)
        non_founder = [record for record in passing if record.get("non_founder_observed") is True]
        if not passing:
            result["status"], result["found"] = "filtered_record", None
        elif non_founder:
            result["status"], result["found"] = "found", True
        else:
            result["status"], result["found"] = "absent_in_non_founder_populations", False
        frequencies = []
        for record in passing:
            value, metric, population = _record_frequency_value(record, scored_ancestries)
            if value is not None:
                frequencies.append((value, metric, population, record))
        if frequencies:
            value, metric, _, best_record = max(frequencies, key=lambda item: item[0])
            result["max_af"], result["frequency_metric"] = value, metric
        else:
            best_record = passing[0] if passing else dataset_records[0]
        callset = _record_to_callset_summary(best_record, scored_ancestries)
        result[config["callset"]] = callset
        result["non_founder_allele_count"] = callset["non_founder_allele_count"]
    else:
        result["status"], result["found"] = "absent", False
        result["dataset"] = config["dataset_names"][0]
    result["coverage"] = lookup_coverage_by_position(
        repository,
        coords,
        config["coverage_dataset_key"],
        config["assembly"],
        coverage_threshold,
        classification_compatible=config["coverage_classification_compatible"],
        compatibility_status=config["coverage_frequency_compatibility"],
        compatibility_reason=config["coverage_compatibility_reason"],
    )
    return result


def get_gnomad_frequencies(
    repository: GnomadRepository,
    gene: str | None,
    grch37: Any | None = None,
    grch38: Any | None = None,
) -> dict[str, Any]:
    result = {
        "status": "not_queried", "gene": gene, "policy_id": None,
        "classification_policy": None,
        "frequency_policy": {}, "found": None, "datasets": {},
        "coverage": {"status": "not_evaluated", "passes_pm2": False, "datasets": {}},
        "max_af": None, "frequency_metric": None, "pm2_absence_established": False,
        "pm2_coverage_ok": False, "errors": [],
        "pm2_coverage_method": dict(PM2_COVERAGE_METHOD_REVIEW),
        "source": str(repository.frequency_path) if repository.frequency_path else None,
        "cache_mode": repository.frequency_status,
    }
    policy = classification_policy_for_gene(gene)
    if policy is None:
        result["status"] = "policy_unavailable"
        result["errors"].append(
            f"no active gene-specific gnomAD classification policy for {gene or 'unspecified gene'}"
        )
        return result
    frequency_policy = policy.get("frequency_criteria") or {}
    pm2_policy = frequency_policy.get("pm2") or {}
    required = list(pm2_policy.get("required_absence_dataset_runtime_keys") or [])
    scored_ancestries = list(frequency_policy.get("scored_non_founder_ancestries") or [])
    threshold = as_float(pm2_policy.get("minimum_mean_depth"))
    if not required or not scored_ancestries or threshold is None or any(
        key not in GNOMAD_LOCAL_DATASET_CONFIG for key in required
    ):
        result["status"] = "policy_unavailable"
        result["errors"].append(f"active gnomAD policy {policy.get('policy_id')} is incomplete")
        return result
    result["policy_id"] = policy["policy_id"]
    result["classification_policy"] = policy
    result["frequency_policy"] = frequency_policy
    coords_by_kind = {"grch37": grch37, "grch38": grch38}
    for dataset_key in required:
        config = GNOMAD_LOCAL_DATASET_CONFIG[dataset_key]
        coords = coords_by_kind.get(config["coord"])
        dataset_result = query_gnomad_dataset(
            repository, variant_id_from_coords(coords), coords, config, threshold, scored_ancestries
        )
        result["datasets"][dataset_key] = dataset_result
        result["errors"].extend(dataset_result.get("errors") or [])
    result["coverage"] = aggregate_coverage(result["datasets"], required, threshold)
    statuses = [result["datasets"].get(key, {}).get("status") for key in required]
    any_found = any(status == "found" for status in statuses)
    all_absent = all(status in {"absent", "absent_in_non_founder_populations"} for status in statuses)
    any_cache_missing = any(
        status in {"cache_missing", "cache_untrusted", "dataset_not_available"}
        for status in statuses
    )
    any_no_coords = any(status == "no_coordinates" for status in statuses)
    any_outside = any(status == "outside_cached_region" for status in statuses)
    frequencies = [
        (as_float(result["datasets"][key].get("max_af")), result["datasets"][key].get("frequency_metric"))
        for key in required
        if as_float(result["datasets"][key].get("max_af")) is not None
        and result["datasets"][key].get("frequency_metric") == "faf95"
    ]
    if frequencies:
        result["max_af"], result["frequency_metric"] = max(frequencies, key=lambda item: item[0])
    result["found"] = any_found
    method_allows_automatic_pm2 = (
        PM2_COVERAGE_METHOD_REVIEW.get("automatic_assignment_allowed") is True
        and pm2_policy.get("coverage_scope")
        == PM2_COVERAGE_METHOD_REVIEW.get("scope")
    )
    result["pm2_coverage_ok"] = (
        result["coverage"]["passes_pm2"] and method_allows_automatic_pm2
    )
    v2_status = result["datasets"].get("v2_1_non_cancer", {}).get("status", "")
    v3_status = result["datasets"].get("v3_1_non_cancer", {}).get("status", "")
    v3_not_in_cache = v3_status in {"outside_cached_region", "dataset_not_available", "cache_missing"}
    if all_absent and result["pm2_coverage_ok"]:
        result["pm2_absence_established"] = True
        result["pm2_datasets_note"] = " + ".join(GNOMAD_LOCAL_DATASET_CONFIG[key]["label"] for key in required)
    else:
        result["pm2_absence_established"] = False
        result["pm2_datasets_note"] = ""
    if any_found:
        result["status"] = "found"
    elif all_absent:
        if not method_allows_automatic_pm2:
            result["status"] = "pm2_coverage_method_unresolved"
        elif not result["coverage"].get("all_sources_classification_compatible"):
            result["status"] = "coverage_release_incompatible"
        else:
            result["status"] = (
                "absent_with_coverage"
                if result["pm2_coverage_ok"]
                else "absent_without_sufficient_coverage"
            )
    elif v2_status == "absent" and v3_not_in_cache:
        result["status"] = "absent_v2_only"
    elif any_cache_missing:
        result["status"] = "cache_missing"
    elif any_no_coords:
        result["status"] = "no_coordinates"
    elif any_outside:
        result["status"] = "outside_cached_region"
    else:
        result["status"] = "partial"
    result["population_frequency_audit"] = population_frequency_audit(result)
    result["founder_context_only_observed"] = result["population_frequency_audit"]["founder_context_only_observed"]
    return result


def gnomad_status_summary(gnomad_data: Mapping[str, Any]) -> str:
    parts = [f"status={gnomad_data.get('status')}"]
    if gnomad_data.get("max_af") is not None:
        parts.append(f"max_{gnomad_data.get('frequency_metric') or 'af'}={gnomad_data['max_af']:.6g}")
    coverage = gnomad_data.get("coverage", {})
    method = gnomad_data.get("pm2_coverage_method") or {}
    parts.extend(
        (
            f"coverage={coverage.get('status')}",
            f"pm2_method={method.get('status') or 'missing'}",
            f"cache={gnomad_data.get('cache_mode')}",
        )
    )
    for key, dataset in gnomad_data.get("datasets", {}).items():
        mean_depth = ((coverage.get("datasets") or {}).get(key, {})).get("mean_depth")
        depth = f",depth={mean_depth:.1f}" if mean_depth is not None else ""
        compatibility = (dataset.get("coverage") or {}).get(
            "compatibility_status"
        )
        parts.append(
            f"{key}:{dataset.get('status')}[{dataset.get('dataset')}]"
            f"{depth},coverage_compatibility={compatibility or 'missing'}"
        )
    return "; ".join(parts)
