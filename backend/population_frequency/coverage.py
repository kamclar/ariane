"""Coverage lookup and quality checks for population-frequency evidence."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from backend.population_frequency.models import GnomadRepository
from backend.population_frequency.policy import GNOMAD_LOCAL_DATASET_CONFIG
from backend.population_frequency.utils import (
    add_chr,
    as_float,
    reference_span_from_coords,
    strip_chr,
)


def lookup_coverage_by_position(
    repository: GnomadRepository,
    coords: Any | None,
    dataset_key: str,
    build: str,
    threshold: float,
    *,
    classification_compatible: bool,
    compatibility_status: str,
    compatibility_reason: str,
) -> dict[str, Any]:
    """Measure depth while keeping source compatibility explicit."""
    compatibility = {
        "classification_compatible": classification_compatible,
        "compatibility_status": compatibility_status,
        "compatibility_reason": compatibility_reason,
    }
    chrom, span_start, span_end = reference_span_from_coords(coords)
    if chrom is None or span_start is None or span_end is None:
        return {
            "mean_depth": None,
            "threshold": threshold,
            "passes": False,
            "source": "no_coordinates",
            "position_key": None,
            "position_keys": [],
            "coverage_scope": "variant_reference_span",
            "span_start": None,
            "span_end": None,
            "positions_expected": 0,
            "positions_available": 0,
            **compatibility,
        }

    found: list[tuple[str, Mapping[str, Any]]] = []
    missing_positions: list[int] = []
    for pos in range(span_start, span_end + 1):
        keys = (
            f"{dataset_key}|{build}|{chrom}|{pos}",
            f"{dataset_key}|{build}|{strip_chr(chrom)}|{pos}",
            f"{dataset_key}|{build}|{add_chr(chrom)}|{pos}",
        )
        match = None
        match_key = None
        for key in dict.fromkeys(keys):
            if key in repository.coverage_by_position:
                match = repository.coverage_by_position[key]
                match_key = key
                break
        if match is None or as_float(match.get("mean_depth")) is None:
            missing_positions.append(pos)
        else:
            found.append((match_key, match))

    expected = span_end - span_start + 1
    if len(found) == expected:
        mean_depth = sum(as_float(cov.get("mean_depth")) for _, cov in found) / expected
        medians = [as_float(cov.get("median_depth")) for _, cov in found]
        over_20_values = [as_float(cov.get("over_20")) for _, cov in found]
        over_25_values = [as_float(cov.get("over_25")) for _, cov in found]
        median_depth = (
            sum(medians) / expected if all(value is not None for value in medians) else None
        )
        over_20 = (
            sum(over_20_values) / expected
            if all(value is not None for value in over_20_values)
            else None
        )
        over_25 = (
            sum(over_25_values) / expected
            if all(value is not None for value in over_25_values)
            else None
        )
        sources = sorted(
            {cov.get("source") or "gnomad_coverage_summary_tsv" for _, cov in found}
        )
        return {
            "mean_depth": mean_depth,
            "median_depth": median_depth,
            "over_20": over_20,
            "over_25": over_25,
            "threshold": threshold,
            "passes": mean_depth >= threshold and classification_compatible,
            "measurement_passes_threshold": mean_depth >= threshold,
            "source": ", ".join(sources),
            "position_key": found[0][1].get("position_key") or found[0][0],
            "position_keys": [cov.get("position_key") or key for key, cov in found],
            "coverage_scope": "variant_reference_span",
            "span_start": span_start,
            "span_end": span_end,
            "positions_expected": expected,
            "positions_available": len(found),
            "missing_positions": [],
            **compatibility,
        }
    return {
        "mean_depth": None,
        "threshold": threshold,
        "passes": False,
        "measurement_passes_threshold": False,
        "source": "not_found_in_coverage_snapshot",
        "position_key": None,
        "position_keys": [cov.get("position_key") or key for key, cov in found],
        "coverage_scope": "variant_reference_span",
        "span_start": span_start,
        "span_end": span_end,
        "positions_expected": expected,
        "positions_available": len(found),
        "missing_positions": missing_positions,
        **compatibility,
    }


def aggregate_coverage(
    dataset_results: Mapping[str, Mapping[str, Any]],
    required_dataset_keys: Sequence[str],
    threshold: float,
) -> dict[str, Any]:
    datasets: dict[str, Any] = {}
    for key, dataset in dataset_results.items():
        coverage = dataset.get("coverage") or {}
        mean_depth = as_float(coverage.get("mean_depth"))
        classification_compatible = (
            coverage.get("classification_compatible") is True
        )
        dataset_threshold = as_float(coverage.get("threshold"))
        measurement_passes = (
            bool(coverage.get("measurement_passes_threshold"))
            if coverage.get("measurement_passes_threshold") is not None
            else (
                mean_depth is not None
                and dataset_threshold is not None
                and mean_depth >= dataset_threshold
            )
        )
        passes = measurement_passes and classification_compatible
        datasets[key] = {
            "available": mean_depth is not None,
            "passes": passes,
            "measurement_passes_threshold": measurement_passes,
            "classification_compatible": classification_compatible,
            "compatibility_status": coverage.get("compatibility_status"),
            "compatibility_reason": coverage.get("compatibility_reason"),
            "mean_depth": mean_depth,
            "threshold": dataset_threshold,
            "source": coverage.get("source"),
            "position_key": coverage.get("position_key"),
            "position_keys": coverage.get("position_keys") or [],
            "coverage_scope": coverage.get("coverage_scope"),
            "span_start": coverage.get("span_start"),
            "span_end": coverage.get("span_end"),
            "positions_expected": coverage.get("positions_expected"),
            "positions_available": coverage.get("positions_available"),
            "missing_positions": coverage.get("missing_positions") or [],
            "callsets": {
                GNOMAD_LOCAL_DATASET_CONFIG[key]["callset"]: {
                    "mean_depth": mean_depth,
                    "passes": passes,
                    "measurement_passes_threshold": measurement_passes,
                    "classification_compatible": classification_compatible,
                    "source": coverage.get("source"),
                }
            },
        }
    required = list(required_dataset_keys)
    all_available = all(datasets.get(key, {}).get("available") for key in required)
    all_compatible = all(
        datasets.get(key, {}).get("classification_compatible") is True
        for key in required
    )
    all_pass = all(datasets.get(key, {}).get("passes") for key in required)
    return {
        "status": (
            "ok"
            if all_pass
            else (
                "missing"
                if not all_available
                else ("incompatible" if not all_compatible else "insufficient")
            )
        ),
        "passes_pm2": all_pass,
        "all_sources_classification_compatible": all_compatible,
        "min_required_mean_depth": threshold,
        "datasets": datasets,
    }


def frequency_depth_ok(
    gnomad_data: Mapping[str, Any], minimum_mean_depth: float
) -> bool:
    """Check policy depth in a dataset contributing the maximum FAF95."""
    max_af = as_float(gnomad_data.get("max_af"))
    for dataset in gnomad_data.get("datasets", {}).values():
        if dataset.get("status") != "found":
            continue
        if max_af is not None and as_float(dataset.get("max_af")) != max_af:
            continue
        mean_depth = as_float((dataset.get("coverage") or {}).get("mean_depth"))
        coverage = dataset.get("coverage") or {}
        if (
            coverage.get("classification_compatible") is True
            and mean_depth is not None
            and mean_depth >= minimum_mean_depth
        ):
            return True
    return False


def frequency_qc_ok(gnomad_data: Mapping[str, Any]) -> bool:
    """Require a PASS record in a dataset contributing the scored FAF95."""
    max_af = as_float(gnomad_data.get("max_af"))
    for dataset in gnomad_data.get("datasets", {}).values():
        if dataset.get("status") != "found":
            continue
        if max_af is not None and as_float(dataset.get("max_af")) != max_af:
            continue
        if dataset.get("quality_filter_passed") is True:
            return True
        for callset_name in ("exomes", "genomes"):
            callset = dataset.get(callset_name) or {}
            if callset.get("available") and not callset.get("filters"):
                return True
    return False
