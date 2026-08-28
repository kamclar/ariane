"""Checksum-validated loading of gnomAD frequency and coverage snapshots."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

from backend.data_health import clear_issue, register_issue
from backend.population_frequency.models import GnomadRepository
from backend.population_frequency.policy import (
    GNOMAD_COVERAGE_SNAPSHOT_PATH,
    GNOMAD_FREQUENCY_SNAPSHOT_PATH,
    GNOMAD_LOCAL_DATASET_CONFIG,
    GNOMAD_PANEL_MANIFEST_JSON,
    approved_manifest,
    approved_manifest_sha256,
    canonical_sha256,
    manifest_scored_ancestries,
    runtime_dataset_binding_error,
)
from backend.population_frequency.utils import add_chr, as_float, as_int, strip_chr


LOGGER = logging.getLogger("ariane.population_frequency.repository")


def validate_frequency_snapshot(payload: Mapping[str, Any]) -> str | None:
    """Validate that every scored record has auditable non-cancer FAF95."""
    mapping = payload.get("variants") or payload.get("by_variant") or {}
    if not isinstance(mapping, dict) or not mapping:
        return "variant mapping is missing or empty"
    metadata = payload.get("metadata") or {}
    if metadata.get("schema_version") != 2:
        return "snapshot schema_version must be 2"
    manifest = approved_manifest()
    if manifest is None:
        return f"approved panel manifest is missing or invalid: {GNOMAD_PANEL_MANIFEST_JSON}"
    binding_error = runtime_dataset_binding_error(manifest)
    if binding_error:
        return binding_error
    if metadata.get("manifest_sha256") != canonical_sha256(manifest):
        return "snapshot was built from a different panel/source manifest"
    if metadata.get("automatic_release_activation") is not False:
        return "automatic release activation must be disabled"
    if metadata.get("records_sha256") != canonical_sha256(mapping):
        return "variant records checksum mismatch"
    if metadata.get("classification_policies") != manifest.get(
        "classification_policies"
    ):
        return "snapshot classification policies differ from the approved manifest"

    expected_datasets = {
        name
        for config in GNOMAD_LOCAL_DATASET_CONFIG.values()
        for name in config["dataset_names"]
    }
    counts = {dataset: 0 for dataset in expected_datasets}
    missing_faf95: list[str] = []
    expected_methods = {
        "gnomad_v2_1_1_exomes_grch37": "official_gnomad_hail_table_non_cancer_faf95",
        "gnomad_v3_1_2_genomes_grch38": (
            "hail.experimental.filtering_allele_frequency_from_official_non_cancer_ac_an"
        ),
    }
    manifest_datasets = {
        item["dataset_key"]: item for item in manifest.get("datasets", [])
    }
    scored_ancestries = manifest_scored_ancestries(manifest)
    if not scored_ancestries:
        return "approved manifest contains no active scored population groups"
    for variant_id, records in mapping.items():
        if not isinstance(records, list):
            return f"record list is invalid for {variant_id}"
        for record in records:
            dataset = record.get("dataset")
            if dataset not in expected_datasets:
                continue
            counts[dataset] += 1
            if (
                as_float(record.get("faf95_max")) is None
                or record.get("faf95_scope") != "non_cancer_non_founder_ancestries"
                or record.get("faf95_method") != expected_methods.get(dataset)
            ):
                missing_faf95.append(variant_id)
                if len(missing_faf95) >= 3:
                    break
                continue
            if set(record.get("non_founder_ac_by_ancestry") or {}) != scored_ancestries:
                return f"non-founder population counts are incomplete for {variant_id}"
            observed = any(
                (as_int(value) or 0) > 0
                for value in record["non_founder_ac_by_ancestry"].values()
            )
            if record.get("non_founder_observed") is not observed:
                return f"non-founder presence flag is inconsistent for {variant_id}"
            expected_context = {
                item["code"]
                for item in manifest_datasets.get(dataset, {}).get(
                    "excluded_population_context", []
                )
            }
            context = record.get("excluded_population_context") or {}
            if set(context) != expected_context or any(
                item.get("used_for_ba1_bs1") is not False
                or item.get("used_for_pm2_presence") is not False
                for item in context.values()
            ):
                return f"excluded population context is invalid for {variant_id}"
        if len(missing_faf95) >= 3:
            break

    absent = sorted(dataset for dataset, count in counts.items() if count == 0)
    if absent:
        return f"required dataset records are missing: {', '.join(absent)}"
    if missing_faf95:
        return (
            "records lack ENIGMA-compatible non-cancer FAF95 provenance: "
            + ", ".join(missing_faf95)
        )
    for release in ("v2_faf95", "v3_faf95"):
        if (metadata.get(release) or {}).get("raw_af_fallback_allowed") is not False:
            return f"metadata {release}.raw_af_fallback_allowed must be false"
    logged = {
        item.get("dataset")
        for item in metadata.get("extraction_log", [])
        if item.get("status") == "ok"
        and (item.get("source_identity") or {}).get("x_goog_hash")
        and (item.get("source_identity") or {}).get("etag")
    }
    if not expected_datasets.issubset(logged):
        return "official source identity or successful extraction log is incomplete"
    return None


def validate_coverage_snapshot(payload: Mapping[str, Any]) -> str | None:
    mapping = payload.get("coverage_by_position") or {}
    if not isinstance(mapping, dict) or not mapping:
        return "coverage mapping is missing or empty"
    metadata = payload.get("metadata") or {}
    if metadata.get("schema_version") != 2:
        return "coverage schema_version must be 2"
    manifest_hash = approved_manifest_sha256()
    if manifest_hash is None:
        return f"approved panel manifest is missing or invalid: {GNOMAD_PANEL_MANIFEST_JSON}"
    if metadata.get("manifest_sha256") != manifest_hash:
        return "coverage was built from a different panel/source manifest"
    manifest = approved_manifest()
    if manifest is not None:
        binding_error = runtime_dataset_binding_error(manifest)
        if binding_error:
            return binding_error
    if manifest is None or metadata.get("classification_policies") != manifest.get(
        "classification_policies"
    ):
        return "coverage classification policies differ from the approved manifest"
    if metadata.get("records") != len(mapping):
        return "coverage record count does not match metadata"
    if metadata.get("records_sha256") != canonical_sha256(mapping):
        return "coverage records checksum mismatch"
    expected_datasets = {
        name
        for config in GNOMAD_LOCAL_DATASET_CONFIG.values()
        for name in config["dataset_names"]
    }
    observed_datasets = {item.get("dataset_key") for item in mapping.values()}
    if observed_datasets != expected_datasets:
        return "coverage does not contain exactly the required gnomAD datasets"
    documented = {
        item.get("dataset")
        for item in metadata.get("datasets", [])
        if (item.get("coverage_source_identity") or {}).get("etag")
        and (item.get("coverage_source_identity") or {}).get("x_goog_hash")
        and item.get("coverage_hail_uri")
    }
    if documented != expected_datasets:
        return "official coverage source identity is incomplete"
    return None


def choose_frequency_snapshot(
    path: Path = GNOMAD_FREQUENCY_SNAPSHOT_PATH,
) -> Path | None:
    return path if path.exists() else None


def _normalize_variant_keys(
    raw_mapping: Mapping[str, Any],
) -> dict[str, Any]:
    normalized: dict[str, Any] = {}
    for key, records in raw_mapping.items():
        normalized[key] = tuple(records)
        parts = str(key).split("-")
        if len(parts) >= 4:
            chrom = parts[0]
            rest = "-".join(parts[1:])
            normalized[f"{strip_chr(chrom)}-{rest}"] = tuple(records)
            normalized[f"{add_chr(chrom)}-{rest}"] = tuple(records)
    return normalized


def _load_json(path: Path) -> tuple[dict[str, Any] | None, str | None]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return None, f"{type(exc).__name__}: {exc}"
    return (value, None) if isinstance(value, dict) else (None, "root is not an object")


def load_gnomad_repository(
    frequency_path: Path = GNOMAD_FREQUENCY_SNAPSHOT_PATH,
    coverage_path: Path = GNOMAD_COVERAGE_SNAPSHOT_PATH,
) -> GnomadRepository:
    """Load both snapshots explicitly and return one immutable runtime binding."""
    selected_frequency = choose_frequency_snapshot(frequency_path)
    variants: Mapping[str, Any] = {}
    metadata: Mapping[str, Any] = {}
    frequency_status = "missing"
    if selected_frequency is None:
        register_issue(
            "gnomAD variant cache",
            f"approved frequency snapshot is missing: {frequency_path}",
        )
    else:
        payload, error = _load_json(selected_frequency)
        if error:
            frequency_status = "load_failed"
            register_issue(
                "gnomAD variant cache",
                f"could not load {selected_frequency}: {error}",
            )
        else:
            assert payload is not None
            metadata = payload.get("metadata", {})
            validation_error = validate_frequency_snapshot(payload)
            if validation_error:
                frequency_status = "invalid_faf95"
                register_issue(
                    "gnomAD variant cache",
                    "cache failed non-cancer FAF95 validation: " + validation_error,
                )
            else:
                raw = payload.get("variants") or payload.get("by_variant") or {}
                variants = _normalize_variant_keys(raw)
                frequency_status = "approved_snapshot"
                clear_issue("gnomAD variant cache")
                LOGGER.info(
                    "Loaded gnomAD frequency snapshot %s with %d variants",
                    selected_frequency,
                    len(raw),
                )

    coverage: Mapping[str, Any] = {}
    coverage_status = "missing"
    if not coverage_path.exists():
        register_issue(
            "gnomAD coverage snapshot",
            f"coverage snapshot is missing: {coverage_path}",
        )
    else:
        payload, error = _load_json(coverage_path)
        if error:
            coverage_status = "load_failed"
            register_issue(
                "gnomAD coverage snapshot",
                f"could not load {coverage_path}: {error}",
            )
        else:
            assert payload is not None
            validation_error = validate_coverage_snapshot(payload)
            if validation_error:
                coverage_status = "invalid"
                register_issue(
                    "gnomAD coverage snapshot",
                    "coverage snapshot validation failed: " + validation_error,
                )
            else:
                coverage = payload.get("coverage_by_position", {}) or {}
                coverage_status = "approved_snapshot"
                clear_issue("gnomAD coverage snapshot")
                LOGGER.info(
                    "Loaded gnomAD coverage snapshot %s with %d positions",
                    coverage_path,
                    len(coverage),
                )

    return GnomadRepository(
        variants=MappingProxyType(dict(variants)),
        metadata=MappingProxyType(dict(metadata)),
        coverage_by_position=MappingProxyType(dict(coverage)),
        frequency_path=selected_frequency,
        frequency_status=frequency_status,
        coverage_path=coverage_path if coverage_path.exists() else None,
        coverage_status=coverage_status,
    )
