"""Versioned gnomAD dataset bindings and gene-specific frequency policies."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


PROJECT_ROOT = Path(__file__).resolve().parent.parent
GNOMAD_DIR = PROJECT_ROOT / "data" / "gnomad"
GNOMAD_FREQUENCY_SNAPSHOT_PATH = GNOMAD_DIR / "gnomad_brca_frequency_snapshot.json"
GNOMAD_COVERAGE_SNAPSHOT_PATH = GNOMAD_DIR / "gnomad_brca_coverage_snapshot.json"
GNOMAD_PANEL_MANIFEST_JSON = GNOMAD_DIR / "gnomad_panel_manifest.json"

GNOMAD_LOCAL_DATASET_CONFIG = {
    "v2_1_non_cancer": {
        "label": "gnomAD v2.1.1 exomes GRCh37",
        "assembly": "GRCh37",
        "coord": "grch37",
        "dataset_names": ["gnomad_v2_1_1_exomes_grch37"],
        "coverage_dataset_key": "gnomad_v2_1_1_exomes_grch37",
        "callset": "exomes",
    },
    "v3_1_non_cancer": {
        "label": "gnomAD v3.1.2 genomes GRCh38",
        "assembly": "GRCh38",
        "coord": "grch38",
        "dataset_names": ["gnomad_v3_1_2_genomes_grch38"],
        "coverage_dataset_key": "gnomad_v3_1_2_genomes_grch38",
        "callset": "genomes",
    },
}


def canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def approved_manifest(
    path: Path = GNOMAD_PANEL_MANIFEST_JSON,
) -> dict[str, Any] | None:
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return manifest if isinstance(manifest, dict) else None


def approved_manifest_sha256(
    path: Path = GNOMAD_PANEL_MANIFEST_JSON,
) -> str | None:
    manifest = approved_manifest(path)
    return canonical_sha256(manifest) if manifest is not None else None


def classification_policy_for_gene(
    gene: str | None,
    manifest: Mapping[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Return the explicitly assigned active policy without inheritance."""
    if not gene:
        return None
    resolved = dict(manifest) if manifest is not None else approved_manifest()
    if not resolved or resolved.get("schema_version") != 2:
        return None
    target = next(
        (
            item
            for item in resolved.get("targets", [])
            if str(item.get("gene", "")).upper() == str(gene).upper()
        ),
        None,
    )
    if not target or target.get("activation_status") != "active":
        return None
    policy_id = target.get("classification_policy_id")
    policy = (resolved.get("classification_policies") or {}).get(policy_id)
    if not policy or policy.get("status") != "active":
        return None
    if str(gene).upper() not in {
        str(item).upper() for item in policy.get("applicable_genes", [])
    }:
        return None
    return {"policy_id": policy_id, **policy}


def manifest_scored_ancestries(manifest: Mapping[str, Any]) -> set[str]:
    return {
        ancestry
        for policy in (manifest.get("classification_policies") or {}).values()
        if policy.get("status") == "active"
        for ancestry in policy.get("frequency_criteria", {}).get(
            "scored_non_founder_ancestries", []
        )
    }


def runtime_dataset_binding_error(
    manifest: Mapping[str, Any],
) -> str | None:
    manifest_by_runtime = {
        item.get("runtime_key"): item for item in manifest.get("datasets", [])
    }
    required = {
        runtime_key
        for policy in (manifest.get("classification_policies") or {}).values()
        if policy.get("status") == "active"
        for runtime_key in policy.get("frequency_criteria", {}).get(
            "required_dataset_runtime_keys", []
        )
    }
    for runtime_key in required:
        item = manifest_by_runtime.get(runtime_key)
        config = GNOMAD_LOCAL_DATASET_CONFIG.get(runtime_key)
        if item is None or config is None:
            return f"runtime dataset {runtime_key!r} is not implemented"
        if (
            item.get("dataset_key") not in config["dataset_names"]
            or item.get("assembly") != config["assembly"]
            or item.get("callset") != config["callset"]
            or item.get("subset") != "non_cancer"
            or item.get("dataset_key") != config["coverage_dataset_key"]
        ):
            return f"runtime dataset {runtime_key!r} differs from the manifest"
    return None
