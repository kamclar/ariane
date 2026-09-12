"""Validated ClinGen ERepo PVS1 RNA records used by the classification DAG.

The live Evidence Repository remains an update and comparison source. Runtime
classification reads only this checksum-bound registry so a remote update
cannot change a cached ARIANE result without a reviewed data release.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from backend.gene_policy import reference_transcript, vcep_specification
from backend.modules.table4 import parse_pvs1_code_strength


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
REGISTRY_PATH = DATA_DIR / "enigma_erepo_pvs1_rna_registry.json"
METADATA_PATH = DATA_DIR / "enigma_erepo_pvs1_rna_registry.metadata.json"


class ErepoPvs1RnaRegistryError(RuntimeError):
    """The pinned curated-evidence registry is missing or inconsistent."""


def _read_object(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise ErepoPvs1RnaRegistryError(f"Required {label} is missing: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ErepoPvs1RnaRegistryError(
            f"Required {label} cannot be loaded: {path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise ErepoPvs1RnaRegistryError(f"Required {label} must be a JSON object")
    return value


def _validate_record(record: Mapping[str, Any], index: int) -> tuple[str, str]:
    required = {
        "gene", "reference_transcript", "c_notation", "code", "source_code",
        "strength", "points", "evidence_mechanism", "guideline_id",
        "guideline_version", "assertion_uuid", "assertion_version",
        "assertion_url", "published_date", "classification", "evidence_summary",
    }
    missing = sorted(required - set(record))
    if missing:
        raise ErepoPvs1RnaRegistryError(
            f"ERepo PVS1 RNA record {index} is incomplete: missing {missing}"
        )
    gene = str(record["gene"])
    c_notation = str(record["c_notation"])
    specification = vcep_specification(gene)
    if str(record["reference_transcript"]) != reference_transcript(gene):
        raise ErepoPvs1RnaRegistryError(
            f"ERepo record {gene} {c_notation} does not use the configured transcript"
        )
    if str(record["guideline_id"]) != specification["id"]:
        raise ErepoPvs1RnaRegistryError(
            f"ERepo record {gene} {c_notation} has the wrong guideline id"
        )
    if str(record["guideline_version"]) != specification["version"]:
        raise ErepoPvs1RnaRegistryError(
            f"ERepo record {gene} {c_notation} has the wrong guideline version"
        )
    if record["code"] != "PVS1_RNA" or record["evidence_mechanism"] != "rna_splicing":
        raise ErepoPvs1RnaRegistryError(
            f"ERepo record {gene} {c_notation} is not explicit PVS1 RNA evidence"
        )
    parsed_strength, parsed_points, requires_rna = parse_pvs1_code_strength(
        str(record["source_code"])
    )
    if not requires_rna or parsed_strength != record["strength"]:
        raise ErepoPvs1RnaRegistryError(
            f"ERepo record {gene} {c_notation} has inconsistent RNA strength"
        )
    if parsed_points != record["points"]:
        raise ErepoPvs1RnaRegistryError(
            f"ERepo record {gene} {c_notation} has inconsistent points"
        )
    if not str(record["assertion_uuid"]).strip() or not str(record["published_date"]).strip():
        raise ErepoPvs1RnaRegistryError(
            f"ERepo record {gene} {c_notation} lacks a published assertion identity"
        )
    return gene, c_notation


def _load_registry() -> tuple[dict[str, Any], dict[tuple[str, str], dict[str, Any]]]:
    registry_bytes = REGISTRY_PATH.read_bytes() if REGISTRY_PATH.is_file() else b""
    registry = _read_object(REGISTRY_PATH, "ERepo PVS1 RNA registry")
    metadata = _read_object(METADATA_PATH, "ERepo PVS1 RNA registry metadata")
    digest = hashlib.sha256(registry_bytes).hexdigest()
    if metadata.get("registry_sha256") != digest:
        raise ErepoPvs1RnaRegistryError(
            "ERepo PVS1 RNA registry checksum does not match metadata"
        )
    if registry.get("schema_version") != 1 or registry.get("status") != "active":
        raise ErepoPvs1RnaRegistryError(
            "ERepo PVS1 RNA registry must be active schema version 1"
        )
    for field in ("registry_id", "registry_version", "coverage_status"):
        if metadata.get(field) != registry.get(field):
            raise ErepoPvs1RnaRegistryError(
                f"ERepo PVS1 RNA registry metadata disagrees on {field}"
            )
    if metadata.get("validation_status") != "approved":
        raise ErepoPvs1RnaRegistryError(
            "ERepo PVS1 RNA registry is not approved for runtime use"
        )
    records = registry.get("records")
    if not isinstance(records, list):
        raise ErepoPvs1RnaRegistryError("ERepo PVS1 RNA records must be a list")
    if metadata.get("record_count") != len(records):
        raise ErepoPvs1RnaRegistryError(
            "ERepo PVS1 RNA record count does not match metadata"
        )
    by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for index, record in enumerate(records, start=1):
        if not isinstance(record, dict):
            raise ErepoPvs1RnaRegistryError(
                f"ERepo PVS1 RNA record {index} must be an object"
            )
        key = _validate_record(record, index)
        if key in by_key:
            raise ErepoPvs1RnaRegistryError(
                f"Duplicate ERepo PVS1 RNA record for {key[0]} {key[1]}"
            )
        by_key[key] = record
    return registry, by_key


REGISTRY, RECORDS_BY_VARIANT = _load_registry()
REGISTRY_SHA256 = hashlib.sha256(REGISTRY_PATH.read_bytes()).hexdigest()


def lookup_erepo_pvs1_rna(gene: str, c_notation: str) -> dict[str, Any]:
    """Return one exact approved record without inferring from nearby variants."""
    record = RECORDS_BY_VARIANT.get((gene, c_notation))
    provenance = {
        "registry_id": REGISTRY["registry_id"],
        "registry_version": REGISTRY["registry_version"],
        "registry_sha256": REGISTRY_SHA256,
        "coverage_status": REGISTRY["coverage_status"],
        "source_url": REGISTRY["source_url"],
    }
    if record is None:
        return {
            "status": "not_in_registry",
            "record": None,
            "reason": REGISTRY["absence_semantics"],
            **provenance,
        }
    return {
        "status": "eligible",
        "record": deepcopy(record),
        "reason": (
            "An exact, published ENIGMA BRCA1/2 VCEP v1.2 PVS1 RNA assertion "
            "is present in the approved local ERepo registry."
        ),
        **provenance,
    }
