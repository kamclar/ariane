"""Validated local ClinGen ERepo assertions for the ENIGMA BRCA1/2 VCEP."""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import hashlib
import json
import re
from pathlib import Path
from typing import Any

from backend.policy.gene import active_genes, reference_transcript


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
REGISTRY_PATH = DATA_DIR / "enigma_erepo_vcep_registry.json"
METADATA_PATH = DATA_DIR / "enigma_erepo_vcep_registry.metadata.json"
PANEL_NAME = "ENIGMA BRCA1 and BRCA2 VCEP"
GUIDELINE_AFFILIATION_ID = "50087"
CURRENT_METHOD_VERSION = "1.2.0"
EXPORT_URL = (
    "https://erepo.clinicalgenome.org/evrepo/api/summary/"
    "classifications/download?type=csv"
)
SERVICE_URL = "https://erepo.clinicalgenome.org/evrepo/api/summary/srvc"
CLASSIFICATIONS = {
    "Pathogenic",
    "Likely Pathogenic",
    "Uncertain Significance",
    "Likely Benign",
    "Benign",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_erepo_vcep_registry(
    registry: dict[str, Any],
    metadata: dict[str, Any],
    *,
    registry_path: Path = REGISTRY_PATH,
) -> None:
    if registry.get("schema_version") != 1 or registry.get("status") != "active":
        raise RuntimeError("ENIGMA ERepo VCEP registry has unsupported metadata")
    if metadata.get("schema_version") != 1 or metadata.get("status") != "active":
        raise RuntimeError("ENIGMA ERepo VCEP registry metadata is unsupported")
    if metadata.get("registry_file") != registry_path.name:
        raise RuntimeError("ENIGMA ERepo VCEP metadata names a different registry file")
    expected_sha = str(metadata.get("registry_sha256") or "")
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha) or _sha256(registry_path) != expected_sha:
        raise RuntimeError("ENIGMA ERepo VCEP registry checksum does not match metadata")
    for payload in (registry, metadata):
        if payload.get("expert_panel") != PANEL_NAME:
            raise RuntimeError("ENIGMA ERepo VCEP registry has an unexpected expert panel")
        if payload.get("guideline_affiliation_id") != GUIDELINE_AFFILIATION_ID:
            raise RuntimeError("ENIGMA ERepo VCEP registry has an unexpected affiliation")
    if registry.get("current_assertion_method_version") != CURRENT_METHOD_VERSION:
        raise RuntimeError("ENIGMA ERepo VCEP registry uses an unsupported current version")
    if metadata.get("source_export_url") != EXPORT_URL or metadata.get("source_service_url") != SERVICE_URL:
        raise RuntimeError("ENIGMA ERepo VCEP metadata names an unexpected source")
    if not re.fullmatch(r"[0-9a-f]{64}", str(metadata.get("source_export_sha256") or "")):
        raise RuntimeError("ENIGMA ERepo VCEP metadata lacks the source export checksum")
    if not re.fullmatch(r"\d+\.\d+\.\d+", str(metadata.get("source_api_version") or "")):
        raise RuntimeError("ENIGMA ERepo VCEP metadata lacks the source API version")
    if not str(metadata.get("retrieved_at") or "").strip():
        raise RuntimeError("ENIGMA ERepo VCEP metadata lacks the retrieval time")
    policy = registry.get("runtime_policy")
    if not isinstance(policy, dict) or policy.get("network_fallback") is not False:
        raise RuntimeError("ENIGMA ERepo VCEP registry must disable runtime network fallback")

    records = registry.get("records")
    if not isinstance(records, list) or not records:
        raise RuntimeError("ENIGMA ERepo VCEP registry contains no records")
    if registry.get("record_count") != len(records) or metadata.get("record_count") != len(records):
        raise RuntimeError("ENIGMA ERepo VCEP registry record count is inconsistent")

    genes = set(active_genes())
    uuids: set[str] = set()
    classification_counts: Counter[str] = Counter()
    version_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    for index, record in enumerate(records):
        prefix = f"ENIGMA ERepo VCEP record #{index}"
        uuid = str(record.get("uuid") or "")
        gene = str(record.get("gene") or "")
        c_notation = str(record.get("c_notation") or "")
        method_version = str(record.get("assertion_method_version") or "")
        source_status = str(record.get("source_status") or "")
        if not re.fullmatch(r"[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}", uuid):
            raise RuntimeError(f"{prefix} has an invalid UUID")
        if uuid in uuids:
            raise RuntimeError(f"{prefix} duplicates UUID {uuid}")
        uuids.add(uuid)
        if gene not in genes or record.get("transcript") != reference_transcript(gene):
            raise RuntimeError(f"{prefix} has an invalid gene or transcript")
        if not c_notation.startswith("c."):
            raise RuntimeError(f"{prefix} has invalid coding HGVS")
        if record.get("classification") not in CLASSIFICATIONS:
            raise RuntimeError(f"{prefix} has an unsupported classification")
        if record.get("expert_panel") != PANEL_NAME:
            raise RuntimeError(f"{prefix} has an unexpected expert panel")
        if GUIDELINE_AFFILIATION_ID not in str(record.get("guideline_url") or ""):
            raise RuntimeError(f"{prefix} has an unexpected guideline")
        if record.get("retracted") is not False:
            raise RuntimeError(f"{prefix} is retracted")
        expected_status = (
            "current_vcep_assertion"
            if method_version == CURRENT_METHOD_VERSION
            else "unversioned_vcep_assertion"
            if not method_version
            else "historical_vcep_assertion"
        )
        if source_status != expected_status:
            raise RuntimeError(f"{prefix} has an inconsistent source status")
        if method_version and not str(record.get("assertion_method_id") or "").startswith("https://cspec."):
            raise RuntimeError(f"{prefix} lacks the assertion method identifier")
        if method_version and gene not in str(record.get("assertion_method_label") or ""):
            raise RuntimeError(f"{prefix} has an inconsistent assertion method label")
        classification_counts[record["classification"]] += 1
        version_counts[method_version] += 1
        status_counts[source_status] += 1

    expected_counts = (
        ("classification_counts", classification_counts),
        ("assertion_method_version_counts", version_counts),
        ("source_status_counts", status_counts),
    )
    for field, counts in expected_counts:
        value = dict(sorted(counts.items()))
        if registry.get(field) != value or metadata.get(field) != value:
            raise RuntimeError(f"ENIGMA ERepo VCEP registry has invalid {field}")


@lru_cache(maxsize=1)
def load_erepo_vcep_registry() -> dict[str, Any]:
    try:
        registry = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
        metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"ENIGMA ERepo VCEP registry cannot be loaded: {exc}") from exc
    validate_erepo_vcep_registry(registry, metadata)
    return registry


def lookup_erepo_vcep_assertion(gene: str, c_notation: str) -> dict[str, Any]:
    """Return the newest current assertion, or a non-current assertion for warning."""
    matches = [
        record
        for record in load_erepo_vcep_registry()["records"]
        if record["gene"] == gene and record["c_notation"] == c_notation
    ]
    if not matches:
        return {"status": "not_found"}
    matches.sort(key=lambda item: (item["published_date"], item["uuid"]), reverse=True)
    current = [item for item in matches if item["source_status"] == "current_vcep_assertion"]
    selected = current[0] if current else matches[0]
    return {"status": selected["source_status"], "record": dict(selected)}
