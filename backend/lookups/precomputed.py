from pathlib import Path
from typing import Any, Optional
import hashlib
import json
from backend.data_health import clear_issue, register_issue


PRECOMPUTED_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "precomputed"
CLASSIFICATION_SNAPSHOT_INDEX = PRECOMPUTED_DIR / "brca_module1_snv_classification_snapshot.index.json"
CLASSIFICATION_SNAPSHOT_METADATA = PRECOMPUTED_DIR / "brca_module1_snv_classification_snapshot.metadata.json"

_CLASSIFICATION_INDEX: Optional[dict[str, dict[str, Any]]] = None
_CLASSIFICATION_METADATA: Optional[dict[str, Any]] = None


def validate_classification_snapshot() -> None:
    """Fail startup if the coding-SNV normalization snapshot is inconsistent."""
    if not CLASSIFICATION_SNAPSHOT_INDEX.is_file():
        raise RuntimeError(
            f"Required classification snapshot index is missing: {CLASSIFICATION_SNAPSHOT_INDEX}"
        )
    if not CLASSIFICATION_SNAPSHOT_METADATA.is_file():
        raise RuntimeError(
            "Required classification snapshot metadata is missing: "
            f"{CLASSIFICATION_SNAPSHOT_METADATA}"
        )
    try:
        index_bytes = CLASSIFICATION_SNAPSHOT_INDEX.read_bytes()
        index = json.loads(index_bytes)
        metadata = json.loads(CLASSIFICATION_SNAPSHOT_METADATA.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Required classification snapshot cannot be loaded: {exc}") from exc
    if not isinstance(index, dict) or len(index) != 47547:
        raise RuntimeError(
            f"Required classification snapshot must contain 47,547 records, found {len(index)}"
        )
    actual_sha = hashlib.sha256(index_bytes).hexdigest()
    if metadata.get("index_sha256") != actual_sha:
        raise RuntimeError(
            "Required classification snapshot checksum mismatch: "
            f"expected {metadata.get('index_sha256')}, found {actual_sha}"
        )
    normalization = metadata.get("protein_consequence_normalization", {})
    if (
        normalization.get("status") != "complete"
        or normalization.get("records_checked") != len(index)
    ):
        raise RuntimeError(
            "Required classification snapshot lacks a complete HGVS protein-consequence audit"
        )


def _normal_key(gene: str, c_notation: str) -> str:
    return f"{gene.strip().upper()}:{c_notation.strip()}"


def load_classification_snapshot_index() -> dict[str, dict[str, Any]]:
    global _CLASSIFICATION_INDEX
    if _CLASSIFICATION_INDEX is not None:
        return _CLASSIFICATION_INDEX
    if not CLASSIFICATION_SNAPSHOT_INDEX.exists():
        _CLASSIFICATION_INDEX = {}
        register_issue("Classification snapshot", f"index is missing: {CLASSIFICATION_SNAPSHOT_INDEX}")
        return _CLASSIFICATION_INDEX
    try:
        with CLASSIFICATION_SNAPSHOT_INDEX.open(encoding="utf-8") as handle:
            _CLASSIFICATION_INDEX = json.load(handle)
        clear_issue("Classification snapshot")
    except (OSError, json.JSONDecodeError) as exc:
        _CLASSIFICATION_INDEX = {}
        register_issue(
            "Classification snapshot",
            f"could not load {CLASSIFICATION_SNAPSHOT_INDEX}: {type(exc).__name__}: {exc}",
        )
    return _CLASSIFICATION_INDEX


def load_classification_snapshot_metadata() -> dict[str, Any]:
    global _CLASSIFICATION_METADATA
    if _CLASSIFICATION_METADATA is not None:
        return _CLASSIFICATION_METADATA
    if not CLASSIFICATION_SNAPSHOT_METADATA.exists():
        _CLASSIFICATION_METADATA = {}
        register_issue("Classification snapshot metadata", f"metadata is missing: {CLASSIFICATION_SNAPSHOT_METADATA}")
        return _CLASSIFICATION_METADATA
    try:
        with CLASSIFICATION_SNAPSHOT_METADATA.open(encoding="utf-8") as handle:
            _CLASSIFICATION_METADATA = json.load(handle)
        clear_issue("Classification snapshot metadata")
    except (OSError, json.JSONDecodeError) as exc:
        _CLASSIFICATION_METADATA = {}
        register_issue(
            "Classification snapshot metadata",
            f"could not load {CLASSIFICATION_SNAPSHOT_METADATA}: {type(exc).__name__}: {exc}",
        )
    return _CLASSIFICATION_METADATA


def lookup_classification_snapshot(gene: str, c_notation: str) -> Optional[dict[str, Any]]:
    entry = load_classification_snapshot_index().get(_normal_key(gene, c_notation))
    if entry is None:
        return None
    metadata = load_classification_snapshot_metadata()
    return {
        "key": _normal_key(gene, c_notation),
        "snapshot_status": metadata.get("status", "snapshot_not_authoritative"),
        "snapshot_created": metadata.get("created", ""),
        "snapshot_version": metadata.get("application_version", ""),
        "snapshot_index_sha256": metadata.get("index_sha256", ""),
        "record": entry,
    }
