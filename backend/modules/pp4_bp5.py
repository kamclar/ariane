"""Automatic PP4/BP5 lookup from the validated local clinical-LR snapshot."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Optional


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_pp4_clinical_lr_snapshot.index.json"
METADATA_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_pp4_clinical_lr_snapshot.metadata.json"
SOURCE_MANIFEST_PATH = REPOSITORY_ROOT / "data" / "sources" / "enigma" / "clinical_lr_sources.manifest.json"
INDEL_SNAPSHOT_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_normalized_indel_snapshot.index.json"
INDEL_METADATA_PATH = REPOSITORY_ROOT / "data" / "precomputed" / "brca_normalized_indel_snapshot.metadata.json"

PP4_POINTS = {"Very Strong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1}
BP5_POINTS = {"Very Strong": -8, "Strong": -4, "Moderate": -2, "Supporting": -1}

_SNAPSHOT: dict[str, dict[str, Any]] | None = None
_ALIASES: dict[str, str] | None = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def lr_to_pp4_strength(lr: float) -> Optional[str]:
    if lr >= 350:
        return "Very Strong"
    if lr >= 18.7:
        return "Strong"
    if lr >= 4.3:
        return "Moderate"
    if lr >= 2.08:
        return "Supporting"
    return None


def lr_to_bp5_strength(lr: float) -> Optional[str]:
    if lr <= 0.00285:
        return "Very Strong"
    if lr <= 0.05:
        return "Strong"
    if lr <= 0.23:
        return "Moderate"
    if lr <= 0.48:
        return "Supporting"
    return None


def load_pp4_bp5_snapshot() -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    """Load and validate the snapshot. Missing or corrupted data is fatal."""
    global _SNAPSHOT, _ALIASES
    if _SNAPSHOT is not None and _ALIASES is not None:
        return _SNAPSHOT, _ALIASES
    if not SNAPSHOT_PATH.is_file() or not METADATA_PATH.is_file():
        raise RuntimeError("PP4/BP5 clinical LR snapshot or its metadata is missing")

    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    if metadata.get("status") != "validated_derived_snapshot":
        raise RuntimeError("PP4/BP5 clinical LR snapshot is not validated")
    if metadata.get("index_sha256") != _sha256(SNAPSHOT_PATH):
        raise RuntimeError("PP4/BP5 clinical LR snapshot checksum does not match metadata")
    if not SOURCE_MANIFEST_PATH.is_file():
        raise RuntimeError("PP4/BP5 clinical LR source manifest is missing")
    if metadata.get("source_manifest_sha256") != _sha256(SOURCE_MANIFEST_PATH):
        raise RuntimeError("PP4/BP5 clinical LR source manifest checksum mismatch")
    normalization = metadata.get("normalization")
    if not isinstance(normalization, dict) or not normalization.get("provenance"):
        raise RuntimeError("PP4/BP5 clinical LR snapshot normalization provenance is missing")
    indel_dependency = normalization.get("normalized_indel_dependency")
    if not isinstance(indel_dependency, dict):
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel dependency is missing")
    if not INDEL_SNAPSHOT_PATH.is_file() or not INDEL_METADATA_PATH.is_file():
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel dependency is unavailable")
    if indel_dependency.get("index_sha256") != _sha256(INDEL_SNAPSHOT_PATH):
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel dependency checksum mismatch")
    if indel_dependency.get("metadata_sha256") != _sha256(INDEL_METADATA_PATH):
        raise RuntimeError("PP4/BP5 clinical LR snapshot indel metadata checksum mismatch")

    snapshot = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    if metadata.get("records") != len(snapshot):
        raise RuntimeError("PP4/BP5 clinical LR snapshot record count does not match metadata")

    aliases: dict[str, str] = {key: key for key in snapshot}
    ambiguous: set[str] = set()
    for canonical_key, record in snapshot.items():
        for notation in record.get("input_c_notations", []):
            alias = f"{record['gene']}:{notation}"
            previous = aliases.get(alias)
            if previous and previous != canonical_key:
                ambiguous.add(alias)
            else:
                aliases[alias] = canonical_key
    for alias in ambiguous:
        aliases.pop(alias, None)
    if ambiguous:
        raise RuntimeError(f"PP4/BP5 clinical LR snapshot has ambiguous aliases: {len(ambiguous)}")

    _SNAPSHOT, _ALIASES = snapshot, aliases
    return snapshot, aliases


def evaluate_pp4_bp5(gene: str, c_notation: str) -> Dict:
    snapshot, aliases = load_pp4_bp5_snapshot()
    query_key = f"{gene}:{c_notation}"
    canonical_key = query_key if query_key in snapshot else aliases.get(query_key)
    entry = snapshot.get(canonical_key) if canonical_key else None
    result = {
        "applies": False, "code": None, "strength": None, "points": 0,
        "reason": "", "likelihood_ratio": None,
        "source_components": [],
    }
    if entry is None:
        result["reason"] = (
            "No informative variant-specific combined clinical LR is available; "
            "PP4/BP5 is not applied under ENIGMA v1.2"
        )
        return result

    lr = entry["combined_lr"]
    result["likelihood_ratio"] = lr
    result["source_components"] = entry.get("source_components", [])
    code = entry.get("criterion")
    if not code:
        result["reason"] = f"Combined clinical LR={lr:.6g} is not informative for PP4 or BP5"
        return result

    result.update({
        "applies": True,
        "code": code,
        "strength": entry["strength"],
        "points": entry["points"],
    })
    pmids = sorted({component["pmid"] for component in result["source_components"]})
    result["reason"] = (
        f"ENIGMA v1.2 combined clinical evidence: "
        f"combined LR={lr:.6g}; {code} {entry['strength']}; "
        f"PMID {', '.join(pmids)}"
    )
    return result
