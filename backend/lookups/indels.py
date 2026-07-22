"""Validated lookup for the versioned normalized BRCA coding-indel snapshot."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Optional


ROOT = Path(__file__).resolve().parents[2]
INDEX_PATH = ROOT / "data/precomputed/brca_normalized_indel_snapshot.index.json"
METADATA_PATH = ROOT / "data/precomputed/brca_normalized_indel_snapshot.metadata.json"
_INDEX: Optional[dict[str, dict[str, Any]]] = None
_ALIASES: Optional[dict[str, str]] = None


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_indel_snapshot() -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    global _INDEX, _ALIASES
    if _INDEX is not None and _ALIASES is not None:
        return _INDEX, _ALIASES
    if not INDEX_PATH.is_file() or not METADATA_PATH.is_file():
        raise RuntimeError("Required normalized BRCA indel snapshot or metadata is missing")
    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    if metadata.get("status") != "validated_reference_snapshot":
        raise RuntimeError("Normalized BRCA indel snapshot is not validated")
    if metadata.get("index_sha256") != _sha256(INDEX_PATH):
        raise RuntimeError("Normalized BRCA indel snapshot checksum mismatch")
    index = json.loads(INDEX_PATH.read_text(encoding="utf-8"))
    if len(index) != metadata.get("records"):
        raise RuntimeError("Normalized BRCA indel snapshot record count mismatch")
    aliases: dict[str, str] = {key: key for key in index}
    for canonical_key, record in index.items():
        gene = record["gene"]
        for notation in record.get("input_c_notations", []):
            alias_key = f"{gene}:{notation}"
            if alias_key in index:
                continue
            previous = aliases.setdefault(alias_key, canonical_key)
            if previous != canonical_key:
                raise RuntimeError(f"Ambiguous normalized BRCA indel alias: {alias_key}")
    _INDEX, _ALIASES = index, aliases
    return _INDEX, _ALIASES


def lookup_indel_snapshot(gene: str, c_notation: str) -> Optional[dict[str, Any]]:
    index, aliases = load_indel_snapshot()
    alias_key = f"{gene.strip().upper()}:{c_notation.strip()}"
    canonical_key = aliases.get(alias_key)
    if not canonical_key:
        return None
    return index[canonical_key]
