#!/usr/bin/env python3
"""Rebuild coding-SNV protein fields with the pinned local HGVS engine."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections import Counter
from datetime import date
from pathlib import Path

from backend.modules.hgvs_engine import derive_protein_consequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INDEX = ROOT / "data/precomputed/brca_module1_snv_classification_snapshot.index.json"
DEFAULT_METADATA = ROOT / "data/precomputed/brca_module1_snv_classification_snapshot.metadata.json"


def atomic_write(path: Path, content: bytes) -> None:
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--created", default=date.today().isoformat())
    args = parser.parse_args()

    records = json.loads(args.index.read_text(encoding="utf-8"))
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    if not isinstance(records, dict) or len(records) != 47547:
        raise RuntimeError(f"Expected 47,547 coding SNVs, found {len(records)}")
    changed = Counter()
    provenance = None
    for key, record in records.items():
        gene, c_notation = key.split(":", 1)
        result = derive_protein_consequence(gene, c_notation)
        provenance = result.provenance
        if result.canonical_c_notation != c_notation:
            raise RuntimeError(
                f"Coding SNV normalization changed key {key} to {result.canonical_c_notation}"
            )
        if record.get("p_notation") != result.p_notation:
            changed[f"{record.get('variant_type', 'unknown')} protein consequence"] += 1
            record["p_notation"] = result.p_notation

    # Preserve record order to keep the generated diff localized to changed
    # fields while making the bytes deterministic for the same input.
    index_content = json.dumps(records, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    index_sha = sha256_bytes(index_content)
    metadata["created"] = args.created
    metadata["index"] = str(args.index.relative_to(ROOT)).replace("\\", "/")
    metadata["index_sha256"] = index_sha
    metadata["protein_consequence_normalization"] = {
        "status": "complete",
        "records_checked": len(records),
        "records_changed": sum(changed.values()),
        "changes": dict(sorted(changed.items())),
        "method": "sequence-derived reference-transcript HGVS",
        "builder": "scripts/refresh_snv_snapshot_protein_consequences.py",
        "provenance": provenance or {},
    }
    metadata_content = (
        json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    atomic_write(args.index, index_content)
    atomic_write(args.metadata, metadata_content)
    print(json.dumps({"index_sha256": index_sha, "changes": dict(changed)}, indent=2))


if __name__ == "__main__":
    main()
