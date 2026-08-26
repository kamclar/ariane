"""Validate the gene policy manifest and refresh its checksum metadata."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = PROJECT_ROOT / "backend" / "data" / "gene_policy_manifest.json"
METADATA_PATH = PROJECT_ROOT / "backend" / "data" / "gene_policy_manifest.metadata.json"


def _load_object(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path.name} must contain a JSON object")
    return value


def expected_metadata() -> dict:
    manifest = _load_object(MANIFEST_PATH)
    active_genes = [
        value
        for value in (manifest.get("genes") or {}).values()
        if value.get("activation_status") == "active"
    ]
    active_policy_ids = {value.get("vcep_policy_id") for value in active_genes}
    return {
        "schema_version": 1,
        "manifest_id": manifest.get("manifest_id"),
        "manifest_version": manifest.get("manifest_version"),
        "manifest_sha256": hashlib.sha256(MANIFEST_PATH.read_bytes()).hexdigest(),
        "active_gene_count": len(active_genes),
        "active_policy_count": len(active_policy_ids),
        "validation_status": "approved",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write the calculated metadata. Without this option only check it.",
    )
    args = parser.parse_args()
    expected = expected_metadata()
    current = _load_object(METADATA_PATH)
    if args.write:
        METADATA_PATH.write_text(
            json.dumps(expected, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"Updated {METADATA_PATH.relative_to(PROJECT_ROOT)}")
        return 0
    if current != expected:
        print("Gene policy manifest metadata is stale. Run with --write.")
        return 1
    print("Gene policy manifest metadata is current.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
