"""Validate or explicitly finalize the curated protein-PS1 reference registry."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from backend.modules.ps1 import (
    compute_approval_basis_checksum,
    validate_ps1_reference_registry,
)


DEFAULT_REGISTRY = Path("backend/data/ps1_protein_reference_registry.json")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("registry", nargs="?", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument(
        "--write-checksums",
        action="store_true",
        help=(
            "Recalculate approval checksums after an intentional curator review. "
            "Without this flag the command is read-only."
        ),
    )
    args = parser.parse_args()
    data = json.loads(args.registry.read_text(encoding="utf-8"))
    if args.write_checksums:
        for record in data.get("references", []):
            record["approval_basis_checksum"] = compute_approval_basis_checksum(record)
        args.registry.write_text(
            json.dumps(data, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    validate_ps1_reference_registry(data)
    print(
        f"PS1 protein registry valid: {len(data.get('references', []))} reference(s), "
        f"statuses {data.get('status_counts', {})}, version {data.get('registry_version')}"
    )


if __name__ == "__main__":
    main()
