"""Export current protein-PS1 candidates for external methodological review."""

from __future__ import annotations

import csv
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = PROJECT_ROOT / "backend" / "data" / "ps1_protein_reference_registry.json"
OUTPUT_PATH = PROJECT_ROOT / "docs" / "ps1_reference_review_candidates.tsv"

FIELDS = [
    "reference_id",
    "gene",
    "transcript",
    "c_notation",
    "p_notation",
    "st7_classification",
    "current_registry_status",
    "current_blocker",
    "protein_mechanism_basis",
    "table9_code",
    "table9_strength",
    "confirmed_splice_status",
    "checked_rna_splice_sources",
    "reference_spliceai_threshold",
    "reference_spliceai_result",
    "vcep_assertion_id",
    "vcep_classification",
    "vcep_ruleset_version",
    "vcep_assertion_url",
    "reference_classification_used_ps1",
    "ps1_dependency_references",
    "reviewer_eligibility_decision",
    "reviewer_note",
]


def export() -> int:
    registry = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
    candidates = [
        record
        for record in registry.get("references", [])
        if record.get("status") == "review_required"
    ]
    candidates.sort(
        key=lambda record: (
            str(record.get("gene") or ""),
            str(record.get("c_notation") or ""),
        )
    )

    with OUTPUT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, dialect="excel-tab")
        writer.writeheader()
        for record in candidates:
            mechanism = record.get("protein_mechanism_evidence") or {}
            splice = record.get("reference_splice_evidence") or {}
            dependency = record.get("classification_ps1_dependency") or {}
            writer.writerow(
                {
                    "reference_id": record.get("reference_id", ""),
                    "gene": record.get("gene", ""),
                    "transcript": record.get("transcript", ""),
                    "c_notation": record.get("c_notation", ""),
                    "p_notation": record.get("p_notation", ""),
                    "st7_classification": record.get("classification", ""),
                    "current_registry_status": record.get("status", ""),
                    "current_blocker": record.get("status_reason", ""),
                    "protein_mechanism_basis": mechanism.get("basis", ""),
                    "table9_code": mechanism.get("table9_code", ""),
                    "table9_strength": mechanism.get("table9_strength", ""),
                    "confirmed_splice_status": splice.get("confirmed_status", ""),
                    "checked_rna_splice_sources": "; ".join(
                        str(value) for value in splice.get("sources_checked", [])
                    ),
                    "reference_spliceai_threshold": splice.get("threshold", ""),
                    "reference_spliceai_result": "runtime lookup required",
                    "vcep_assertion_id": "",
                    "vcep_classification": "",
                    "vcep_ruleset_version": "",
                    "vcep_assertion_url": "",
                    "reference_classification_used_ps1": (
                        "" if dependency.get("used") == "unknown" else dependency.get("used", "")
                    ),
                    "ps1_dependency_references": "; ".join(
                        str(value) for value in dependency.get("reference_ids", [])
                    ),
                    "reviewer_eligibility_decision": "",
                    "reviewer_note": "",
                }
            )
    return len(candidates)


def main() -> None:
    count = export()
    print(f"Wrote {count} review candidates to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
