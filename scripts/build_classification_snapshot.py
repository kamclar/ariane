import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
DEFAULT_SUMMARY = REPO_ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.summary.json"
DEFAULT_INDEX = REPO_ROOT / "data" / "precomputed" / "brca_module1_snv_classification_snapshot.index.json"
DEFAULT_METADATA = REPO_ROOT / "data" / "precomputed" / "brca_module1_snv_classification_snapshot.metadata.json"
DEFAULT_SPLICEAI_METADATA = REPO_ROOT / "data" / "spliceai" / "spliceai_brca_snv_reference_cache.metadata.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def as_int(value: Any) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def as_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def criteria_codes(value: str) -> list[str]:
    codes = []
    for item in (value or "").split(";"):
        if not item:
            continue
        codes.append(item.split(":", 1)[0])
    return codes


def compact_record(row: dict[str, str]) -> dict[str, Any]:
    return {
        "gene": row["gene"],
        "c_notation": row["c_notation"],
        "p_notation": row.get("p_notation", ""),
        "variant_type": row.get("variant_type", ""),
        "predicted_class": as_int(row.get("predicted_class")),
        "predicted_label": row.get("predicted_label", ""),
        "total_points": as_int(row.get("total_points")),
        "criteria": row.get("criteria", ""),
        "criteria_codes": criteria_codes(row.get("criteria", "")),
        "classification_note": row.get("classification_note", ""),
        "spliceai_score": as_float(row.get("spliceai_score")),
        "bayesdel_score": as_float(row.get("bayesdel_score")),
        "alphamissense_score": as_float(row.get("alphamissense_score")),
        "alphamissense_class": row.get("alphamissense_class", ""),
        "coordinate_status": row.get("coordinate_status", ""),
        "grch37": row.get("grch37", ""),
        "grch38": row.get("grch38", ""),
        "gnomad_status": row.get("gnomad_status", ""),
        "warnings": row.get("warnings", ""),
    }


def build_index(input_csv: Path) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    index = {}
    class_counter: Counter[str] = Counter()
    type_counter: Counter[str] = Counter()
    criteria_counter: Counter[str] = Counter()
    duplicate_keys = []

    with input_csv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = f"{row['gene']}:{row['c_notation']}"
            if key in index:
                duplicate_keys.append(key)
            record = compact_record(row)
            index[key] = record
            class_counter[str(record["predicted_class"])] += 1
            type_counter[record["variant_type"]] += 1
            criteria_counter.update(record["criteria_codes"])

    summary = {
        "n_records": len(index),
        "duplicate_keys": duplicate_keys,
        "classes": dict(class_counter),
        "variant_types": dict(type_counter),
        "criteria": dict(criteria_counter),
        "with_spliceai_score": sum(1 for row in index.values() if row["spliceai_score"] is not None),
        "coordinate_status": dict(Counter(row["coordinate_status"] for row in index.values())),
    }
    return index, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build ARIANE precomputed classification snapshot index.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--spliceai-metadata", type=Path, default=DEFAULT_SPLICEAI_METADATA)
    args = parser.parse_args()

    index, computed_summary = build_index(args.input)
    source_summary = json.loads(args.summary.read_text(encoding="utf-8")) if args.summary.exists() else {}
    spliceai_metadata = (
        json.loads(args.spliceai_metadata.read_text(encoding="utf-8"))
        if args.spliceai_metadata.exists()
        else {}
    )

    args.index.parent.mkdir(parents=True, exist_ok=True)
    args.index.write_text(json.dumps(index, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")

    metadata = {
        "dataset": "BRCA1/BRCA2 coding SNV ARIANE Module 1 classification snapshot",
        "created": str(date.today()),
        "status": "snapshot_not_authoritative",
        "application": "ARIANE",
        "application_version": "1.8.0",
        "scope": "BRCA1/BRCA2 coding SNVs from variant_space_scan manifest",
        "input_csv": str(args.input),
        "input_csv_sha256": sha256_file(args.input),
        "source_summary": str(args.summary),
        "source_summary_sha256": sha256_file(args.summary) if args.summary.exists() else "",
        "index": str(args.index),
        "index_sha256": sha256_file(args.index),
        "n_records": computed_summary["n_records"],
        "duplicate_keys_count": len(computed_summary["duplicate_keys"]),
        "duplicate_keys_first_20": computed_summary["duplicate_keys"][:20],
        "classes": computed_summary["classes"],
        "variant_types": computed_summary["variant_types"],
        "criteria": computed_summary["criteria"],
        "with_spliceai_score": computed_summary["with_spliceai_score"],
        "coordinate_status": computed_summary["coordinate_status"],
        "source_summary_counts": source_summary,
        "dependencies": {
            "spliceai_cache": spliceai_metadata,
            "classification_overlay_report": "variant_space_scan/docs/full_snv_classification_overlay_report.md",
            "spliceai_precompute_report": "variant_space_scan/docs/final_spliceai_precompute_report.md",
        },
        "intended_use": (
            "Fast lookup and landscape analysis snapshot. Do not use as the sole "
            "authoritative classification without checking that rule and data versions match."
        ),
    }
    args.metadata.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {args.index}")
    print(f"wrote {args.metadata}")
    print(json.dumps({
        "n_records": metadata["n_records"],
        "classes": metadata["classes"],
        "index_sha256": metadata["index_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
