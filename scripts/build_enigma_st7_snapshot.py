"""Build the lossless ENIGMA Supplementary Table 7 reference snapshot."""
from __future__ import annotations

import argparse
import json
import re
from datetime import date
from pathlib import Path

import openpyxl


SOURCE_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "cb4a09fe-30f4-4aa8-9d76-d7ea407c9754/data"
)

FIELDS = [
    "gene", "c_notation", "p_notation", "reference_set",
    "prior_probability", "posterior_probability", "iarc_class", "source",
    "frequency_category", "maximum_maf", "total_alleles_gnomad",
    "total_alleles_outbred", "homozygote_count", "allele_count_african",
    "allele_number_african", "maf_african", "allele_count_latino",
    "allele_number_latino", "maf_latino", "allele_count_east_asian",
    "allele_number_east_asian", "maf_east_asian",
    "allele_count_non_finnish_european",
    "allele_number_non_finnish_european", "maf_non_finnish_european",
    "allele_count_south_asian", "allele_number_south_asian",
    "maf_south_asian",
]


def _clean(value):
    if isinstance(value, str):
        value = value.strip()
        return None if value in {"", "N/A"} else value
    return value


def _probability(value):
    value = _clean(value)
    return float(value) if isinstance(value, (int, float)) else None


def _iarc_class(value) -> int:
    match = re.match(r"\s*([1245])(?:\D|$)", str(value))
    if not match:
        raise RuntimeError(f"Invalid ST7 IARC class: {value!r}")
    return int(match.group(1))


def build_snapshot(source: Path) -> dict:
    sheet = openpyxl.load_workbook(
        source, read_only=True, data_only=True
    )["ST7 Reference Set Vts"]
    variants = []
    seen = set()
    for row in sheet.iter_rows(min_row=4, values_only=True):
        if row[0] not in {"BRCA1", "BRCA2"}:
            continue
        values = [_clean(value) for value in row[: len(FIELDS)]]
        entry = dict(zip(FIELDS, values))
        entry["prior_probability"] = _probability(row[4])
        entry["posterior_probability"] = _probability(row[5])
        entry["iarc_class"] = _iarc_class(row[6])
        key = (entry["gene"], entry["c_notation"])
        if key in seen:
            raise RuntimeError(f"Duplicate ST7 variant: {key}")
        seen.add(key)
        variants.append(entry)

    return {
        "schema_version": 2,
        "version": "1.2.0",
        "released": "2025-01-09",
        "source": "ClinGen ENIGMA BRCA1/2 VCEP Supplementary Table 7",
        "source_url": SOURCE_URL,
        "generated": date.today().isoformat(),
        "source_columns": len(FIELDS),
        "total_variants": len(variants),
        "variants": variants,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path("backend/data/st7_reference_set.json")
    )
    args = parser.parse_args()
    snapshot = build_snapshot(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {snapshot['total_variants']} lossless ST7 rows to {args.output}")


if __name__ == "__main__":
    main()
