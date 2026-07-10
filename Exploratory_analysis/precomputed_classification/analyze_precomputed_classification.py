"""Exploratory analysis for the BRCA Module 1 precomputed SNV snapshot.

The analysis is intentionally kept outside the ARIANE server code. It reads the
already generated full SNV classification CSV and writes small derived tables
plus a markdown report for manual review.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
INPUT_SUMMARY = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.summary.json"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = OUT_DIR / "tables"
REPORT = OUT_DIR / "precomputed_classification_analysis.md"

CLASS_LABELS = {
    "1": "Benign",
    "2": "Likely Benign",
    "3": "VUS",
    "4": "Likely Pathogenic",
    "5": "Pathogenic",
}
GROUP_LABELS = {
    "benign": "Benign/Likely Benign",
    "vus": "VUS",
    "pathogenic": "Likely Pathogenic/Pathogenic",
}
GROUP_ORDER = ["benign", "vus", "pathogenic"]


def class_group(predicted_class: str | object) -> str:
    cls = str(predicted_class)
    if cls in {"1", "2"}:
        return "benign"
    if cls in {"4", "5"}:
        return "pathogenic"
    return "vus"


def pct(part: int, whole: int) -> str:
    if whole == 0:
        return "0.00"
    return f"{(part / whole) * 100:.2f}"


def parse_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def spliceai_bin(score: float | None) -> str:
    if score is None:
        return "missing"
    if score == 0:
        return "0"
    if score < 0.10:
        return "0-0.099"
    if score < 0.20:
        return "0.10-0.199"
    if score < 0.50:
        return "0.20-0.499"
    if score < 0.80:
        return "0.50-0.799"
    return "0.80-1.00"


def gnomad_coarse_status(status: str | None) -> str:
    if not status:
        return "missing"
    if status.startswith("status=found"):
        return "found"
    if status.startswith("status=absent_v2_only"):
        return "absent_v2_only"
    if "outside_cached_region" in status:
        return "outside_cached_region"
    return status.split(";", 1)[0].replace("status=", "")


def parse_criteria(criteria: str) -> list[dict[str, str | int]]:
    parsed = []
    if not criteria:
        return parsed
    for item in criteria.split(";"):
        parts = item.split(":")
        if len(parts) != 3:
            parsed.append({"code": item, "strength": "", "points": 0})
            continue
        code, strength, points = parts
        try:
            point_value = int(points)
        except ValueError:
            point_value = 0
        parsed.append({"code": code, "strength": strength, "points": point_value})
    return parsed


def write_counter(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_dict_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Missing input CSV: {INPUT_CSV}")

    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["_spliceai"] = parse_float(row.get("spliceai_score"))
            row["_criteria"] = parse_criteria(row.get("criteria", ""))
            row["_criteria_codes"] = {item["code"] for item in row["_criteria"]}
            row["_gnomad_coarse"] = gnomad_coarse_status(row.get("gnomad_status"))
            rows.append(row)

    total = len(rows)
    by_class = Counter(row["predicted_class"] for row in rows)
    by_gene = Counter(row["gene"] for row in rows)
    by_type = Counter(row["variant_type"] for row in rows)
    by_gnomad = Counter(row["_gnomad_coarse"] for row in rows)
    by_class_gene = Counter((row["gene"], row["predicted_class"]) for row in rows)
    by_class_type = Counter((row["variant_type"], row["predicted_class"]) for row in rows)
    by_spliceai_bin = Counter((spliceai_bin(row["_spliceai"]), row["predicted_class"]) for row in rows)
    by_gene_group = Counter((row["gene"], class_group(row["predicted_class"])) for row in rows)
    by_gene_type = Counter((row["gene"], row["variant_type"]) for row in rows)
    by_gene_spliceai_bin = Counter((row["gene"], spliceai_bin(row["_spliceai"])) for row in rows)
    by_gene_gnomad = Counter((row["gene"], row["_gnomad_coarse"]) for row in rows)

    criteria_counts = Counter()
    criteria_by_class = Counter()
    criteria_points = defaultdict(Counter)
    for row in rows:
        cls = row["predicted_class"]
        for criterion in row["_criteria"]:
            code = criterion["code"]
            criteria_counts[code] += 1
            criteria_by_class[(code, cls)] += 1
            criteria_points[code][criterion["points"]] += 1

    write_counter(
        TABLE_DIR / "class_distribution.csv",
        ["predicted_class", "label", "count", "percent"],
        [
            [cls, CLASS_LABELS.get(cls, ""), by_class[cls], pct(by_class[cls], total)]
            for cls in sorted(CLASS_LABELS, key=int)
        ],
    )

    write_counter(
        TABLE_DIR / "gene_by_class.csv",
        ["gene", "predicted_class", "label", "count", "percent_within_gene"],
        [
            [gene, cls, CLASS_LABELS.get(cls, ""), by_class_gene[(gene, cls)], pct(by_class_gene[(gene, cls)], by_gene[gene])]
            for gene in sorted(by_gene)
            for cls in sorted(CLASS_LABELS, key=int)
        ],
    )

    write_counter(
        TABLE_DIR / "gene_normalized_class_distribution.csv",
        ["gene", "predicted_class", "label", "count", "percent_within_gene", "per_1000_gene_variants"],
        [
            [
                gene,
                cls,
                CLASS_LABELS.get(cls, ""),
                by_class_gene[(gene, cls)],
                pct(by_class_gene[(gene, cls)], by_gene[gene]),
                f"{(by_class_gene[(gene, cls)] / by_gene[gene] * 1000):.2f}" if by_gene[gene] else "0.00",
            ]
            for gene in sorted(by_gene)
            for cls in sorted(CLASS_LABELS, key=int)
        ],
    )

    write_counter(
        TABLE_DIR / "gene_normalized_grouped_class_distribution.csv",
        ["gene", "grouped_class", "label", "count", "percent_within_gene", "per_1000_gene_variants"],
        [
            [
                gene,
                group,
                GROUP_LABELS[group],
                by_gene_group[(gene, group)],
                pct(by_gene_group[(gene, group)], by_gene[gene]),
                f"{(by_gene_group[(gene, group)] / by_gene[gene] * 1000):.2f}" if by_gene[gene] else "0.00",
            ]
            for gene in sorted(by_gene)
            for group in GROUP_ORDER
        ],
    )

    write_counter(
        TABLE_DIR / "variant_type_by_class.csv",
        ["variant_type", "predicted_class", "label", "count", "percent_within_type"],
        [
            [vtype, cls, CLASS_LABELS.get(cls, ""), by_class_type[(vtype, cls)], pct(by_class_type[(vtype, cls)], by_type[vtype])]
            for vtype in sorted(by_type)
            for cls in sorted(CLASS_LABELS, key=int)
        ],
    )

    write_counter(
        TABLE_DIR / "gene_normalized_variant_type_distribution.csv",
        ["gene", "variant_type", "count", "percent_within_gene", "per_1000_gene_variants"],
        [
            [
                gene,
                vtype,
                by_gene_type[(gene, vtype)],
                pct(by_gene_type[(gene, vtype)], by_gene[gene]),
                f"{(by_gene_type[(gene, vtype)] / by_gene[gene] * 1000):.2f}" if by_gene[gene] else "0.00",
            ]
            for gene in sorted(by_gene)
            for vtype in sorted(by_type)
        ],
    )

    write_counter(
        TABLE_DIR / "criteria_counts.csv",
        ["criterion", "count", "percent_of_all_variants", "point_values"],
        [
            [
                code,
                count,
                pct(count, total),
                ";".join(f"{points}:{n}" for points, n in sorted(criteria_points[code].items())),
            ]
            for code, count in criteria_counts.most_common()
        ],
    )

    write_counter(
        TABLE_DIR / "criteria_by_class.csv",
        ["criterion", "predicted_class", "label", "count"],
        [
            [code, cls, CLASS_LABELS.get(cls, ""), criteria_by_class[(code, cls)]]
            for code in sorted(criteria_counts)
            for cls in sorted(CLASS_LABELS, key=int)
            if criteria_by_class[(code, cls)]
        ],
    )

    criteria_by_gene = Counter()
    for row in rows:
        for criterion in row["_criteria"]:
            criteria_by_gene[(row["gene"], criterion["code"])] += 1

    write_counter(
        TABLE_DIR / "gene_normalized_criteria_counts.csv",
        ["gene", "criterion", "count", "percent_within_gene", "per_1000_gene_variants"],
        [
            [
                gene,
                code,
                criteria_by_gene[(gene, code)],
                pct(criteria_by_gene[(gene, code)], by_gene[gene]),
                f"{(criteria_by_gene[(gene, code)] / by_gene[gene] * 1000):.2f}" if by_gene[gene] else "0.00",
            ]
            for gene in sorted(by_gene)
            for code in sorted(criteria_counts)
        ],
    )

    splice_bins = ["0", "0-0.099", "0.10-0.199", "0.20-0.499", "0.50-0.799", "0.80-1.00", "missing"]
    write_counter(
        TABLE_DIR / "spliceai_bins_by_class.csv",
        ["spliceai_bin", "predicted_class", "label", "count"],
        [
            [bin_name, cls, CLASS_LABELS.get(cls, ""), by_spliceai_bin[(bin_name, cls)]]
            for bin_name in splice_bins
            for cls in sorted(CLASS_LABELS, key=int)
            if by_spliceai_bin[(bin_name, cls)]
        ],
    )

    write_counter(
        TABLE_DIR / "gene_normalized_spliceai_bins.csv",
        ["gene", "spliceai_bin", "count", "percent_within_gene", "per_1000_gene_variants"],
        [
            [
                gene,
                bin_name,
                by_gene_spliceai_bin[(gene, bin_name)],
                pct(by_gene_spliceai_bin[(gene, bin_name)], by_gene[gene]),
                f"{(by_gene_spliceai_bin[(gene, bin_name)] / by_gene[gene] * 1000):.2f}" if by_gene[gene] else "0.00",
            ]
            for gene in sorted(by_gene)
            for bin_name in splice_bins
        ],
    )

    write_counter(
        TABLE_DIR / "gnomad_coarse_status.csv",
        ["gnomad_status", "count", "percent"],
        [[status, count, pct(count, total)] for status, count in by_gnomad.most_common()],
    )

    write_counter(
        TABLE_DIR / "gene_normalized_gnomad_coarse_status.csv",
        ["gene", "gnomad_status", "count", "percent_within_gene", "per_1000_gene_variants"],
        [
            [
                gene,
                status,
                by_gene_gnomad[(gene, status)],
                pct(by_gene_gnomad[(gene, status)], by_gene[gene]),
                f"{(by_gene_gnomad[(gene, status)] / by_gene[gene] * 1000):.2f}" if by_gene[gene] else "0.00",
            ]
            for gene in sorted(by_gene)
            for status in sorted(by_gnomad)
        ],
    )

    high_spliceai = sorted(
        [row for row in rows if row["_spliceai"] is not None and row["_spliceai"] >= 0.20],
        key=lambda row: (row["_spliceai"], row["gene"], row["c_notation"]),
        reverse=True,
    )
    key_fields = [
        "gene",
        "c_notation",
        "p_notation",
        "variant_type",
        "spliceai_score",
        "criteria",
        "total_points",
        "predicted_class",
        "predicted_label",
        "gnomad_status",
        "warnings",
    ]
    write_dict_rows(TABLE_DIR / "high_spliceai_variants_ge_0_20.csv", high_spliceai, key_fields)

    attention = []
    for row in rows:
        if row["predicted_class"] != "3":
            continue
        codes = row["_criteria_codes"]
        reasons = []
        if row["_spliceai"] is not None and row["_spliceai"] >= 0.20:
            reasons.append("high_spliceai")
        if "PVS1" in codes:
            reasons.append("pvs1_but_vus")
        if "PS3" in codes or "BS3" in codes:
            reasons.append("functional_evidence_but_vus")
        if "PM5_PTC" in codes:
            reasons.append("ptc_same_codon_but_vus")
        if reasons:
            copied = {field: row.get(field, "") for field in key_fields}
            copied["attention_reason"] = ";".join(reasons)
            attention.append(copied)

    write_dict_rows(
        TABLE_DIR / "vus_attention_list.csv",
        attention,
        ["attention_reason"] + key_fields,
    )

    lp_p = [row for row in rows if row["predicted_class"] in {"4", "5"}]
    write_dict_rows(TABLE_DIR / "likely_pathogenic_pathogenic_variants.csv", lp_p, key_fields)

    summary = {}
    if INPUT_SUMMARY.exists():
        summary = json.loads(INPUT_SUMMARY.read_text(encoding="utf-8"))

    report = build_report(
        total=total,
        by_class=by_class,
        by_gene=by_gene,
        by_type=by_type,
        by_gnomad=by_gnomad,
        by_class_gene=by_class_gene,
        by_gene_group=by_gene_group,
        criteria_counts=criteria_counts,
        high_spliceai=high_spliceai,
        attention=attention,
        summary=summary,
    )
    REPORT.write_text(report, encoding="utf-8")
    print(f"Wrote analysis to {OUT_DIR}")


def build_report(
    *,
    total: int,
    by_class: Counter,
    by_gene: Counter,
    by_type: Counter,
    by_gnomad: Counter,
    by_class_gene: Counter,
    by_gene_group: Counter,
    criteria_counts: Counter,
    high_spliceai: list[dict[str, object]],
    attention: list[dict[str, object]],
    summary: dict,
) -> str:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    class_lines = "\n".join(
        f"| {cls} | {CLASS_LABELS[cls]} | {by_class[cls]} | {pct(by_class[cls], total)}% |"
        for cls in sorted(CLASS_LABELS, key=int)
    )
    gene_lines = "\n".join(f"| {gene} | {count} | {pct(count, total)}% |" for gene, count in sorted(by_gene.items()))
    gene_class_lines = "\n".join(
        f"| {gene} | {cls} | {CLASS_LABELS[cls]} | {by_class_gene[(gene, cls)]} | {pct(by_class_gene[(gene, cls)], by_gene[gene])}% |"
        for gene in sorted(by_gene)
        for cls in sorted(CLASS_LABELS, key=int)
    )
    gene_group_lines = "\n".join(
        f"| {gene} | {GROUP_LABELS[group]} | {by_gene_group[(gene, group)]} | {pct(by_gene_group[(gene, group)], by_gene[gene])}% |"
        for gene in sorted(by_gene)
        for group in GROUP_ORDER
    )
    type_lines = "\n".join(f"| {vtype} | {count} | {pct(count, total)}% |" for vtype, count in sorted(by_type.items()))
    gnomad_lines = "\n".join(f"| {status} | {count} | {pct(count, total)}% |" for status, count in by_gnomad.most_common())
    criterion_lines = "\n".join(
        f"| {code} | {count} | {pct(count, total)}% |"
        for code, count in criteria_counts.most_common(15)
    )

    high_spliceai_by_class = Counter(row["predicted_class"] for row in high_spliceai)
    high_spliceai_lines = "\n".join(
        f"| {cls} | {CLASS_LABELS[cls]} | {high_spliceai_by_class[cls]} |"
        for cls in sorted(CLASS_LABELS, key=int)
        if high_spliceai_by_class[cls]
    )

    return f"""# Exploratory Analysis: Precomputed BRCA Module 1 SNV Classification

Generated: {generated}

Input snapshot: `{INPUT_CSV.relative_to(ROOT)}`

This is an exploratory analysis of the automated Module 1 overlay. It is not an expert final classification and it does not add non-automated ACMG/ENIGMA criteria such as PS4, PM3, PP1, BS2, BS4, or curated RNA conclusions.

## Dataset

- Total variants: {total}
- Coordinate status: {summary.get("coordinate_status", {})}
- Variants with SpliceAI score: {summary.get("with_spliceai_score", "unknown")}

## Class Distribution

| Class | Label | Count | Percent |
|---|---:|---:|---:|
{class_lines}

## Gene Coverage

| Gene | Count | Percent |
|---|---:|---:|
{gene_lines}

## Gene-Normalized Class Distribution

Because BRCA1 and BRCA2 contribute different numbers of coding SNVs to the
snapshot, cross-gene comparisons should use within-gene percentages or rates
per 1000 generated variants, not raw combined counts.

| Gene | Class | Label | Count | Percent within gene |
|---|---:|---|---:|---:|
{gene_class_lines}

### Grouped Classes By Gene

| Gene | Grouped class | Count | Percent within gene |
|---|---|---:|---:|
{gene_group_lines}

## Variant Types

| Variant type | Count | Percent |
|---|---:|---:|
{type_lines}

## Most Frequent Criteria

| Criterion | Count | Percent of all variants |
|---|---:|---:|
{criterion_lines}

## SpliceAI Signal

Variants with reference-transcript SpliceAI score >= 0.20: {len(high_spliceai)}

| Class | Label | Count |
|---|---:|---:|
{high_spliceai_lines}

The full list is in `tables/high_spliceai_variants_ge_0_20.csv`.

## VUS Attention List

VUS variants flagged for manual review: {len(attention)}

The list is in `tables/vus_attention_list.csv`. Reasons include high SpliceAI, PVS1 still ending as VUS, functional evidence still ending as VUS, or PM5_PTC still ending as VUS.

## gnomAD Status

| gnomAD status | Count | Percent |
|---|---:|---:|
{gnomad_lines}

Interpret this cautiously: the current overlay uses the available local frequency lookup state. It is useful for spotting patterns, but it should not be treated as a complete population-frequency audit.

## Output Tables

- `tables/class_distribution.csv`
- `tables/gene_by_class.csv`
- `tables/gene_normalized_class_distribution.csv`
- `tables/gene_normalized_grouped_class_distribution.csv`
- `tables/gene_normalized_variant_type_distribution.csv`
- `tables/gene_normalized_criteria_counts.csv`
- `tables/gene_normalized_spliceai_bins.csv`
- `tables/gene_normalized_gnomad_coarse_status.csv`
- `tables/variant_type_by_class.csv`
- `tables/criteria_counts.csv`
- `tables/criteria_by_class.csv`
- `tables/spliceai_bins_by_class.csv`
- `tables/gnomad_coarse_status.csv`
- `tables/high_spliceai_variants_ge_0_20.csv`
- `tables/vus_attention_list.csv`
- `tables/likely_pathogenic_pathogenic_variants.csv`

## Next Analytical Questions

1. Compare VUS with high SpliceAI against exon/intron boundary context and known ENIGMA splice rules.
2. Review PVS1-only or PVS1-dominated VUS variants to see whether Table 3 combinations are behaving as expected.
3. Separate biological patterns from current evidence availability, especially for PM2 and sparse functional or curated evidence.
4. Add plots after the tabular checks are stable, for example class by CDS position, criteria heatmaps, and SpliceAI score density.
"""


if __name__ == "__main__":
    main()
