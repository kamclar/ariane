"""Structure/function mapping for benign and likely benign generated variants."""

from __future__ import annotations

import csv
import html
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Optional

from structure_function_mapping import (
    INPUT,
    ROOT,
    STRUCTURE_FEATURES,
    aa_position,
    criteria_combo,
    curated_active_site_status,
    feature_hits,
    load_curated_active_site_annotations,
    load_curated_interface_regions,
    load_pathogenic_positions,
    parse_criteria,
    structure_domain_context,
)


ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = ANALYSIS_DIR / "tables" / "benign_structure_function"
PLOT_DIR = ANALYSIS_DIR / "plots" / "19_benign_structure_function"
REPORT = ANALYSIS_DIR / "benign_structure_function_report.md"


def load_rows() -> list[dict]:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def class_group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


def annotate_row(
    row: dict,
    pathogenic_positions: dict,
    active_site_annotations: dict,
    interface_regions: list[dict],
) -> dict:
    criteria = parse_criteria(row["criteria"])
    pos = aa_position(row["p_notation"])
    hits = feature_hits(row["gene"], pos)
    feature_names = [hit["feature"] for hit in hits] or ["outside_mapped_structure_function_features"]
    known_context = ";".join(pathogenic_positions.get(row["gene"], {}).get(pos or -1, []))
    return {
        "gene": row["gene"],
        "c_notation": row["c_notation"],
        "p_notation": row["p_notation"],
        "aa_position": "" if pos is None else pos,
        "variant_type": row["variant_type"],
        "predicted_class": row["predicted_class"],
        "predicted_label": row["predicted_label"],
        "class_group": class_group(row["predicted_class"]),
        "total_points": row["total_points"],
        "spliceai_score": row["spliceai_score"],
        "structure_features": ";".join(feature_names),
        "structure_domain_context": structure_domain_context(hits, known_context),
        "curated_active_site_status": curated_active_site_status(
            row["gene"], pos, active_site_annotations, interface_regions
        ),
        "known_pathogenic_residue_context": known_context,
        "criteria_combo": criteria_combo(criteria),
        "criteria": row["criteria"],
    }


def analyze() -> tuple[list[dict], dict[str, list[dict]]]:
    pathogenic_positions = load_pathogenic_positions()
    active_site_annotations = load_curated_active_site_annotations()
    interface_regions = load_curated_interface_regions()
    annotated_all = [
        annotate_row(row, pathogenic_positions, active_site_annotations, interface_regions)
        for row in load_rows()
    ]
    benign = [row for row in annotated_all if row["class_group"] == "benign"]
    pathogenic = [row for row in annotated_all if row["class_group"] == "pathogenic"]

    tables = {
        "benign_feature_summary": counter_rows(benign, "structure_features", "feature"),
        "benign_gene_feature_summary": gene_feature_rows(benign),
        "benign_active_site_summary": counter_rows(benign, "curated_active_site_status", "curated_active_site_status"),
        "benign_variant_type_summary": counter_rows(benign, "variant_type", "variant_type"),
        "benign_criteria_combo_summary": counter_rows(benign, "criteria_combo", "criteria_combo"),
        "benign_known_pathogenic_residue_summary": known_residue_rows(benign),
        "benign_vs_pathogenic_feature_comparison": feature_comparison(benign, pathogenic),
        "benign_vs_pathogenic_active_site_comparison": active_site_comparison(benign, pathogenic),
    }
    return benign, tables


def counter_rows(rows: list[dict], key: str, output_key: str) -> list[dict]:
    counts = Counter()
    for row in rows:
        for value in str(row[key]).split(";"):
            counts[value] += 1
    return [{output_key: key, "count": value} for key, value in counts.most_common()]


def gene_feature_rows(rows: list[dict]) -> list[dict]:
    counts = Counter()
    for row in rows:
        for feature in str(row["structure_features"]).split(";"):
            counts[(row["gene"], feature)] += 1
    return [
        {"gene": key[0], "feature": key[1], "count": value}
        for key, value in counts.most_common()
    ]


def known_residue_rows(rows: list[dict]) -> list[dict]:
    counts = Counter(
        "known_pathogenic_residue_position" if row["known_pathogenic_residue_context"] else "no_known_pathogenic_residue_position"
        for row in rows
    )
    return [{"known_pathogenic_residue_context": key, "count": value} for key, value in counts.most_common()]


def feature_set(row: dict) -> set[str]:
    return set(str(row["structure_features"]).split(";"))


def feature_comparison(benign: list[dict], pathogenic: list[dict]) -> list[dict]:
    features = sorted(
        {
            feature
            for row in benign + pathogenic
            for feature in feature_set(row)
        }
    )
    benign_total = len(benign)
    pathogenic_total = len(pathogenic)
    rows = []
    for feature in features:
        b = sum(1 for row in benign if feature in feature_set(row))
        p = sum(1 for row in pathogenic if feature in feature_set(row))
        b_rate = b / benign_total if benign_total else 0.0
        p_rate = p / pathogenic_total if pathogenic_total else 0.0
        rows.append(
            {
                "feature": feature,
                "benign_count": b,
                "pathogenic_count": p,
                "benign_percent": f"{100 * b_rate:.3f}",
                "pathogenic_percent": f"{100 * p_rate:.3f}",
                "pathogenic_to_benign_rate_ratio": f"{(p_rate / b_rate):.3f}" if b_rate else "NA",
            }
        )
    return sorted(rows, key=lambda row: float(row["pathogenic_percent"]), reverse=True)


def active_site_comparison(benign: list[dict], pathogenic: list[dict]) -> list[dict]:
    statuses = sorted(
        set(row["curated_active_site_status"] for row in benign + pathogenic)
    )
    benign_total = len(benign)
    pathogenic_total = len(pathogenic)
    rows = []
    for status in statuses:
        b = sum(1 for row in benign if row["curated_active_site_status"] == status)
        p = sum(1 for row in pathogenic if row["curated_active_site_status"] == status)
        if b == 0 and p == 0:
            continue
        b_rate = b / benign_total if benign_total else 0.0
        p_rate = p / pathogenic_total if pathogenic_total else 0.0
        rows.append(
            {
                "curated_active_site_status": status,
                "benign_count": b,
                "pathogenic_count": p,
                "benign_percent": f"{100 * b_rate:.3f}",
                "pathogenic_percent": f"{100 * p_rate:.3f}",
                "pathogenic_to_benign_rate_ratio": f"{(p_rate / b_rate):.3f}" if b_rate else "NA",
            }
        )
    return sorted(rows, key=lambda row: float(row["pathogenic_percent"]), reverse=True)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def bar_svg(path: Path, rows: list[dict], label_key: str, value_key: str, title: str, max_items: int = 20) -> None:
    rows = rows[:max_items]
    width = 1120
    left = 470
    top = 55
    row_h = 31
    chart_w = 560
    height = top + row_h * len(rows) + 55
    max_value = max((int(row[value_key]) for row in rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:13px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".small{font-size:11px;fill:#555}",
        ".bar{fill:#2f7d62}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = int(row[value_key])
        bar_w = int(chart_w * value / max_value) if max_value else 0
        parts.append(f'<text x="24" y="{y + 20}">{html.escape(str(row[label_key]))}</text>')
        parts.append(f'<rect x="{left}" y="{y + 6}" width="{bar_w}" height="19" class="bar"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 21}" class="small">{value}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def comparison_svg(path: Path, rows: list[dict], label_key: str, title: str, max_items: int = 14) -> None:
    rows = rows[:max_items]
    width = 1200
    left = 470
    top = 58
    row_h = 34
    chart_w = 600
    height = top + row_h * len(rows) + 65
    max_percent = max(
        [float(row["benign_percent"]) for row in rows] + [float(row["pathogenic_percent"]) for row in rows] + [1.0]
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:12px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".benign{fill:#2f7d62}",
        ".pathogenic{fill:#b23a48}",
        ".small{font-size:10px;fill:#555}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
        f'<text x="{left}" y="48" class="small">green: benign 1+2, red: pathogenic 4+5, percent within group</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        benign_w = chart_w * float(row["benign_percent"]) / max_percent
        path_w = chart_w * float(row["pathogenic_percent"]) / max_percent
        parts.append(f'<text x="24" y="{y + 19}">{html.escape(str(row[label_key])[:72])}</text>')
        parts.append(f'<rect x="{left}" y="{y + 3}" width="{benign_w:.1f}" height="12" class="benign"/>')
        parts.append(f'<rect x="{left}" y="{y + 17}" width="{path_w:.1f}" height="12" class="pathogenic"/>')
        parts.append(
            f'<text x="{left + max(benign_w, path_w) + 8}" y="{y + 20}" class="small">'
            f'B {row["benign_percent"]}%, P {row["pathogenic_percent"]}%</text>'
        )
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def lollipop_svg(path: Path, rows: list[dict], gene: str) -> None:
    gene_rows = [row for row in rows if row["gene"] == gene and row["aa_position"]]
    if not gene_rows:
        return
    length = 1863 if gene == "BRCA1" else 3418
    width = 1200
    height = 230
    left = 80
    right = 40
    axis_y = 120
    scale = (width - left - right) / length
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:12px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".axis{stroke:#333;stroke-width:3}",
        ".feature{fill:#e2f2ea;stroke:#5a9b7c}",
        ".b1{fill:#1f6f50;stroke:#222;opacity:.7}",
        ".b2{fill:#55a37f;stroke:#222;opacity:.7}",
        ".small{font-size:10px;fill:#555}",
        "</style>",
        f'<text x="24" y="30" class="title">{gene} benign/likely benign variants in protein features</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width - right}" y2="{axis_y}" class="axis"/>',
    ]
    for feature in STRUCTURE_FEATURES.get(gene, []):
        x = left + int(feature["start"]) * scale
        w = max(3, (int(feature["end"]) - int(feature["start"]) + 1) * scale)
        parts.append(f'<rect x="{x}" y="{axis_y - 28}" width="{w}" height="18" class="feature"/>')
        parts.append(f'<text x="{x}" y="{axis_y - 35}" class="small">{html.escape(feature["feature"].split("/")[0])}</text>')
    for row in gene_rows:
        pos = int(row["aa_position"])
        x = left + pos * scale
        y2 = axis_y - 44 if row["predicted_class"] == "1" else axis_y + 44
        klass = "b1" if row["predicted_class"] == "1" else "b2"
        parts.append(f'<line x1="{x}" y1="{axis_y}" x2="{x}" y2="{y2}" stroke="#777" opacity=".4"/>')
        parts.append(f'<circle cx="{x}" cy="{y2}" r="3" class="{klass}"/>')
    parts.append(f'<text x="{left}" y="{axis_y + 38}" class="small">0</text>')
    parts.append(f'<text x="{width - right - 35}" y="{axis_y + 38}" class="small">{length} aa</text>')
    parts.append('<text x="24" y="208" class="small">Upper dots: class 1. Lower dots: class 2. Green boxes: broad structural/function features.</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def markdown_table(rows: list[dict], columns: list[str], limit: int = 30) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(html.escape(str(row.get(col, ""))) for col in columns) + " |")
    return "\n".join(lines)


def build_report(benign: list[dict], tables: dict[str, list[dict]]) -> None:
    total = len(benign)
    feature_hits_count = sum(
        1
        for row in benign
        if row["structure_features"] != "outside_mapped_structure_function_features"
    )
    active_hits = sum(
        1
        for row in benign
        if row["curated_active_site_status"] != "not_curated_in_current_dataset"
    )
    text = f"""# Benign Structure Function Mapping

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Purpose

This is the benign counterpart to the pathogenic structure/function mapping.
It asks where generated benign and likely benign variants, classes 1 and 2, lie
relative to the same structural domains and UniProt active/interface annotations.

This does not mean that a structure annotation is benign evidence by itself. It
is a landscape check: if benign variants also occur in a region, then location
alone is not sufficient for interpretation.

## Dataset

- Benign/likely benign variants analyzed: {total}
- Variants inside mapped structure/function features: {feature_hits_count}
- Variants with exact or region-level UniProt active/interface annotation: {active_hits}

## Feature Summary

{markdown_table(tables["benign_feature_summary"], ["feature", "count"], 20)}

## Gene Feature Summary

{markdown_table(tables["benign_gene_feature_summary"], ["gene", "feature", "count"], 30)}

## Curated Active Site Or Interface Status

{markdown_table(tables["benign_active_site_summary"], ["curated_active_site_status", "count"], 20)}

## Variant Types

{markdown_table(tables["benign_variant_type_summary"], ["variant_type", "count"], 20)}

## Criteria Combinations

{markdown_table(tables["benign_criteria_combo_summary"], ["criteria_combo", "count"], 20)}

## Benign Versus Pathogenic Feature Comparison

The percentages are within each group. This helps avoid reading raw counts as
enrichment when benign variants greatly outnumber pathogenic variants in the
generated snapshot.

{markdown_table(tables["benign_vs_pathogenic_feature_comparison"], ["feature", "benign_count", "pathogenic_count", "benign_percent", "pathogenic_percent", "pathogenic_to_benign_rate_ratio"], 30)}

## Benign Versus Pathogenic Active/Interface Comparison

{markdown_table(tables["benign_vs_pathogenic_active_site_comparison"], ["curated_active_site_status", "benign_count", "pathogenic_count", "benign_percent", "pathogenic_percent", "pathogenic_to_benign_rate_ratio"], 30)}

## Top Benign Variants In Mapped Features

{markdown_table([row for row in benign if row["structure_features"] != "outside_mapped_structure_function_features"], ["gene", "c_notation", "p_notation", "variant_type", "predicted_class", "structure_features", "curated_active_site_status", "criteria_combo"], 40)}

## Interpretation

- Benign/likely benign variants are expected to be numerous because BP1, BP4,
  BP7, BS3 and related benign-direction rules classify many generated SNVs.
- Structural location alone is not enough: benign variants can occur inside
  broad domains such as BRCT, RING, BRC/RAD51 context, or BRCA2 DBD.
- Exact residue-level functional annotations are more useful than broad domain
  membership, but they still require expert interpretation and should not be
  treated as an automatic criterion outside ACMG/ENIGMA rules.

## Outputs

- `tables/benign_structure_function/benign_structure_function_variants.csv`
- `tables/benign_structure_function/benign_feature_summary.csv`
- `tables/benign_structure_function/benign_gene_feature_summary.csv`
- `tables/benign_structure_function/benign_active_site_summary.csv`
- `tables/benign_structure_function/benign_vs_pathogenic_feature_comparison.csv`
- `tables/benign_structure_function/benign_vs_pathogenic_active_site_comparison.csv`
- `plots/19_benign_structure_function/benign_feature_summary.svg`
- `plots/19_benign_structure_function/benign_active_site_summary.svg`
- `plots/19_benign_structure_function/benign_vs_pathogenic_feature_comparison.svg`
- `plots/19_benign_structure_function/benign_vs_pathogenic_active_site_comparison.svg`
- `plots/19_benign_structure_function/brca1_benign_lollipop.svg`
- `plots/19_benign_structure_function/brca2_benign_lollipop.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    benign, tables = analyze()
    fields = [
        "gene",
        "c_notation",
        "p_notation",
        "aa_position",
        "variant_type",
        "predicted_class",
        "predicted_label",
        "total_points",
        "spliceai_score",
        "structure_features",
        "structure_domain_context",
        "curated_active_site_status",
        "known_pathogenic_residue_context",
        "criteria_combo",
        "criteria",
    ]
    write_csv(TABLE_DIR / "benign_structure_function_variants.csv", benign, fields)
    write_csv(TABLE_DIR / "benign_feature_summary.csv", tables["benign_feature_summary"], ["feature", "count"])
    write_csv(TABLE_DIR / "benign_gene_feature_summary.csv", tables["benign_gene_feature_summary"], ["gene", "feature", "count"])
    write_csv(TABLE_DIR / "benign_active_site_summary.csv", tables["benign_active_site_summary"], ["curated_active_site_status", "count"])
    write_csv(TABLE_DIR / "benign_variant_type_summary.csv", tables["benign_variant_type_summary"], ["variant_type", "count"])
    write_csv(TABLE_DIR / "benign_criteria_combo_summary.csv", tables["benign_criteria_combo_summary"], ["criteria_combo", "count"])
    write_csv(TABLE_DIR / "benign_known_pathogenic_residue_summary.csv", tables["benign_known_pathogenic_residue_summary"], ["known_pathogenic_residue_context", "count"])
    write_csv(
        TABLE_DIR / "benign_vs_pathogenic_feature_comparison.csv",
        tables["benign_vs_pathogenic_feature_comparison"],
        ["feature", "benign_count", "pathogenic_count", "benign_percent", "pathogenic_percent", "pathogenic_to_benign_rate_ratio"],
    )
    write_csv(
        TABLE_DIR / "benign_vs_pathogenic_active_site_comparison.csv",
        tables["benign_vs_pathogenic_active_site_comparison"],
        ["curated_active_site_status", "benign_count", "pathogenic_count", "benign_percent", "pathogenic_percent", "pathogenic_to_benign_rate_ratio"],
    )
    bar_svg(PLOT_DIR / "benign_feature_summary.svg", tables["benign_feature_summary"], "feature", "count", "Benign/Likely Benign Variants By Protein Feature")
    bar_svg(PLOT_DIR / "benign_active_site_summary.svg", tables["benign_active_site_summary"], "curated_active_site_status", "count", "Benign/Likely Benign Active Site Or Interface Status")
    comparison_svg(PLOT_DIR / "benign_vs_pathogenic_feature_comparison.svg", tables["benign_vs_pathogenic_feature_comparison"], "feature", "Benign vs Pathogenic Feature Rates")
    comparison_svg(PLOT_DIR / "benign_vs_pathogenic_active_site_comparison.svg", tables["benign_vs_pathogenic_active_site_comparison"], "curated_active_site_status", "Benign vs Pathogenic Active/Interface Rates")
    lollipop_svg(PLOT_DIR / "brca1_benign_lollipop.svg", benign, "BRCA1")
    lollipop_svg(PLOT_DIR / "brca2_benign_lollipop.svg", benign, "BRCA2")
    build_report(benign, tables)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
