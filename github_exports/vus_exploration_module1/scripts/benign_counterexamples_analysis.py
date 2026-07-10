"""Benign counterexamples inside structurally important or pathogenic-rich regions."""

from __future__ import annotations

import csv
import html
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
BENIGN_INPUT = ANALYSIS_DIR / "tables" / "benign_structure_function" / "benign_structure_function_variants.csv"
FEATURE_COMPARISON = ANALYSIS_DIR / "tables" / "benign_structure_function" / "benign_vs_pathogenic_feature_comparison.csv"
TABLE_DIR = ANALYSIS_DIR / "tables" / "benign_counterexamples"
PLOT_DIR = ANALYSIS_DIR / "plots" / "20_benign_counterexamples"
REPORT = ANALYSIS_DIR / "benign_counterexamples_report.md"


IMPORTANT_FEATURES = {
    "RING zinc-binding/E3 ligase region",
    "BRCT phosphopeptide-binding region",
    "DNA-binding domain / helical-OB-DSS1 region",
    "PALB2 binding N-terminal region",
    "coiled-coil PALB2 interaction region",
    "BRC repeats / RAD51-binding region",
    "C-terminal RAD51-binding/nuclear localization context",
}

TIER_RANK = {
    "A_exact_or_known_residue_context": 4,
    "B_pathogenic_enriched_domain": 3,
    "C_region_level_interface_context": 2,
    "D_broad_domain_context": 1,
}


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def feature_set(row: dict) -> set[str]:
    return set(str(row["structure_features"]).split(";"))


def is_counterexample(row: dict) -> bool:
    features = feature_set(row)
    in_important_feature = bool(features & IMPORTANT_FEATURES)
    has_active_context = row["curated_active_site_status"] != "not_curated_in_current_dataset"
    known_pathogenic_residue = bool(row["known_pathogenic_residue_context"])
    return in_important_feature or has_active_context or known_pathogenic_residue


def counterexample_tier(row: dict, pathogenic_rate_ratios: dict[str, float]) -> str:
    features = feature_set(row)
    ratios = [pathogenic_rate_ratios.get(feature, 0.0) for feature in features]
    max_ratio = max(ratios) if ratios else 0.0
    exact_or_source_active = (
        row["curated_active_site_status"] != "not_curated_in_current_dataset"
        and not row["curated_active_site_status"].startswith("region_level:")
    )
    if exact_or_source_active or bool(row["known_pathogenic_residue_context"]):
        return "A_exact_or_known_residue_context"
    if max_ratio >= 10:
        return "B_pathogenic_enriched_domain"
    if row["curated_active_site_status"].startswith("region_level:"):
        return "C_region_level_interface_context"
    return "D_broad_domain_context"


def score_row(row: dict, pathogenic_rate_ratios: dict[str, float]) -> int:
    features = feature_set(row)
    max_ratio = max((pathogenic_rate_ratios.get(feature, 0.0) for feature in features), default=0.0)
    score = int(max_ratio * 10)
    if row["curated_active_site_status"] != "not_curated_in_current_dataset":
        score += 25
    if row["curated_active_site_status"].startswith("region_level:"):
        score -= 8
    if row["known_pathogenic_residue_context"]:
        score += 40
    if "BS3" in row["criteria_combo"]:
        score += 10
    if row["predicted_class"] == "1":
        score += 5
    return score


def analyze() -> tuple[list[dict], dict[str, list[dict]]]:
    benign_rows = read_csv(BENIGN_INPUT)
    feature_comparison = read_csv(FEATURE_COMPARISON)
    pathogenic_rate_ratios = {}
    for row in feature_comparison:
        value = row.get("pathogenic_to_benign_rate_ratio", "0")
        try:
            pathogenic_rate_ratios[row["feature"]] = float(value)
        except ValueError:
            pathogenic_rate_ratios[row["feature"]] = 0.0

    rows = []
    for row in benign_rows:
        if not is_counterexample(row):
            continue
        out = dict(row)
        out["counterexample_tier"] = counterexample_tier(row, pathogenic_rate_ratios)
        out["counterexample_score"] = score_row(row, pathogenic_rate_ratios)
        out["max_pathogenic_to_benign_feature_ratio"] = f"{max((pathogenic_rate_ratios.get(feature, 0.0) for feature in feature_set(row)), default=0.0):.3f}"
        rows.append(out)

    rows.sort(
        key=lambda row: (
            TIER_RANK.get(row["counterexample_tier"], 0),
            int(row["counterexample_score"]),
            float(row["max_pathogenic_to_benign_feature_ratio"]),
        ),
        reverse=True,
    )
    tables = {
        "tier_summary": counter_rows(rows, "counterexample_tier", "counterexample_tier"),
        "feature_summary": feature_summary(rows),
        "active_site_summary": counter_rows(rows, "curated_active_site_status", "curated_active_site_status"),
        "criteria_summary": counter_rows(rows, "criteria_combo", "criteria_combo"),
        "gene_feature_summary": gene_feature_summary(rows),
    }
    return rows, tables


def counter_rows(rows: list[dict], key: str, output_key: str) -> list[dict]:
    counts = Counter(row[key] for row in rows)
    return [{output_key: key, "count": value} for key, value in counts.most_common()]


def feature_summary(rows: list[dict]) -> list[dict]:
    counts = Counter()
    for row in rows:
        for feature in feature_set(row):
            if feature != "outside_mapped_structure_function_features":
                counts[feature] += 1
    return [{"feature": key, "count": value} for key, value in counts.most_common()]


def gene_feature_summary(rows: list[dict]) -> list[dict]:
    counts = Counter()
    for row in rows:
        for feature in feature_set(row):
            if feature != "outside_mapped_structure_function_features":
                counts[(row["gene"], feature)] += 1
    return [
        {"gene": key[0], "feature": key[1], "count": value}
        for key, value in counts.most_common()
    ]


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def bar_svg(path: Path, rows: list[dict], label_key: str, value_key: str, title: str, max_items: int = 18) -> None:
    rows = rows[:max_items]
    width = 1120
    left = 500
    top = 54
    row_h = 31
    chart_w = 520
    height = top + len(rows) * row_h + 58
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
        parts.append(f'<text x="24" y="{y + 20}">{html.escape(str(row[label_key])[:72])}</text>')
        parts.append(f'<rect x="{left}" y="{y + 6}" width="{bar_w}" height="19" class="bar"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 21}" class="small">{value}</text>')
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


def write_report(rows: list[dict], tables: dict[str, list[dict]]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# Benign Counterexamples In Pathogenic Domains

Generated: {generated}

Input: `{BENIGN_INPUT.relative_to(ROOT)}`

## Purpose

This analysis finds benign/likely benign generated variants that fall inside
structurally or functionally important protein regions. These are
counterexamples to the tempting but unsafe shortcut: "variant is in a domain,
therefore it is pathogenic".

They are not errors by default. They are useful reminders that region/domain
context is not an ACMG/ENIGMA criterion by itself and must be interpreted with
the actual evidence, such as BS3, BP1, BP4, BP7, PS3, PP3, PS1, PVS1, or PM5.

## Counterexample Tiers

- `A_exact_or_known_residue_context`: benign variant overlaps a source-extracted
  functional annotation or a known pathogenic-residue context
- `B_pathogenic_enriched_domain`: benign variant lies in a domain that is much
  more frequent in pathogenic than benign variants after group normalization
- `C_region_level_interface_context`: benign variant lies in a UniProt
  region-level interaction annotation
- `D_broad_domain_context`: benign variant lies in a broad important domain

## Summary

Total benign counterexamples: {len(rows)}

{markdown_table(tables["tier_summary"], ["counterexample_tier", "count"], 20)}

## Feature Summary

{markdown_table(tables["feature_summary"], ["feature", "count"], 20)}

## Gene Feature Summary

{markdown_table(tables["gene_feature_summary"], ["gene", "feature", "count"], 30)}

## Active/Interface Summary

{markdown_table(tables["active_site_summary"], ["curated_active_site_status", "count"], 20)}

## Criteria Summary

{markdown_table(tables["criteria_summary"], ["criteria_combo", "count"], 20)}

## Top Counterexamples

{markdown_table(rows, ["counterexample_score", "counterexample_tier", "gene", "c_notation", "p_notation", "variant_type", "predicted_class", "structure_features", "curated_active_site_status", "criteria_combo"], 40)}

## Interpretation

The most important use of this table is negative control. If a future narrative
or model says that all variants in BRCT, RING, BRCA2 DBD, PALB2, RAD51/BRC, or
SEM1/DSS1 contexts should be treated as suspicious, these benign counterexamples
show why that is too simple.

The counterexamples should be read together with the pathogenic-region and VUS
prioritization analyses. A VUS in a pathogenic region becomes more interesting
when its local neighborhood, criteria, and residue-level evidence point in the
same direction. Domain membership alone is not enough.

## Outputs

- `tables/benign_counterexamples/benign_counterexamples.csv`
- `tables/benign_counterexamples/tier_summary.csv`
- `tables/benign_counterexamples/feature_summary.csv`
- `tables/benign_counterexamples/gene_feature_summary.csv`
- `tables/benign_counterexamples/active_site_summary.csv`
- `tables/benign_counterexamples/criteria_summary.csv`
- `plots/20_benign_counterexamples/tier_summary.svg`
- `plots/20_benign_counterexamples/feature_summary.svg`
- `plots/20_benign_counterexamples/active_site_summary.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rows, tables = analyze()
    fields = [
        "counterexample_score",
        "counterexample_tier",
        "max_pathogenic_to_benign_feature_ratio",
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
    write_csv(TABLE_DIR / "benign_counterexamples.csv", rows, fields)
    write_csv(TABLE_DIR / "tier_summary.csv", tables["tier_summary"], ["counterexample_tier", "count"])
    write_csv(TABLE_DIR / "feature_summary.csv", tables["feature_summary"], ["feature", "count"])
    write_csv(TABLE_DIR / "gene_feature_summary.csv", tables["gene_feature_summary"], ["gene", "feature", "count"])
    write_csv(TABLE_DIR / "active_site_summary.csv", tables["active_site_summary"], ["curated_active_site_status", "count"])
    write_csv(TABLE_DIR / "criteria_summary.csv", tables["criteria_summary"], ["criteria_combo", "count"])
    bar_svg(PLOT_DIR / "tier_summary.svg", tables["tier_summary"], "counterexample_tier", "count", "Benign Counterexample Tiers")
    bar_svg(PLOT_DIR / "feature_summary.svg", tables["feature_summary"], "feature", "count", "Benign Counterexamples By Protein Feature")
    bar_svg(PLOT_DIR / "active_site_summary.svg", tables["active_site_summary"], "curated_active_site_status", "count", "Benign Counterexamples By Active/Interface Annotation")
    write_report(rows, tables)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
