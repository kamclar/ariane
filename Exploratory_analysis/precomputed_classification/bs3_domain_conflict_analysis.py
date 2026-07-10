"""BS3 and benign evidence inside functional protein domains."""

from __future__ import annotations

import csv
import html
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ANALYSIS_DIR / "tables" / "benign_structure_function" / "benign_structure_function_variants.csv"
TABLE_DIR = ANALYSIS_DIR / "tables" / "bs3_domain_conflicts"
PLOT_DIR = ANALYSIS_DIR / "plots" / "21_bs3_domain_conflicts"
REPORT = ANALYSIS_DIR / "bs3_domain_conflict_report.md"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_float(value: str) -> float:
    try:
        return float(value or 0)
    except ValueError:
        return 0.0


def criteria_codes(criteria_combo: str) -> set[str]:
    return {item for item in criteria_combo.split("+") if item and item != "none"}


def feature_values(row: dict) -> list[str]:
    return [item for item in row["structure_features"].split(";") if item]


def has_domain_or_interface(row: dict) -> bool:
    return (
        row["structure_features"] != "outside_mapped_structure_function_features"
        or row["curated_active_site_status"] != "not_curated_in_current_dataset"
        or bool(row["known_pathogenic_residue_context"])
    )


def has_exact_or_source_active(row: dict) -> bool:
    status = row["curated_active_site_status"]
    return status != "not_curated_in_current_dataset" and not status.startswith("region_level:")


def conflict_flags(row: dict) -> list[str]:
    codes = criteria_codes(row["criteria_combo"])
    spliceai = parse_float(row["spliceai_score"])
    flags = []
    if "BS3" in codes:
        flags.append("BS3")
    if "PP3" in codes:
        flags.append("PP3")
    if spliceai >= 0.20:
        flags.append("SpliceAI_ge_0_20")
    if "PS3" in codes:
        flags.append("PS3")
    if "PVS1" in codes or "PM5_PTC" in codes:
        flags.append("truncating_pathogenic_criterion")
    if row["structure_features"] != "outside_mapped_structure_function_features":
        flags.append("functional_domain")
    if row["curated_active_site_status"] != "not_curated_in_current_dataset":
        flags.append("active_or_interface_context")
    if has_exact_or_source_active(row):
        flags.append("exact_or_source_active_context")
    return flags


def review_bucket(row: dict) -> str:
    flags = set(conflict_flags(row))
    if "BS3" not in flags:
        return "non_BS3_benign_domain_context"
    if "PP3" in flags and "SpliceAI_ge_0_20" in flags:
        return "BS3_with_PP3_and_high_SpliceAI"
    if "PP3" in flags:
        return "BS3_with_PP3"
    if "SpliceAI_ge_0_20" in flags:
        return "BS3_with_high_SpliceAI"
    if "exact_or_source_active_context" in flags:
        return "BS3_in_exact_or_source_active_context"
    if "active_or_interface_context" in flags:
        return "BS3_in_region_level_interface_context"
    return "BS3_in_broad_domain_context"


def analyze() -> tuple[list[dict], dict[str, list[dict]]]:
    rows = []
    for row in read_csv(INPUT):
        if not has_domain_or_interface(row):
            continue
        out = dict(row)
        out["spliceai_float"] = parse_float(row["spliceai_score"])
        out["has_BS3"] = "yes" if "BS3" in criteria_codes(row["criteria_combo"]) else "no"
        out["has_PP3"] = "yes" if "PP3" in criteria_codes(row["criteria_combo"]) else "no"
        out["high_spliceai"] = "yes" if out["spliceai_float"] >= 0.20 else "no"
        out["variant_mechanism_group"] = mechanism_group(row)
        out["conflict_flags"] = ";".join(conflict_flags(row))
        out["review_bucket"] = review_bucket(row)
        out["bs3_domain_conflict_score"] = score(out)
        rows.append(out)

    rows.sort(
        key=lambda row: (
            row["has_BS3"] == "yes",
            row["review_bucket"] in {"BS3_with_PP3_and_high_SpliceAI", "BS3_with_PP3", "BS3_with_high_SpliceAI"},
            int(row["bs3_domain_conflict_score"]),
            float(row["spliceai_float"]),
        ),
        reverse=True,
    )
    tables = {
        "review_bucket_summary": counter_rows(rows, "review_bucket", "review_bucket"),
        "feature_summary": feature_summary([row for row in rows if row["has_BS3"] == "yes"]),
        "variant_mechanism_summary": counter_rows([row for row in rows if row["has_BS3"] == "yes"], "variant_mechanism_group", "variant_mechanism_group"),
        "criteria_summary": counter_rows([row for row in rows if row["has_BS3"] == "yes"], "criteria_combo", "criteria_combo"),
        "domain_conflict_summary": domain_conflict_summary(rows),
    }
    return rows, tables


def mechanism_group(row: dict) -> str:
    variant_type = row["variant_type"]
    if variant_type == "missense":
        return "missense"
    if variant_type == "synonymous":
        return "synonymous"
    if variant_type in {"splice_site", "splice_region"}:
        return "splice_region"
    if variant_type in {"nonsense", "frameshift"}:
        return "truncating"
    return variant_type or "other"


def score(row: dict) -> int:
    value = 0
    if row["has_BS3"] == "yes":
        value += 20
    if row["has_PP3"] == "yes":
        value += 20
    if row["high_spliceai"] == "yes":
        value += 20
    if row["curated_active_site_status"] != "not_curated_in_current_dataset":
        value += 15
    if has_exact_or_source_active(row):
        value += 15
    if row["known_pathogenic_residue_context"]:
        value += 15
    if row["variant_mechanism_group"] == "missense":
        value += 10
    if row["variant_mechanism_group"] == "synonymous":
        value -= 5
    return value


def counter_rows(rows: list[dict], key: str, output_key: str) -> list[dict]:
    counts = Counter(row[key] for row in rows)
    return [{output_key: key, "count": value} for key, value in counts.most_common()]


def feature_summary(rows: list[dict]) -> list[dict]:
    counts = Counter()
    for row in rows:
        for feature in feature_values(row):
            if feature != "outside_mapped_structure_function_features":
                counts[feature] += 1
    return [{"feature": key, "count": value} for key, value in counts.most_common()]


def domain_conflict_summary(rows: list[dict]) -> list[dict]:
    output = []
    features = sorted(
        {
            feature
            for row in rows
            for feature in feature_values(row)
            if feature != "outside_mapped_structure_function_features"
        }
    )
    for feature in features:
        feature_rows = [row for row in rows if feature in feature_values(row)]
        bs3_rows = [row for row in feature_rows if row["has_BS3"] == "yes"]
        output.append(
            {
                "feature": feature,
                "benign_domain_count": len(feature_rows),
                "bs3_count": len(bs3_rows),
                "bs3_percent": f"{100 * len(bs3_rows) / len(feature_rows):.2f}" if feature_rows else "0.00",
                "bs3_pp3_count": sum(1 for row in bs3_rows if row["has_PP3"] == "yes"),
                "bs3_high_spliceai_count": sum(1 for row in bs3_rows if row["high_spliceai"] == "yes"),
                "bs3_missense_count": sum(1 for row in bs3_rows if row["variant_mechanism_group"] == "missense"),
                "bs3_synonymous_count": sum(1 for row in bs3_rows if row["variant_mechanism_group"] == "synonymous"),
            }
        )
    return sorted(output, key=lambda row: int(row["bs3_count"]), reverse=True)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def bar_svg(path: Path, rows: list[dict], label_key: str, value_key: str, title: str, max_items: int = 18) -> None:
    rows = rows[:max_items]
    width = 1120
    left = 450
    top = 54
    row_h = 31
    chart_w = 560
    height = top + len(rows) * row_h + 60
    max_value = max((int(row[value_key]) for row in rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:13px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".small{font-size:11px;fill:#555}",
        ".bar{fill:#b45309}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = int(row[value_key])
        bar_w = int(chart_w * value / max_value) if max_value else 0
        parts.append(f'<text x="24" y="{y + 20}">{html.escape(str(row[label_key])[:64])}</text>')
        parts.append(f'<rect x="{left}" y="{y + 6}" width="{bar_w}" height="19" class="bar"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 21}" class="small">{value}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def stacked_domain_svg(path: Path, rows: list[dict]) -> None:
    rows = rows[:12]
    width = 1160
    left = 410
    top = 58
    row_h = 34
    chart_w = 620
    height = top + len(rows) * row_h + 65
    max_count = max((int(row["bs3_count"]) for row in rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:12px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".missense{fill:#b45309}",
        ".syn{fill:#64748b}",
        ".other{fill:#94a3b8}",
        ".small{font-size:10px;fill:#555}",
        "</style>",
        '<text x="24" y="30" class="title">BS3 benign variants by domain and mechanism</text>',
        f'<text x="{left}" y="48" class="small">brown: missense, gray: synonymous, pale: other</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        missense = int(row["bs3_missense_count"])
        syn = int(row["bs3_synonymous_count"])
        total = int(row["bs3_count"])
        other = max(0, total - missense - syn)
        x = left
        parts.append(f'<text x="24" y="{y + 19}">{html.escape(row["feature"][:58])}</text>')
        for klass, value in [("missense", missense), ("syn", syn), ("other", other)]:
            w = chart_w * value / max_count if max_count else 0
            parts.append(f'<rect x="{x}" y="{y + 5}" width="{w:.1f}" height="20" class="{klass}"/>')
            x += w
        parts.append(f'<text x="{x + 8}" y="{y + 20}" class="small">{total}</text>')
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
    bs3_rows = [row for row in rows if row["has_BS3"] == "yes"]
    conflict_rows = [
        row for row in bs3_rows if row["has_PP3"] == "yes" or row["high_spliceai"] == "yes"
    ]
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# BS3 In Functional Domains

Generated: {generated}

Input: `{INPUT.relative_to(ROOT)}`

## Purpose

This analysis focuses on benign/likely benign variants in functional domains or
active/interface contexts where `BS3` is present. It asks when functional
benign evidence is acting as the main counterweight to domain or interface
context, and where conflicts such as `BS3+PP3` or high SpliceAI appear.

This is a review aid only. It does not change classifications and it does not
create new ACMG/ENIGMA criteria.

## Summary

- Benign domain/interface variants analyzed: {len(rows)}
- With `BS3`: {len(bs3_rows)}
- With `BS3` plus `PP3` or SpliceAI >= 0.20: {len(conflict_rows)}

## Review Buckets

{markdown_table(tables["review_bucket_summary"], ["review_bucket", "count"], 20)}

## BS3 By Domain

{markdown_table(tables["domain_conflict_summary"], ["feature", "benign_domain_count", "bs3_count", "bs3_percent", "bs3_pp3_count", "bs3_high_spliceai_count", "bs3_missense_count", "bs3_synonymous_count"], 20)}

## BS3 Variant Mechanisms

{markdown_table(tables["variant_mechanism_summary"], ["variant_mechanism_group", "count"], 20)}

## BS3 Criteria Combinations

{markdown_table(tables["criteria_summary"], ["criteria_combo", "count"], 20)}

## Top BS3 Domain Conflicts

{markdown_table([row for row in rows if row["has_BS3"] == "yes"], ["bs3_domain_conflict_score", "review_bucket", "gene", "c_notation", "p_notation", "variant_type", "spliceai_score", "structure_features", "curated_active_site_status", "criteria_combo"], 40)}

## Interpretation

The main question is not whether domain context exists, but whether the BS3
assay is measuring the relevant mechanism for that region. A benign functional
result is strongest when the assay captures the disease-relevant mechanism for
that domain. It is more ambiguous when the same variant also has a PP3/splice
signal or when the assay may not detect splicing, NMD, or another mechanism.

Missense and synonymous counterexamples should be interpreted separately.
Missense variants are more directly relevant to protein-function hypotheses.
Synonymous variants in protein domains are often more about splice/RNA context
than protein active-site biology.

## Outputs

- `tables/bs3_domain_conflicts/bs3_domain_conflict_variants.csv`
- `tables/bs3_domain_conflicts/review_bucket_summary.csv`
- `tables/bs3_domain_conflicts/domain_conflict_summary.csv`
- `tables/bs3_domain_conflicts/variant_mechanism_summary.csv`
- `tables/bs3_domain_conflicts/criteria_summary.csv`
- `plots/21_bs3_domain_conflicts/review_bucket_summary.svg`
- `plots/21_bs3_domain_conflicts/domain_conflict_summary.svg`
- `plots/21_bs3_domain_conflicts/domain_mechanism_stacked.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rows, tables = analyze()
    fields = [
        "bs3_domain_conflict_score",
        "review_bucket",
        "conflict_flags",
        "has_BS3",
        "has_PP3",
        "high_spliceai",
        "variant_mechanism_group",
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
    write_csv(TABLE_DIR / "bs3_domain_conflict_variants.csv", rows, fields)
    write_csv(TABLE_DIR / "review_bucket_summary.csv", tables["review_bucket_summary"], ["review_bucket", "count"])
    write_csv(TABLE_DIR / "feature_summary.csv", tables["feature_summary"], ["feature", "count"])
    write_csv(TABLE_DIR / "domain_conflict_summary.csv", tables["domain_conflict_summary"], ["feature", "benign_domain_count", "bs3_count", "bs3_percent", "bs3_pp3_count", "bs3_high_spliceai_count", "bs3_missense_count", "bs3_synonymous_count"])
    write_csv(TABLE_DIR / "variant_mechanism_summary.csv", tables["variant_mechanism_summary"], ["variant_mechanism_group", "count"])
    write_csv(TABLE_DIR / "criteria_summary.csv", tables["criteria_summary"], ["criteria_combo", "count"])
    bar_svg(PLOT_DIR / "review_bucket_summary.svg", tables["review_bucket_summary"], "review_bucket", "count", "Benign domain variants by BS3 review bucket")
    bar_svg(PLOT_DIR / "domain_conflict_summary.svg", tables["domain_conflict_summary"], "feature", "bs3_count", "BS3 benign variants by domain")
    stacked_domain_svg(PLOT_DIR / "domain_mechanism_stacked.svg", tables["domain_conflict_summary"])
    write_report(rows, tables)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
