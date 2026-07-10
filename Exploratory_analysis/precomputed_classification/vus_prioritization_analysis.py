"""Focused VUS prioritization summaries for manual review."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT_CSV = OUT_DIR / "tables" / "exon_vus_conflict" / "vus_priority_list.csv"
TABLE_DIR = OUT_DIR / "tables" / "vus_prioritization"
PLOT_DIR = OUT_DIR / "plots" / "08_vus_prioritization"
REPORT = OUT_DIR / "vus_prioritization_report.md"


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_svg(path: Path, body: str, width: int, height: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
text {{ font-family: Arial, Helvetica, sans-serif; fill: #172033; }}
.small {{ font-size: 11px; }}
.label {{ font-size: 12px; }}
.title {{ font-size: 18px; font-weight: 700; }}
.grid {{ stroke: #e2e8f0; stroke-width: 1; }}
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def load_rows() -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["priority_score_int"] = int(row["priority_score"])
            row["spliceai_float"] = float(row["spliceai_score"] or 0)
            row["reasons"] = [item for item in row["priority_reasons"].split(";") if item]
            row["priority_tier"] = priority_tier(row["priority_score_int"])
            row["review_category"] = review_category(row)
            rows.append(row)
    return sorted(rows, key=lambda item: (item["priority_score_int"], item["spliceai_float"]), reverse=True)


def priority_tier(score: int) -> str:
    if score >= 80:
        return "tier1_urgent"
    if score >= 50:
        return "tier2_high"
    if score >= 25:
        return "tier3_medium"
    return "tier4_low"


def review_category(row: dict[str, object]) -> str:
    reasons = set(row["reasons"])
    boundary = int(row["boundary_distance"]) if str(row["boundary_distance"]).isdigit() else 9999
    if "PVS1" in reasons:
        return "vus_with_pvs1"
    if "PS3" in reasons and "PP3" in reasons:
        return "functional_plus_splice_prediction"
    if "SpliceAI>=0.20" in reasons and boundary <= 2:
        return "splice_boundary"
    if "SpliceAI>=0.20" in reasons and boundary > 50:
        return "splice_far_from_boundary"
    if "near_pathogenic_dense" in reasons:
        return "pathogenic_neighborhood"
    if "SpliceAI>=0.20" in reasons:
        return "splice_other"
    if "PS1" in reasons or "PS3" in reasons:
        return "functional_or_same_aa"
    return "other_vus"


def top_by_category(rows: list[dict[str, object]], limit: int = 50) -> list[dict[str, object]]:
    output = []
    by_category: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_category[str(row["review_category"])].append(row)
    for category in sorted(by_category):
        output.extend(by_category[category][:limit])
    return output


def exon_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_exon: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_exon[(row["gene"], row["cds_exon"])].append(row)
    output = []
    for (gene, exon), items in sorted(by_exon.items(), key=lambda item: (item[0][0], int(item[0][1] or 0))):
        scores = [int(row["priority_score"]) for row in items]
        tier1 = sum(1 for row in items if row["priority_tier"] == "tier1_urgent")
        tier2 = sum(1 for row in items if row["priority_tier"] == "tier2_high")
        high_splice = sum(1 for row in items if row["spliceai_float"] >= 0.20)
        output.append(
            {
                "gene": gene,
                "cds_exon": exon,
                "vus_count": len(items),
                "tier1_urgent": tier1,
                "tier2_high": tier2,
                "max_priority_score": max(scores),
                "mean_priority_score": f"{sum(scores) / len(scores):.2f}",
                "high_spliceai_count": high_splice,
                "high_spliceai_percent": f"{high_splice / len(items) * 100:.2f}",
            }
        )
    return sorted(output, key=lambda item: (int(item["tier1_urgent"]), int(item["tier2_high"]), float(item["mean_priority_score"])), reverse=True)


def reason_counts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts = Counter()
    for row in rows:
        counts.update(row["reasons"])
    return [{"reason": reason, "count": count} for reason, count in counts.most_common()]


def category_counts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts = Counter(row["review_category"] for row in rows)
    return [{"review_category": category, "count": count} for category, count in counts.most_common()]


def tier_counts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts = Counter(row["priority_tier"] for row in rows)
    order = ["tier1_urgent", "tier2_high", "tier3_medium", "tier4_low"]
    return [{"priority_tier": tier, "count": counts[tier]} for tier in order]


def draw_score_histogram(rows: list[dict[str, object]]) -> None:
    bins = [(0, 24), (25, 49), (50, 64), (65, 79), (80, 100)]
    counts = []
    for start, end in bins:
        counts.append(sum(1 for row in rows if start <= row["priority_score_int"] <= end))
    width = 760
    height = 360
    left = 80
    top = 58
    plot_w = 590
    plot_h = 230
    max_count = max(counts) if counts else 1
    body = ['<text x="28" y="34" class="title">VUS priority score distribution</text>']
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    bar_w = plot_w / len(bins) - 18
    for i, ((start, end), count) in enumerate(zip(bins, counts)):
        h = plot_h * count / max_count if max_count else 0
        x = left + i * (plot_w / len(bins)) + 9
        y = top + plot_h - h
        color = "#7c3aed" if start >= 50 else "#94a3b8"
        body.append(f'<rect x="{x}" y="{y}" width="{bar_w}" height="{h}" fill="{color}"/>')
        body.append(f'<text x="{x + 4}" y="{y - 6}" class="small">{count}</text>')
        body.append(f'<text x="{x + 2}" y="{top + plot_h + 20}" class="small">{start}-{end}</text>')
    write_svg(PLOT_DIR / "vus_priority_score_histogram.svg", "\n".join(body), width, height)


def draw_reason_barplot(reason_rows: list[dict[str, object]]) -> None:
    rows = reason_rows[:12]
    width = 950
    height = 470
    left = 250
    top = 60
    bar_h = 24
    gap = 9
    plot_w = 590
    max_count = max((int(row["count"]) for row in rows), default=1)
    body = ['<text x="28" y="34" class="title">Most common VUS priority reasons</text>']
    for i, row in enumerate(rows):
        y = top + i * (bar_h + gap)
        count = int(row["count"])
        w = plot_w * count / max_count
        body.append(f'<text x="28" y="{y + 17}" class="small">{esc(row["reason"])}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="#2563eb"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 17}" class="small">{count}</text>')
    write_svg(PLOT_DIR / "vus_priority_reason_counts.svg", "\n".join(body), width, height)


def draw_exon_barplot(exon_rows: list[dict[str, object]]) -> None:
    rows = exon_rows[:20]
    width = 980
    height = 650
    left = 270
    top = 60
    bar_h = 21
    gap = 8
    plot_w = 600
    max_count = max((int(row["tier1_urgent"]) + int(row["tier2_high"]) for row in rows), default=1)
    body = ['<text x="28" y="34" class="title">Exons enriched for high-priority VUS</text>']
    for i, row in enumerate(rows):
        y = top + i * (bar_h + gap)
        high = int(row["tier1_urgent"]) + int(row["tier2_high"])
        w = plot_w * high / max_count
        label = f'{row["gene"]} exon {row["cds_exon"]}'
        body.append(f'<text x="28" y="{y + 15}" class="small">{esc(label)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="#dc2626"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 15}" class="small">high={high}, all VUS={row["vus_count"]}</text>')
    write_svg(PLOT_DIR / "vus_priority_exon_barplot.svg", "\n".join(body), width, height)


def report_table(rows: list[dict[str, object]], columns: list[str], limit: int = 12) -> str:
    if not rows:
        return "| none |\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows[:limit]:
        body.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join([header, sep] + body)


def write_report(rows: list[dict[str, object]], reasons: list[dict[str, object]], categories: list[dict[str, object]], exons: list[dict[str, object]]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    tier_summary = tier_counts(rows)
    text = f"""# VUS Prioritization Report

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

This is a focused manual-review layer over the automated Module 1 VUS set. It does not change classifications.

## What This Is For

This report helps prioritize VUS variants for manual curator review. It is not a classifier and it does not upgrade or downgrade any variant. The goal is to turn thousands of VUS variants into a smaller review queue by highlighting variants with stronger biological or rule-based reasons to look again.

The most useful output is `vus_priority_annotated.csv`, where every VUS gets a priority score, a tier, a review category, and explicit reasons. The `tier1_urgent` and `tier2_high` variants are the practical first-pass review list.

`Pathogenic density` means that many variants close to the VUS in coding-sequence position are already classified as likely pathogenic/pathogenic by the automated Module 1 snapshot. In this analysis, `near_pathogenic_dense` is assigned when at least 5 pathogenic-group variants occur within +/-20 coding bases of the VUS. This does not mean the VUS is pathogenic by itself. It means the local region has many pathogenic signals and deserves earlier curator attention.

Typical high-priority reasons include:

- strong or moderate SpliceAI signal
- proximity to a CDS exon boundary
- PP3, PS3, PS1, or PVS1 evidence still ending as VUS
- location in a dense pathogenic neighborhood
- high SpliceAI in a region otherwise rich in benign evidence, which is useful as a sanity check

## How This Was Derived

The starting point was the full precomputed BRCA1/BRCA2 coding SNV Module 1 snapshot. First, only variants classified as class 3 VUS by the automated snapshot were selected. Then each VUS was annotated with its CDS position, CDS exon, distance to the nearest CDS exon boundary, SpliceAI score, applied ACMG/ENIGMA Module 1 criteria, and local neighborhood context.

The priority score is a heuristic review score. Points are added for signals that make a VUS more useful to inspect manually:

- SpliceAI >= 0.20 adds strong splice-priority weight
- SpliceAI >= 0.10 adds moderate splice-priority weight
- distance <= 2 bp or <= 10 bp from a CDS exon boundary adds splice-site proximity weight
- PVS1, PS3, PS1, PP3, PP4, or BS3 add rule/evidence-based review weight
- `near_pathogenic_dense` adds weight when the local +/-20 coding-base neighborhood contains at least 5 pathogenic-group variants
- `high_spliceai_in_benign_dense_region` adds weight as a sanity-check flag when high SpliceAI appears in a locally benign-rich region

The tiers are then assigned from the total score: `tier1_urgent` for score >= 80, `tier2_high` for score >= 50, `tier3_medium` for score >= 25, and `tier4_low` below 25. Review categories are assigned from the dominant reason pattern, for example splice-boundary, splice-far-from-boundary, functional-plus-splice-prediction, or pathogenic-neighborhood.

This is intentionally conservative as a triage tool. It does not introduce new ACMG/ENIGMA criteria and it does not claim that a high-priority VUS should be reclassified. It only identifies variants that are more informative to open first.

## Outputs

- `tables/vus_prioritization/vus_priority_annotated.csv`
- `tables/vus_prioritization/vus_priority_top_by_category.csv`
- `tables/vus_prioritization/vus_priority_reason_counts.csv`
- `tables/vus_prioritization/vus_priority_category_counts.csv`
- `tables/vus_prioritization/vus_priority_exon_summary.csv`
- `plots/08_vus_prioritization/vus_priority_score_histogram.svg`
- `plots/08_vus_prioritization/vus_priority_reason_counts.svg`
- `plots/08_vus_prioritization/vus_priority_exon_barplot.svg`

## Tier Counts

{report_table(tier_summary, ["priority_tier", "count"], limit=10)}

## Review Categories

{report_table(categories, ["review_category", "count"], limit=20)}

## Most Common Reasons

{report_table(reasons, ["reason", "count"], limit=20)}

## Top Exons for High-priority VUS

{report_table(exons, ["gene", "cds_exon", "vus_count", "tier1_urgent", "tier2_high", "max_priority_score", "mean_priority_score", "high_spliceai_percent"], limit=15)}

## Top VUS Variants

{report_table(rows, ["priority_score", "priority_tier", "review_category", "priority_reasons", "gene", "c_notation", "p_notation", "cds_exon", "boundary_distance", "spliceai_score"], limit=20)}
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    reasons = reason_counts(rows)
    categories = category_counts(rows)
    exons = exon_summary(rows)
    fields = [
        "priority_score",
        "priority_tier",
        "review_category",
        "priority_reasons",
        "gene",
        "c_notation",
        "p_notation",
        "cds_pos",
        "cds_exon",
        "boundary_distance",
        "variant_type",
        "spliceai_score",
        "criteria",
        "total_points",
        "predicted_class",
        "near_pathogenic_20bp",
        "near_benign_20bp",
        "near_vus_20bp",
    ]
    write_csv(TABLE_DIR / "vus_priority_annotated.csv", rows, fields)
    write_csv(TABLE_DIR / "vus_priority_top_by_category.csv", top_by_category(rows), fields)
    write_csv(TABLE_DIR / "vus_priority_reason_counts.csv", reasons, ["reason", "count"])
    write_csv(TABLE_DIR / "vus_priority_category_counts.csv", categories, ["review_category", "count"])
    write_csv(
        TABLE_DIR / "vus_priority_exon_summary.csv",
        exons,
        [
            "gene",
            "cds_exon",
            "vus_count",
            "tier1_urgent",
            "tier2_high",
            "max_priority_score",
            "mean_priority_score",
            "high_spliceai_count",
            "high_spliceai_percent",
        ],
    )
    draw_score_histogram(rows)
    draw_reason_barplot(reasons)
    draw_exon_barplot(exons)
    write_report(rows, reasons, categories, exons)
    print(f"Wrote VUS prioritization report to {REPORT}")


if __name__ == "__main__":
    main()
