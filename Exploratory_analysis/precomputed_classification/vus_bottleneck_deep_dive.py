"""Deep dive for the three main VUS bottleneck groups."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
BOTTLENECK_CSV = OUT_DIR / "tables" / "vus_bottleneck" / "vus_bottleneck_variant_examples.csv"
PRIORITY_CSV = OUT_DIR / "tables" / "vus_prioritization" / "vus_priority_annotated.csv"
TABLE_DIR = OUT_DIR / "tables" / "vus_bottleneck_deep_dive"
PLOT_DIR = OUT_DIR / "plots" / "13_vus_bottleneck_deep_dive"
REPORT = OUT_DIR / "vus_bottleneck_deep_dive_report.md"

TARGET_CATEGORIES = {
    "strong_pathogenic_combo_one_step_short",
    "PM2_only",
    "mixed_splice_benign_and_pathogenic",
}


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def as_float(value: object, default: float = 0.0) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def key(row: dict[str, str]) -> tuple[str, str]:
    return (row["gene"], row["c_notation"])


def load_rows() -> list[dict[str, object]]:
    priority_by_key = {key(row): row for row in read_csv(PRIORITY_CSV)}
    rows = []
    for row in read_csv(BOTTLENECK_CSV):
        if row["bottleneck_category"] not in TARGET_CATEGORIES:
            continue
        priority = priority_by_key.get(key(row), {})
        merged: dict[str, object] = dict(row)
        for name in [
            "priority_score",
            "priority_tier",
            "review_category",
            "priority_reasons",
            "cds_pos",
            "cds_exon",
            "boundary_distance",
            "near_pathogenic_20bp",
            "near_benign_20bp",
            "near_vus_20bp",
        ]:
            merged[name] = priority.get(name, "")
        merged["priority_score_int"] = as_int(merged["priority_score"])
        merged["spliceai_float"] = as_float(merged["spliceai_score"])
        merged["near_pathogenic_int"] = as_int(merged["near_pathogenic_20bp"])
        merged["boundary_distance_int"] = as_int(merged["boundary_distance"], 9999)
        merged["deep_bucket"] = deep_bucket(merged)
        merged["what_it_means"] = what_it_means(str(merged["deep_bucket"]))
        merged["practical_next_step"] = practical_next_step(str(merged["deep_bucket"]))
        rows.append(merged)
    return rows


def deep_bucket(row: dict[str, object]) -> str:
    category = str(row["bottleneck_category"])
    combo = str(row["criteria_combo"])
    priority_score = int(row["priority_score_int"])
    near_path = int(row["near_pathogenic_int"])
    boundary = int(row["boundary_distance_int"])
    spliceai = float(row["spliceai_float"])

    if category == "strong_pathogenic_combo_one_step_short":
        if combo == "PM2_Supporting+PS3":
            return "PS3_PM2_needs_one_more_pathogenic_support"
        if combo == "PS3":
            return "PS3_only_needs_independent_support"
        if "PP3" in combo and "PS3" in combo:
            return "PS3_PP3_splice_or_functional_followup"
        return "complex_strong_pathogenic_combo_review"

    if category == "PM2_only":
        if priority_score >= 50:
            return "PM2_only_high_priority_context"
        if near_path >= 5:
            return "PM2_only_pathogenic_neighborhood"
        if spliceai >= 0.20 or boundary <= 2:
            return "PM2_only_splice_or_boundary_context"
        return "PM2_only_low_information_background"

    if category == "mixed_splice_benign_and_pathogenic":
        if near_path >= 5:
            return "BP4_BP7_PM2_near_pathogenic_density"
        if boundary <= 2:
            return "BP4_BP7_PM2_boundary_splice_benign"
        return "BP4_BP7_PM2_low_splice_benign_prediction"

    return "other"


def what_it_means(bucket: str) -> str:
    meanings = {
        "PS3_PM2_needs_one_more_pathogenic_support": "Strong functional pathogenic evidence plus PM2 is close, but still lacks the exact ENIGMA combination.",
        "PS3_only_needs_independent_support": "Functional pathogenic evidence alone is not enough for likely pathogenic.",
        "PS3_PP3_splice_or_functional_followup": "Functional evidence and computational splice/pathogenic signal point together, but still need another independent evidence item.",
        "complex_strong_pathogenic_combo_review": "Several pathogenic signals are present, but the exact combination still misses the rule threshold.",
        "PM2_only_high_priority_context": "PM2 is the only applied criterion, but local context makes the variant worth prioritizing.",
        "PM2_only_pathogenic_neighborhood": "The variant is near a dense pathogenic region, but has no direct automated evidence beyond PM2; in full-SNV scans this is a broad context signal, not a decision signal.",
        "PM2_only_splice_or_boundary_context": "The variant is PM2-only but sits in splice/boundary context that may deserve checking.",
        "PM2_only_low_information_background": "The variant has only absence from population data and no strong local prioritization signal.",
        "BP4_BP7_PM2_near_pathogenic_density": "Benign splice prediction is present, but the local region contains many pathogenic variants.",
        "BP4_BP7_PM2_boundary_splice_benign": "The variant is at or near a splice boundary, but SpliceAI/BP7 logic points benign.",
        "BP4_BP7_PM2_low_splice_benign_prediction": "Benign splice prediction plus PM2 remains VUS because evidence is weak and mixed.",
    }
    return meanings.get(bucket, "")


def practical_next_step(bucket: str) -> str:
    steps = {
        "PS3_PM2_needs_one_more_pathogenic_support": "Look for PP1, PM3, PS4, PP4, curated PS1, or another accepted supporting pathogenic item.",
        "PS3_only_needs_independent_support": "Look first for population absence, segregation, case-control, phenotype specificity, or curated functional interpretation.",
        "PS3_PP3_splice_or_functional_followup": "Prioritize RNA or independent functional/clinical evidence because signals already point in the same direction.",
        "complex_strong_pathogenic_combo_review": "Manual rule-combination review; check if non-automated criteria already exist in external sources.",
        "PM2_only_high_priority_context": "Review these before generic PM2-only variants; use neighborhood, exon/domain, and external databases.",
        "PM2_only_pathogenic_neighborhood": "Do not prioritize on density alone; use only if another independent signal or external evidence is available.",
        "PM2_only_splice_or_boundary_context": "Check transcript, splice boundary annotation, RNA evidence, and whether SpliceAI remains low on reference transcript.",
        "PM2_only_low_information_background": "Lowest priority; needs external data before classification can move.",
        "BP4_BP7_PM2_near_pathogenic_density": "Do not over-trust benign splice prediction; review local pathogenic mechanism and possible non-splice effect.",
        "BP4_BP7_PM2_boundary_splice_benign": "RNA evidence would be most useful; one additional benign supporting item could move these toward likely benign.",
        "BP4_BP7_PM2_low_splice_benign_prediction": "Consider as benign-leaning VUS; needs additional benign evidence such as BS2, BS4, or functional/RNA confirmation.",
    }
    return steps.get(bucket, "")


def summarize_buckets(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_bucket: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_bucket[str(row["deep_bucket"])].append(row)
    output = []
    for bucket, items in by_bucket.items():
        category_counts = Counter(row["bottleneck_category"] for row in items)
        high_priority = sum(1 for row in items if int(row["priority_score_int"]) >= 50)
        near_path = sum(1 for row in items if int(row["near_pathogenic_int"]) >= 5)
        high_splice = sum(1 for row in items if float(row["spliceai_float"]) >= 0.20)
        output.append(
            {
                "deep_bucket": bucket,
                "parent_bottleneck": category_counts.most_common(1)[0][0],
                "count": len(items),
                "high_priority_count": high_priority,
                "near_pathogenic_dense_count": near_path,
                "high_spliceai_count": high_splice,
                "top_combo": Counter(row["criteria_combo"] for row in items).most_common(1)[0][0],
                "what_it_means": what_it_means(bucket),
                "practical_next_step": practical_next_step(bucket),
            }
        )
    return sorted(output, key=lambda row: int(row["count"]), reverse=True)


def summarize_exons(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_exon: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_exon[(str(row["bottleneck_category"]), str(row["gene"]), str(row["cds_exon"]))].append(row)
    output = []
    for (category, gene, exon), items in by_exon.items():
        if not exon:
            continue
        output.append(
            {
                "bottleneck_category": category,
                "gene": gene,
                "cds_exon": exon,
                "count": len(items),
                "high_priority_count": sum(1 for row in items if int(row["priority_score_int"]) >= 50),
                "near_pathogenic_dense_count": sum(1 for row in items if int(row["near_pathogenic_int"]) >= 5),
                "max_spliceai": f"{max(float(row['spliceai_float']) for row in items):.2f}",
                "top_deep_bucket": Counter(row["deep_bucket"] for row in items).most_common(1)[0][0],
            }
        )
    return sorted(output, key=lambda row: (int(row["high_priority_count"]), int(row["count"])), reverse=True)


def action_queue(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    selected = sorted(
        rows,
        key=lambda row: (
            action_rank(str(row["deep_bucket"])),
            -int(row["priority_score_int"]),
            -float(row["spliceai_float"]),
            -int(row["near_pathogenic_int"]),
            str(row["gene"]),
            str(row["c_notation"]),
        ),
    )
    output = []
    for row in selected:
        output.append(
            {
                "deep_bucket": row["deep_bucket"],
                "bottleneck_category": row["bottleneck_category"],
                "gene": row["gene"],
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "variant_type": row["variant_type"],
                "criteria_combo": row["criteria_combo"],
                "total_points": row["total_points"],
                "priority_score": row["priority_score"],
                "priority_reasons": row["priority_reasons"],
                "cds_exon": row["cds_exon"],
                "boundary_distance": row["boundary_distance"],
                "near_pathogenic_20bp": row["near_pathogenic_20bp"],
                "spliceai_score": row["spliceai_score"],
                "practical_next_step": row["practical_next_step"],
            }
        )
    return output


def action_rank(bucket: str) -> int:
    order = {
        "PS3_PP3_splice_or_functional_followup": 0,
        "PS3_PM2_needs_one_more_pathogenic_support": 1,
        "complex_strong_pathogenic_combo_review": 2,
        "PS3_only_needs_independent_support": 3,
        "BP4_BP7_PM2_near_pathogenic_density": 4,
        "BP4_BP7_PM2_boundary_splice_benign": 5,
        "PM2_only_high_priority_context": 6,
        "PM2_only_pathogenic_neighborhood": 7,
        "PM2_only_splice_or_boundary_context": 8,
        "BP4_BP7_PM2_low_splice_benign_prediction": 9,
        "PM2_only_low_information_background": 10,
    }
    return order.get(bucket, 99)


def subset(rows: list[dict[str, object]], category: str) -> list[dict[str, object]]:
    return [row for row in rows if row["bottleneck_category"] == category]


def draw_barplot(path: Path, title: str, rows: list[dict[str, object]], label_col: str, value_col: str) -> None:
    selected = rows[:18]
    width = 1260
    height = max(250, 78 + len(selected) * 34)
    left = 460
    top = 58
    plot_w = 650
    max_value = max((int(row[value_col]) for row in selected), default=1)
    body = [f'<text x="28" y="34" class="title">{esc(title)}</text>']
    for index, row in enumerate(selected):
        y = top + index * 34
        value = int(row[value_col])
        bar_w = plot_w * value / max_value
        body.append(f'<text x="24" y="{y + 18}" class="small">{esc(row[label_col])}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{bar_w:.1f}" height="22" fill="#4f46e5"/>')
        body.append(f'<text x="{left + bar_w + 8:.1f}" y="{y + 17}" class="small">{value}</text>')
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
text {{ font-family: Arial, Helvetica, sans-serif; fill: #172033; }}
.small {{ font-size: 12px; }}
.title {{ font-size: 18px; font-weight: 700; }}
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{chr(10).join(body)}
</svg>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(svg, encoding="utf-8")


def report_table(rows: list[dict[str, object]], columns: list[str], limit: int = 12) -> str:
    if not rows:
        return "| none |\n| --- |\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows[:limit]:
        body.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join([header, sep] + body)


def write_report(
    bucket_rows: list[dict[str, object]],
    exon_rows: list[dict[str, object]],
    queue_rows: list[dict[str, object]],
    rows: list[dict[str, object]],
) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    one_step = subset(rows, "strong_pathogenic_combo_one_step_short")
    pm2_only = subset(rows, "PM2_only")
    splice_benign = subset(rows, "mixed_splice_benign_and_pathogenic")
    one_step_queue = [row for row in queue_rows if row["bottleneck_category"] == "strong_pathogenic_combo_one_step_short"]
    pm2_queue = [row for row in queue_rows if row["bottleneck_category"] == "PM2_only"]
    splice_queue = [row for row in queue_rows if row["bottleneck_category"] == "mixed_splice_benign_and_pathogenic"]

    text = f"""# VUS Bottleneck Deep Dive

Generated: {generated}

Inputs:

- `tables/vus_bottleneck/vus_bottleneck_variant_examples.csv`
- `tables/vus_prioritization/vus_priority_annotated.csv`

## Why We Did This

The VUS bottleneck analysis showed three groups worth looking at more closely. The goal here is not to reclassify variants. The goal is to understand what kind of evidence is missing and which subsets are worth curator time.

We focused on:

1. `strong_pathogenic_combo_one_step_short`: VUS variants with strong pathogenic evidence that miss the exact ENIGMA combination threshold.
2. `PM2_only`: the largest VUS block, where absence from population data is the only automated evidence.
3. `BP4+BP7+PM2_Supporting`: benign splice prediction mixed with PM2 background, still ending as VUS.

## Overall Deep Buckets

{report_table(bucket_rows, ["deep_bucket", "parent_bottleneck", "count", "high_priority_count", "near_pathogenic_dense_count", "high_spliceai_count", "top_combo"], limit=20)}

## 1. One-Step-Short Pathogenic VUS

Rows reviewed: `{len(one_step)}`

These are the most actionable pathogenic-side VUS. They often already have `PS3` and sometimes `PP3`, but they still do not satisfy the exact ENIGMA combination. This means the application is behaving conservatively: strong evidence alone, or strong plus one weak signal, is not enough.

What can be concluded:

- these are not automatic upgrades
- they are a good manual review queue
- the most useful missing evidence types are PP1, PM3, PS4, PP4, curated PS1, or another accepted supporting pathogenic item
- variants with both `PS3` and `PP3` are especially interesting because functional and computational/splice evidence point in the same direction

Top examples:

{report_table(one_step_queue, ["deep_bucket", "gene", "c_notation", "p_notation", "criteria_combo", "total_points", "priority_score", "spliceai_score", "near_pathogenic_20bp"], limit=18)}

## 2. PM2-Only VUS

Rows reviewed: `{len(pm2_only)}`

This is the largest group, but also the least informative by itself. `PM2_Supporting` means the variant is absent or very rare in the available population data. That does not tell us whether the variant is damaging.

What can be concluded:

- most PM2-only variants should stay low priority until external evidence appears
- pathogenic density alone is too broad in this full-SNV synthetic landscape and should not be treated as an action trigger
- a smaller subset becomes interesting if PM2-only also has splice/boundary context or another independent signal
- PM2-only is a data-gap group, not a biological mechanism

Top PM2-only examples selected by local context. These are not proposed upgrades; they are examples of how broad the PM2-only background remains even after adding local context:

{report_table(pm2_queue, ["deep_bucket", "gene", "c_notation", "p_notation", "variant_type", "priority_score", "priority_reasons", "cds_exon", "near_pathogenic_20bp"], limit=18)}

## 3. BP4/BP7 Plus PM2

Rows reviewed: `{len(splice_benign)}`

This group has benign splice prediction (`BP4+BP7`) but also PM2 background, so the final point total usually remains VUS. This is not a contradiction in the biological sense; it is a weak benign signal mixed with weak pathogenic background.

What can be concluded:

- these are benign-leaning VUS, not likely benign by current automated rules
- one additional benign supporting item could often move them toward likely benign in the contradictory-evidence point path
- RNA evidence would be the cleanest way to strengthen the benign splice branch
- cases in pathogenic-dense regions deserve manual caution, because benign splice prediction does not exclude non-splice pathogenic mechanisms

Top BP4/BP7 plus PM2 examples:

{report_table(splice_queue, ["deep_bucket", "gene", "c_notation", "p_notation", "variant_type", "priority_score", "boundary_distance", "near_pathogenic_20bp", "spliceai_score"], limit=18)}

## Exon-Level Concentrations

{report_table(exon_rows, ["bottleneck_category", "gene", "cds_exon", "count", "high_priority_count", "near_pathogenic_dense_count", "max_spliceai", "top_deep_bucket"], limit=25)}

## Practical Interpretation

For curation, the best next queue is not the largest group. It is the group where one credible additional evidence item can change the interpretation:

1. start with `PS3_PP3_splice_or_functional_followup` and `PS3_PM2_needs_one_more_pathogenic_support`
2. then review BP4/BP7 plus PM2 variants only when there is an additional reason beyond pathogenic density
3. keep most PM2-only variants as low-priority background unless another evidence source appears

For the application, this suggests that the UI could eventually expose a VUS explanation like: "VUS because evidence is one supporting pathogenic criterion short", "VUS because PM2 is the only evidence", or "benign-leaning VUS, needs additional benign confirmation."

## Outputs

- `tables/vus_bottleneck_deep_dive/deep_bucket_summary.csv`
- `tables/vus_bottleneck_deep_dive/deep_dive_action_queue.csv`
- `tables/vus_bottleneck_deep_dive/deep_dive_exon_summary.csv`
- `tables/vus_bottleneck_deep_dive/one_step_short_variants.csv`
- `tables/vus_bottleneck_deep_dive/pm2_only_variants.csv`
- `tables/vus_bottleneck_deep_dive/bp4_bp7_pm2_variants.csv`
- `plots/13_vus_bottleneck_deep_dive/deep_bucket_counts.svg`
- `plots/13_vus_bottleneck_deep_dive/action_queue_top_buckets.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    bucket_rows = summarize_buckets(rows)
    exon_rows = summarize_exons(rows)
    queue_rows = action_queue(rows)

    base_fields = [
        "deep_bucket",
        "bottleneck_category",
        "gene",
        "c_notation",
        "p_notation",
        "variant_type",
        "criteria_combo",
        "total_points",
        "priority_score",
        "priority_reasons",
        "cds_exon",
        "boundary_distance",
        "near_pathogenic_20bp",
        "spliceai_score",
        "practical_next_step",
    ]

    write_csv(
        TABLE_DIR / "deep_bucket_summary.csv",
        bucket_rows,
        [
            "deep_bucket",
            "parent_bottleneck",
            "count",
            "high_priority_count",
            "near_pathogenic_dense_count",
            "high_spliceai_count",
            "top_combo",
            "what_it_means",
            "practical_next_step",
        ],
    )
    write_csv(TABLE_DIR / "deep_dive_action_queue.csv", queue_rows, base_fields)
    write_csv(
        TABLE_DIR / "deep_dive_exon_summary.csv",
        exon_rows,
        [
            "bottleneck_category",
            "gene",
            "cds_exon",
            "count",
            "high_priority_count",
            "near_pathogenic_dense_count",
            "max_spliceai",
            "top_deep_bucket",
        ],
    )
    write_csv(TABLE_DIR / "one_step_short_variants.csv", [row for row in queue_rows if row["bottleneck_category"] == "strong_pathogenic_combo_one_step_short"], base_fields)
    write_csv(TABLE_DIR / "pm2_only_variants.csv", [row for row in queue_rows if row["bottleneck_category"] == "PM2_only"], base_fields)
    write_csv(TABLE_DIR / "bp4_bp7_pm2_variants.csv", [row for row in queue_rows if row["bottleneck_category"] == "mixed_splice_benign_and_pathogenic"], base_fields)

    action_bucket_rows = [
        {"deep_bucket": bucket, "count": count}
        for bucket, count in Counter(row["deep_bucket"] for row in queue_rows[:500]).most_common()
    ]
    draw_barplot(PLOT_DIR / "deep_bucket_counts.svg", "Deep bucket counts", bucket_rows, "deep_bucket", "count")
    draw_barplot(PLOT_DIR / "action_queue_top_buckets.svg", "Top buckets among first 500 action-queue rows", action_bucket_rows, "deep_bucket", "count")
    write_report(bucket_rows, exon_rows, queue_rows, rows)
    print(f"Wrote VUS bottleneck deep dive to {REPORT}")


if __name__ == "__main__":
    main()
