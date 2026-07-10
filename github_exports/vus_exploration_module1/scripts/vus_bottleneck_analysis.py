"""Decision bottleneck analysis for VUS variants in the BRCA Module 1 snapshot."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = OUT_DIR / "tables" / "vus_bottleneck"
PLOT_DIR = OUT_DIR / "plots" / "12_vus_bottleneck"
REPORT = OUT_DIR / "vus_bottleneck_report.md"

PATH_CODES = {"PVS1", "PM5_PTC", "PM2_Supporting", "PP3", "PP4", "PS1", "PS3"}
BENIGN_CODES = {"BA1", "BS1_Strong", "BS1_Supporting", "BS3", "BP1", "BP4", "BP5", "BP7"}


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_criteria(criteria: str) -> dict[str, dict[str, object]]:
    parsed: dict[str, dict[str, object]] = {}
    if not criteria:
        return parsed
    for item in criteria.split(";"):
        parts = item.split(":")
        if len(parts) < 3:
            continue
        code, strength, points = parts[0], parts[1], parts[2]
        try:
            point_value = int(points)
        except ValueError:
            point_value = 0
        parsed[code] = {"strength": strength, "points": point_value}
    return parsed


def criteria_combo(criteria: dict[str, dict[str, object]]) -> str:
    return "+".join(sorted(criteria)) if criteria else "none"


def spliceai_float(row: dict[str, str]) -> float:
    try:
        return float(row["spliceai_score"])
    except (KeyError, TypeError, ValueError):
        return 0.0


def has_pathogenic(criteria: dict[str, dict[str, object]]) -> bool:
    return any(item["points"] > 0 for item in criteria.values())


def has_benign(criteria: dict[str, dict[str, object]]) -> bool:
    return any(item["points"] < 0 for item in criteria.values())


def strength_counts(criteria: dict[str, dict[str, object]]) -> dict[str, int]:
    counts = {"very_strong": 0, "strong": 0, "moderate": 0, "supporting": 0}
    for item in criteria.values():
        points = int(item["points"])
        if points <= 0:
            continue
        strength = str(item["strength"]).lower().replace("-", " ")
        if points >= 8 or "very strong" in strength:
            counts["very_strong"] += 1
        elif points >= 4 or "strong" in strength:
            counts["strong"] += 1
        elif points >= 2 or "moderate" in strength:
            counts["moderate"] += 1
        else:
            counts["supporting"] += 1
    return counts


def bottleneck_category(row: dict[str, object]) -> str:
    criteria = row["criteria_parsed"]
    codes = set(criteria)
    path = has_pathogenic(criteria)
    benign = has_benign(criteria)
    strengths = strength_counts(criteria)

    if path and benign:
        if "PVS1" in codes and "BS3" in codes:
            return "conflicting_PVS1_BS3"
        if "PP3" in codes and "BS3" in codes:
            return "conflicting_PP3_BS3"
        if {"BP4", "BP7"} & codes and PATH_CODES & codes:
            return "mixed_splice_benign_and_pathogenic"
        return "other_conflicting_evidence"

    if not path and not benign:
        return "no_automated_evidence"

    if benign and not path:
        if "BS3" in codes:
            return "benign_functional_evidence_not_enough"
        if {"BP4", "BP7"} & codes:
            return "benign_splice_prediction_not_enough"
        return "benign_evidence_not_enough"

    if "PVS1" in codes:
        return "PVS1_without_required_support"
    if "PS3" in codes:
        if strengths["strong"] >= 1 and strengths["moderate"] == 0 and strengths["supporting"] < 2:
            return "strong_pathogenic_combo_one_step_short"
        return "functional_pathogenic_evidence_not_enough"
    if "PP3" in codes:
        return "computational_pathogenic_evidence_not_enough"
    if codes == {"PM2_Supporting"}:
        return "PM2_only"
    if "PM2_Supporting" in codes:
        return "PM2_plus_weak_context"
    return "pathogenic_evidence_not_enough"


def suggested_next_evidence(category: str) -> str:
    suggestions = {
        "conflicting_PVS1_BS3": "Expert adjudication of functional BS3 versus PVS1 boundary; consider curated ENIGMA interpretation.",
        "conflicting_PP3_BS3": "Manual review of functional assay evidence and splice prediction; RNA evidence may clarify if splice signal is real.",
        "mixed_splice_benign_and_pathogenic": "Review whether benign splice/RNA branch and pathogenic evidence should both stand for this variant.",
        "other_conflicting_evidence": "Expert review of mixed benign and pathogenic evidence before any automatic class change.",
        "no_automated_evidence": "Needs external evidence not in Module 1, such as case-control, segregation, curated clinical assertions, or functional data.",
        "benign_functional_evidence_not_enough": "Additional benign evidence such as BS2, BS4, population evidence, or curated functional strength could resolve toward benign.",
        "benign_splice_prediction_not_enough": "RNA/functional confirmation or stronger benign evidence would be needed to move from VUS to likely benign.",
        "benign_evidence_not_enough": "Needs an additional benign criterion or stronger benign evidence.",
    "PVS1_without_required_support": "Additional supporting pathogenic evidence would be needed; PVS1 strength and context should be checked before interpretation.",
        "strong_pathogenic_combo_one_step_short": "One additional supporting pathogenic criterion, or one moderate criterion, could move many cases toward likely pathogenic.",
        "functional_pathogenic_evidence_not_enough": "Needs additional pathogenic evidence beyond functional evidence, such as PP1, PM3, PP4, PS4, or curated PS1 context.",
        "computational_pathogenic_evidence_not_enough": "Computational evidence alone is insufficient; needs functional, segregation, case-control, trans, or curated clinical evidence.",
        "PM2_only": "Absent population frequency alone is insufficient; needs independent pathogenic or benign evidence.",
        "PM2_plus_weak_context": "Needs stronger independent evidence; PM2 is background support, not a decision driver.",
        "pathogenic_evidence_not_enough": "Needs additional pathogenic evidence to meet ENIGMA combination rules.",
    }
    return suggestions.get(category, "Additional variant-level evidence would be needed.")


def load_vus_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["predicted_class"] != "3":
                continue
            criteria = parse_criteria(row["criteria"])
            row["criteria_parsed"] = criteria
            row["criteria_combo"] = criteria_combo(criteria)
            row["spliceai_float"] = spliceai_float(row)
            row["total_points_int"] = int(row["total_points"])
            row["bottleneck_category"] = bottleneck_category(row)
            row["suggested_next_evidence"] = suggested_next_evidence(str(row["bottleneck_category"]))
            rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize_categories(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    total = len(rows)
    output = []
    by_category: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_category[str(row["bottleneck_category"])].append(row)
    for category, items in by_category.items():
        genes = Counter(row["gene"] for row in items)
        high_splice = sum(1 for row in items if float(row["spliceai_float"]) >= 0.20)
        combos = Counter(row["criteria_combo"] for row in items)
        output.append(
            {
                "bottleneck_category": category,
                "count": len(items),
                "percent_of_vus": f"{len(items) / total * 100:.2f}" if total else "0.00",
                "brca1_count": genes["BRCA1"],
                "brca2_count": genes["BRCA2"],
                "high_spliceai_count": high_splice,
                "high_spliceai_percent": f"{high_splice / len(items) * 100:.2f}" if items else "0.00",
                "top_combo": combos.most_common(1)[0][0] if combos else "",
                "suggested_next_evidence": suggested_next_evidence(category),
            }
        )
    return sorted(output, key=lambda row: int(row["count"]), reverse=True)


def summarize_combos(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts = Counter()
    genes = Counter()
    points = defaultdict(list)
    high_splice = Counter()
    for row in rows:
        key = (row["bottleneck_category"], row["criteria_combo"])
        counts[key] += 1
        genes[(key, row["gene"])] += 1
        points[key].append(int(row["total_points_int"]))
        if float(row["spliceai_float"]) >= 0.20:
            high_splice[key] += 1
    output = []
    for (category, combo), count in counts.most_common():
        point_values = points[(category, combo)]
        output.append(
            {
                "bottleneck_category": category,
                "criteria_combo": combo,
                "count": count,
                "brca1_count": genes[((category, combo), "BRCA1")],
                "brca2_count": genes[((category, combo), "BRCA2")],
                "median_points": median(point_values),
                "high_spliceai_count": high_splice[(category, combo)],
                "suggested_next_evidence": suggested_next_evidence(category),
            }
        )
    return output


def detail_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    selected = sorted(
        rows,
        key=lambda row: (
            category_rank(str(row["bottleneck_category"])),
            -float(row["spliceai_float"]),
            row["gene"],
            row["c_notation"],
        ),
    )
    output = []
    for row in selected:
        output.append(
            {
                "bottleneck_category": row["bottleneck_category"],
                "gene": row["gene"],
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "variant_type": row["variant_type"],
                "criteria_combo": row["criteria_combo"],
                "total_points": row["total_points"],
                "spliceai_score": row["spliceai_score"],
                "classification_note": row["classification_note"],
                "suggested_next_evidence": row["suggested_next_evidence"],
            }
        )
    return output


def category_rank(category: str) -> int:
    order = {
        "conflicting_PVS1_BS3": 0,
        "conflicting_PP3_BS3": 1,
        "other_conflicting_evidence": 2,
        "mixed_splice_benign_and_pathogenic": 3,
        "strong_pathogenic_combo_one_step_short": 4,
        "functional_pathogenic_evidence_not_enough": 5,
        "computational_pathogenic_evidence_not_enough": 6,
        "benign_splice_prediction_not_enough": 7,
        "PM2_only": 8,
        "no_automated_evidence": 9,
    }
    return order.get(category, 99)


def median(values: list[int]) -> str:
    ordered = sorted(values)
    if not ordered:
        return ""
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return str(ordered[mid])
    return f"{(ordered[mid - 1] + ordered[mid]) / 2:.1f}"


def draw_barplot(
    path: Path,
    title: str,
    rows: list[dict[str, object]],
    label_col: str,
    value_col: str,
    subtitle: str = "",
) -> None:
    selected = rows[:18]
    width = 1180
    height = max(250, 102 + len(selected) * 34)
    left = 410
    top = 86
    plot_w = 610
    max_value = max((int(row[value_col]) for row in selected), default=1)
    body = [f'<text x="28" y="34" class="title">{esc(title)}</text>']
    if subtitle:
        body.append(f'<text x="28" y="56" class="small">{esc(subtitle)}</text>')
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
    rows: list[dict[str, object]],
    category_rows: list[dict[str, object]],
    combo_rows: list[dict[str, object]],
    detail: list[dict[str, object]],
) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    total = len(rows)
    top_actionable = [
        row for row in detail
        if row["bottleneck_category"] in {
            "conflicting_PVS1_BS3",
            "conflicting_PP3_BS3",
            "strong_pathogenic_combo_one_step_short",
            "functional_pathogenic_evidence_not_enough",
        }
    ]
    text = f"""# VUS Decision Bottleneck Analysis

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

## Purpose

This analysis asks why automated Module 1 leaves variants as VUS and what kind
of evidence type would theoretically be needed next. It is not another
classifier and it does not change any variant class. It should be read as an
evidence-pattern audit, not as a validated list of variants that must be
manually reviewed.

The main idea is to group VUS variants by decision mechanism:

- not enough pathogenic evidence to meet ENIGMA combinations
- not enough benign evidence to meet benign combinations
- contradictory benign and pathogenic evidence
- no automated evidence beyond the input variant itself

This is more useful than a simple BRCA1 versus BRCA2 comparison because it points to what is missing.

## Overall VUS Count

Total VUS rows analyzed: `{total}`

## Bottleneck Summary

{report_table(category_rows, ["bottleneck_category", "count", "percent_of_vus", "brca1_count", "brca2_count", "high_spliceai_count", "top_combo"], limit=20)}

## Most Common Blocking Combinations

{report_table(combo_rows, ["bottleneck_category", "criteria_combo", "count", "brca1_count", "brca2_count", "median_points", "high_spliceai_count"], limit=25)}

## Example Variants From Selected Bottleneck Groups

{report_table(top_actionable, ["bottleneck_category", "gene", "c_notation", "p_notation", "variant_type", "criteria_combo", "total_points", "spliceai_score"], limit=25)}

## Interpretation

The largest bottleneck is expected: variants with only `PM2_Supporting` remain VUS because absence from population data is not independent pathogenic evidence. These variants mostly need outside evidence, not more arithmetic.

The second major bottleneck is benign splice prediction plus PM2 background, commonly `BP4+BP7+PM2_Supporting`. This is a VUS because the benign evidence is not strong enough and PM2 pulls in the opposite direction. RNA/functional confirmation or stronger benign evidence would be needed to resolve these.

The most actionable pathogenic-side group is `strong_pathogenic_combo_one_step_short`, usually `PM2_Supporting+PS3`. These are close to likely pathogenic: one additional supporting pathogenic criterion, or one moderate pathogenic criterion, could move many cases. This group also catches rarer high-point combinations that still miss the exact ENIGMA Table 3 combination pattern.

The conflict groups are smaller but more important for quality control. `PP3+BS3` and `PVS1+BS3` should be reviewed as evidence-adjudication cases rather than treated as automatic upgrades or downgrades.

## Outputs

- `tables/vus_bottleneck/vus_bottleneck_summary.csv`
- `tables/vus_bottleneck/vus_bottleneck_combinations.csv`
- `tables/vus_bottleneck/vus_bottleneck_variant_examples.csv`
- `plots/12_vus_bottleneck/vus_bottleneck_categories.svg`
- `plots/12_vus_bottleneck/vus_bottleneck_combinations.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_vus_rows()
    category_rows = summarize_categories(rows)
    combo_rows = summarize_combos(rows)
    detail = detail_rows(rows)

    write_csv(
        TABLE_DIR / "vus_bottleneck_summary.csv",
        category_rows,
        [
            "bottleneck_category",
            "count",
            "percent_of_vus",
            "brca1_count",
            "brca2_count",
            "high_spliceai_count",
            "high_spliceai_percent",
            "top_combo",
            "suggested_next_evidence",
        ],
    )
    write_csv(
        TABLE_DIR / "vus_bottleneck_combinations.csv",
        combo_rows,
        [
            "bottleneck_category",
            "criteria_combo",
            "count",
            "brca1_count",
            "brca2_count",
            "median_points",
            "high_spliceai_count",
            "suggested_next_evidence",
        ],
    )
    write_csv(
        TABLE_DIR / "vus_bottleneck_variant_examples.csv",
        detail,
        [
            "bottleneck_category",
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "criteria_combo",
            "total_points",
            "spliceai_score",
            "classification_note",
            "suggested_next_evidence",
        ],
    )
    draw_barplot(
        PLOT_DIR / "vus_bottleneck_categories.svg",
        "Why Module 1 VUS remain unresolved",
        category_rows,
        "bottleneck_category",
        "count",
        "Evidence-pattern groups only; not a validated manual-review priority criterion.",
    )
    draw_barplot(
        PLOT_DIR / "vus_bottleneck_combinations.svg",
        "Most common Module 1 VUS evidence combinations",
        combo_rows,
        "criteria_combo",
        "count",
        "Shows applied automatable criteria among VUS; not an independent classification rule.",
    )
    write_report(rows, category_rows, combo_rows, detail)
    print(f"Wrote VUS bottleneck analysis to {REPORT}")


if __name__ == "__main__":
    main()
