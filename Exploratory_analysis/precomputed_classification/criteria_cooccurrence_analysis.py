"""Criteria co-occurrence analysis for the precomputed BRCA Module 1 snapshot."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from datetime import datetime
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = OUT_DIR / "tables" / "criteria_cooccurrence"
PLOT_DIR = OUT_DIR / "plots" / "10_criteria_cooccurrence"
REPORT = OUT_DIR / "criteria_cooccurrence_report.md"

GROUP_ORDER = ["benign", "vus", "pathogenic"]
GROUP_LABELS = {
    "benign": "Benign/Likely Benign",
    "vus": "VUS",
    "pathogenic": "Likely Pathogenic/Pathogenic",
}
GROUP_COLORS = {
    "benign": "#2b8cbe",
    "vus": "#8a8f9c",
    "pathogenic": "#d7301f",
}
CRITERIA_ORDER = [
    "BA1",
    "BS1_Strong",
    "BS1_Supporting",
    "BS3",
    "BP1",
    "BP4",
    "BP5",
    "BP7",
    "PM2_Supporting",
    "PM5_PTC",
    "PP3",
    "PP4",
    "PS1",
    "PS3",
    "PVS1",
]
PATHOGENIC_LEANING = {"PVS1", "PM5_PTC", "PP3", "PP4", "PS1", "PS3", "PM2_Supporting"}
CONFLICT_PATHOGENIC_LEANING = {"PVS1", "PM5_PTC", "PP3", "PP4", "PS1", "PS3"}
BENIGN_LEANING = {"BA1", "BS1_Strong", "BS1_Supporting", "BS3", "BP1", "BP4", "BP5", "BP7"}


def class_group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


def parse_criteria(criteria: str) -> list[str]:
    if not criteria:
        return []
    codes = []
    for item in criteria.split(";"):
        code = item.split(":", 1)[0]
        if code:
            codes.append(code)
    return sorted(set(codes), key=lambda code: CRITERIA_ORDER.index(code) if code in CRITERIA_ORDER else 999)


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def heat_color(value: float, max_value: float, base: tuple[int, int, int]) -> str:
    intensity = 0.0 if max_value <= 0 else min(1.0, math.sqrt(max(value, 0.0) / max_value))
    bg = (247, 250, 252)
    rgb = tuple(round(bg[i] * (1 - intensity) + base[i] * intensity) for i in range(3))
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


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
            codes = parse_criteria(row["criteria"])
            row["criteria_codes"] = codes
            row["criteria_combo"] = "+".join(codes) if codes else "none"
            row["group"] = class_group(row["predicted_class"])
            rows.append(row)
    return rows


def criteria_counts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    total = len(rows)
    counts = Counter(code for row in rows for code in row["criteria_codes"])
    output = []
    for code in CRITERIA_ORDER:
        count = counts[code]
        if not count:
            continue
        record = {"criterion": code, "count": count, "percent": f"{count / total * 100:.2f}"}
        for group in GROUP_ORDER:
            group_rows = [row for row in rows if row["group"] == group]
            group_count = sum(1 for row in group_rows if code in row["criteria_codes"])
            record[f"{group}_count"] = group_count
            record[f"{group}_percent"] = f"{group_count / len(group_rows) * 100:.2f}" if group_rows else "0.00"
        output.append(record)
    return output


def pair_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    total_by_group = Counter(row["group"] for row in rows)
    pair_counts = Counter()
    pair_group_counts = Counter()
    for row in rows:
        codes = row["criteria_codes"]
        for left, right in combinations(codes, 2):
            pair = tuple(sorted((left, right), key=lambda code: CRITERIA_ORDER.index(code) if code in CRITERIA_ORDER else 999))
            pair_counts[pair] += 1
            pair_group_counts[(pair, row["group"])] += 1
    output = []
    for pair, count in pair_counts.most_common():
        record = {"criterion_a": pair[0], "criterion_b": pair[1], "count": count}
        for group in GROUP_ORDER:
            group_count = pair_group_counts[(pair, group)]
            record[f"{group}_count"] = group_count
            record[f"{group}_percent_of_pair"] = f"{group_count / count * 100:.2f}" if count else "0.00"
            record[f"{group}_percent_of_group"] = f"{group_count / total_by_group[group] * 100:.2f}" if total_by_group[group] else "0.00"
        output.append(record)
    return output


def combo_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    combo_counts = Counter(row["criteria_combo"] for row in rows)
    combo_group_counts = Counter((row["criteria_combo"], row["group"]) for row in rows)
    output = []
    for combo, count in combo_counts.most_common():
        record = {"criteria_combo": combo, "count": count}
        for group in GROUP_ORDER:
            group_count = combo_group_counts[(combo, group)]
            record[f"{group}_count"] = group_count
            record[f"{group}_percent_of_combo"] = f"{group_count / count * 100:.2f}" if count else "0.00"
        output.append(record)
    return output


def vus_holding_combos(combo_output: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = [
        row
        for row in combo_output
        if int(row["count"]) >= 5 and float(row["vus_percent_of_combo"]) >= 80.0
    ]
    return sorted(rows, key=lambda row: (float(row["vus_percent_of_combo"]), int(row["count"])), reverse=True)


def pathogenic_combos(combo_output: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = [
        row
        for row in combo_output
        if int(row["count"]) >= 5 and float(row["pathogenic_percent_of_combo"]) >= 80.0
    ]
    return sorted(rows, key=lambda row: (float(row["pathogenic_percent_of_combo"]), int(row["count"])), reverse=True)


def sanity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for row in rows:
        codes = set(row["criteria_codes"])
        reasons = []
        if codes & CONFLICT_PATHOGENIC_LEANING and codes & BENIGN_LEANING:
            reasons.append("pathogenic_and_benign_evidence")
        if "PP3" in codes and "BS3" in codes:
            reasons.append("PP3_with_BS3")
        if "PP3" in codes and "BP7" in codes:
            reasons.append("PP3_with_BP7")
        if "PVS1" in codes and row["group"] == "vus":
            reasons.append("PVS1_but_VUS")
        if "PS3" in codes and row["group"] == "vus":
            reasons.append("PS3_but_VUS")
        if "BS3" in codes and row["group"] == "pathogenic":
            reasons.append("BS3_but_pathogenic")
        if not reasons:
            continue
        output.append(
            {
                "sanity_reason": ";".join(reasons),
                "gene": row["gene"],
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "variant_type": row["variant_type"],
                "criteria_combo": row["criteria_combo"],
                "total_points": row["total_points"],
                "predicted_class": row["predicted_class"],
                "predicted_label": row["predicted_label"],
                "spliceai_score": row["spliceai_score"],
            }
        )
    return output


def draw_pair_heatmap(pair_output: list[dict[str, object]], group: str | None = None) -> None:
    criteria = [code for code in CRITERIA_ORDER if any(row["criterion_a"] == code or row["criterion_b"] == code for row in pair_output)]
    values = {}
    for row in pair_output:
        key = (row["criterion_a"], row["criterion_b"])
        value = int(row[f"{group}_count"]) if group else int(row["count"])
        values[key] = value
        values[(key[1], key[0])] = value
    cell = 34
    left = 145
    top = 95
    width = left + len(criteria) * cell + 40
    height = top + len(criteria) * cell + 50
    max_value = max(values.values()) if values else 1
    title = f"Criteria co-occurrence heatmap {GROUP_LABELS[group]}" if group else "Criteria co-occurrence heatmap all variants"
    body = [f'<text x="28" y="34" class="title">{esc(title)}</text>']
    for xi, code in enumerate(criteria):
        x = left + xi * cell
        body.append(f'<text x="{x + 4}" y="{top - 10}" class="small" transform="rotate(-45 {x + 4} {top - 10})">{esc(code)}</text>')
    for yi, row_code in enumerate(criteria):
        y = top + yi * cell
        body.append(f'<text x="20" y="{y + 21}" class="small">{esc(row_code)}</text>')
        for xi, col_code in enumerate(criteria):
            x = left + xi * cell
            if row_code == col_code:
                color = "#e5e7eb"
                value = 0
            else:
                value = values.get((row_code, col_code), 0)
                color = heat_color(value, max_value, (37, 99, 235))
            body.append(f'<rect x="{x}" y="{y}" width="{cell - 2}" height="{cell - 2}" fill="{color}"/>')
            if value:
                body.append(f'<text x="{x + 4}" y="{y + 20}" class="small">{value}</text>')
    suffix = group if group else "all"
    write_svg(PLOT_DIR / f"criteria_pair_heatmap_{suffix}.svg", "\n".join(body), width, height)


def draw_combo_barplot(combo_output: list[dict[str, object]], filename: str, title: str, rows: list[dict[str, object]]) -> None:
    selected = rows[:18]
    width = 1180
    height = 650
    left = 430
    top = 60
    bar_h = 22
    gap = 9
    plot_w = 620
    max_count = max((int(row["count"]) for row in selected), default=1)
    body = [f'<text x="28" y="34" class="title">{esc(title)}</text>']
    for i, row in enumerate(selected):
        y = top + i * (bar_h + gap)
        count = int(row["count"])
        w = plot_w * count / max_count
        body.append(f'<text x="24" y="{y + 16}" class="small">{esc(row["criteria_combo"][:65])}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="#7c3aed"/>')
        body.append(
            f'<text x="{left + w + 8}" y="{y + 16}" class="small">n={count}, VUS={row["vus_percent_of_combo"]}%, P={row["pathogenic_percent_of_combo"]}%</text>'
        )
    write_svg(PLOT_DIR / filename, "\n".join(body), width, height)


def report_table(rows: list[dict[str, object]], columns: list[str], limit: int = 12) -> str:
    if not rows:
        return "| none |\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows[:limit]:
        body.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join([header, sep] + body)


def write_report(
    counts: list[dict[str, object]],
    pairs: list[dict[str, object]],
    combos: list[dict[str, object]],
    vus_combos: list[dict[str, object]],
    path_combos: list[dict[str, object]],
    sanity: list[dict[str, object]],
) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    sanity_counts = Counter()
    for row in sanity:
        for reason in str(row["sanity_reason"]).split(";"):
            sanity_counts[reason] += 1
    sanity_lines = "\n".join(f"- {reason}: {count}" for reason, count in sanity_counts.most_common())
    text = f"""# Criteria Co-occurrence Analysis

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

This report describes which automated Module 1 criteria occur together and how their combinations relate to benign, VUS, and pathogenic grouped outcomes. It is descriptive and does not add new ACMG/ENIGMA criteria.

## Outputs

- `tables/criteria_cooccurrence/criteria_counts_by_group.csv`
- `tables/criteria_cooccurrence/criteria_pair_cooccurrence.csv`
- `tables/criteria_cooccurrence/criteria_combo_outcomes.csv`
- `tables/criteria_cooccurrence/vus_holding_combinations.csv`
- `tables/criteria_cooccurrence/pathogenic_dominant_combinations.csv`
- `tables/criteria_cooccurrence/criteria_sanity_checks.csv`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_all.svg`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_benign.svg`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_vus.svg`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_pathogenic.svg`
- `plots/10_criteria_cooccurrence/top_vus_holding_combinations.svg`
- `plots/10_criteria_cooccurrence/top_pathogenic_combinations.svg`

## Most Frequent Criteria

{report_table(counts, ["criterion", "count", "percent", "benign_percent", "vus_percent", "pathogenic_percent"], limit=15)}

## Top Criterion Pairs

{report_table(pairs, ["criterion_a", "criterion_b", "count", "benign_count", "vus_count", "pathogenic_count"], limit=15)}

## VUS-holding Combinations

Combinations with at least 5 variants and at least 80% VUS.

{report_table(vus_combos, ["criteria_combo", "count", "vus_count", "vus_percent_of_combo", "benign_count", "pathogenic_count"], limit=15)}

## Pathogenic-dominant Combinations

Combinations with at least 5 variants and at least 80% pathogenic grouped outcome.

{report_table(path_combos, ["criteria_combo", "count", "pathogenic_count", "pathogenic_percent_of_combo", "vus_count", "benign_count"], limit=15)}

## Sanity-check Patterns

{sanity_lines if sanity_lines else "- none"}

Top sanity-check rows:

{report_table(sanity, ["sanity_reason", "gene", "c_notation", "p_notation", "criteria_combo", "predicted_class", "spliceai_score"], limit=20)}

## Reading Notes

`PM2_Supporting` is common because much of this synthetic coding SNV landscape is absent from the available local population-frequency lookup. Pair counts involving `PM2_Supporting` should therefore be read as background context, not as independent biological structure. For the sanity-check table, `PM2_Supporting` is not counted as conflicting pathogenic evidence by itself; this avoids flagging expected combinations such as `BP1+PM2_Supporting` as false conflicts.
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    counts = criteria_counts(rows)
    pairs = pair_rows(rows)
    combos = combo_rows(rows)
    vus_combos = vus_holding_combos(combos)
    path_combos = pathogenic_combos(combos)
    sanity = sanity_rows(rows)
    write_csv(
        TABLE_DIR / "criteria_counts_by_group.csv",
        counts,
        ["criterion", "count", "percent"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent"]],
    )
    write_csv(
        TABLE_DIR / "criteria_pair_cooccurrence.csv",
        pairs,
        ["criterion_a", "criterion_b", "count"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent_of_pair", "percent_of_group"]],
    )
    write_csv(
        TABLE_DIR / "criteria_combo_outcomes.csv",
        combos,
        ["criteria_combo", "count"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent_of_combo"]],
    )
    write_csv(
        TABLE_DIR / "vus_holding_combinations.csv",
        vus_combos,
        ["criteria_combo", "count"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent_of_combo"]],
    )
    write_csv(
        TABLE_DIR / "pathogenic_dominant_combinations.csv",
        path_combos,
        ["criteria_combo", "count"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent_of_combo"]],
    )
    write_csv(
        TABLE_DIR / "criteria_sanity_checks.csv",
        sanity,
        [
            "sanity_reason",
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "criteria_combo",
            "total_points",
            "predicted_class",
            "predicted_label",
            "spliceai_score",
        ],
    )
    draw_pair_heatmap(pairs)
    for group in GROUP_ORDER:
        draw_pair_heatmap(pairs, group)
    draw_combo_barplot(combos, "top_vus_holding_combinations.svg", "Top VUS-holding criterion combinations", vus_combos)
    draw_combo_barplot(combos, "top_pathogenic_combinations.svg", "Top pathogenic-dominant criterion combinations", path_combos)
    write_report(counts, pairs, combos, vus_combos, path_combos, sanity)
    print(f"Wrote criteria co-occurrence analysis to {REPORT}")


if __name__ == "__main__":
    main()
