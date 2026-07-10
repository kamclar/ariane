"""Compare BRCA1 RING/BRCT VUS tiers with continuous Findlay SGE scores."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
VUS_PRIORITY = ANALYSIS_DIR / "tables" / "vus_prioritization" / "vus_priority_annotated.csv"
FINDLAY = (
    ANALYSIS_DIR
    / "external_sources"
    / "brca1_bs3_functional_evidence_audit_sources"
    / "findlay_2018_brca1_sge"
    / "Findlay_2018_Nature_Supplementary_Table_1_SGE_scores.normalized.csv"
)
ACTION_PLAN = ANALYSIS_DIR / "tables" / "vus_evidence_action_plan" / "vus_evidence_action_plan_variants.csv"
TABLE_DIR = ANALYSIS_DIR / "tables" / "findlay_sge_vus_tier"
PLOT_DIR = ANALYSIS_DIR / "plots" / "28_findlay_sge_vus_tier"
REPORT = ANALYSIS_DIR / "findlay_sge_vus_tier_report.md"


NEIGHBORHOOD_REASONS = {"near_pathogenic_dense", "high_spliceai_in_benign_dense_region"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_float(value: str) -> float | None:
    if value in {"", "NA", "nan", "None"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def fmt_float(value: float | None) -> str:
    return "" if value is None else f"{value:.3f}"


def brca1_region(cds_pos: str) -> str:
    try:
        pos = int(cds_pos)
    except ValueError:
        return "unknown"
    aa = (pos + 2) // 3
    if 1 <= aa <= 109:
        return "N_terminal_RING_region"
    if 1642 <= aa <= 1855:
        return "BRCT_region"
    return "outside_Findlay_Ring_BRCT_core"


def tier_order(tier: str) -> int:
    order = {"tier1_urgent": 0, "tier2_high": 1, "tier3_medium": 2, "tier4_low": 3}
    return order.get(tier, 99)


def load_findlay_by_variant() -> dict[str, dict[str, str]]:
    rows = read_csv(FINDLAY)
    return {row["transcript_variant"]: row for row in rows if row.get("gene") == "BRCA1"}


def annotate_vus() -> list[dict[str, object]]:
    findlay = load_findlay_by_variant()
    rows = []
    for row in read_csv(VUS_PRIORITY):
        if row["gene"] != "BRCA1":
            continue
        source = findlay.get(row["c_notation"])
        if not source:
            continue
        reasons = {item for item in row["priority_reasons"].split(";") if item}
        function_score = parse_float(source.get("function.score.mean", ""))
        rna_score = parse_float(source.get("mean.rna.score", ""))
        rows.append(
            {
                **row,
                "findlay_function_score_mean": fmt_float(function_score),
                "findlay_mean_rna_score": fmt_float(rna_score),
                "findlay_func_class": source.get("func.class", ""),
                "findlay_consequence": source.get("consequence", ""),
                "findlay_region": brca1_region(row["cds_pos"]),
                "neighborhood_driven": "yes" if reasons & NEIGHBORHOOD_REASONS else "no",
                "neighborhood_only": "yes" if reasons and reasons <= NEIGHBORHOOD_REASONS else "no",
            }
        )
    return rows


def summarize(rows: list[dict[str, object]], group_key: str) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[group_key])].append(row)
    output = []
    for group, items in grouped.items():
        scores = [
            float(row["findlay_function_score_mean"])
            for row in items
            if row["findlay_function_score_mean"] != ""
        ]
        classes = Counter(str(row["findlay_func_class"] or "missing") for row in items)
        output.append(
            {
                group_key: group,
                "count": len(items),
                "median_function_score": fmt_float(median(scores)),
                "mean_function_score": fmt_float(sum(scores) / len(scores) if scores else None),
                "lof_count": classes["LOF"],
                "int_count": classes["INT"],
                "func_count": classes["FUNC"],
                "missing_func_class_count": classes["missing"],
            }
        )
    return sorted(output, key=lambda row: (tier_order(str(row.get(group_key, ""))), str(row.get(group_key, ""))))


def summarize_lof_gap(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    lof_rows = [row for row in rows if row["findlay_func_class"] == "LOF"]
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in lof_rows:
        grouped[str(row["priority_tier"])].append(row)
    output = []
    for tier, items in grouped.items():
        with_ps3 = [row for row in items if "PS3" in str(row["criteria"])]
        with_bs3 = [row for row in items if "BS3" in str(row["criteria"])]
        with_splice = [row for row in items if "SpliceAI>=" in str(row["priority_reasons"])]
        with_boundary = [row for row in items if "boundary" in str(row["priority_reasons"])]
        with_neighborhood = [row for row in items if "near_pathogenic_dense" in str(row["priority_reasons"])]
        output.append(
            {
                "priority_tier": tier,
                "lof_count": len(items),
                "with_ps3": len(with_ps3),
                "with_bs3": len(with_bs3),
                "with_splice_priority_reason": len(with_splice),
                "with_boundary_priority_reason": len(with_boundary),
                "with_pathogenic_neighborhood_reason": len(with_neighborhood),
                "top_criteria_combo": most_common(str(row["criteria"]) for row in items),
                "top_priority_reasons": most_common(str(row["priority_reasons"]) for row in items),
            }
        )
    return sorted(output, key=lambda row: tier_order(str(row["priority_tier"])))


def summarize_lof_action_overlap(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    if not ACTION_PLAN.exists():
        return []
    action_index = {
        (row["gene"], row["c_notation"]): row
        for row in read_csv(ACTION_PLAN)
    }
    tier3_lof = [
        row for row in rows
        if row["priority_tier"] == "tier3_medium" and row["findlay_func_class"] == "LOF"
    ]
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in tier3_lof:
        action = action_index.get((str(row["gene"]), str(row["c_notation"])), {})
        grouped[
            (
                action.get("action_group", "missing_action_plan"),
                action.get("bottleneck_category", "missing_action_plan"),
            )
        ].append(row)
    output = []
    for (action_group, bottleneck_category), items in grouped.items():
        output.append(
            {
                "action_group": action_group,
                "bottleneck_category": bottleneck_category,
                "count": len(items),
                "top_criteria_combo": most_common(str(row["criteria"]) for row in items),
                "top_priority_reasons": most_common(str(row["priority_reasons"]) for row in items),
            }
        )
    return sorted(output, key=lambda row: -int(row["count"]))


def most_common(values) -> str:
    counts = Counter(value for value in values if value)
    return counts.most_common(1)[0][0] if counts else ""


def markdown_table(rows: list[dict[str, object]], columns: list[str], limit: int | None = None) -> str:
    selected = rows if limit is None else rows[:limit]
    if not selected:
        return "| none |\n"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in selected:
        lines.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def write_svg(path: Path, body: str, width: int, height: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
text {{ font-family: Arial, Helvetica, sans-serif; fill: #172033; }}
.title {{ font-size: 18px; font-weight: 700; }}
.label {{ font-size: 12px; }}
.small {{ font-size: 11px; }}
.grid {{ stroke: #e2e8f0; stroke-width: 1; }}
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def draw_tier_score_plot(rows: list[dict[str, object]]) -> None:
    width = 900
    height = 420
    left = 90
    top = 60
    plot_w = 650
    plot_h = 270
    tiers = ["tier1_urgent", "tier2_high", "tier3_medium", "tier4_low"]
    colors = {"LOF": "#dc2626", "INT": "#f59e0b", "FUNC": "#16a34a", "": "#64748b"}
    scores = [float(row["findlay_function_score_mean"]) for row in rows if row["findlay_function_score_mean"] != ""]
    min_y = min(scores + [-4.0])
    max_y = max(scores + [1.0])
    pad = 0.15 * (max_y - min_y)
    min_y -= pad
    max_y += pad

    def y_pos(score: float) -> float:
        return top + plot_h - (score - min_y) / (max_y - min_y) * plot_h

    body = [
        '<text x="28" y="34" class="title">Findlay SGE function score by VUS priority tier</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>',
    ]
    for value in [-3, -2, -1, 0, 1]:
        if min_y <= value <= max_y:
            y = y_pos(value)
            body.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>')
            body.append(f'<text x="{left - 34}" y="{y + 4:.1f}" class="small">{value}</text>')
    for idx, tier in enumerate(tiers):
        x_center = left + 85 + idx * 155
        tier_rows = [row for row in rows if row["priority_tier"] == tier]
        for j, row in enumerate(tier_rows):
            score_text = row["findlay_function_score_mean"]
            if score_text == "":
                continue
            score = float(score_text)
            jitter = ((j % 17) - 8) * 2.0
            color = colors.get(str(row["findlay_func_class"]), "#64748b")
            stroke = "#111827" if row["neighborhood_driven"] == "yes" else "none"
            body.append(
                f'<circle cx="{x_center + jitter:.1f}" cy="{y_pos(score):.1f}" r="3.2" fill="{color}" opacity="0.62" stroke="{stroke}" stroke-width="0.6"/>'
            )
        body.append(f'<text x="{x_center - 43}" y="{top + plot_h + 24}" class="small">{tier}</text>')
    legend_x = left + plot_w + 28
    legend_y = top + 24
    body.append(f'<text x="{legend_x}" y="{legend_y}" class="label">Findlay class</text>')
    for i, label in enumerate(["LOF", "INT", "FUNC"]):
        y = legend_y + 28 + i * 24
        body.append(f'<circle cx="{legend_x + 8}" cy="{y}" r="5" fill="{colors[label]}"/>')
        body.append(f'<text x="{legend_x + 22}" y="{y + 4}" class="small">{label}</text>')
    body.append(f'<circle cx="{legend_x + 8}" cy="{legend_y + 112}" r="5" fill="#ffffff" stroke="#111827"/>')
    body.append(f'<text x="{legend_x + 22}" y="{legend_y + 116}" class="small">neighborhood reason</text>')
    body.append(f'<text x="{left + 210}" y="{height - 28}" class="label">Each dot is one BRCA1 VUS matched to Findlay 2018 SGE.</text>')
    body.append(f'<text x="{left - 68}" y="{top + 120}" class="label" transform="rotate(-90 {left - 68} {top + 120})">function.score.mean</text>')
    write_svg(PLOT_DIR / "findlay_sge_score_by_priority_tier.svg", "\n".join(body), width, height)


def main() -> None:
    rows = annotate_vus()
    tier_summary = summarize(rows, "priority_tier")
    region_summary = summarize(rows, "findlay_region")
    neighborhood_summary = summarize(rows, "neighborhood_driven")
    lof_gap_summary = summarize_lof_gap(rows)
    lof_action_overlap = summarize_lof_action_overlap(rows)
    top_lof = sorted(
        [row for row in rows if row["findlay_func_class"] == "LOF"],
        key=lambda row: (tier_order(str(row["priority_tier"])), float(row["findlay_function_score_mean"] or 0)),
    )

    fields = [
        "priority_score",
        "priority_tier",
        "review_category",
        "priority_reasons",
        "neighborhood_driven",
        "neighborhood_only",
        "gene",
        "c_notation",
        "p_notation",
        "cds_pos",
        "findlay_region",
        "variant_type",
        "spliceai_score",
        "findlay_func_class",
        "findlay_function_score_mean",
        "findlay_mean_rna_score",
        "findlay_consequence",
        "criteria",
        "total_points",
    ]
    write_csv(TABLE_DIR / "brca1_vus_with_findlay_sge_scores.csv", rows, fields)
    write_csv(
        TABLE_DIR / "findlay_sge_by_priority_tier.csv",
        tier_summary,
        ["priority_tier", "count", "median_function_score", "mean_function_score", "lof_count", "int_count", "func_count", "missing_func_class_count"],
    )
    write_csv(
        TABLE_DIR / "findlay_sge_by_region.csv",
        region_summary,
        ["findlay_region", "count", "median_function_score", "mean_function_score", "lof_count", "int_count", "func_count", "missing_func_class_count"],
    )
    write_csv(
        TABLE_DIR / "findlay_sge_by_neighborhood_driven.csv",
        neighborhood_summary,
        ["neighborhood_driven", "count", "median_function_score", "mean_function_score", "lof_count", "int_count", "func_count", "missing_func_class_count"],
    )
    write_csv(
        TABLE_DIR / "findlay_sge_lof_gap_by_tier.csv",
        lof_gap_summary,
        [
            "priority_tier",
            "lof_count",
            "with_ps3",
            "with_bs3",
            "with_splice_priority_reason",
            "with_boundary_priority_reason",
            "with_pathogenic_neighborhood_reason",
            "top_criteria_combo",
            "top_priority_reasons",
        ],
    )
    write_csv(
        TABLE_DIR / "findlay_sge_tier3_lof_action_overlap.csv",
        lof_action_overlap,
        ["action_group", "bottleneck_category", "count", "top_criteria_combo", "top_priority_reasons"],
    )
    write_csv(TABLE_DIR / "findlay_sge_lof_vus_candidates.csv", top_lof, fields)
    draw_tier_score_plot(rows)

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# Findlay SGE Continuous Score Vs VUS Tier Analysis

Generated: {generated}

## Purpose

This analysis compares ARIANE Module 1 VUS priority tiers with the continuous
Findlay et al. 2018 BRCA1 saturation genome editing `function.score.mean`.
This is intentionally different from comparing against the derived `BS3` label:
the continuous score is used as an external functional axis where available.

Scope is limited to BRCA1 regions covered by Findlay 2018, mainly RING and
BRCT. Lack of comparable dense BRCA2 functional data is itself a limitation.

## Matched VUS Count

Matched BRCA1 VUS with Findlay SGE rows: `{len(rows)}`

## By Priority Tier

{markdown_table(tier_summary, ["priority_tier", "count", "median_function_score", "mean_function_score", "lof_count", "int_count", "func_count"])}

## LOF Distribution And Priority Gap

{markdown_table(lof_gap_summary, ["priority_tier", "lof_count", "with_ps3", "with_bs3", "with_splice_priority_reason", "with_boundary_priority_reason", "with_pathogenic_neighborhood_reason", "top_criteria_combo", "top_priority_reasons"])}

The Findlay SGE axis does not support a monotonic interpretation of the current
priority tiers. Tier4 is clearly more function-preserved or benign-leaning, but
tier3 contains more LOF-class variants than tier1 and tier2 combined. In this
snapshot, most tier3 LOF variants already have `PS3`; they remain tier3 because
`PS3` alone, or `PS3+PM2_Supporting` with weak context, does not reach the
tier1/tier2 review threshold unless splice, boundary, PS1, PP4, or stronger
neighborhood signals are also present.

This is a useful critique of the current triage heuristic: it prioritizes
splice/boundary/context combinations, and may under-prioritize functionally
damaging missense VUS that have PS3 but lack additional review signals.

## Tier3 LOF Overlap With Evidence Action Plan

{markdown_table(lof_action_overlap, ["action_group", "bottleneck_category", "count", "top_criteria_combo", "top_priority_reasons"])}

This connects the Findlay SGE axis to the action-plan analysis. The tier3 LOF
group is not a separate story: almost all of it overlaps with the
`near_pathogenic_threshold` / `strong_pathogenic_combo_one_step_short` queue.
This means review priority and distance-to-reclassification are related but not
identical. A variant can be close to likely pathogenic under the evidence model
while remaining only tier3 in the heuristic priority score if its signal is
mainly `PS3+PM2_Supporting`.

## By Region

{markdown_table(region_summary, ["findlay_region", "count", "median_function_score", "mean_function_score", "lof_count", "int_count", "func_count"])}

## By Neighborhood-Driven Status

{markdown_table(neighborhood_summary, ["neighborhood_driven", "count", "median_function_score", "mean_function_score", "lof_count", "int_count", "func_count"])}

## LOF-Class VUS Candidates

{markdown_table(top_lof, ["priority_score", "priority_tier", "review_category", "priority_reasons", "gene", "c_notation", "p_notation", "findlay_region", "findlay_function_score_mean", "findlay_mean_rna_score", "spliceai_score", "criteria"], limit=30)}

## Interpretation Boundary

Findlay SGE is not used as a clinical reference set here. It is a functional
measurement axis for a restricted BRCA1 subset. It is also not fully
independent of the generated Module 1 map, because BRCA1 RING/BRCT `PS3` and
`BS3` evidence in the local ENIGMA Table 9 lookup is partly derived from the
same Findlay SGE source. The least circular subset is therefore VUS with a
Findlay score but without applied `PS3` or `BS3`.

The main supported interpretation is coarse: tier4 is function-preserved or
benign-leaning compared with the other tiers. The data do not validate the fine
ordering tier1 over tier2 over tier3, and small counts in tier1 and tier2 make
their medians unstable.

## Outputs

- `tables/findlay_sge_vus_tier/brca1_vus_with_findlay_sge_scores.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_by_priority_tier.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_by_region.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_by_neighborhood_driven.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_lof_gap_by_tier.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_tier3_lof_action_overlap.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_lof_vus_candidates.csv`
- `plots/28_findlay_sge_vus_tier/findlay_sge_score_by_priority_tier.svg`
"""
    REPORT.write_text(text, encoding="utf-8")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
