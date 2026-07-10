"""Analyses responding to internal critique of the VUS mini manuscript."""

from __future__ import annotations

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ANALYSIS_DIR / "tables" / "vus_prioritization" / "vus_priority_annotated.csv"
BOTTLENECK = ANALYSIS_DIR / "tables" / "vus_bottleneck" / "vus_bottleneck_summary.csv"
OUT_DIR = ANALYSIS_DIR / "tables" / "vus_manuscript_critique_response"
PLOT_DIR = ANALYSIS_DIR / "plots" / "27_vus_manuscript_critique_response"
REPORT = ANALYSIS_DIR / "vus_manuscript_critique_response_report.md"

NEIGHBORHOOD_WEIGHTS = {
    "near_pathogenic_dense": 12,
    "high_spliceai_in_benign_dense_region": 8,
}


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def tier(score: int) -> str:
    if score >= 80:
        return "tier1_urgent"
    if score >= 50:
        return "tier2_high"
    if score >= 25:
        return "tier3_medium"
    return "tier4_low"


def high_tier(value: str) -> bool:
    return value in {"tier1_urgent", "tier2_high"}


def load_priority_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with INPUT.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            reasons = [item for item in row["priority_reasons"].split(";") if item]
            score = int(row["priority_score"])
            neighborhood_points = sum(NEIGHBORHOOD_WEIGHTS.get(reason, 0) for reason in reasons)
            score_without_neighborhood = score - neighborhood_points
            row["reasons_list"] = reasons
            row["score_int"] = score
            row["neighborhood_points"] = neighborhood_points
            row["score_without_neighborhood"] = score_without_neighborhood
            row["tier_without_neighborhood"] = tier(score_without_neighborhood)
            row["only_neighborhood_reason"] = all(reason in NEIGHBORHOOD_WEIGHTS for reason in reasons)
            row["has_neighborhood_reason"] = neighborhood_points > 0
            row["lost_high_priority_without_neighborhood"] = (
                high_tier(row["priority_tier"]) and not high_tier(str(row["tier_without_neighborhood"]))
            )
            rows.append(row)
    return rows


def load_bottleneck_rows() -> list[dict[str, str]]:
    if not BOTTLENECK.exists():
        return []
    with BOTTLENECK.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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


def markdown_table(rows: list[dict[str, object]], columns: list[str], limit: int | None = None) -> str:
    selected = rows if limit is None else rows[:limit]
    if not selected:
        return "| none |\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = ["| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |" for row in selected]
    return "\n".join([header, sep] + body)


def summarize_neighborhood_increment(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    total = len(rows)
    high_original = sum(1 for row in rows if high_tier(str(row["priority_tier"])))
    high_without = sum(1 for row in rows if high_tier(str(row["tier_without_neighborhood"])))
    with_neighborhood = sum(1 for row in rows if row["has_neighborhood_reason"])
    only_neighborhood = sum(1 for row in rows if row["only_neighborhood_reason"])
    lost_high = sum(1 for row in rows if row["lost_high_priority_without_neighborhood"])
    return [
        {"metric": "total_vus", "count": total, "percent_of_vus": "100.00"},
        {
            "metric": "vus_with_any_neighborhood_reason",
            "count": with_neighborhood,
            "percent_of_vus": f"{with_neighborhood / total * 100:.2f}",
        },
        {
            "metric": "vus_with_only_neighborhood_reasons",
            "count": only_neighborhood,
            "percent_of_vus": f"{only_neighborhood / total * 100:.2f}",
        },
        {
            "metric": "tier1_or_tier2_original",
            "count": high_original,
            "percent_of_vus": f"{high_original / total * 100:.2f}",
        },
        {
            "metric": "tier1_or_tier2_without_neighborhood_points",
            "count": high_without,
            "percent_of_vus": f"{high_without / total * 100:.2f}",
        },
        {
            "metric": "lost_tier1_or_tier2_without_neighborhood_points",
            "count": lost_high,
            "percent_of_vus": f"{lost_high / total * 100:.2f}",
        },
    ]


def tier_transition_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts = Counter((row["priority_tier"], row["tier_without_neighborhood"]) for row in rows)
    order = ["tier1_urgent", "tier2_high", "tier3_medium", "tier4_low"]
    output = []
    for original in order:
        for recalculated in order:
            count = counts[(original, recalculated)]
            if count:
                output.append(
                    {
                        "original_tier": original,
                        "tier_without_neighborhood": recalculated,
                        "count": count,
                    }
                )
    return output


def neighborhood_only_examples(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    selected = [row for row in rows if row["only_neighborhood_reason"] or row["lost_high_priority_without_neighborhood"]]
    return sorted(
        selected,
        key=lambda row: (
            -int(row["score_int"]),
            -int(row["neighborhood_points"]),
            row["gene"],
            row["c_notation"],
        ),
    )


def draw_tier_transition(rows: list[dict[str, object]]) -> None:
    order = ["tier1_urgent", "tier2_high", "tier3_medium", "tier4_low"]
    counts = Counter((row["priority_tier"], row["tier_without_neighborhood"]) for row in rows)
    width = 820
    height = 420
    left = 170
    top = 90
    cell = 62
    max_count = max(counts.values(), default=1)
    body = [
        '<text x="26" y="34" class="title">Priority tier change after removing neighborhood points</text>',
        '<text x="26" y="56" class="small">Rows: original tier. Columns: recalculated tier without local-neighborhood points.</text>',
    ]
    for i, original in enumerate(order):
        body.append(f'<text x="24" y="{top + i * cell + 38}" class="small">{original}</text>')
    for j, recalculated in enumerate(order):
        body.append(
            f'<text x="{left + j * cell + 6}" y="{top - 16}" class="small" transform="rotate(-35 {left + j * cell + 6} {top - 16})">{recalculated}</text>'
        )
    for i, original in enumerate(order):
        for j, recalculated in enumerate(order):
            count = counts[(original, recalculated)]
            shade = int(245 - 175 * count / max_count) if max_count else 245
            color = f"rgb(220,{shade},{shade})" if original != recalculated else f"rgb({shade},{shade},245)"
            x = left + j * cell
            y = top + i * cell
            body.append(f'<rect x="{x}" y="{y}" width="{cell - 4}" height="{cell - 4}" fill="{color}" stroke="#cbd5e1"/>')
            if count:
                body.append(f'<text x="{x + 12}" y="{y + 35}" class="label">{count}</text>')
    write_svg(PLOT_DIR / "tier_transition_without_neighborhood.svg", "\n".join(body), width, height)


def draw_neighborhood_increment(summary: list[dict[str, object]]) -> None:
    rows = [row for row in summary if row["metric"] != "total_vus"]
    width = 900
    height = 360
    left = 310
    top = 70
    bar_h = 28
    gap = 13
    plot_w = 500
    max_count = max(int(row["count"]) for row in rows) if rows else 1
    body = ['<text x="26" y="34" class="title">Neighborhood contribution to VUS priority</text>']
    for i, row in enumerate(rows):
        y = top + i * (bar_h + gap)
        count = int(row["count"])
        w = plot_w * count / max_count if max_count else 0
        body.append(f'<text x="26" y="{y + 19}" class="small">{esc(row["metric"])}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="#7c3aed"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 19}" class="small">{count} ({row["percent_of_vus"]}%)</text>')
    write_svg(PLOT_DIR / "neighborhood_increment_summary.svg", "\n".join(body), width, height)


def bottleneck_note_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    total = sum(int(row["count"]) for row in rows) if rows else 0
    manuscript_major = {
        "PM2_only",
        "mixed_splice_benign_and_pathogenic",
        "computational_pathogenic_evidence_not_enough",
        "strong_pathogenic_combo_one_step_short",
        "no_automated_evidence",
        "benign_functional_evidence_not_enough",
    }
    shown = sum(int(row["count"]) for row in rows if row["bottleneck_category"] in manuscript_major)
    return [
        {"metric": "all_bottleneck_categories_total", "count": total},
        {"metric": "categories_shown_in_manuscript_table", "count": shown},
        {"metric": "other_categories_not_shown_in_short_table", "count": total - shown},
    ]


def write_report(
    summary: list[dict[str, object]],
    transitions: list[dict[str, object]],
    examples: list[dict[str, object]],
    bottleneck_notes: list[dict[str, object]],
) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# VUS Manuscript Critique Response Analysis

Generated: {generated}

Input: `{INPUT.relative_to(ROOT)}`

## Purpose

This analysis responds to manuscript critique about possible circularity and
redundancy of the local-neighborhood signal. It asks what happens if the two
local-neighborhood priority components are removed from the VUS priority score:

- `near_pathogenic_dense`: 12 points
- `high_spliceai_in_benign_dense_region`: 8 points

This does not validate the priority score externally. It only quantifies how
much the current tiering depends on local-neighborhood context.

## Neighborhood Contribution

{markdown_table(summary, ["metric", "count", "percent_of_vus"])}

## Tier Transitions

{markdown_table(transitions, ["original_tier", "tier_without_neighborhood", "count"])}

## Examples Affected By Neighborhood Removal

{markdown_table(examples, ["priority_score", "score_without_neighborhood", "priority_tier", "tier_without_neighborhood", "neighborhood_points", "priority_reasons", "gene", "c_notation", "p_notation"], limit=25)}

## Bottleneck Table Note

The manuscript bottleneck table intentionally shows the largest bottleneck
groups, not all mutually exclusive categories. The complete bottleneck table is
available in `tables/vus_bottleneck/vus_bottleneck_summary.csv`.

{markdown_table(bottleneck_notes, ["metric", "count"])}

## Interpretation

- The neighborhood signal should be described as a secondary triage and
  explanation signal, not as an independent classifier.
- Any VUS that remains tier1/tier2 after removing neighborhood points is
  prioritized mainly by variant-level criteria such as SpliceAI, boundary
  proximity, PS3, PP3, PS1, PVS1, PP4, or BS3.
- VUS that drop out of tier1/tier2 after removing neighborhood points should be
  labelled as neighborhood-driven candidates and should be interpreted more
  cautiously.
- External validation against ENIGMA, ClinVar high-review variants, EQA cases,
  or local expert review remains required before the score can be presented as
  a validated prioritization method.

## Outputs

- `tables/vus_manuscript_critique_response/neighborhood_increment_summary.csv`
- `tables/vus_manuscript_critique_response/tier_transition_without_neighborhood.csv`
- `tables/vus_manuscript_critique_response/neighborhood_affected_examples.csv`
- `tables/vus_manuscript_critique_response/bottleneck_short_table_note.csv`
- `plots/27_vus_manuscript_critique_response/neighborhood_increment_summary.svg`
- `plots/27_vus_manuscript_critique_response/tier_transition_without_neighborhood.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_priority_rows()
    bottleneck_rows = load_bottleneck_rows()
    summary = summarize_neighborhood_increment(rows)
    transitions = tier_transition_rows(rows)
    examples = neighborhood_only_examples(rows)
    bottleneck_notes = bottleneck_note_rows(bottleneck_rows)

    write_csv(
        OUT_DIR / "neighborhood_increment_summary.csv",
        summary,
        ["metric", "count", "percent_of_vus"],
    )
    write_csv(
        OUT_DIR / "tier_transition_without_neighborhood.csv",
        transitions,
        ["original_tier", "tier_without_neighborhood", "count"],
    )
    write_csv(
        OUT_DIR / "neighborhood_affected_examples.csv",
        examples,
        [
            "priority_score",
            "score_without_neighborhood",
            "priority_tier",
            "tier_without_neighborhood",
            "neighborhood_points",
            "priority_reasons",
            "gene",
            "c_notation",
            "p_notation",
            "spliceai_score",
            "criteria",
            "near_pathogenic_20bp",
            "near_benign_20bp",
        ],
    )
    write_csv(OUT_DIR / "bottleneck_short_table_note.csv", bottleneck_notes, ["metric", "count"])

    draw_neighborhood_increment(summary)
    draw_tier_transition(rows)
    write_report(summary, transitions, examples, bottleneck_notes)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
