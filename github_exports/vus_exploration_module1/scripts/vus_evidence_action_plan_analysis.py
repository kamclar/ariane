"""Build an evidence action plan for unresolved Module 1 VUS."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
VUS_PRIORITY = ANALYSIS_DIR / "tables" / "vus_prioritization" / "vus_priority_annotated.csv"
BOTTLENECK_EXAMPLES = ANALYSIS_DIR / "tables" / "vus_bottleneck" / "vus_bottleneck_variant_examples.csv"
OUT_DIR = ANALYSIS_DIR / "tables" / "vus_evidence_action_plan"
PLOT_DIR = ANALYSIS_DIR / "plots" / "29_vus_evidence_action_plan"
REPORT = ANALYSIS_DIR / "vus_evidence_action_plan_report.md"


ACTION_BY_BOTTLENECK = {
    "strong_pathogenic_combo_one_step_short": {
        "action_group": "near_pathogenic_threshold",
        "review_question": "Which single independent pathogenic evidence item would move this VUS out of VUS?",
        "candidate_evidence": "supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion",
    },
    "benign_functional_evidence_not_enough": {
        "action_group": "near_benign_threshold",
        "review_question": "Is there additional benign evidence or stronger functional benign evidence?",
        "candidate_evidence": "BS2, BS4, population evidence, or curated functional strength review",
    },
    "benign_evidence_not_enough": {
        "action_group": "near_benign_threshold",
        "review_question": "Is another benign criterion available?",
        "candidate_evidence": "BS2, BS4, population evidence, or ClinGen/ENIGMA assertion context",
    },
    "computational_pathogenic_evidence_not_enough": {
        "action_group": "computational_signal_needs_noncomputational_support",
        "review_question": "Does the computational signal have RNA, case-control, segregation, trans, or functional support?",
        "candidate_evidence": "RNA, functional, PS4, PP1, PM3, or curated assertion evidence",
    },
    "mixed_splice_benign_and_pathogenic": {
        "action_group": "conflicting_computational_splice_context",
        "review_question": "Do benign RNA branch and pathogenic splice prediction both truly apply?",
        "candidate_evidence": "manual RNA/splicing review and transcript-specific interpretation",
    },
    "PM2_only": {
        "action_group": "background_absence_only",
        "review_question": "Is there any independent evidence beyond absence from population data?",
        "candidate_evidence": "any independent variant-level evidence; otherwise low priority",
    },
    "no_automated_evidence": {
        "action_group": "module1_no_signal",
        "review_question": "Is there external evidence outside automated Module 1?",
        "candidate_evidence": "ClinVar/ClinGen/ENIGMA assertion context, literature, segregation, functional data",
    },
    "other_conflicting_evidence": {
        "action_group": "manual_conflict_adjudication",
        "review_question": "Which applied criteria are valid and which are conditional or overridden?",
        "candidate_evidence": "manual review of mixed benign and pathogenic evidence",
    },
    "conflicting_PVS1_BS3": {
        "action_group": "manual_conflict_adjudication",
        "review_question": "How should PVS1 and BS3 be reconciled for this exact molecular mechanism?",
        "candidate_evidence": "functional assay metadata, NMD/splicing review, ClinGen/ENIGMA assertion context",
    },
    "PM2_plus_weak_context": {
        "action_group": "background_absence_only",
        "review_question": "Does the weak context signal correspond to actual variant-level evidence?",
        "candidate_evidence": "stronger independent evidence; do not classify from context alone",
    },
}


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


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except ValueError:
        return default


def parse_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except ValueError:
        return default


def criterion_names(combo: str) -> list[str]:
    if combo in {"", "none"}:
        return []
    return [item.split(":")[0] for item in combo.split("+") if item]


def exon_key(row: dict[str, str]) -> str:
    if row.get("cds_exon"):
        return row["cds_exon"]
    try:
        pos = int(row.get("c_notation", "").split(".")[1].split(">")[0])
    except Exception:
        return ""
    return str((pos + 149) // 150)


def load_priority_index() -> dict[tuple[str, str], dict[str, str]]:
    return {(row["gene"], row["c_notation"]): row for row in read_csv(VUS_PRIORITY)}


def build_action_rows() -> list[dict[str, object]]:
    priority = load_priority_index()
    rows = []
    for row in read_csv(BOTTLENECK_EXAMPLES):
        key = (row["gene"], row["c_notation"])
        pr = priority.get(key, {})
        action = ACTION_BY_BOTTLENECK.get(
            row["bottleneck_category"],
            {
                "action_group": "other",
                "review_question": "Manual review needed.",
                "candidate_evidence": "case-specific evidence",
            },
        )
        criteria = criterion_names(row.get("criteria_combo", ""))
        rows.append(
            {
                **row,
                "priority_score": pr.get("priority_score", ""),
                "priority_tier": pr.get("priority_tier", ""),
                "review_category": pr.get("review_category", ""),
                "priority_reasons": pr.get("priority_reasons", ""),
                "cds_exon": pr.get("cds_exon", ""),
                "boundary_distance": pr.get("boundary_distance", ""),
                "near_pathogenic_20bp": pr.get("near_pathogenic_20bp", ""),
                "near_benign_20bp": pr.get("near_benign_20bp", ""),
                "near_vus_20bp": pr.get("near_vus_20bp", ""),
                "action_group": action["action_group"],
                "review_question": action["review_question"],
                "candidate_evidence": action["candidate_evidence"],
                "criteria_count": len(criteria),
                "has_ps3": "yes" if "PS3" in criteria else "no",
                "has_pp3": "yes" if "PP3" in criteria else "no",
                "has_bp7": "yes" if "BP7" in criteria else "no",
                "has_bp5": "yes" if "BP5" in criteria else "no",
                "has_pvs1": "yes" if "PVS1" in criteria else "no",
            }
        )
    return rows


def summarize(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(str(row.get(key, "")) for key in keys)].append(row)
    output = []
    for key_values, items in grouped.items():
        high_priority = [
            row for row in items if str(row.get("priority_tier", "")) in {"tier1_urgent", "tier2_high"}
        ]
        splice_high = [row for row in items if parse_float(str(row.get("spliceai_score", ""))) >= 0.20]
        out = {key: value for key, value in zip(keys, key_values)}
        out.update(
            {
                "count": len(items),
                "high_priority_count": len(high_priority),
                "high_spliceai_count": len(splice_high),
                "median_points": median_int([parse_int(str(row.get("total_points", ""))) for row in items]),
                "top_candidate_evidence": most_common(str(row.get("candidate_evidence", "")) for row in items),
            }
        )
        output.append(out)
    return sorted(output, key=lambda row: (-int(row["count"]), str(row)))


def most_common(values) -> str:
    counts = Counter(value for value in values if value)
    return counts.most_common(1)[0][0] if counts else ""


def median_int(values: list[int]) -> str:
    if not values:
        return ""
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return str(ordered[mid])
    return f"{(ordered[mid - 1] + ordered[mid]) / 2:.1f}"


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
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def draw_action_group_bar(summary: list[dict[str, object]]) -> None:
    width = 980
    height = 520
    left = 310
    top = 60
    bar_h = 34
    gap = 12
    max_count = max(int(row["count"]) for row in summary)
    colors = {
        "background_absence_only": "#94a3b8",
        "conflicting_computational_splice_context": "#f59e0b",
        "computational_signal_needs_noncomputational_support": "#6366f1",
        "near_pathogenic_threshold": "#dc2626",
        "near_benign_threshold": "#16a34a",
        "module1_no_signal": "#64748b",
        "manual_conflict_adjudication": "#7c3aed",
    }
    body = ['<text x="28" y="34" class="title">VUS evidence action groups</text>']
    for idx, row in enumerate(summary):
        y = top + idx * (bar_h + gap)
        group = str(row["action_group"])
        count = int(row["count"])
        w = 560 * count / max_count
        color = colors.get(group, "#334155")
        body.append(f'<text x="28" y="{y + 22}" class="small">{esc(group)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w:.1f}" height="{bar_h}" rx="5" fill="{color}" opacity="0.82"/>')
        body.append(f'<text x="{left + w + 10:.1f}" y="{y + 22}" class="small">{count}</text>')
    body.append('<text x="28" y="495" class="label">Counts are unresolved Module 1 VUS, grouped by the likely next evidence question.</text>')
    write_svg(PLOT_DIR / "vus_evidence_action_groups.svg", "\n".join(body), width, height)


def main() -> None:
    rows = build_action_rows()
    by_action = summarize(rows, ["action_group"])
    by_bottleneck = summarize(rows, ["bottleneck_category", "action_group"])
    by_gene_action = summarize(rows, ["gene", "action_group"])
    by_gene_exon_action = summarize(rows, ["gene", "cds_exon", "action_group"])
    high_priority = sorted(
        rows,
        key=lambda row: (
            {"tier1_urgent": 0, "tier2_high": 1, "tier3_medium": 2, "tier4_low": 3}.get(str(row.get("priority_tier", "")), 9),
            -parse_int(str(row.get("priority_score", ""))),
            str(row.get("gene", "")),
            str(row.get("c_notation", "")),
        ),
    )

    fields = [
        "action_group",
        "bottleneck_category",
        "candidate_evidence",
        "review_question",
        "priority_score",
        "priority_tier",
        "review_category",
        "priority_reasons",
        "gene",
        "c_notation",
        "p_notation",
        "variant_type",
        "criteria_combo",
        "total_points",
        "spliceai_score",
        "cds_exon",
        "boundary_distance",
        "near_pathogenic_20bp",
        "near_benign_20bp",
        "classification_note",
        "suggested_next_evidence",
    ]
    write_csv(OUT_DIR / "vus_evidence_action_plan_variants.csv", rows, fields)
    write_csv(
        OUT_DIR / "vus_evidence_action_plan_by_action_group.csv",
        by_action,
        ["action_group", "count", "high_priority_count", "high_spliceai_count", "median_points", "top_candidate_evidence"],
    )
    write_csv(
        OUT_DIR / "vus_evidence_action_plan_by_bottleneck.csv",
        by_bottleneck,
        ["bottleneck_category", "action_group", "count", "high_priority_count", "high_spliceai_count", "median_points", "top_candidate_evidence"],
    )
    write_csv(
        OUT_DIR / "vus_evidence_action_plan_by_gene.csv",
        by_gene_action,
        ["gene", "action_group", "count", "high_priority_count", "high_spliceai_count", "median_points", "top_candidate_evidence"],
    )
    write_csv(
        OUT_DIR / "vus_evidence_action_plan_by_gene_exon.csv",
        by_gene_exon_action,
        ["gene", "cds_exon", "action_group", "count", "high_priority_count", "high_spliceai_count", "median_points", "top_candidate_evidence"],
    )
    write_csv(OUT_DIR / "vus_evidence_action_plan_top_candidates.csv", high_priority[:250], fields)
    draw_action_group_bar(by_action)

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# VUS Evidence Action Plan

Generated: {generated}

## Purpose

This analysis converts unresolved ARIANE Module 1 VUS bottlenecks into a
practical evidence worklist. It does not add ACMG/ENIGMA points and it does
not treat local context as evidence. The aim is to ask: for this VUS group,
which kind of additional variant-level evidence would be most useful?

## Action Groups

{markdown_table(by_action, ["action_group", "count", "high_priority_count", "high_spliceai_count", "median_points", "top_candidate_evidence"])}

## Bottleneck To Action Mapping

{markdown_table(by_bottleneck, ["bottleneck_category", "action_group", "count", "high_priority_count", "high_spliceai_count", "top_candidate_evidence"])}

## Per Gene

{markdown_table(by_gene_action, ["gene", "action_group", "count", "high_priority_count", "high_spliceai_count", "top_candidate_evidence"])}

## Top Candidates

{markdown_table(high_priority, ["priority_score", "priority_tier", "action_group", "bottleneck_category", "gene", "c_notation", "p_notation", "criteria_combo", "spliceai_score", "candidate_evidence"], limit=30)}

## Interpretation

The largest group is usually `background_absence_only`, meaning PM2 alone or
PM2 with weak context. These variants are not the most efficient first manual
review target unless external public assertions, literature, or laboratory data
already exist.

The most actionable groups are:

1. `near_pathogenic_threshold`, where one independent pathogenic evidence item
   may be enough to move the case from VUS toward likely pathogenic.
2. `near_benign_threshold`, where stronger benign evidence or an additional
   benign criterion may resolve the case.
3. `computational_signal_needs_noncomputational_support`, where SpliceAI/PP3
   points to a hypothesis but cannot classify the variant by itself.

## Outputs

- `tables/vus_evidence_action_plan/vus_evidence_action_plan_variants.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_action_group.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_bottleneck.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_gene.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_gene_exon.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_top_candidates.csv`
- `plots/29_vus_evidence_action_plan/vus_evidence_action_groups.svg`
"""
    REPORT.write_text(text, encoding="utf-8")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
