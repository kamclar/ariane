"""Summarize VUS prioritization as overlapping review queues."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
PRIORITY = ANALYSIS_DIR / "tables" / "vus_prioritization" / "vus_priority_annotated.csv"
ACTION_PLAN = ANALYSIS_DIR / "tables" / "vus_evidence_action_plan" / "vus_evidence_action_plan_variants.csv"
PUBLIC = ANALYSIS_DIR / "tables" / "public_classification_snapshot" / "public_classification_snapshot_variants.csv"
FINDLAY = ANALYSIS_DIR / "tables" / "findlay_sge_vus_tier" / "brca1_vus_with_findlay_sge_scores.csv"
OUT_DIR = ANALYSIS_DIR / "tables" / "priority_queue_synthesis"
PLOT_DIR = ANALYSIS_DIR / "plots" / "30_priority_queue_synthesis"
REPORT = ANALYSIS_DIR / "priority_queue_synthesis_report.md"


QUEUE_COLUMNS = [
    "near_threshold_queue",
    "splice_rna_queue",
    "functional_queue",
    "evidence_conflict_queue",
    "public_assertion_queue",
    "neighborhood_context_queue",
    "low_information_queue",
]

QUEUE_LABELS = {
    "near_threshold_queue": "near threshold",
    "splice_rna_queue": "splice/RNA",
    "functional_queue": "functional",
    "evidence_conflict_queue": "evidence conflict",
    "public_assertion_queue": "public assertion",
    "neighborhood_context_queue": "neighborhood context",
    "low_information_queue": "low information",
}

PATHOGENIC_CRITERIA = {"PVS1", "PS1", "PS3", "PP3", "PP4", "PM5_PTC"}
BENIGN_CRITERIA = {"BA1", "BS1_Supporting", "BS1_Strong", "BS3", "BP1", "BP4", "BP5", "BP7"}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def key(row: dict[str, str]) -> tuple[str, str]:
    return row.get("gene", ""), row.get("c_notation", "")


def criterion_names(criteria: str) -> set[str]:
    names: set[str] = set()
    for item in criteria.split(";"):
        if item:
            names.add(item.split(":")[0])
    return names


def yes(value: bool) -> str:
    return "yes" if value else "no"


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def draw_bar_svg(path: Path, rows: list[dict[str, object]]) -> None:
    width, height = 980, 440
    margin_left, margin_right, margin_top, margin_bottom = 210, 40, 50, 80
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    max_value = max(int(row["count"]) for row in rows) if rows else 1
    bar_h = plot_h / max(len(rows), 1) * 0.62
    gap = plot_h / max(len(rows), 1)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033} .small{font-size:12px;fill:#43506a} .axis{stroke:#d7deea} .bar{fill:#5578c2}</style>',
        f'<text x="{margin_left}" y="28" font-size="22" font-weight="700">VUS review queues are overlapping</text>',
        f'<text x="{margin_left}" y="48" class="small">A single VUS may belong to several manual-review queues.</text>',
    ]
    for tick in range(0, max_value + 1, max(1, round(max_value / 5))):
        x = margin_left + plot_w * tick / max_value
        parts.append(f'<line x1="{x:.1f}" y1="{margin_top}" x2="{x:.1f}" y2="{height - margin_bottom}" class="axis"/>')
        parts.append(f'<text x="{x:.1f}" y="{height - margin_bottom + 24}" text-anchor="middle" class="small">{tick}</text>')

    for idx, row in enumerate(rows):
        y = margin_top + idx * gap + (gap - bar_h) / 2
        value = int(row["count"])
        bar_w = plot_w * value / max_value
        label = esc(row["queue_label"])
        parts.append(f'<text x="{margin_left - 12}" y="{y + bar_h * 0.68:.1f}" text-anchor="end" font-size="14">{label}</text>')
        parts.append(f'<rect x="{margin_left}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" rx="4" class="bar"/>')
        parts.append(f'<text x="{margin_left + bar_w + 8:.1f}" y="{y + bar_h * 0.68:.1f}" font-size="13">{value}</text>')
    parts.append(f'<text x="{margin_left + plot_w / 2}" y="{height - 20}" text-anchor="middle" class="small">number of VUS in queue</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def draw_overlap_svg(path: Path, rows: list[dict[str, object]]) -> None:
    labels = [QUEUE_LABELS[col] for col in QUEUE_COLUMNS]
    values = {(row["queue_a"], row["queue_b"]): int(row["count"]) for row in rows}
    width, height = 900, 780
    left, top, cell = 230, 125, 72
    max_value = max(values.values()) if values else 1
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.small{font-size:12px;fill:#43506a}.grid{stroke:#d7deea;stroke-width:1}</style>',
        f'<text x="{left}" y="36" font-size="22" font-weight="700">Priority queue overlap</text>',
        f'<text x="{left}" y="58" class="small">Diagonal shows queue size; off-diagonal cells show VUS shared by two queues.</text>',
    ]
    for i, label in enumerate(labels):
        x = left + i * cell + cell / 2
        parts.append(f'<text x="{x:.1f}" y="112" text-anchor="middle" class="small" transform="rotate(-45 {x:.1f} 112)">{esc(label)}</text>')
        y = top + i * cell + cell / 2 + 4
        parts.append(f'<text x="{left - 10}" y="{y:.1f}" text-anchor="end" class="small">{esc(label)}</text>')
    for i, col_a in enumerate(QUEUE_COLUMNS):
        for j, col_b in enumerate(QUEUE_COLUMNS):
            value = values.get((col_a, col_b), 0)
            intensity = value / max_value
            fill = f"rgb({int(245 - 160 * intensity)}, {int(248 - 120 * intensity)}, {int(255 - 20 * intensity)})"
            x = left + j * cell
            y = top + i * cell
            parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" class="grid"/>')
            if value:
                parts.append(f'<text x="{x + cell / 2:.1f}" y="{y + cell / 2 + 5:.1f}" text-anchor="middle" font-size="13">{value}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    priority_rows = read_csv(PRIORITY)
    action_rows = {key(row): row for row in read_csv(ACTION_PLAN)}
    public_rows = {key(row): row for row in read_csv(PUBLIC)}
    findlay_rows = {key(row): row for row in read_csv(FINDLAY)}

    synthesized: list[dict[str, object]] = []
    for row in priority_rows:
        k = key(row)
        action = action_rows.get(k, {})
        public = public_rows.get(k, {})
        findlay = findlay_rows.get(k, {})
        criteria = criterion_names(row.get("criteria", ""))
        action_group = action.get("action_group", "")
        bottleneck = action.get("bottleneck_category", "")
        reasons = row.get("priority_reasons", "")
        spliceai = parse_float(row.get("spliceai_score", "0"))
        boundary = parse_int(row.get("boundary_distance", "9999"), 9999)
        near_path = parse_int(row.get("near_pathogenic_20bp", "0"))
        near_benign = parse_int(row.get("near_benign_20bp", "0"))

        has_path = bool(criteria & PATHOGENIC_CRITERIA)
        has_benign = bool(criteria & BENIGN_CRITERIA)
        has_public = bool(public) and public.get("discordance_label") not in {"", "no_public_classification"}

        flags = {
            "near_threshold_queue": action_group in {"near_pathogenic_threshold", "near_benign_threshold"}
            or "one_step_short" in bottleneck,
            "splice_rna_queue": spliceai >= 0.10
            or boundary <= 2
            or action_group in {"computational_signal_needs_noncomputational_support", "conflicting_computational_splice_context"},
            "functional_queue": bool(criteria & {"PS3", "BS3"}) or bool(findlay),
            "evidence_conflict_queue": action_group == "manual_conflict_adjudication"
            or (has_path and has_benign)
            or public.get("discordance_label") == "public_conflict",
            "public_assertion_queue": has_public,
            "neighborhood_context_queue": near_path >= 5
            or near_benign >= 20
            or "near_pathogenic_dense" in reasons
            or "benign_dense" in reasons,
            "low_information_queue": action_group in {"background_absence_only", "module1_no_signal"}
            or row.get("criteria", "") in {"PM2_Supporting:Supporting:1", ""},
        }

        queue_names = [QUEUE_LABELS[col] for col in QUEUE_COLUMNS if flags[col]]
        synthesized.append(
            {
                **{name: yes(value) for name, value in flags.items()},
                "queue_count": sum(flags.values()),
                "queue_labels": ";".join(queue_names),
                "action_group": action_group,
                "bottleneck_category": bottleneck,
                "public_snapshot_category": public.get("public_snapshot_category", ""),
                "discordance_label": public.get("discordance_label", ""),
                "findlay_func_class": findlay.get("findlay_func_class", ""),
                "findlay_function_score_mean": findlay.get("findlay_function_score_mean", ""),
                "priority_score": row.get("priority_score", ""),
                "priority_tier": row.get("priority_tier", ""),
                "review_category": row.get("review_category", ""),
                "priority_reasons": reasons,
                "gene": row.get("gene", ""),
                "c_notation": row.get("c_notation", ""),
                "p_notation": row.get("p_notation", ""),
                "cds_pos": row.get("cds_pos", ""),
                "cds_exon": row.get("cds_exon", ""),
                "variant_type": row.get("variant_type", ""),
                "criteria": row.get("criteria", ""),
                "total_points": row.get("total_points", ""),
                "spliceai_score": row.get("spliceai_score", ""),
                "boundary_distance": row.get("boundary_distance", ""),
                "near_pathogenic_20bp": row.get("near_pathogenic_20bp", ""),
                "near_benign_20bp": row.get("near_benign_20bp", ""),
            }
        )

    queue_counts = []
    for col in QUEUE_COLUMNS:
        count = sum(1 for row in synthesized if row[col] == "yes")
        queue_counts.append({"queue": col, "queue_label": QUEUE_LABELS[col], "count": count})
    queue_counts.sort(key=lambda item: int(item["count"]), reverse=True)

    overlap_rows = []
    for col_a in QUEUE_COLUMNS:
        for col_b in QUEUE_COLUMNS:
            count = sum(1 for row in synthesized if row[col_a] == "yes" and row[col_b] == "yes")
            overlap_rows.append(
                {
                    "queue_a": col_a,
                    "queue_b": col_b,
                    "queue_a_label": QUEUE_LABELS[col_a],
                    "queue_b_label": QUEUE_LABELS[col_b],
                    "count": count,
                }
            )

    by_action: dict[tuple[str, str], int] = defaultdict(int)
    for row in synthesized:
        action_group = str(row["action_group"]) or "unassigned"
        for col in QUEUE_COLUMNS:
            if row[col] == "yes":
                by_action[(action_group, col)] += 1
    by_action_rows = [
        {
            "action_group": action_group,
            "queue": col,
            "queue_label": QUEUE_LABELS[col],
            "count": count,
        }
        for (action_group, col), count in sorted(by_action.items())
    ]

    top_examples = sorted(
        synthesized,
        key=lambda row: (
            -parse_int(str(row["queue_count"])),
            -parse_int(str(row["priority_score"])),
            str(row["gene"]),
            str(row["c_notation"]),
        ),
    )[:200]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(
        OUT_DIR / "priority_queue_synthesis_variants.csv",
        synthesized,
        QUEUE_COLUMNS
        + [
            "queue_count",
            "queue_labels",
            "action_group",
            "bottleneck_category",
            "public_snapshot_category",
            "discordance_label",
            "findlay_func_class",
            "findlay_function_score_mean",
            "priority_score",
            "priority_tier",
            "review_category",
            "priority_reasons",
            "gene",
            "c_notation",
            "p_notation",
            "cds_pos",
            "cds_exon",
            "variant_type",
            "criteria",
            "total_points",
            "spliceai_score",
            "boundary_distance",
            "near_pathogenic_20bp",
            "near_benign_20bp",
        ],
    )
    write_csv(OUT_DIR / "priority_queue_counts.csv", queue_counts, ["queue", "queue_label", "count"])
    write_csv(OUT_DIR / "priority_queue_overlap.csv", overlap_rows, ["queue_a", "queue_b", "queue_a_label", "queue_b_label", "count"])
    write_csv(OUT_DIR / "priority_queue_by_action_group.csv", by_action_rows, ["action_group", "queue", "queue_label", "count"])
    write_csv(
        OUT_DIR / "priority_queue_top_multiqueue_examples.csv",
        top_examples,
        [
            "queue_count",
            "queue_labels",
            "priority_score",
            "priority_tier",
            "action_group",
            "public_snapshot_category",
            "discordance_label",
            "findlay_func_class",
            "gene",
            "c_notation",
            "p_notation",
            "criteria",
            "spliceai_score",
            "boundary_distance",
            "near_pathogenic_20bp",
            "near_benign_20bp",
        ],
    )

    draw_bar_svg(PLOT_DIR / "priority_queue_counts.svg", queue_counts)
    draw_overlap_svg(PLOT_DIR / "priority_queue_overlap.svg", overlap_rows)

    queue_count_distribution = Counter(int(row["queue_count"]) for row in synthesized)
    report_lines = [
        "# Priority Queue Synthesis",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Purpose",
        "",
        "This analysis operationalizes VUS prioritization as overlapping manual-review queues. "
        "The queues do not add ACMG/ENIGMA evidence and do not reclassify variants. They only "
        "label why a VUS may be worth opening earlier.",
        "",
        "## Queue Counts",
        "",
        "| Queue | VUS count |",
        "| --- | ---: |",
    ]
    for row in queue_counts:
        report_lines.append(f"| {row['queue_label']} | {row['count']} |")
    report_lines.extend(
        [
            "",
            "## Number Of Queues Per VUS",
            "",
            "| Queue count | VUS count |",
            "| ---: | ---: |",
        ]
    )
    for queue_count, count in sorted(queue_count_distribution.items()):
        report_lines.append(f"| {queue_count} | {count} |")
    report_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "A high-priority VUS can be high priority for different reasons. Some variants are "
            "near an ACMG/ENIGMA threshold, some are splice/RNA candidates, some are useful "
            "because public assertions disagree with the local Module 1 output, and many are "
            "low-information variants where absence from population data is the main signal.",
            "",
            "The most useful manual-review queue is therefore not a universal top-N list. It "
            "depends on the curation question: resolving conflicts, finding RNA candidates, "
            "reviewing functional assay metadata, or selecting variants close to an accepted "
            "classification threshold.",
            "",
            "## Output Files",
            "",
            "- `tables/priority_queue_synthesis/priority_queue_synthesis_variants.csv`",
            "- `tables/priority_queue_synthesis/priority_queue_counts.csv`",
            "- `tables/priority_queue_synthesis/priority_queue_overlap.csv`",
            "- `tables/priority_queue_synthesis/priority_queue_by_action_group.csv`",
            "- `tables/priority_queue_synthesis/priority_queue_top_multiqueue_examples.csv`",
            "- `plots/30_priority_queue_synthesis/priority_queue_counts.svg`",
            "- `plots/30_priority_queue_synthesis/priority_queue_overlap.svg`",
            "",
        ]
    )
    REPORT.write_text("\n".join(report_lines), encoding="utf-8")


if __name__ == "__main__":
    main()
