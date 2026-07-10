from __future__ import annotations

import csv
import html
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from backend.modules.classifier import classify_by_enigma_combination


ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
TABLE_DIR = ANALYSIS_DIR / "tables" / "frequency_baseline"
PLOT_DIR = ANALYSIS_DIR / "plots" / "15_frequency_baseline"
REPORT = ANALYSIS_DIR / "frequency_baseline_audit_report.md"


LABELS = {
    1: "Benign",
    2: "Likely Benign",
    3: "VUS",
    4: "Likely Pathogenic",
    5: "Pathogenic",
}


FREQ_CODES = {"BA1", "BS1_Strong", "BS1_Supporting", "PM2_Supporting"}


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_criteria(value: str) -> Dict[str, dict]:
    criteria: Dict[str, dict] = {}
    if not value:
        return criteria
    for item in value.split(";"):
        parts = item.split(":")
        if len(parts) < 3:
            continue
        code = parts[0].strip()
        if not code:
            continue
        criteria[code] = {
            "applies": True,
            "strength": parts[1].strip(),
            "points": parse_int(parts[2].strip()),
            "reason": "parsed from precomputed snapshot",
        }
    return criteria


def criteria_combo(criteria: Dict[str, dict]) -> str:
    if not criteria:
        return "none"
    return "+".join(sorted(criteria))


def group_label(cls: int) -> str:
    if cls in (1, 2):
        return "benign_group"
    if cls == 3:
        return "vus"
    if cls in (4, 5):
        return "pathogenic_group"
    return "unknown"


def class_label(cls: int) -> str:
    return f"{cls} {LABELS.get(cls, 'Unknown')}"


def write_csv(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_rows() -> List[dict]:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def bar_svg(path: Path, rows: List[dict], label_key: str, value_key: str, title: str, max_items: int = 18) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = rows[:max_items]
    width = 980
    left = 330
    top = 55
    row_h = 30
    chart_w = 560
    height = top + row_h * len(rows) + 45
    max_value = max((int(row[value_key]) for row in rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:13px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".small{font-size:11px;fill:#555}",
        ".bar{fill:#3f7f5f}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = int(row[value_key])
        bar_w = 0 if max_value == 0 else int(chart_w * value / max_value)
        parts.append(f'<text x="24" y="{y + 18}">{html.escape(str(row[label_key]))}</text>')
        parts.append(f'<rect x="{left}" y="{y + 5}" width="{bar_w}" height="18" class="bar"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 19}" class="small">{value}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def status_group(status: str) -> str:
    if "status=found" in status:
        return "found"
    if "absent_v2_only" in status:
        return "absent_v2_only"
    if "absent" in status:
        return "absent"
    if not status:
        return "missing"
    return "other"


def summarize_baseline(rows: List[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    by_status_class = Counter()
    by_freq_code = Counter()
    by_combo = Counter()
    examples = []

    for row in rows:
        cls = parse_int(row["predicted_class"])
        criteria = parse_criteria(row["criteria"])
        codes = set(criteria)
        freq_codes = sorted(codes & FREQ_CODES)
        status = status_group(row.get("gnomad_status", ""))

        by_status_class[(status, class_label(cls))] += 1
        for code in freq_codes:
            by_freq_code[(code, class_label(cls))] += 1
        if freq_codes:
            by_combo[("+".join(freq_codes), group_label(cls))] += 1
            if len(examples) < 500:
                examples.append(
                    {
                        "gene": row["gene"],
                        "c_notation": row["c_notation"],
                        "p_notation": row["p_notation"],
                        "variant_type": row["variant_type"],
                        "predicted_class": cls,
                        "predicted_label": row["predicted_label"],
                        "total_points": row["total_points"],
                        "frequency_codes": "+".join(freq_codes),
                        "criteria": row["criteria"],
                        "gnomad_status_group": status,
                        "gnomad_status": row["gnomad_status"],
                    }
                )

    status_rows = [
        {"gnomad_status_group": status, "class": cls, "count": count}
        for (status, cls), count in by_status_class.items()
    ]
    status_rows.sort(key=lambda r: (r["gnomad_status_group"], r["class"]))

    freq_rows = [
        {"frequency_criterion": code, "class": cls, "count": count}
        for (code, cls), count in by_freq_code.items()
    ]
    freq_rows.sort(key=lambda r: (r["frequency_criterion"], r["class"]))

    combo_rows = [
        {"frequency_combo": combo, "group": group, "count": count}
        for (combo, group), count in by_combo.items()
    ]
    combo_rows.sort(key=lambda r: int(r["count"]), reverse=True)
    return status_rows, freq_rows, combo_rows, examples


def simulate_pm2_loss(rows: List[dict]) -> tuple[list[dict], list[dict]]:
    changed = []
    transitions = Counter()
    for row in rows:
        criteria = parse_criteria(row["criteria"])
        if "PM2_Supporting" not in criteria:
            continue
        original_class = parse_int(row["predicted_class"])
        original_points = parse_int(row["total_points"])
        reduced = {k: v for k, v in criteria.items() if k != "PM2_Supporting"}
        new_points = original_points - parse_int(criteria["PM2_Supporting"]["points"])
        new_class, new_label, note = classify_by_enigma_combination(reduced, new_points)
        transition = f"{class_label(original_class)} -> {class_label(new_class)}"
        transitions[transition] += 1
        if new_class != original_class:
            changed.append(
                {
                    "gene": row["gene"],
                    "c_notation": row["c_notation"],
                    "p_notation": row["p_notation"],
                    "variant_type": row["variant_type"],
                    "original_class": original_class,
                    "new_class_without_pm2": new_class,
                    "original_points": original_points,
                    "new_points_without_pm2": new_points,
                    "original_combo": criteria_combo(criteria),
                    "new_combo_without_pm2": criteria_combo(reduced),
                    "gnomad_status": row["gnomad_status"],
                    "note": note,
                }
            )

    transition_rows = [
        {"transition": transition, "count": count}
        for transition, count in transitions.items()
    ]
    transition_rows.sort(key=lambda r: int(r["count"]), reverse=True)
    changed.sort(key=lambda r: (int(r["original_class"]), int(r["new_class_without_pm2"]), r["gene"], r["c_notation"]))
    return transition_rows, changed


def markdown_table(rows: List[dict], columns: List[str], limit: int = 20) -> str:
    if not rows:
        return "_No rows._"
    shown = rows[:limit]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in shown:
        lines.append("| " + " | ".join(html.escape(str(row.get(col, ""))) for col in columns) + " |")
    return "\n".join(lines)


def build_report(
    rows: List[dict],
    status_rows: List[dict],
    freq_rows: List[dict],
    combo_rows: List[dict],
    transition_rows: List[dict],
    changed_rows: List[dict],
) -> None:
    total = len(rows)
    pm2_count = sum(1 for row in rows if "PM2_Supporting" in parse_criteria(row["criteria"]))
    bs1_count = sum(1 for row in rows if {"BS1_Strong", "BS1_Supporting"} & set(parse_criteria(row["criteria"])))
    ba1_count = sum(1 for row in rows if "BA1" in parse_criteria(row["criteria"]))
    absent_v2_only = sum(1 for row in rows if status_group(row.get("gnomad_status", "")) == "absent_v2_only")
    found = sum(1 for row in rows if status_group(row.get("gnomad_status", "")) == "found")

    text = f"""# Frequency Baseline Audit

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Why We Did This

Budapest 2026 notes indicate that frequency calibration is moving toward
gnomAD v4.1. Before changing data sources, this audit documents how much the
current precomputed snapshot depends on the existing local frequency evidence.

This is a baseline analysis, not a gnomAD v4.1 recalculation. The local cache
currently represents a BRCA regional gnomAD v2.1.1 exome source with real
coverage, while the snapshot often marks variants as `absent_v2_only`.

## Dataset

- Total variants: {total}
- Variants with `PM2_Supporting`: {pm2_count}
- Variants with `BS1_Strong` or `BS1_Supporting`: {bs1_count}
- Variants with `BA1`: {ba1_count}
- gnomAD status `absent_v2_only`: {absent_v2_only}
- gnomAD status `found`: {found}

## gnomAD Status By Class

{markdown_table(status_rows, ["gnomad_status_group", "class", "count"], 20)}

## Frequency Criteria By Class

{markdown_table(freq_rows, ["frequency_criterion", "class", "count"], 30)}

## Frequency Criterion Combos By Group

{markdown_table(combo_rows, ["frequency_combo", "group", "count"], 20)}

## PM2 Loss Simulation

This simulation removes `PM2_Supporting` from each variant that has it and
re-runs the same ENIGMA combination logic. It approximates what could happen if
future gnomAD v4.1 calibration shows that PM2 should not be applied for some
variants.

{markdown_table(transition_rows, ["transition", "count"], 20)}

Variants whose class changes when PM2 is removed: {len(changed_rows)}

Examples:

{markdown_table(changed_rows, ["gene", "c_notation", "p_notation", "original_class", "new_class_without_pm2", "original_points", "new_points_without_pm2", "original_combo"], 25)}

## Interpretation

- The current snapshot is highly dependent on PM2 because most synthetic coding
  SNVs are absent from the local population cache.
- PM2-only VUS variants remain the largest low-information group.
- Loss of PM2 is most important when PM2 completes a pathogenic combination or
  participates in contradictory evidence.
- The `VUS -> Likely Benign` transitions after removing PM2 should be read in
  terms of what PM2 means. PM2 is not direct proof that a variant is damaging;
  it means the variant is absent or very rare in population data and therefore
  lacks population-based benign reassurance. In combinations such as
  `BP4+BP7+PM2_Supporting`, PM2 acts as a cautionary counterweight to weak
  benign computational/splice-neutral evidence. If PM2 is not established, only
  the benign criteria remain, and the same combination logic can move the
  variant to Likely Benign.
- A real gnomAD v4.1 analysis should not simply replace this simulation. It
  should recompute BA1, BS1, and PM2 from v4.1/v4.1.1 frequency and allele
  number data, including exome/genome discordance flags and coverage/mappability
  context.

## Outputs

- `tables/frequency_baseline/gnomad_status_by_class.csv`
- `tables/frequency_baseline/frequency_criteria_by_class.csv`
- `tables/frequency_baseline/frequency_combo_by_group.csv`
- `tables/frequency_baseline/frequency_examples.csv`
- `tables/frequency_baseline/pm2_loss_transitions.csv`
- `tables/frequency_baseline/pm2_loss_changed_variants.csv`
- `plots/15_frequency_baseline/frequency_combo_by_group.svg`
- `plots/15_frequency_baseline/pm2_loss_transitions.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    status_rows, freq_rows, combo_rows, examples = summarize_baseline(rows)
    transition_rows, changed_rows = simulate_pm2_loss(rows)

    write_csv(TABLE_DIR / "gnomad_status_by_class.csv", status_rows, ["gnomad_status_group", "class", "count"])
    write_csv(TABLE_DIR / "frequency_criteria_by_class.csv", freq_rows, ["frequency_criterion", "class", "count"])
    write_csv(TABLE_DIR / "frequency_combo_by_group.csv", combo_rows, ["frequency_combo", "group", "count"])
    write_csv(
        TABLE_DIR / "frequency_examples.csv",
        examples,
        [
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "predicted_class",
            "predicted_label",
            "total_points",
            "frequency_codes",
            "criteria",
            "gnomad_status_group",
            "gnomad_status",
        ],
    )
    write_csv(TABLE_DIR / "pm2_loss_transitions.csv", transition_rows, ["transition", "count"])
    write_csv(
        TABLE_DIR / "pm2_loss_changed_variants.csv",
        changed_rows,
        [
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "original_class",
            "new_class_without_pm2",
            "original_points",
            "new_points_without_pm2",
            "original_combo",
            "new_combo_without_pm2",
            "gnomad_status",
            "note",
        ],
    )

    plot_combo_rows = [
        {"label": f"{row['frequency_combo']} / {row['group']}", "count": row["count"]}
        for row in combo_rows
    ]
    bar_svg(PLOT_DIR / "frequency_combo_by_group.svg", plot_combo_rows, "label", "count", "Frequency Criteria By Group")
    bar_svg(PLOT_DIR / "pm2_loss_transitions.svg", transition_rows, "transition", "count", "PM2 Loss Simulation Transitions")

    build_report(rows, status_rows, freq_rows, combo_rows, transition_rows, changed_rows)


if __name__ == "__main__":
    main()
