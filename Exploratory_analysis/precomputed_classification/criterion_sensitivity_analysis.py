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
TABLE_DIR = ANALYSIS_DIR / "tables" / "criterion_sensitivity"
PLOT_DIR = ANALYSIS_DIR / "plots" / "14_criterion_sensitivity"
REPORT = ANALYSIS_DIR / "criterion_sensitivity_report.md"


LABELS = {
    1: "Benign",
    2: "Likely Benign",
    3: "VUS",
    4: "Likely Pathogenic",
    5: "Pathogenic",
}


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_float(value: str) -> float | None:
    try:
        if value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_criteria(value: str) -> Dict[str, dict]:
    criteria: Dict[str, dict] = {}
    if not value:
        return criteria
    for item in value.split(";"):
        parts = item.split(":")
        if len(parts) < 3:
            continue
        code = parts[0].strip()
        strength = parts[1].strip()
        points = parse_int(parts[2].strip())
        if code:
            criteria[code] = {
                "applies": True,
                "strength": strength,
                "points": points,
                "reason": "parsed from precomputed snapshot",
            }
    return criteria


def criteria_combo(criteria: Dict[str, dict]) -> str:
    if not criteria:
        return "none"
    return "+".join(sorted(criteria))


def class_label(cls: int) -> str:
    return f"{cls} {LABELS.get(cls, 'Unknown')}"


def write_csv(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def bar_svg(path: Path, rows: List[dict], label_key: str, value_key: str, title: str, max_items: int = 18) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = rows[:max_items]
    width = 980
    left = 300
    top = 55
    row_h = 30
    chart_w = 600
    height = top + row_h * len(rows) + 45
    max_value = max((int(row[value_key]) for row in rows), default=1)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:13px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".small{font-size:11px;fill:#555}",
        ".bar{fill:#4f6bed}",
        ".grid{stroke:#ddd;stroke-width:1}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = int(row[value_key])
        bar_w = 0 if max_value == 0 else int(chart_w * value / max_value)
        label = str(row[label_key])
        parts.append(f'<text x="24" y="{y + 18}">{html.escape(label)}</text>')
        parts.append(f'<rect x="{left}" y="{y + 5}" width="{bar_w}" height="18" class="bar"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 19}" class="small">{value}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def load_rows() -> List[dict]:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def removal_sensitivity(rows: List[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    variant_rows = []
    criterion_summary: Counter[tuple] = Counter()
    transition_summary: Counter[tuple] = Counter()

    for row in rows:
        criteria = parse_criteria(row["criteria"])
        if not criteria:
            continue
        original_class = parse_int(row["predicted_class"])
        original_points = parse_int(row["total_points"])
        original_combo = criteria_combo(criteria)

        for removed_code, removed_criterion in criteria.items():
            reduced = {k: v for k, v in criteria.items() if k != removed_code}
            reduced_points = original_points - parse_int(removed_criterion.get("points", 0))
            new_class, new_label, note = classify_by_enigma_combination(reduced, reduced_points)
            changed = new_class != original_class
            direction = "unchanged"
            if new_class > original_class:
                direction = "higher_class_after_removal"
            elif new_class < original_class:
                direction = "lower_class_after_removal"

            transition = f"{class_label(original_class)} -> {class_label(new_class)}"
            criterion_summary[(removed_code, changed, direction)] += 1
            transition_summary[(removed_code, transition)] += 1

            if changed:
                variant_rows.append(
                    {
                        "gene": row["gene"],
                        "c_notation": row["c_notation"],
                        "p_notation": row["p_notation"],
                        "variant_type": row["variant_type"],
                        "removed_criterion": removed_code,
                        "removed_points": removed_criterion.get("points", 0),
                        "original_class": original_class,
                        "new_class": new_class,
                        "original_points": original_points,
                        "new_points": reduced_points,
                        "direction": direction,
                        "original_combo": original_combo,
                        "new_combo": criteria_combo(reduced),
                        "spliceai_score": row["spliceai_score"],
                        "note": note,
                    }
                )

    criterion_rows = []
    totals_by_code = Counter()
    changed_by_code = Counter()
    for (code, changed, direction), count in criterion_summary.items():
        totals_by_code[code] += count
        if changed:
            changed_by_code[code] += count
    for code in sorted(totals_by_code):
        criterion_rows.append(
            {
                "criterion": code,
                "applications": totals_by_code[code],
                "class_changed_when_removed": changed_by_code[code],
                "percent_changed": f"{100 * changed_by_code[code] / totals_by_code[code]:.2f}",
            }
        )
    criterion_rows.sort(key=lambda r: (int(r["class_changed_when_removed"]), int(r["applications"])), reverse=True)

    transition_rows = [
        {"removed_criterion": code, "transition": transition, "count": count}
        for (code, transition), count in transition_summary.items()
    ]
    transition_rows.sort(key=lambda r: int(r["count"]), reverse=True)
    variant_rows.sort(
        key=lambda r: (
            abs(int(r["original_class"]) - int(r["new_class"])),
            abs(int(r["removed_points"])),
            parse_float(r["spliceai_score"]) or 0,
        ),
        reverse=True,
    )
    return criterion_rows, transition_rows, variant_rows


def hypothetical_vus_support(rows: List[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    variant_rows = []
    summary_counter: Counter[tuple] = Counter()
    combo_counter: Counter[tuple] = Counter()

    hypothetical = {
        "pathogenic_supporting": {
            "HYP_PATH_SUPPORT": {
                "applies": True,
                "strength": "Supporting",
                "points": 1,
                "reason": "hypothetical additional pathogenic supporting evidence",
            }
        },
        "benign_supporting": {
            "HYP_BEN_SUPPORT": {
                "applies": True,
                "strength": "Supporting",
                "points": -1,
                "reason": "hypothetical additional benign supporting evidence",
            }
        },
        "pathogenic_strong": {
            "HYP_PATH_STRONG": {
                "applies": True,
                "strength": "Strong",
                "points": 4,
                "reason": "hypothetical additional pathogenic strong evidence",
            }
        },
        "benign_strong": {
            "HYP_BEN_STRONG": {
                "applies": True,
                "strength": "Strong",
                "points": -4,
                "reason": "hypothetical additional benign strong evidence",
            }
        },
    }

    for row in rows:
        original_class = parse_int(row["predicted_class"])
        if original_class != 3:
            continue
        criteria = parse_criteria(row["criteria"])
        original_points = parse_int(row["total_points"])
        original_combo = criteria_combo(criteria)

        for scenario, added in hypothetical.items():
            added_code, added_criterion = next(iter(added.items()))
            new_criteria = dict(criteria)
            new_criteria[added_code] = added_criterion
            new_points = original_points + parse_int(added_criterion["points"])
            new_class, new_label, note = classify_by_enigma_combination(new_criteria, new_points)
            transition = f"3 VUS -> {class_label(new_class)}"
            summary_counter[(scenario, transition)] += 1
            combo_counter[(scenario, original_combo, transition)] += 1
            if new_class != 3:
                variant_rows.append(
                    {
                        "scenario": scenario,
                        "gene": row["gene"],
                        "c_notation": row["c_notation"],
                        "p_notation": row["p_notation"],
                        "variant_type": row["variant_type"],
                        "original_points": original_points,
                        "new_points": new_points,
                        "new_class": new_class,
                        "original_combo": original_combo,
                        "spliceai_score": row["spliceai_score"],
                        "note": note,
                    }
                )

    summary_rows = [
        {"scenario": scenario, "transition": transition, "count": count}
        for (scenario, transition), count in summary_counter.items()
    ]
    summary_rows.sort(key=lambda r: (r["scenario"], int(r["count"])), reverse=True)
    combo_rows = [
        {"scenario": scenario, "original_combo": combo, "transition": transition, "count": count}
        for (scenario, combo, transition), count in combo_counter.items()
    ]
    combo_rows.sort(key=lambda r: int(r["count"]), reverse=True)
    variant_rows.sort(
        key=lambda r: (
            r["scenario"],
            int(r["new_class"]),
            int(r["original_points"]),
            parse_float(r["spliceai_score"]) or 0,
        ),
        reverse=True,
    )
    return summary_rows, combo_rows, variant_rows


def markdown_table(rows: List[dict], columns: List[str], limit: int = 20) -> str:
    if not rows:
        return "_No rows._"
    shown = rows[:limit]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in shown:
        values = [html.escape(str(row.get(col, ""))) for col in columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_report(
    rows: List[dict],
    criterion_rows: List[dict],
    transition_rows: List[dict],
    changed_rows: List[dict],
    hypothetical_summary: List[dict],
    hypothetical_combo: List[dict],
    hypothetical_variants: List[dict],
) -> None:
    total = len(rows)
    vus_count = sum(1 for row in rows if parse_int(row["predicted_class"]) == 3)
    changed_total = len(changed_rows)
    text = f"""# Criterion Sensitivity Analysis

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Why We Did This

This analysis asks which automated criteria actually hold the generated
classification in place. It is not a new classifier and it does not change any
variant class.

For each variant, the analysis removes one applied criterion at a time and
re-runs the same ENIGMA combination logic on the remaining criteria. If the
class changes, that criterion is classification-sensitive for that variant.

For VUS variants, it also asks a separate hypothetical question: what would
happen if one additional accepted supporting or strong evidence item appeared?
This is a triage model only. The hypothetical evidence is not applied to the
real classification.

## Dataset

- Total variants: {total}
- VUS variants: {vus_count}
- Single-criterion removals that changed class: {changed_total}

## Criteria Most Often Changing Class When Removed

{markdown_table(criterion_rows, ["criterion", "applications", "class_changed_when_removed", "percent_changed"], 20)}

## Most Common Removal Transitions

{markdown_table(transition_rows, ["removed_criterion", "transition", "count"], 20)}

## VUS Hypothetical Evidence Impact

{markdown_table(hypothetical_summary, ["scenario", "transition", "count"], 20)}

## VUS Combos Most Affected By Hypothetical Evidence

{markdown_table(hypothetical_combo, ["scenario", "original_combo", "transition", "count"], 25)}

## Example Variants Whose Class Changes After Criterion Removal

{markdown_table(changed_rows, ["gene", "c_notation", "p_notation", "removed_criterion", "original_class", "new_class", "original_points", "new_points", "original_combo"], 25)}

## Example VUS Variants Moved By Hypothetical Evidence

{markdown_table(hypothetical_variants, ["scenario", "gene", "c_notation", "p_notation", "original_points", "new_points", "new_class", "original_combo"], 25)}

## Interpretation

- Criteria with many class changes are not necessarily more important
  biologically; they are more important for the current automated decision
  boundary.
- `BP1`, `PVS1`, `PM5_PTC`, `BS3`, `PS3`, and frequency criteria are expected to
  be highly sensitive because they carry strong point weight or satisfy Table 3
  combinations.
- The VUS hypothetical section is useful for curation planning: it identifies
  which VUS groups are one accepted evidence item away from a different class.
- This analysis should not be used to add evidence that is not present in the
  real data.

## Outputs

- `tables/criterion_sensitivity/criterion_removal_summary.csv`
- `tables/criterion_sensitivity/criterion_removal_transitions.csv`
- `tables/criterion_sensitivity/criterion_removal_changed_variants.csv`
- `tables/criterion_sensitivity/vus_hypothetical_evidence_summary.csv`
- `tables/criterion_sensitivity/vus_hypothetical_combo_impact.csv`
- `tables/criterion_sensitivity/vus_hypothetical_changed_variants.csv`
- `plots/14_criterion_sensitivity/criterion_removal_changed_counts.svg`
- `plots/14_criterion_sensitivity/vus_hypothetical_evidence_impact.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()

    criterion_rows, transition_rows, changed_rows = removal_sensitivity(rows)
    hypothetical_summary, hypothetical_combo, hypothetical_variants = hypothetical_vus_support(rows)

    write_csv(
        TABLE_DIR / "criterion_removal_summary.csv",
        criterion_rows,
        ["criterion", "applications", "class_changed_when_removed", "percent_changed"],
    )
    write_csv(
        TABLE_DIR / "criterion_removal_transitions.csv",
        transition_rows,
        ["removed_criterion", "transition", "count"],
    )
    write_csv(
        TABLE_DIR / "criterion_removal_changed_variants.csv",
        changed_rows,
        [
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "removed_criterion",
            "removed_points",
            "original_class",
            "new_class",
            "original_points",
            "new_points",
            "direction",
            "original_combo",
            "new_combo",
            "spliceai_score",
            "note",
        ],
    )
    write_csv(
        TABLE_DIR / "vus_hypothetical_evidence_summary.csv",
        hypothetical_summary,
        ["scenario", "transition", "count"],
    )
    write_csv(
        TABLE_DIR / "vus_hypothetical_combo_impact.csv",
        hypothetical_combo,
        ["scenario", "original_combo", "transition", "count"],
    )
    write_csv(
        TABLE_DIR / "vus_hypothetical_changed_variants.csv",
        hypothetical_variants,
        [
            "scenario",
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "original_points",
            "new_points",
            "new_class",
            "original_combo",
            "spliceai_score",
            "note",
        ],
    )

    bar_svg(
        PLOT_DIR / "criterion_removal_changed_counts.svg",
        criterion_rows,
        "criterion",
        "class_changed_when_removed",
        "Class Changes After Removing One Criterion",
    )
    hyp_plot_rows = [
        {
            "label": f"{row['scenario']}: {row['transition']}",
            "count": row["count"],
        }
        for row in hypothetical_summary
        if not str(row["transition"]).endswith("3 VUS")
    ]
    hyp_plot_rows.sort(key=lambda r: int(r["count"]), reverse=True)
    bar_svg(
        PLOT_DIR / "vus_hypothetical_evidence_impact.svg",
        hyp_plot_rows,
        "label",
        "count",
        "VUS Moved By Hypothetical Evidence",
    )

    build_report(
        rows,
        criterion_rows,
        transition_rows,
        changed_rows,
        hypothetical_summary,
        hypothetical_combo,
        hypothetical_variants,
    )


if __name__ == "__main__":
    main()
