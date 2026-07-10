"""Detailed audit of criteria co-occurrence sanity-check patterns."""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = OUT_DIR / "tables" / "criteria_sanity_audit"
PLOT_DIR = OUT_DIR / "plots" / "11_criteria_sanity_audit"
REPORT = OUT_DIR / "criteria_sanity_audit_report.md"

CONFLICT_PATHOGENIC_LEANING = {"PVS1", "PM5_PTC", "PP3", "PP4", "PS1", "PS3"}
BENIGN_LEANING = {"BA1", "BS1_Strong", "BS1_Supporting", "BS3", "BP1", "BP4", "BP5", "BP7"}


def parse_criteria(criteria: str) -> list[str]:
    if not criteria:
        return []
    return sorted({part.split(":", 1)[0] for part in criteria.split(";") if part})


def class_group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


def coding_position(c_notation: str) -> int | None:
    match = re.search(r"c\.(\d+)", c_notation or "")
    return int(match.group(1)) if match else None


def score_bin(value: str) -> str:
    try:
        score = float(value)
    except (TypeError, ValueError):
        return "missing"
    if score >= 0.80:
        return ">=0.80"
    if score >= 0.50:
        return "0.50-0.79"
    if score >= 0.20:
        return "0.20-0.49"
    return "<0.20"


def reasons_for(row: dict[str, str]) -> list[str]:
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
    return reasons


def load_rows() -> list[dict[str, str]]:
    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            codes = parse_criteria(row["criteria"])
            row["criteria_codes"] = codes
            row["criteria_combo"] = "+".join(codes) if codes else "none"
            row["group"] = class_group(row["predicted_class"])
            row["coding_position"] = coding_position(row["c_notation"])
            row["spliceai_bin"] = score_bin(row["spliceai_score"])
            reasons = reasons_for(row)
            row["sanity_reasons"] = reasons
            if reasons:
                rows.append(row)
    return rows


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


def draw_barplot(path: Path, title: str, rows: list[dict[str, object]], label_col: str, value_col: str) -> None:
    selected = rows[:18]
    width = 980
    height = max(230, 80 + len(selected) * 32)
    left = 330
    top = 55
    plot_w = 520
    max_value = max((int(row[value_col]) for row in selected), default=1)
    body = [
        f'<text x="28" y="34" class="title">{esc(title)}</text>',
    ]
    for index, row in enumerate(selected):
        y = top + index * 32
        value = int(row[value_col])
        bar_w = plot_w * value / max_value
        body.append(f'<text x="24" y="{y + 18}" class="small">{esc(row[label_col])}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{bar_w:.1f}" height="21" fill="#4f46e5"/>')
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


def summarize_by_reason(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    reason_rows = []
    reasons = sorted({reason for row in rows for reason in row["sanity_reasons"]})
    for reason in reasons:
        subset = [row for row in rows if reason in row["sanity_reasons"]]
        genes = Counter(row["gene"] for row in subset)
        groups = Counter(row["group"] for row in subset)
        types = Counter(row["variant_type"] for row in subset)
        scores = [float(row["spliceai_score"]) for row in subset if row["spliceai_score"] not in {"", "NA"}]
        reason_rows.append(
            {
                "sanity_reason": reason,
                "count": len(subset),
                "brca1_count": genes["BRCA1"],
                "brca2_count": genes["BRCA2"],
                "benign_count": groups["benign"],
                "vus_count": groups["vus"],
                "pathogenic_count": groups["pathogenic"],
                "top_variant_type": types.most_common(1)[0][0] if types else "",
                "median_spliceai": f"{median(scores):.2f}" if scores else "",
                "max_spliceai": f"{max(scores):.2f}" if scores else "",
            }
        )
    return sorted(reason_rows, key=lambda row: int(row["count"]), reverse=True)


def summarize_by_reason_gene_type(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    counts = Counter()
    for row in rows:
        for reason in row["sanity_reasons"]:
            counts[(reason, row["gene"], row["variant_type"], row["criteria_combo"], row["group"])] += 1
    output = []
    for (reason, gene, variant_type, combo, group), count in counts.most_common():
        output.append(
            {
                "sanity_reason": reason,
                "gene": gene,
                "variant_type": variant_type,
                "criteria_combo": combo,
                "group": group,
                "count": count,
            }
        )
    return output


def summarize_splice_bins(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    counts = Counter()
    for row in rows:
        for reason in row["sanity_reasons"]:
            counts[(reason, row["spliceai_bin"], row["gene"], row["group"])] += 1
    return [
        {
            "sanity_reason": reason,
            "spliceai_bin": score_bin_name,
            "gene": gene,
            "group": group,
            "count": count,
        }
        for (reason, score_bin_name, gene, group), count in counts.most_common()
    ]


def summarize_positions(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    output = []
    for row in rows:
        output.append(
            {
                "sanity_reasons": ";".join(row["sanity_reasons"]),
                "gene": row["gene"],
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "coding_position": row["coding_position"] or "",
                "variant_type": row["variant_type"],
                "criteria_combo": row["criteria_combo"],
                "total_points": row["total_points"],
                "predicted_class": row["predicted_class"],
                "spliceai_score": row["spliceai_score"],
                "interpretation_bucket": interpretation_bucket(row),
            }
        )
    return sorted(output, key=lambda row: (str(row["gene"]), int(row["coding_position"] or 0), str(row["c_notation"])))


def interpretation_bucket(row: dict[str, str]) -> str:
    reasons = set(row["sanity_reasons"])
    combo = set(row["criteria_codes"])
    if "PVS1_but_VUS" in reasons:
        return "manual_review_single_pvs1_case"
    if "PS3_but_VUS" in reasons and row["predicted_class"] == "3":
        if combo <= {"PM2_Supporting", "PS3"}:
            return "expected_vus_single_strong_evidence"
        return "vus_with_ps3_and_extra_context"
    if "PP3_with_BS3" in reasons:
        return "benign_functional_overrides_pp3_signal"
    if "pathogenic_and_benign_evidence" in reasons:
        return "mixed_direction_evidence"
    return "other"


def median(values: list[float]) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def top_examples(rows: list[dict[str, str]], reason: str, limit: int = 12) -> list[dict[str, object]]:
    subset = [row for row in rows if reason in row["sanity_reasons"]]
    subset = sorted(
        subset,
        key=lambda row: (
            -float(row["spliceai_score"] or 0),
            row["gene"],
            int(row["coding_position"] or 0),
            row["c_notation"],
        ),
    )
    return [
        {
            "gene": row["gene"],
            "c_notation": row["c_notation"],
            "p_notation": row["p_notation"],
            "variant_type": row["variant_type"],
            "criteria_combo": row["criteria_combo"],
            "total_points": row["total_points"],
            "predicted_class": row["predicted_class"],
            "spliceai_score": row["spliceai_score"],
            "interpretation_bucket": interpretation_bucket(row),
        }
        for row in subset[:limit]
    ]


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
    reason_summary: list[dict[str, object]],
    reason_gene_type: list[dict[str, object]],
    position_rows: list[dict[str, object]],
    example_tables: dict[str, list[dict[str, object]]],
) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    buckets = Counter(row["interpretation_bucket"] for row in position_rows)
    bucket_rows = [{"bucket": bucket, "count": count} for bucket, count in buckets.most_common()]
    text = f"""# Criteria Sanity Audit

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

This report expands the criteria co-occurrence sanity checks. It asks whether flagged combinations are likely to be real rule conflicts, expected borderline outcomes, or useful candidates for manual review. It does not add or change ACMG/ENIGMA criteria.

## What Was Checked

The audit focuses on automatically detectable patterns:

- pathogenic and benign evidence in the same variant, excluding `PM2_Supporting` as a conflict driver
- `PP3` together with `BS3`
- `PVS1` variants that remain VUS
- `PS3` variants that remain VUS
- any `BS3` variant that would still classify as pathogenic

`PM2_Supporting` is intentionally treated as background context for conflict detection because it is very common in this synthetic SNV landscape. Counting it as a conflict driver would falsely flag many expected combinations such as `BP1+PM2_Supporting`.

## Main Result

Most sanity hits are not software errors. The largest group is `PS3_but_VUS`: variants with strong pathogenic evidence from the automated rule set, but without enough additional evidence to cross the ENIGMA Table 3 combination threshold. This is an expected Module 1 limitation and a useful triage group.

The most review-worthy groups are:

- the single `PVS1_but_VUS` row, because PVS1 normally contributes very strong evidence but can still remain VUS if the Table 3 combination is not met
- `PP3_with_BS3`, because a computational pathogenic signal and functional benign evidence point in opposite directions
- `pathogenic_and_benign_evidence`, after excluding PM2 background, because these are genuine mixed-direction evidence cases

## Summary By Pattern

{report_table(reason_summary, ["sanity_reason", "count", "brca1_count", "brca2_count", "benign_count", "vus_count", "pathogenic_count", "top_variant_type", "median_spliceai", "max_spliceai"], limit=20)}

## Interpretation Buckets

{report_table(bucket_rows, ["bucket", "count"], limit=20)}

## Top Detailed Pattern Splits

{report_table(reason_gene_type, ["sanity_reason", "gene", "variant_type", "criteria_combo", "group", "count"], limit=25)}

## PS3 But VUS Examples

These are mostly expected VUS outcomes: `PS3` is strong evidence, but alone or with only `PM2_Supporting` it does not satisfy the pathogenic combination rule.

{report_table(example_tables["PS3_but_VUS"], ["gene", "c_notation", "p_notation", "variant_type", "criteria_combo", "total_points", "predicted_class", "spliceai_score", "interpretation_bucket"], limit=12)}

## PP3 With BS3 Examples

These are useful manual-review candidates because the splice/computational signal (`PP3`) and functional benign evidence (`BS3`) point in opposite directions. In the current point system the benign evidence usually dominates.

{report_table(example_tables["PP3_with_BS3"], ["gene", "c_notation", "p_notation", "variant_type", "criteria_combo", "total_points", "predicted_class", "spliceai_score", "interpretation_bucket"], limit=12)}

## PVS1 But VUS

This is the single highest-priority sanity row to inspect manually.

{report_table(example_tables["PVS1_but_VUS"], ["gene", "c_notation", "p_notation", "variant_type", "criteria_combo", "total_points", "predicted_class", "spliceai_score", "interpretation_bucket"], limit=12)}

## Practical Takeaways

1. `PS3_but_VUS` should be treated as a triage label, not as an error label. It identifies variants where Module 1 found meaningful evidence but the rule combination is still incomplete.
2. `PP3_with_BS3` is the best conflict set for curator review, especially where SpliceAI is high.
3. `PVS1_but_VUS` should be manually checked because it tests whether the PVS1 rule strength and final combination logic are behaving as intended.
4. No `BS3_but_pathogenic` pattern was observed, which is reassuring for the current automated evidence combination.

## Outputs

- `tables/criteria_sanity_audit/sanity_reason_summary.csv`
- `tables/criteria_sanity_audit/sanity_reason_gene_type_combo.csv`
- `tables/criteria_sanity_audit/sanity_spliceai_bins.csv`
- `tables/criteria_sanity_audit/sanity_variant_level_audit.csv`
- `tables/criteria_sanity_audit/sanity_examples_PS3_but_VUS.csv`
- `tables/criteria_sanity_audit/sanity_examples_PP3_with_BS3.csv`
- `tables/criteria_sanity_audit/sanity_examples_PVS1_but_VUS.csv`
- `plots/11_criteria_sanity_audit/sanity_reason_counts.svg`
- `plots/11_criteria_sanity_audit/sanity_interpretation_buckets.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    reason_summary = summarize_by_reason(rows)
    reason_gene_type = summarize_by_reason_gene_type(rows)
    splice_bins = summarize_splice_bins(rows)
    position_rows = summarize_positions(rows)
    example_tables = {
        "PS3_but_VUS": top_examples(rows, "PS3_but_VUS"),
        "PP3_with_BS3": top_examples(rows, "PP3_with_BS3"),
        "PVS1_but_VUS": top_examples(rows, "PVS1_but_VUS"),
    }

    write_csv(
        TABLE_DIR / "sanity_reason_summary.csv",
        reason_summary,
        [
            "sanity_reason",
            "count",
            "brca1_count",
            "brca2_count",
            "benign_count",
            "vus_count",
            "pathogenic_count",
            "top_variant_type",
            "median_spliceai",
            "max_spliceai",
        ],
    )
    write_csv(
        TABLE_DIR / "sanity_reason_gene_type_combo.csv",
        reason_gene_type,
        ["sanity_reason", "gene", "variant_type", "criteria_combo", "group", "count"],
    )
    write_csv(
        TABLE_DIR / "sanity_spliceai_bins.csv",
        splice_bins,
        ["sanity_reason", "spliceai_bin", "gene", "group", "count"],
    )
    write_csv(
        TABLE_DIR / "sanity_variant_level_audit.csv",
        position_rows,
        [
            "sanity_reasons",
            "gene",
            "c_notation",
            "p_notation",
            "coding_position",
            "variant_type",
            "criteria_combo",
            "total_points",
            "predicted_class",
            "spliceai_score",
            "interpretation_bucket",
        ],
    )
    for reason, examples in example_tables.items():
        write_csv(
            TABLE_DIR / f"sanity_examples_{reason}.csv",
            examples,
            [
                "gene",
                "c_notation",
                "p_notation",
                "variant_type",
                "criteria_combo",
                "total_points",
                "predicted_class",
                "spliceai_score",
                "interpretation_bucket",
            ],
        )

    bucket_rows = [{"bucket": bucket, "count": count} for bucket, count in Counter(row["interpretation_bucket"] for row in position_rows).most_common()]
    draw_barplot(PLOT_DIR / "sanity_reason_counts.svg", "Sanity-check pattern counts", reason_summary, "sanity_reason", "count")
    draw_barplot(PLOT_DIR / "sanity_interpretation_buckets.svg", "Sanity-check interpretation buckets", bucket_rows, "bucket", "count")
    write_report(reason_summary, reason_gene_type, position_rows, example_tables)
    print(f"Wrote criteria sanity audit to {REPORT}")


if __name__ == "__main__":
    main()
