"""Class diversity at the same coding position and codon."""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
TABLE_DIR = ANALYSIS_DIR / "tables" / "position_class_conflicts"
PLOT_DIR = ANALYSIS_DIR / "plots" / "22_position_class_conflicts"
REPORT = ANALYSIS_DIR / "position_class_conflict_report.md"

C_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")


def parse_cds_pos(c_notation: str) -> int | None:
    match = C_RE.match(c_notation or "")
    if not match:
        return None
    return int(match.group(1))


def parse_float(value: str | None) -> float:
    if not value:
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def class_group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


def parse_criteria(criteria: str) -> set[str]:
    return {item.split(":", 1)[0] for item in criteria.split(";") if item}


def load_rows() -> list[dict]:
    rows = []
    with INPUT.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pos = parse_cds_pos(row["c_notation"])
            if pos is None:
                continue
            row["cds_pos_int"] = pos
            row["codon"] = (pos - 1) // 3 + 1
            row["class_group"] = class_group(row["predicted_class"])
            row["spliceai_float"] = parse_float(row.get("spliceai_score"))
            row["criteria_codes"] = parse_criteria(row.get("criteria", ""))
            rows.append(row)
    return rows


def pattern_from_groups(groups: set[str]) -> str:
    if groups == {"benign"}:
        return "benign_only"
    if groups == {"pathogenic"}:
        return "pathogenic_only"
    if groups == {"vus"}:
        return "vus_only"
    if groups == {"benign", "pathogenic"}:
        return "benign_and_pathogenic"
    if groups == {"benign", "vus"}:
        return "benign_and_vus"
    if groups == {"pathogenic", "vus"}:
        return "pathogenic_and_vus"
    if groups == {"benign", "pathogenic", "vus"}:
        return "all_three_groups"
    return "other"


def summarize_bucket(gene: str, key_value: int, rows: list[dict], key_name: str) -> dict:
    groups = Counter(row["class_group"] for row in rows)
    classes = Counter(row["predicted_class"] for row in rows)
    criteria = sorted({code for row in rows for code in row["criteria_codes"]})
    variant_types = sorted({row["variant_type"] for row in rows})
    variants = sorted(f"{row['c_notation']} {row['p_notation']} [{row['predicted_class']}]" for row in rows)
    codons = sorted({str(row["codon"]) for row in rows})
    positions = sorted({str(row["cds_pos_int"]) for row in rows})
    out = {
        "gene": gene,
        key_name: key_value,
        "variant_count": len(rows),
        "benign_count": groups["benign"],
        "vus_count": groups["vus"],
        "pathogenic_count": groups["pathogenic"],
        "class_1_count": classes["1"],
        "class_2_count": classes["2"],
        "class_3_count": classes["3"],
        "class_4_count": classes["4"],
        "class_5_count": classes["5"],
        "group_pattern": pattern_from_groups(set(groups)),
        "class_set": "+".join(sorted(classes)),
        "variant_types": "+".join(variant_types),
        "criteria_codes": "+".join(criteria),
        "max_spliceai": f"{max((row['spliceai_float'] for row in rows), default=0.0):.3f}",
        "min_points": min(int(row["total_points"]) for row in rows),
        "max_points": max(int(row["total_points"]) for row in rows),
        "example_variants": "; ".join(variants[:8]),
    }
    if key_name != "codon":
        out["codon"] = codons[0] if codons else ""
    if key_name != "cds_pos":
        out["cds_positions"] = "+".join(positions)
    return out


def build_tables(rows: list[dict]) -> dict[str, list[dict]]:
    by_position: dict[tuple[str, int], list[dict]] = defaultdict(list)
    by_codon: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in rows:
        by_position[(row["gene"], row["cds_pos_int"])].append(row)
        by_codon[(row["gene"], row["codon"])].append(row)

    position_summary = [
        summarize_bucket(gene, pos, bucket_rows, "cds_pos")
        for (gene, pos), bucket_rows in sorted(by_position.items())
    ]
    codon_summary = [
        summarize_bucket(gene, codon, bucket_rows, "codon")
        for (gene, codon), bucket_rows in sorted(by_codon.items())
    ]
    position_mixed = [
        row for row in position_summary
        if row["group_pattern"] in {"benign_and_pathogenic", "all_three_groups"}
    ]
    codon_mixed = [
        row for row in codon_summary
        if row["group_pattern"] in {"benign_and_pathogenic", "all_three_groups"}
    ]
    vus_with_mixed_position_context = [
        row for row in position_summary
        if row["vus_count"] and row["benign_count"] and row["pathogenic_count"]
    ]
    vus_with_mixed_codon_context = [
        row for row in codon_summary
        if row["vus_count"] and row["benign_count"] and row["pathogenic_count"]
    ]
    return {
        "position_summary": position_summary,
        "codon_summary": codon_summary,
        "position_mixed": position_mixed,
        "codon_mixed": codon_mixed,
        "vus_with_mixed_position_context": vus_with_mixed_position_context,
        "vus_with_mixed_codon_context": vus_with_mixed_codon_context,
        "position_mixed_variant_type_summary": summarize_field(position_mixed, "variant_types"),
        "codon_mixed_variant_type_summary": summarize_field(codon_mixed, "variant_types"),
        "position_mixed_criteria_summary": summarize_field(position_mixed, "criteria_codes"),
        "codon_mixed_criteria_summary": summarize_field(codon_mixed, "criteria_codes"),
    }


def summarize_field(rows: list[dict], field: str) -> list[dict]:
    counts = Counter(row[field] for row in rows)
    return [
        {"field": field, "value": value, "count": count}
        for value, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
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


def bar_summary_svg(path: Path, title: str, rows: list[dict]) -> None:
    counts = Counter(row["group_pattern"] for row in rows)
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    width = 980
    height = 72 + len(ordered) * 34
    left = 250
    top = 50
    plot_w = 600
    max_count = max(counts.values(), default=1)
    body = [f'<text x="28" y="30" class="title">{esc(title)}</text>']
    for i, (label, count) in enumerate(ordered):
        y = top + i * 34
        w = plot_w * count / max_count
        color = "#b91c1c" if "pathogenic" in label else "#2563eb" if "vus" in label else "#16a34a"
        if label == "all_three_groups":
            color = "#7c3aed"
        body.append(f'<text x="28" y="{y + 17}" class="small">{esc(label)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w:.1f}" height="22" fill="{color}" opacity="0.85"/>')
        body.append(f'<text x="{left + w + 8:.1f}" y="{y + 17}" class="small">{count}</text>')
    write_svg(path, "\n".join(body), width, height)


def mixed_examples_svg(path: Path, title: str, rows: list[dict], key_name: str) -> None:
    selected = sorted(
        rows,
        key=lambda row: (
            row["group_pattern"] != "all_three_groups",
            -int(row["pathogenic_count"]),
            -int(row["benign_count"]),
            row["gene"],
            int(row[key_name]),
        ),
    )[:35]
    width = 1180
    height = 90 + len(selected) * 24
    left = 330
    top = 52
    scale = 34
    body = [f'<text x="28" y="30" class="title">{esc(title)}</text>']
    body.append(f'<text x="{left}" y="30" class="small">green benign, gray VUS, red pathogenic</text>')
    for i, row in enumerate(selected):
        y = top + i * 24
        label = f"{row['gene']} {key_name} {row[key_name]} ({row['group_pattern']})"
        body.append(f'<text x="28" y="{y + 15}" class="small">{esc(label)}</text>')
        parts = [
            ("benign", int(row["benign_count"]), "#16a34a"),
            ("vus", int(row["vus_count"]), "#64748b"),
            ("pathogenic", int(row["pathogenic_count"]), "#b91c1c"),
        ]
        x = left
        for name, count, color in parts:
            if count:
                w = count * scale
                body.append(f'<rect x="{x}" y="{y}" width="{w}" height="16" fill="{color}" opacity="0.82"/>')
                body.append(f'<text x="{x + 4}" y="{y + 12}" class="small" fill="#ffffff">{count}</text>')
                x += w + 4
        body.append(f'<text x="{left + 240}" y="{y + 15}" class="small">{esc(row["class_set"])}; max SpliceAI {row["max_spliceai"]}</text>')
    write_svg(path, "\n".join(body), width, height)


def markdown_table(rows: list[dict], columns: list[str], limit: int = 20) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def write_report(tables: dict[str, list[dict]]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    position_counts = Counter(row["group_pattern"] for row in tables["position_summary"])
    codon_counts = Counter(row["group_pattern"] for row in tables["codon_summary"])
    position_summary_rows = [
        {"level": "cds_position", "pattern": pattern, "count": count}
        for pattern, count in sorted(position_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    codon_summary_rows = [
        {"level": "codon", "pattern": pattern, "count": count}
        for pattern, count in sorted(codon_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    text = f"""# Position And Codon Class Conflict Analysis

Generated: {generated}

Input: `{INPUT.relative_to(ROOT)}`

## Purpose

This analysis separates a coding position from a concrete variant. A single
coding position can have up to three possible SNV substitutions, and those
substitutions can receive different generated ARIANE classes. The same is even
more true at codon level, where up to nine SNVs can affect one amino-acid
position.

This matters for VUS interpretation because a local pathogenic-looking region
does not mean every possible substitution at that position is pathogenic. The
opposite is also true: a benign variant at a position does not prove that a
different substitution at the same position is benign.

This is an exploratory sanity-check of the generated Module 1 snapshot only. It
does not reclassify variants.

## Angles Tested

1. Same coding nucleotide position: do different SNVs at the same `c.` position
   fall into different grouped classes?
2. Same codon: do different SNVs in the same codon fall into different grouped
   classes?
3. VUS inside mixed contexts: are there VUS at positions or codons where both
   benign and pathogenic generated variants also exist?
4. Mechanistic hints: do mixed positions differ by criteria, variant type,
   points, or SpliceAI score?

## Position-Level Summary

{markdown_table(position_summary_rows, ["level", "pattern", "count"], 20)}

## Codon-Level Summary

{markdown_table(codon_summary_rows, ["level", "pattern", "count"], 20)}

## Same Coding Position With Both Benign And Pathogenic Generated Variants

{markdown_table(tables["position_mixed"], ["gene", "cds_pos", "codon", "group_pattern", "variant_count", "benign_count", "vus_count", "pathogenic_count", "class_set", "variant_types", "criteria_codes", "max_spliceai", "example_variants"], 25)}

## Same Codon With Both Benign And Pathogenic Generated Variants

{markdown_table(tables["codon_mixed"], ["gene", "codon", "cds_positions", "group_pattern", "variant_count", "benign_count", "vus_count", "pathogenic_count", "class_set", "variant_types", "criteria_codes", "max_spliceai", "example_variants"], 25)}

## VUS In Fully Mixed Contexts

Same-position contexts with VUS plus benign plus pathogenic generated variants:
{len(tables["vus_with_mixed_position_context"])}

Same-codon contexts with VUS plus benign plus pathogenic generated variants:
{len(tables["vus_with_mixed_codon_context"])}

These are useful manual-review candidates because the local context is not
uniform. For such variants, the exact nucleotide change, protein consequence,
splice signal, and criterion combination matter more than the position label.

## What Drives The Mixed Contexts?

Same-position mixed contexts by variant-type combination:

{markdown_table(tables["position_mixed_variant_type_summary"], ["field", "value", "count"], 20)}

Same-codon mixed contexts by variant-type combination:

{markdown_table(tables["codon_mixed_variant_type_summary"], ["field", "value", "count"], 20)}

Same-position mixed contexts by criterion-code combination:

{markdown_table(tables["position_mixed_criteria_summary"], ["field", "value", "count"], 20)}

Same-codon mixed contexts by criterion-code combination:

{markdown_table(tables["codon_mixed_criteria_summary"], ["field", "value", "count"], 20)}

The dominant pattern should be read biologically and technically: many mixed
positions are expected because one possible nucleotide substitution creates a
premature stop codon while another creates a missense or synonymous consequence.
The more interesting subset is where the mixed context is driven by splice
signal, functional evidence, PS1/PS3, or other non-truncating mechanisms.

## Outputs

- `tables/position_class_conflicts/position_class_summary.csv`
- `tables/position_class_conflicts/codon_class_summary.csv`
- `tables/position_class_conflicts/same_position_benign_pathogenic.csv`
- `tables/position_class_conflicts/same_codon_benign_pathogenic.csv`
- `tables/position_class_conflicts/vus_with_mixed_position_context.csv`
- `tables/position_class_conflicts/vus_with_mixed_codon_context.csv`
- `tables/position_class_conflicts/position_mixed_variant_type_summary.csv`
- `tables/position_class_conflicts/codon_mixed_variant_type_summary.csv`
- `tables/position_class_conflicts/position_mixed_criteria_summary.csv`
- `tables/position_class_conflicts/codon_mixed_criteria_summary.csv`
- `plots/22_position_class_conflicts/position_group_patterns.svg`
- `plots/22_position_class_conflicts/codon_group_patterns.svg`
- `plots/22_position_class_conflicts/same_position_mixed_examples.svg`
- `plots/22_position_class_conflicts/same_codon_mixed_examples.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    tables = build_tables(rows)
    position_fields = [
        "gene",
        "cds_pos",
        "codon",
        "variant_count",
        "benign_count",
        "vus_count",
        "pathogenic_count",
        "class_1_count",
        "class_2_count",
        "class_3_count",
        "class_4_count",
        "class_5_count",
        "group_pattern",
        "class_set",
        "variant_types",
        "criteria_codes",
        "max_spliceai",
        "min_points",
        "max_points",
        "example_variants",
    ]
    codon_fields = [
        "gene",
        "codon",
        "cds_positions",
        "variant_count",
        "benign_count",
        "vus_count",
        "pathogenic_count",
        "class_1_count",
        "class_2_count",
        "class_3_count",
        "class_4_count",
        "class_5_count",
        "group_pattern",
        "class_set",
        "variant_types",
        "criteria_codes",
        "max_spliceai",
        "min_points",
        "max_points",
        "example_variants",
    ]
    write_csv(TABLE_DIR / "position_class_summary.csv", tables["position_summary"], position_fields)
    write_csv(TABLE_DIR / "codon_class_summary.csv", tables["codon_summary"], codon_fields)
    write_csv(TABLE_DIR / "same_position_benign_pathogenic.csv", tables["position_mixed"], position_fields)
    write_csv(TABLE_DIR / "same_codon_benign_pathogenic.csv", tables["codon_mixed"], codon_fields)
    write_csv(TABLE_DIR / "vus_with_mixed_position_context.csv", tables["vus_with_mixed_position_context"], position_fields)
    write_csv(TABLE_DIR / "vus_with_mixed_codon_context.csv", tables["vus_with_mixed_codon_context"], codon_fields)
    summary_fields = ["field", "value", "count"]
    write_csv(TABLE_DIR / "position_mixed_variant_type_summary.csv", tables["position_mixed_variant_type_summary"], summary_fields)
    write_csv(TABLE_DIR / "codon_mixed_variant_type_summary.csv", tables["codon_mixed_variant_type_summary"], summary_fields)
    write_csv(TABLE_DIR / "position_mixed_criteria_summary.csv", tables["position_mixed_criteria_summary"], summary_fields)
    write_csv(TABLE_DIR / "codon_mixed_criteria_summary.csv", tables["codon_mixed_criteria_summary"], summary_fields)
    bar_summary_svg(PLOT_DIR / "position_group_patterns.svg", "Grouped class patterns per CDS position", tables["position_summary"])
    bar_summary_svg(PLOT_DIR / "codon_group_patterns.svg", "Grouped class patterns per codon", tables["codon_summary"])
    mixed_examples_svg(PLOT_DIR / "same_position_mixed_examples.svg", "Same CDS position with mixed generated classes", tables["position_mixed"], "cds_pos")
    mixed_examples_svg(PLOT_DIR / "same_codon_mixed_examples.svg", "Same codon with mixed generated classes", tables["codon_mixed"], "codon")
    write_report(tables)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
