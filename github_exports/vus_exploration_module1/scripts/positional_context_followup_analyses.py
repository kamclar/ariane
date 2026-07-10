"""Follow-up analyses for interpreting positional context."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
TRANSCRIPT_MAPS = ROOT / "variant_space_scan" / "coordinate_mapping" / "data" / "transcript_maps.json"
TABLE_DIR = ANALYSIS_DIR / "tables" / "positional_context_followup"
PLOT_DIR = ANALYSIS_DIR / "plots" / "26_positional_context_followup"
REPORT = ANALYSIS_DIR / "positional_context_followup_report.md"

CDS_LENGTHS = {"BRCA1": 5592, "BRCA2": 10257}
BIN_SIZE = 100
C_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")
P_RE = re.compile(r"[A-Z][a-z]{2}(\d+)")

GROUP_ORDER = ["benign", "vus", "pathogenic"]
GROUP_COLORS = {"benign": "#16a34a", "vus": "#64748b", "pathogenic": "#b91c1c"}
VARIANT_TYPES = ["missense", "synonymous", "nonsense", "initiation_codon"]


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_cds_pos(c_notation: str) -> int | None:
    match = C_RE.match(c_notation or "")
    if not match:
        return None
    return int(match.group(1))


def aa_position(p_notation: str) -> int | None:
    match = P_RE.search(p_notation or "")
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


def parse_criteria(criteria: str) -> set[str]:
    return {item.split(":", 1)[0] for item in criteria.split(";") if item}


def group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


def cds_boundaries() -> dict[str, list[int]]:
    data = json.loads(TRANSCRIPT_MAPS.read_text(encoding="utf-8"))
    output = {}
    for gene, info in data["assemblies"]["GRCh38"].items():
        features = list(info["cds"])
        ordered = sorted(features, key=lambda item: item["start"], reverse=info["strand"] == "-")
        pos = 1
        boundaries = []
        for feature in ordered:
            length = feature["end"] - feature["start"] + 1
            start = pos
            end = pos + length - 1
            boundaries.extend([start, end])
            pos = end + 1
        output[gene] = sorted(set(boundaries))
    return output


def nearest_boundary_distance(gene: str, pos: int, boundaries: dict[str, list[int]]) -> int:
    return min(abs(pos - boundary) for boundary in boundaries[gene])


def load_rows() -> list[dict]:
    boundaries = cds_boundaries()
    rows = []
    with INPUT.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pos = parse_cds_pos(row["c_notation"])
            if pos is None:
                continue
            criteria = parse_criteria(row.get("criteria", ""))
            spliceai = parse_float(row.get("spliceai_score"))
            row["cds_pos"] = pos
            row["aa_pos"] = aa_position(row.get("p_notation", ""))
            row["group"] = group(row["predicted_class"])
            row["criteria_codes"] = criteria
            row["spliceai_float"] = spliceai
            row["boundary_distance"] = nearest_boundary_distance(row["gene"], pos, boundaries)
            row["bin"] = (pos - 1) // BIN_SIZE + 1
            row["truncation_driver"] = row["variant_type"] == "nonsense" or bool(criteria & {"PVS1", "PM5_PTC"})
            row["splice_driver"] = "PP3" in criteria or spliceai >= 0.20
            row["functional_pathogenic_driver"] = bool(criteria & {"PS3", "PS1", "PP4"})
            row["functional_benign_driver"] = "BS3" in criteria
            rows.append(row)
    return rows


def count_by_group(rows: list[dict]) -> Counter:
    return Counter(row["group"] for row in rows)


def pathogenic_fraction(rows: list[dict]) -> float:
    return sum(1 for row in rows if row["group"] == "pathogenic") / len(rows) if rows else 0.0


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
.grid {{ stroke: #e2e8f0; stroke-width: 1; }}
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def variant_type_summary(rows: list[dict]) -> list[dict]:
    output = []
    for variant_type in VARIANT_TYPES:
        items = [row for row in rows if row["variant_type"] == variant_type]
        counts = count_by_group(items)
        output.append(
            {
                "variant_type": variant_type,
                "total_count": len(items),
                "benign_count": counts["benign"],
                "vus_count": counts["vus"],
                "pathogenic_count": counts["pathogenic"],
                "pathogenic_fraction": f"{pathogenic_fraction(items):.4f}",
                "truncation_driver_count": sum(1 for row in items if row["truncation_driver"]),
                "splice_driver_count": sum(1 for row in items if row["splice_driver"]),
                "functional_pathogenic_driver_count": sum(1 for row in items if row["functional_pathogenic_driver"]),
                "functional_benign_driver_count": sum(1 for row in items if row["functional_benign_driver"]),
            }
        )
    return output


def regional_ablation_by_type(rows: list[dict]) -> list[dict]:
    output = []
    by_bin: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in rows:
        by_bin[(row["gene"], row["bin"])].append(row)
    for (gene, idx), bucket in sorted(by_bin.items()):
        no_trunc = [row for row in bucket if not row["truncation_driver"]]
        no_splice = [row for row in bucket if not row["splice_driver"]]
        no_trunc_splice = [row for row in bucket if not row["truncation_driver"] and not row["splice_driver"]]
        no_functional = [row for row in bucket if not row["functional_pathogenic_driver"] and not row["functional_benign_driver"]]
        start = (idx - 1) * BIN_SIZE + 1
        end = min(idx * BIN_SIZE, CDS_LENGTHS[gene])
        output.append(
            {
                "gene": gene,
                "bin": idx,
                "cds_start": start,
                "cds_end": end,
                "all_pathogenic_fraction": f"{pathogenic_fraction(bucket):.4f}",
                "no_truncation_pathogenic_fraction": f"{pathogenic_fraction(no_trunc):.4f}",
                "no_splice_pathogenic_fraction": f"{pathogenic_fraction(no_splice):.4f}",
                "no_truncation_no_splice_pathogenic_fraction": f"{pathogenic_fraction(no_trunc_splice):.4f}",
                "no_functional_pathogenic_fraction": f"{pathogenic_fraction(no_functional):.4f}",
                "total_count": len(bucket),
                "pathogenic_count": sum(1 for row in bucket if row["group"] == "pathogenic"),
                "truncation_count": sum(1 for row in bucket if row["truncation_driver"]),
                "splice_count": sum(1 for row in bucket if row["splice_driver"]),
                "functional_count": sum(1 for row in bucket if row["functional_pathogenic_driver"] or row["functional_benign_driver"]),
            }
        )
    return output


def splice_driver_regions(rows: list[dict]) -> list[dict]:
    output = []
    by_bin: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in rows:
        by_bin[(row["gene"], row["bin"])].append(row)
    for (gene, idx), bucket in sorted(by_bin.items()):
        start = (idx - 1) * BIN_SIZE + 1
        end = min(idx * BIN_SIZE, CDS_LENGTHS[gene])
        pp3_count = sum(1 for row in bucket if "PP3" in row["criteria_codes"])
        high_splice = sum(1 for row in bucket if row["spliceai_float"] >= 0.20)
        near_boundary = sum(1 for row in bucket if row["boundary_distance"] <= 10)
        output.append(
            {
                "gene": gene,
                "bin": idx,
                "cds_start": start,
                "cds_end": end,
                "pathogenic_fraction": f"{pathogenic_fraction(bucket):.4f}",
                "splice_driver_fraction": f"{sum(1 for row in bucket if row['splice_driver']) / len(bucket):.4f}",
                "pp3_count": pp3_count,
                "high_spliceai_count": high_splice,
                "near_boundary_count": near_boundary,
                "total_count": len(bucket),
            }
        )
    return sorted(output, key=lambda row: (float(row["splice_driver_fraction"]), float(row["pathogenic_fraction"])), reverse=True)


def functional_driver_regions(rows: list[dict]) -> list[dict]:
    output = []
    by_bin: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in rows:
        by_bin[(row["gene"], row["bin"])].append(row)
    for (gene, idx), bucket in sorted(by_bin.items()):
        start = (idx - 1) * BIN_SIZE + 1
        end = min(idx * BIN_SIZE, CDS_LENGTHS[gene])
        output.append(
            {
                "gene": gene,
                "bin": idx,
                "cds_start": start,
                "cds_end": end,
                "pathogenic_fraction": f"{pathogenic_fraction(bucket):.4f}",
                "ps3_ps1_pp4_count": sum(1 for row in bucket if row["functional_pathogenic_driver"]),
                "bs3_count": sum(1 for row in bucket if row["functional_benign_driver"]),
                "bs3_pp3_count": sum(1 for row in bucket if row["functional_benign_driver"] and "PP3" in row["criteria_codes"]),
                "ps3_pp3_count": sum(1 for row in bucket if row["functional_pathogenic_driver"] and "PP3" in row["criteria_codes"]),
                "total_count": len(bucket),
            }
        )
    return sorted(output, key=lambda row: (int(row["ps3_ps1_pp4_count"]) + int(row["bs3_count"]), float(row["pathogenic_fraction"])), reverse=True)


def mixed_context_ablation(rows: list[dict]) -> list[dict]:
    filters = {
        "all_variants": lambda row: True,
        "no_truncation": lambda row: not row["truncation_driver"],
        "no_splice": lambda row: not row["splice_driver"],
        "no_truncation_no_splice": lambda row: not row["truncation_driver"] and not row["splice_driver"],
        "no_truncation_no_splice_no_functional": lambda row: (
            not row["truncation_driver"]
            and not row["splice_driver"]
            and not row["functional_pathogenic_driver"]
            and not row["functional_benign_driver"]
        ),
    }
    output = []
    for label, keep in filters.items():
        selected = [row for row in rows if keep(row)]
        pos_groups: dict[tuple[str, int], set[str]] = defaultdict(set)
        codon_groups: dict[tuple[str, int], set[str]] = defaultdict(set)
        for row in selected:
            pos_groups[(row["gene"], row["cds_pos"])].add(row["group"])
            if row["aa_pos"]:
                codon_groups[(row["gene"], row["aa_pos"])].add(row["group"])
        output.append(
            {
                "filter": label,
                "variant_count": len(selected),
                "positions_total": len(pos_groups),
                "positions_mixed": sum(1 for groups in pos_groups.values() if len(groups) > 1),
                "positions_with_benign_and_pathogenic": sum(1 for groups in pos_groups.values() if {"benign", "pathogenic"} <= groups),
                "codons_total": len(codon_groups),
                "codons_mixed": sum(1 for groups in codon_groups.values() if len(groups) > 1),
                "codons_with_benign_and_pathogenic": sum(1 for groups in codon_groups.values() if {"benign", "pathogenic"} <= groups),
            }
        )
    return output


def stacked_variant_type_svg(rows: list[dict]) -> None:
    width = 920
    height = 360
    left = 190
    top = 58
    bar_h = 36
    gap = 30
    plot_w = 560
    body = ['<text x="28" y="32" class="title">Grouped classes by variant type</text>']
    for i, row in enumerate(rows):
        y = top + i * (bar_h + gap)
        total = int(row["total_count"])
        x = left
        body.append(f'<text x="28" y="{y + 23}" class="label">{esc(row["variant_type"])}</text>')
        for group_name, key in [("benign", "benign_count"), ("vus", "vus_count"), ("pathogenic", "pathogenic_count")]:
            count = int(row[key])
            w = plot_w * count / total if total else 0
            body.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{GROUP_COLORS[group_name]}" opacity="0.84"/>')
            if w > 32:
                body.append(f'<text x="{x + 6:.1f}" y="{y + 23}" class="small" fill="#ffffff">{count}</text>')
            x += w
        body.append(f'<text x="{left + plot_w + 12}" y="{y + 23}" class="small">path frac {row["pathogenic_fraction"]}</text>')
    legend_x = left
    legend_y = height - 34
    for i, group_name in enumerate(GROUP_ORDER):
        x = legend_x + i * 130
        body.append(f'<rect x="{x}" y="{legend_y - 12}" width="14" height="14" fill="{GROUP_COLORS[group_name]}"/>')
        body.append(f'<text x="{x + 20}" y="{legend_y}" class="small">{group_name}</text>')
    write_svg(PLOT_DIR / "variant_type_grouped_class_stacked.svg", "\n".join(body), width, height)


def scatter_svg(rows: list[dict], x_key: str, y_key: str, path: Path, title: str, x_label: str, y_label: str) -> None:
    width = 720
    height = 540
    left = 74
    top = 56
    plot_w = 520
    plot_h = 360
    xs = [float(row[x_key]) for row in rows]
    ys = [float(row[y_key]) for row in rows]
    max_x = max(xs + [0.01])
    max_y = max(ys + [0.01])
    body = [f'<text x="28" y="30" class="title">{esc(title)}</text>']
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        x = left + tick * plot_w
        y = top + plot_h - tick * plot_h
        body.append(f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + plot_h}" class="grid"/>')
        body.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_w}" y2="{y}" class="grid"/>')
        body.append(f'<text x="{x - 8}" y="{top + plot_h + 18}" class="small">{tick:.2f}</text>')
        body.append(f'<text x="38" y="{y + 4}" class="small">{tick:.2f}</text>')
    for row in rows:
        x = left + (float(row[x_key]) / max_x) * plot_w
        y = top + plot_h - (float(row[y_key]) / max_y) * plot_h
        color = "#b91c1c" if row["gene"] == "BRCA1" else "#2563eb"
        body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="{color}" opacity="0.58"/>')
    body.append(f'<text x="{left + 130}" y="{height - 36}" class="label">{esc(x_label)}</text>')
    body.append(f'<text x="18" y="{top + 230}" class="label" transform="rotate(-90 18 {top + 230})">{esc(y_label)}</text>')
    body.append(f'<text x="{left + plot_w + 24}" y="{top + 26}" class="small" fill="#b91c1c">BRCA1</text>')
    body.append(f'<text x="{left + plot_w + 24}" y="{top + 46}" class="small" fill="#2563eb">BRCA2</text>')
    write_svg(path, "\n".join(body), width, height)


def mixed_ablation_svg(rows: list[dict]) -> None:
    width = 980
    height = 420
    left = 90
    top = 70
    plot_w = 760
    plot_h = 250
    metrics = ["positions_mixed", "positions_with_benign_and_pathogenic", "codons_mixed", "codons_with_benign_and_pathogenic"]
    colors = ["#64748b", "#b91c1c", "#7c3aed", "#f97316"]
    max_value = max(int(row[m]) for row in rows for m in metrics)
    group_w = plot_w / len(rows)
    bar_w = group_w / len(metrics) - 3
    body = ['<text x="28" y="32" class="title">Mixed position/codon contexts after driver exclusions</text>']
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    for i, row in enumerate(rows):
        x0 = left + i * group_w
        for mi, metric in enumerate(metrics):
            value = int(row[metric])
            h = plot_h * value / max_value if max_value else 0
            x = x0 + mi * (bar_w + 3) + 4
            y = top + plot_h - h
            body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{colors[mi]}" opacity="0.82"/>')
        body.append(f'<text x="{x0 + 4:.1f}" y="{top + plot_h + 22}" class="small">{esc(row["filter"]).replace("_", " ")}</text>')
    legend_x = left + plot_w - 310
    legend_y = top + 20
    for mi, metric in enumerate(metrics):
        body.append(f'<rect x="{legend_x}" y="{legend_y + mi * 20 - 10}" width="12" height="12" fill="{colors[mi]}"/>')
        body.append(f'<text x="{legend_x + 18}" y="{legend_y + mi * 20}" class="small">{esc(metric)}</text>')
    write_svg(PLOT_DIR / "mixed_position_codon_driver_ablation.svg", "\n".join(body), width, height)


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


def write_report(variant_summary: list[dict], splice_regions: list[dict], functional_regions: list[dict], mixed: list[dict]) -> None:
    text = f"""# Positional Context Follow-up Analyses

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Purpose

These follow-up analyses ask why positional context can look informative in the
precomputed Module 1 map. The goal is to separate regional signals driven by
variant type, splice prediction, functional evidence, and mixed positions or
codons.

This remains a context analysis only. It does not create ACMG/ENIGMA evidence.

## 1. Variant Type Effect

{markdown_table(variant_summary, ["variant_type", "total_count", "benign_count", "vus_count", "pathogenic_count", "pathogenic_fraction", "truncation_driver_count", "splice_driver_count", "functional_pathogenic_driver_count", "functional_benign_driver_count"], 10)}

## 2. Splice-driven Regional Signal

Top bins by splice-driver fraction:

{markdown_table(splice_regions, ["gene", "bin", "cds_start", "cds_end", "pathogenic_fraction", "splice_driver_fraction", "pp3_count", "high_spliceai_count", "near_boundary_count", "total_count"], 20)}

## 3. Functional-evidence Regional Signal

Top bins by functional evidence count:

{markdown_table(functional_regions, ["gene", "bin", "cds_start", "cds_end", "pathogenic_fraction", "ps3_ps1_pp4_count", "bs3_count", "bs3_pp3_count", "ps3_pp3_count", "total_count"], 20)}

## 4. Mixed Position/Codon Context After Driver Exclusions

{markdown_table(mixed, ["filter", "variant_count", "positions_total", "positions_mixed", "positions_with_benign_and_pathogenic", "codons_total", "codons_mixed", "codons_with_benign_and_pathogenic"], 10)}

## Interpretation

- Nonsense/PTC variants strongly explain many pathogenic-enriched regions.
- Splice driver signal identifies bins where local context is tied to `PP3`,
  high SpliceAI, or proximity to CDS exon boundaries.
- Functional evidence explains BRCA1 RING and BRCT mixed regions, especially
  where `PS3` and `BS3` are both common.
- Mixed positions and codons remain useful for review interpretation because a
  coordinate can host multiple specific variant consequences.

The practical conclusion is unchanged: positional context is useful for triage
and explanation, but classification must remain variant-level and criterion
based.

## Outputs

- `tables/positional_context_followup/variant_type_driver_summary.csv`
- `tables/positional_context_followup/regional_driver_ablation_by_type.csv`
- `tables/positional_context_followup/splice_driver_regions.csv`
- `tables/positional_context_followup/functional_driver_regions.csv`
- `tables/positional_context_followup/mixed_context_driver_ablation.csv`
- `plots/26_positional_context_followup/variant_type_grouped_class_stacked.svg`
- `plots/26_positional_context_followup/splice_driver_vs_pathogenic_fraction.svg`
- `plots/26_positional_context_followup/functional_driver_vs_pathogenic_fraction.svg`
- `plots/26_positional_context_followup/mixed_position_codon_driver_ablation.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    variant_summary = variant_type_summary(rows)
    regional_ablation = regional_ablation_by_type(rows)
    splice_regions = splice_driver_regions(rows)
    functional_regions = functional_driver_regions(rows)
    mixed = mixed_context_ablation(rows)

    write_csv(
        TABLE_DIR / "variant_type_driver_summary.csv",
        variant_summary,
        [
            "variant_type",
            "total_count",
            "benign_count",
            "vus_count",
            "pathogenic_count",
            "pathogenic_fraction",
            "truncation_driver_count",
            "splice_driver_count",
            "functional_pathogenic_driver_count",
            "functional_benign_driver_count",
        ],
    )
    write_csv(
        TABLE_DIR / "regional_driver_ablation_by_type.csv",
        regional_ablation,
        [
            "gene",
            "bin",
            "cds_start",
            "cds_end",
            "all_pathogenic_fraction",
            "no_truncation_pathogenic_fraction",
            "no_splice_pathogenic_fraction",
            "no_truncation_no_splice_pathogenic_fraction",
            "no_functional_pathogenic_fraction",
            "total_count",
            "pathogenic_count",
            "truncation_count",
            "splice_count",
            "functional_count",
        ],
    )
    write_csv(TABLE_DIR / "splice_driver_regions.csv", splice_regions, ["gene", "bin", "cds_start", "cds_end", "pathogenic_fraction", "splice_driver_fraction", "pp3_count", "high_spliceai_count", "near_boundary_count", "total_count"])
    write_csv(TABLE_DIR / "functional_driver_regions.csv", functional_regions, ["gene", "bin", "cds_start", "cds_end", "pathogenic_fraction", "ps3_ps1_pp4_count", "bs3_count", "bs3_pp3_count", "ps3_pp3_count", "total_count"])
    write_csv(TABLE_DIR / "mixed_context_driver_ablation.csv", mixed, ["filter", "variant_count", "positions_total", "positions_mixed", "positions_with_benign_and_pathogenic", "codons_total", "codons_mixed", "codons_with_benign_and_pathogenic"])

    stacked_variant_type_svg(variant_summary)
    scatter_svg(
        splice_regions,
        "splice_driver_fraction",
        "pathogenic_fraction",
        PLOT_DIR / "splice_driver_vs_pathogenic_fraction.svg",
        "Splice driver fraction versus regional pathogenic fraction",
        "Splice driver fraction",
        "Pathogenic fraction",
    )
    functional_scatter_rows = [
        {
            **row,
            "functional_driver_fraction": f"{(int(row['ps3_ps1_pp4_count']) + int(row['bs3_count'])) / int(row['total_count']):.4f}",
        }
        for row in functional_regions
    ]
    scatter_svg(
        functional_scatter_rows,
        "functional_driver_fraction",
        "pathogenic_fraction",
        PLOT_DIR / "functional_driver_vs_pathogenic_fraction.svg",
        "Functional evidence fraction versus regional pathogenic fraction",
        "Functional evidence fraction",
        "Pathogenic fraction",
    )
    mixed_ablation_svg(mixed)
    write_report(variant_summary, splice_regions, functional_regions, mixed)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
