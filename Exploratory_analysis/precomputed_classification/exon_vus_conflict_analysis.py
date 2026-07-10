"""Exon-level summary, VUS prioritization, and conflict sanity checks."""

from __future__ import annotations

import csv
import json
import math
import re
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
TRANSCRIPT_MAPS = ROOT / "variant_space_scan" / "coordinate_mapping" / "data" / "transcript_maps.json"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = OUT_DIR / "tables" / "exon_vus_conflict"
PLOT_DIR = OUT_DIR / "plots" / "06_exon_vus_conflict"
REPORT = OUT_DIR / "exon_vus_conflict_report.md"

CDS_LENGTHS = {"BRCA1": 5592, "BRCA2": 10257}
C_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")
GROUP_ORDER = ["benign", "vus", "pathogenic"]
GROUP_COLORS = {
    "benign": "#2b8cbe",
    "vus": "#8a8f9c",
    "pathogenic": "#d7301f",
}
GROUP_LABELS = {
    "benign": "benign 1/2",
    "vus": "VUS 3",
    "pathogenic": "pathogenic 4/5",
}


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


def parse_criteria(criteria: str) -> set[str]:
    if not criteria:
        return set()
    return {item.split(":", 1)[0] for item in criteria.split(";") if item}


def class_group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


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


def write_svg(path: Path, body: str, width: int, height: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
text {{ font-family: Arial, Helvetica, sans-serif; fill: #172033; }}
.small {{ font-size: 11px; }}
.label {{ font-size: 12px; }}
.title {{ font-size: 18px; font-weight: 700; }}
.grid {{ stroke: #e2e8f0; stroke-width: 1; }}
.axis {{ stroke: #9aa4b2; stroke-width: 1; }}
.panel {{ fill: #f3f4f6; }}
.tile {{ stroke: #ffffff; stroke-width: 1; }}
.legend-title {{ font-size: 11px; font-weight: 700; }}
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_exon_maps() -> dict[str, dict[str, object]]:
    data = json.loads(TRANSCRIPT_MAPS.read_text(encoding="utf-8"))
    maps = {}
    for gene, info in data["assemblies"]["GRCh38"].items():
        cds_features = list(info["cds"])
        if info["strand"] == "+":
            ordered = sorted(cds_features, key=lambda item: item["start"])
        else:
            ordered = sorted(cds_features, key=lambda item: item["start"], reverse=True)
        exons = []
        pos = 1
        boundaries = []
        for idx, feature in enumerate(ordered, start=1):
            length = feature["end"] - feature["start"] + 1
            start = pos
            end = pos + length - 1
            exons.append(
                {
                    "gene": gene,
                    "cds_exon": idx,
                    "cds_start": start,
                    "cds_end": end,
                    "length": length,
                    "genomic_start": feature["start"],
                    "genomic_end": feature["end"],
                    "strand": feature["strand"],
                }
            )
            boundaries.extend([start, end])
            pos = end + 1
        maps[gene] = {"exons": exons, "boundaries": sorted(set(boundaries))}
    return maps


def annotate_exon(row: dict[str, object], exon_maps: dict[str, dict[str, object]]) -> None:
    gene = str(row["gene"])
    pos = row.get("cds_pos")
    if not isinstance(pos, int):
        row["cds_exon"] = ""
        row["exon_cds_start"] = ""
        row["exon_cds_end"] = ""
        row["boundary_distance"] = ""
        return
    exons = exon_maps[gene]["exons"]
    exon = next((item for item in exons if item["cds_start"] <= pos <= item["cds_end"]), None)
    row["cds_exon"] = exon["cds_exon"] if exon else ""
    row["exon_cds_start"] = exon["cds_start"] if exon else ""
    row["exon_cds_end"] = exon["cds_end"] if exon else ""
    row["boundary_distance"] = min(abs(pos - boundary) for boundary in exon_maps[gene]["boundaries"])


def load_rows(exon_maps: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["cds_pos"] = parse_cds_pos(row["c_notation"])
            row["spliceai"] = parse_float(row["spliceai_score"])
            row["criteria_codes"] = parse_criteria(row.get("criteria", ""))
            row["group"] = class_group(row["predicted_class"])
            annotate_exon(row, exon_maps)
            rows.append(row)
    return rows


def exon_summary(rows: list[dict[str, object]], exon_maps: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    summary = []
    for gene in ["BRCA1", "BRCA2"]:
        for exon in exon_maps[gene]["exons"]:
            items = [row for row in rows if row["gene"] == gene and row["cds_exon"] == exon["cds_exon"]]
            total = len(items)
            if total == 0:
                continue
            record = dict(exon)
            record["n_variants"] = total
            for group in GROUP_ORDER:
                count = sum(1 for row in items if row["group"] == group)
                record[f"{group}_count"] = count
                record[f"{group}_percent"] = f"{count / total * 100:.2f}"
            record["high_spliceai_count"] = sum(1 for row in items if row["spliceai"] >= 0.20)
            record["high_spliceai_percent"] = f"{record['high_spliceai_count'] / total * 100:.2f}"
            for criterion in ["PVS1", "PM5_PTC", "PP3", "PS3", "BS3", "BP1", "BP7"]:
                count = sum(1 for row in items if criterion in row["criteria_codes"])
                record[f"{criterion.lower()}_count"] = count
                record[f"{criterion.lower()}_percent"] = f"{count / total * 100:.2f}"
            summary.append(record)
    return summary


def build_context_index(rows: list[dict[str, object]]) -> dict[str, dict[str, list[int]]]:
    index = {}
    for gene in ["BRCA1", "BRCA2"]:
        gene_rows = [row for row in rows if row["gene"] == gene and isinstance(row.get("cds_pos"), int)]
        index[gene] = {
            "all": sorted(int(row["cds_pos"]) for row in gene_rows),
            "pathogenic": sorted(int(row["cds_pos"]) for row in gene_rows if row["group"] == "pathogenic"),
            "benign": sorted(int(row["cds_pos"]) for row in gene_rows if row["group"] == "benign"),
            "vus": sorted(int(row["cds_pos"]) for row in gene_rows if row["group"] == "vus"),
        }
    return index


def count_in_window(values: list[int], pos: int, window: int) -> int:
    return bisect_right(values, pos + window) - bisect_left(values, pos - window)


def local_context(index: dict[str, dict[str, list[int]]], row: dict[str, object], window: int = 20) -> dict[str, int]:
    gene = row["gene"]
    pos = row["cds_pos"]
    if not isinstance(pos, int):
        return {"near_total": 0, "near_pathogenic": 0, "near_benign": 0, "near_vus": 0}
    gene_index = index[str(gene)]
    return {
        "near_total": count_in_window(gene_index["all"], pos, window),
        "near_pathogenic": count_in_window(gene_index["pathogenic"], pos, window),
        "near_benign": count_in_window(gene_index["benign"], pos, window),
        "near_vus": count_in_window(gene_index["vus"], pos, window),
    }


def vus_priority(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    context_index = build_context_index(rows)
    prioritized = []
    for row in rows:
        if row["group"] != "vus":
            continue
        codes = row["criteria_codes"]
        reasons = []
        score = 0
        if row["spliceai"] >= 0.20:
            score += 30
            reasons.append("SpliceAI>=0.20")
        elif row["spliceai"] >= 0.10:
            score += 12
            reasons.append("SpliceAI>=0.10")
        if isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 2:
            score += 20
            reasons.append("within_2bp_boundary")
        elif isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 10:
            score += 8
            reasons.append("within_10bp_boundary")
        for criterion, points in [("PVS1", 35), ("PS3", 25), ("PS1", 18), ("PP3", 15), ("PP4", 15), ("BS3", 8)]:
            if criterion in codes:
                score += points
                reasons.append(criterion)
        context = local_context(context_index, row)
        if context["near_pathogenic"] >= 5:
            score += 12
            reasons.append("near_pathogenic_dense")
        if context["near_benign"] >= 20 and row["spliceai"] >= 0.20:
            score += 8
            reasons.append("high_spliceai_in_benign_dense_region")
        output = {
            "priority_score": score,
            "priority_reasons": ";".join(reasons),
            "gene": row["gene"],
            "c_notation": row["c_notation"],
            "p_notation": row["p_notation"],
            "cds_pos": row["cds_pos"],
            "cds_exon": row["cds_exon"],
            "boundary_distance": row["boundary_distance"],
            "variant_type": row["variant_type"],
            "spliceai_score": f"{row['spliceai']:.3f}",
            "criteria": row["criteria"],
            "total_points": row["total_points"],
            "predicted_class": row["predicted_class"],
            "near_pathogenic_20bp": context["near_pathogenic"],
            "near_benign_20bp": context["near_benign"],
            "near_vus_20bp": context["near_vus"],
        }
        prioritized.append(output)
    return sorted(prioritized, key=lambda item: (int(item["priority_score"]), item["spliceai_score"]), reverse=True)


def conflict_checks(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    conflicts = []
    for row in rows:
        codes = row["criteria_codes"]
        reasons = []
        if row["group"] == "benign" and row["spliceai"] >= 0.20:
            reasons.append("benign_with_high_spliceai")
        if row["group"] == "benign" and ("PVS1" in codes or "PS3" in codes or "PP3" in codes):
            reasons.append("benign_with_pathogenic_evidence")
        if row["group"] == "pathogenic" and ("BS3" in codes or "BP1" in codes or "BP7" in codes):
            reasons.append("pathogenic_with_benign_evidence")
        if row["group"] == "vus" and "PVS1" in codes:
            reasons.append("vus_with_pvs1")
        if row["group"] == "vus" and row["spliceai"] >= 0.20 and isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] > 100:
            reasons.append("vus_high_spliceai_far_from_boundary")
        if row["group"] == "benign" and row["spliceai"] >= 0.20 and isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 2:
            reasons.append("benign_high_spliceai_at_boundary")
        if not reasons:
            continue
        conflicts.append(
            {
                "conflict_reasons": ";".join(reasons),
                "gene": row["gene"],
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "cds_pos": row["cds_pos"],
                "cds_exon": row["cds_exon"],
                "boundary_distance": row["boundary_distance"],
                "variant_type": row["variant_type"],
                "spliceai_score": f"{row['spliceai']:.3f}",
                "criteria": row["criteria"],
                "total_points": row["total_points"],
                "predicted_class": row["predicted_class"],
                "predicted_label": row["predicted_label"],
            }
        )
    return sorted(conflicts, key=lambda item: (item["conflict_reasons"], -float(item["spliceai_score"])))


def draw_exon_group_heatmap(summary: list[dict[str, object]]) -> None:
    for gene in ["BRCA1", "BRCA2"]:
        rows = [row for row in summary if row["gene"] == gene]
        cell_w = 32
        row_h = 34
        left = 122
        top = 104
        panel_w = len(rows) * cell_w
        panel_h = len(GROUP_ORDER) * row_h
        width = left + panel_w + 130
        height = top + panel_h + 92
        body = [f'<text x="28" y="34" class="title">{gene} exon-level generated class mix</text>']
        body.append(f'<text x="28" y="56" class="small">Cells show percent within exon. This is a descriptive Module 1 map, not evidence for reclassification.</text>')
        body.append(f'<text x="{left}" y="82" class="small">CDS exon</text>')
        body.append(f'<rect x="{left - 2}" y="{top - 2}" width="{panel_w + 2}" height="{panel_h + 2}" class="panel"/>')
        for xi, row in enumerate(rows):
            x = left + xi * cell_w
            body.append(f'<text x="{x + 4}" y="{top - 14}" class="small">E{row["cds_exon"]}</text>')
        for yi, group in enumerate(GROUP_ORDER):
            y = top + yi * row_h
            body.append(f'<text x="20" y="{y + 22}" class="label">{GROUP_LABELS[group]}</text>')
            for xi, row in enumerate(rows):
                value = float(row[f"{group}_percent"])
                x = left + xi * cell_w
                base = tuple(int(GROUP_COLORS[group].lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
                body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h}" fill="{heat_color(value, 100, base)}" class="tile"/>')
                if value >= 20:
                    body.append(f'<text x="{x + 4}" y="{y + 22}" class="small">{value:.0f}</text>')
        legend_x = left + panel_w + 34
        legend_y = top
        body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">percent</text>')
        for step in range(6):
            frac = 1 - step / 5
            color = heat_color(frac * 100, 100, (37, 99, 235))
            y = legend_y + step * 18
            body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{color}" class="tile"/>')
        body.append(f'<text x="{legend_x + 24}" y="{legend_y + 12}" class="small">100</text>')
        body.append(f'<text x="{legend_x + 24}" y="{legend_y + 102}" class="small">0</text>')
        write_svg(PLOT_DIR / f"{gene.lower()}_exon_group_heatmap.svg", "\n".join(body), width, height)


def draw_exon_signal_heatmap(summary: list[dict[str, object]]) -> None:
    metrics = [
        ("VUS", "vus_percent"),
        ("pathogenic 4/5", "pathogenic_percent"),
        ("SpAI >=0.20", "high_spliceai_percent"),
        ("PVS1", "pvs1_percent"),
        ("PP3", "pp3_percent"),
        ("PS3", "ps3_percent"),
        ("BS3", "bs3_percent"),
        ("BP1", "bp1_percent"),
    ]
    for gene in ["BRCA1", "BRCA2"]:
        rows = [row for row in summary if row["gene"] == gene]
        cell_w = 32
        row_h = 30
        left = 126
        top = 112
        panel_w = len(rows) * cell_w
        panel_h = len(metrics) * row_h
        width = left + panel_w + 132
        height = top + panel_h + 88
        max_value = max(float(row[col]) for row in rows for _, col in metrics)
        body = [f'<text x="28" y="34" class="title">{gene} exon-level Module 1 signals</text>']
        body.append(f'<text x="28" y="56" class="small">Cells show percent within exon for selected generated signals and criteria.</text>')
        body.append(f'<text x="{left}" y="84" class="small">CDS exon</text>')
        body.append(f'<rect x="{left - 2}" y="{top - 2}" width="{panel_w + 2}" height="{panel_h + 2}" class="panel"/>')
        for xi, row in enumerate(rows):
            x = left + xi * cell_w
            body.append(f'<text x="{x + 4}" y="{top - 14}" class="small">E{row["cds_exon"]}</text>')
        for yi, (label, col) in enumerate(metrics):
            y = top + yi * row_h
            body.append(f'<text x="18" y="{y + 20}" class="small">{esc(label)}</text>')
            for xi, row in enumerate(rows):
                x = left + xi * cell_w
                value = float(row[col])
                body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h}" fill="{heat_color(value, max_value, (37, 99, 235))}" class="tile"/>')
                if value >= 20:
                    body.append(f'<text x="{x + 4}" y="{y + 19}" class="small">{value:.0f}</text>')
        legend_x = left + panel_w + 34
        legend_y = top
        body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">percent</text>')
        for step in range(6):
            frac = 1 - step / 5
            color = heat_color(frac * max_value, max_value, (37, 99, 235))
            y = legend_y + step * 18
            body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{color}" class="tile"/>')
        body.append(f'<text x="{legend_x + 24}" y="{legend_y + 12}" class="small">{max_value:.1f}</text>')
        body.append(f'<text x="{legend_x + 24}" y="{legend_y + 102}" class="small">0</text>')
        write_svg(PLOT_DIR / f"{gene.lower()}_exon_signal_heatmap.svg", "\n".join(body), width, height)


def draw_vus_score_barplot(priority_rows: list[dict[str, object]]) -> None:
    rows = priority_rows[:25]
    width = 1100
    left = 310
    top = 88
    bar_h = 20
    gap = 8
    plot_w = 680
    height = top + len(rows) * (bar_h + gap) + 60
    max_score = max((int(row["priority_score"]) for row in rows), default=1)
    body = ['<text x="28" y="34" class="title">Highest heuristic VUS scores</text>']
    body.append('<text x="28" y="56" class="small">Heuristic score examples only. Not a validated manual-review priority criterion.</text>')
    for i, row in enumerate(rows):
        y = top + i * (bar_h + gap)
        label = f"{row['gene']} {row['c_notation']} exon {row['cds_exon']}"
        score = int(row["priority_score"])
        w = plot_w * score / max_score if max_score else 0
        body.append(f'<text x="24" y="{y + 15}" class="small">{esc(label)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="#7c3aed"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 15}" class="small">{score}: {esc(row["priority_reasons"][:80])}</text>')
    svg_body = "\n".join(body)
    write_svg(PLOT_DIR / "heuristic_vus_score_examples.svg", svg_body, width, height)
    write_svg(PLOT_DIR / "top_vus_priority_barplot.svg", svg_body, width, height)


def draw_conflict_reason_barplot(conflicts: list[dict[str, object]]) -> None:
    reason_gene_counts = Counter()
    reason_counts = Counter()
    for row in conflicts:
        for reason in str(row["conflict_reasons"]).split(";"):
            reason_counts[reason] += 1
            reason_gene_counts[(reason, row["gene"])] += 1
    reasons = [reason for reason, _ in reason_counts.most_common()]
    width = 980
    left = 300
    top = 86
    bar_h = 24
    gap = 12
    plot_w = 560
    height = top + len(reasons) * (bar_h + gap) + 82
    max_value = max(reason_counts.values(), default=1)
    body = ['<text x="28" y="34" class="title">Conflict sanity-check patterns</text>']
    body.append('<text x="28" y="56" class="small">Counts of generated patterns that should be interpreted cautiously, not classification errors by default.</text>')
    for idx, reason in enumerate(reasons):
        y = top + idx * (bar_h + gap)
        brca1 = reason_gene_counts[(reason, "BRCA1")]
        brca2 = reason_gene_counts[(reason, "BRCA2")]
        total = brca1 + brca2
        w1 = plot_w * brca1 / max_value if max_value else 0
        w2 = plot_w * brca2 / max_value if max_value else 0
        body.append(f'<text x="24" y="{y + 17}" class="small">{esc(reason)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w1:.1f}" height="{bar_h}" fill="#5578c2"/>')
        body.append(f'<rect x="{left + w1:.1f}" y="{y}" width="{w2:.1f}" height="{bar_h}" fill="#c26a55"/>')
        body.append(f'<text x="{left + w1 + w2 + 8:.1f}" y="{y + 17}" class="small">{total}</text>')
    legend_y = height - 42
    body.append(f'<rect x="{left}" y="{legend_y}" width="14" height="14" fill="#5578c2"/>')
    body.append(f'<text x="{left + 20}" y="{legend_y + 12}" class="small">BRCA1</text>')
    body.append(f'<rect x="{left + 90}" y="{legend_y}" width="14" height="14" fill="#c26a55"/>')
    body.append(f'<text x="{left + 110}" y="{legend_y + 12}" class="small">BRCA2</text>')
    write_svg(PLOT_DIR / "conflict_sanity_check_patterns.svg", "\n".join(body), width, height)


def report_table(rows: list[dict[str, object]], columns: list[str], limit: int = 10) -> str:
    if not rows:
        return "| none |\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows[:limit]:
        body.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join([header, sep] + body)


def write_report(summary: list[dict[str, object]], priority_rows: list[dict[str, object]], conflicts: list[dict[str, object]]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    high_vus = [row for row in priority_rows if int(row["priority_score"]) >= 50]
    conflict_counts = Counter()
    for row in conflicts:
        for reason in str(row["conflict_reasons"]).split(";"):
            conflict_counts[reason] += 1
    top_exons_vus = sorted(summary, key=lambda row: float(row["vus_percent"]), reverse=True)[:10]
    top_exons_splice = sorted(summary, key=lambda row: float(row["high_spliceai_percent"]), reverse=True)[:10]
    conflict_lines = "\n".join(f"- {reason}: {count}" for reason, count in conflict_counts.most_common())
    text = f"""# Exon, VUS Priority, and Conflict Analysis

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

This report is based on the automated Module 1 coding SNV snapshot. It is a
descriptive exon-level audit of generated class mix, selected Module 1 signals,
and sanity-check patterns. It is not final expert classification and it is not
a validated manual-review prioritization criterion.

## Outputs

- `tables/exon_vus_conflict/exon_level_summary.csv`
- `tables/exon_vus_conflict/vus_priority_list.csv`
- `tables/exon_vus_conflict/conflict_sanity_checks.csv`
- `plots/06_exon_vus_conflict/brca1_exon_group_heatmap.svg`
- `plots/06_exon_vus_conflict/brca2_exon_group_heatmap.svg`
- `plots/06_exon_vus_conflict/brca1_exon_signal_heatmap.svg`
- `plots/06_exon_vus_conflict/brca2_exon_signal_heatmap.svg`
- `plots/06_exon_vus_conflict/heuristic_vus_score_examples.svg`
- `plots/06_exon_vus_conflict/conflict_sanity_check_patterns.svg`

## Summary

- Exon summary rows: {len(summary)}
- VUS variants with heuristic score: {len(priority_rows)}
- VUS variants with priority score >= 50: {len(high_vus)}
- Conflict sanity-check rows: {len(conflicts)}

## Conflict Types

{conflict_lines if conflict_lines else "- none"}

## Top Exons by VUS Percent

{report_table(top_exons_vus, ["gene", "cds_exon", "cds_start", "cds_end", "n_variants", "vus_count", "vus_percent"])}

## Top Exons by High SpliceAI Percent

{report_table(top_exons_splice, ["gene", "cds_exon", "cds_start", "cds_end", "n_variants", "high_spliceai_count", "high_spliceai_percent"])}

## Highest Heuristic VUS Score Examples

{report_table(priority_rows, ["priority_score", "priority_reasons", "gene", "c_notation", "p_notation", "cds_exon", "boundary_distance", "spliceai_score"], limit=15)}

## Top Conflict Sanity Checks

{report_table(conflicts, ["conflict_reasons", "gene", "c_notation", "p_notation", "cds_exon", "boundary_distance", "spliceai_score", "predicted_class"], limit=15)}
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    exon_maps = build_exon_maps()
    rows = load_rows(exon_maps)
    summary = exon_summary(rows, exon_maps)
    priority_rows = vus_priority(rows)
    conflicts = conflict_checks(rows)
    write_csv(
        TABLE_DIR / "exon_level_summary.csv",
        summary,
        [
            "gene",
            "cds_exon",
            "cds_start",
            "cds_end",
            "length",
            "genomic_start",
            "genomic_end",
            "strand",
            "n_variants",
            "benign_count",
            "benign_percent",
            "vus_count",
            "vus_percent",
            "pathogenic_count",
            "pathogenic_percent",
            "high_spliceai_count",
            "high_spliceai_percent",
            "pvs1_count",
            "pvs1_percent",
            "pm5_ptc_count",
            "pm5_ptc_percent",
            "pp3_count",
            "pp3_percent",
            "ps3_count",
            "ps3_percent",
            "bs3_count",
            "bs3_percent",
            "bp1_count",
            "bp1_percent",
            "bp7_count",
            "bp7_percent",
        ],
    )
    write_csv(
        TABLE_DIR / "vus_priority_list.csv",
        priority_rows,
        [
            "priority_score",
            "priority_reasons",
            "gene",
            "c_notation",
            "p_notation",
            "cds_pos",
            "cds_exon",
            "boundary_distance",
            "variant_type",
            "spliceai_score",
            "criteria",
            "total_points",
            "predicted_class",
            "near_pathogenic_20bp",
            "near_benign_20bp",
            "near_vus_20bp",
        ],
    )
    write_csv(
        TABLE_DIR / "conflict_sanity_checks.csv",
        conflicts,
        [
            "conflict_reasons",
            "gene",
            "c_notation",
            "p_notation",
            "cds_pos",
            "cds_exon",
            "boundary_distance",
            "variant_type",
            "spliceai_score",
            "criteria",
            "total_points",
            "predicted_class",
            "predicted_label",
        ],
    )
    draw_exon_group_heatmap(summary)
    draw_exon_signal_heatmap(summary)
    draw_vus_score_barplot(priority_rows)
    draw_conflict_reason_barplot(conflicts)
    write_report(summary, priority_rows, conflicts)
    print(f"Wrote exon/VUS/conflict analysis to {REPORT}")


if __name__ == "__main__":
    main()
