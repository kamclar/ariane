"""Detailed splice-boundary enrichment analysis for coding SNVs."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
TRANSCRIPT_MAPS = ROOT / "variant_space_scan" / "coordinate_mapping" / "data" / "transcript_maps.json"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = OUT_DIR / "tables" / "boundary_enrichment"
PLOT_DIR = OUT_DIR / "plots" / "09_boundary_enrichment"
REPORT = OUT_DIR / "boundary_enrichment_report.md"

CDS_LENGTHS = {"BRCA1": 5592, "BRCA2": 10257}
C_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")
GROUP_ORDER = ["benign", "vus", "pathogenic"]
GROUP_LABELS = {
    "benign": "Benign/Likely Benign",
    "vus": "VUS",
    "pathogenic": "Likely Pathogenic/Pathogenic",
}
GROUP_COLORS = {
    "benign": "#2b8cbe",
    "vus": "#8a8f9c",
    "pathogenic": "#d7301f",
}
DISTANCE_BINS = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"]
SITE_TYPES = ["acceptor_like", "donor_like"]
SPLICEAI_THRESHOLDS = [0.10, 0.20, 0.50, 0.80]
RELATIVE_POSITIONS = list(range(-20, 21))
VARIANT_TYPES = ["missense", "nonsense", "synonymous", "initiation_codon"]


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


def class_group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


def distance_bin(distance: int | None) -> str:
    if distance is None:
        return "missing"
    if distance == 0:
        return "0"
    if distance <= 2:
        return "1-2"
    if distance <= 5:
        return "3-5"
    if distance <= 10:
        return "6-10"
    if distance <= 20:
        return "11-20"
    if distance <= 50:
        return "21-50"
    if distance <= 100:
        return "51-100"
    return ">100"


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


def build_splice_boundaries() -> dict[str, list[dict[str, object]]]:
    data = json.loads(TRANSCRIPT_MAPS.read_text(encoding="utf-8"))
    output = {}
    for gene, info in data["assemblies"]["GRCh38"].items():
        cds_features = list(info["cds"])
        if info["strand"] == "+":
            ordered = sorted(cds_features, key=lambda item: item["start"])
        else:
            ordered = sorted(cds_features, key=lambda item: item["start"], reverse=True)
        pos = 1
        exons = []
        for idx, feature in enumerate(ordered, start=1):
            length = feature["end"] - feature["start"] + 1
            start = pos
            end = pos + length - 1
            exons.append({"cds_exon": idx, "cds_start": start, "cds_end": end})
            pos = end + 1
        boundaries = []
        n = len(exons)
        for exon in exons:
            idx = int(exon["cds_exon"])
            if idx > 1:
                boundaries.append(
                    {
                        "gene": gene,
                        "cds_exon": idx,
                        "site_type": "acceptor_like",
                        "cds_pos": exon["cds_start"],
                    }
                )
            if idx < n:
                boundaries.append(
                    {
                        "gene": gene,
                        "cds_exon": idx,
                        "site_type": "donor_like",
                        "cds_pos": exon["cds_end"],
                    }
                )
        output[gene] = boundaries
    return output


def nearest_boundary(gene: str, pos: int, boundaries: dict[str, list[dict[str, object]]]) -> dict[str, object]:
    best = min(boundaries[gene], key=lambda item: abs(int(item["cds_pos"]) - pos))
    distance = abs(int(best["cds_pos"]) - pos)
    return {
        "nearest_site_type": best["site_type"],
        "nearest_boundary_cds_pos": best["cds_pos"],
        "nearest_boundary_exon": best["cds_exon"],
        "relative_boundary_position": pos - int(best["cds_pos"]),
        "splice_boundary_distance": distance,
        "splice_boundary_distance_bin": distance_bin(distance),
    }


def load_rows(boundaries: dict[str, list[dict[str, object]]]) -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["cds_pos"] = parse_cds_pos(row["c_notation"])
            row["spliceai"] = parse_float(row["spliceai_score"])
            row["group"] = class_group(row["predicted_class"])
            if isinstance(row["cds_pos"], int):
                row.update(nearest_boundary(row["gene"], int(row["cds_pos"]), boundaries))
            rows.append(row)
    return rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for gene in ["BRCA1", "BRCA2"]:
        for site_type in SITE_TYPES:
            for group in GROUP_ORDER:
                for bin_name in DISTANCE_BINS:
                    items = [
                        row
                        for row in rows
                        if row["gene"] == gene
                        and row["group"] == group
                        and row.get("nearest_site_type") == site_type
                        and row.get("splice_boundary_distance_bin") == bin_name
                    ]
                    total = len(items)
                    high = sum(1 for row in items if float(row["spliceai"]) >= 0.20)
                    output.append(
                        {
                            "gene": gene,
                            "site_type": site_type,
                            "group": group,
                            "group_label": GROUP_LABELS[group],
                            "distance_bin": bin_name,
                            "total_variants": total,
                            "spliceai_ge_0_20": high,
                            "percent": f"{high / total * 100:.2f}" if total else "0.00",
                        }
                    )
    return output


def site_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for gene in ["BRCA1", "BRCA2"]:
        for site_type in SITE_TYPES:
            items = [row for row in rows if row["gene"] == gene and row.get("nearest_site_type") == site_type]
            total = len(items)
            high = sum(1 for row in items if float(row["spliceai"]) >= 0.20)
            record = {
                "gene": gene,
                "site_type": site_type,
                "total_variants": total,
                "spliceai_ge_0_20": high,
                "percent": f"{high / total * 100:.2f}" if total else "0.00",
            }
            for group in GROUP_ORDER:
                group_items = [row for row in items if row["group"] == group]
                group_high = sum(1 for row in group_items if float(row["spliceai"]) >= 0.20)
                record[f"{group}_total"] = len(group_items)
                record[f"{group}_high_spliceai"] = group_high
                record[f"{group}_percent"] = f"{group_high / len(group_items) * 100:.2f}" if group_items else "0.00"
            output.append(record)
    return output


def threshold_sensitivity(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for gene in ["BRCA1", "BRCA2"]:
        for site_type in SITE_TYPES:
            for group in GROUP_ORDER:
                for bin_name in DISTANCE_BINS:
                    items = [
                        row
                        for row in rows
                        if row["gene"] == gene
                        and row["group"] == group
                        and row.get("nearest_site_type") == site_type
                        and row.get("splice_boundary_distance_bin") == bin_name
                    ]
                    total = len(items)
                    record = {
                        "gene": gene,
                        "site_type": site_type,
                        "group": group,
                        "distance_bin": bin_name,
                        "total_variants": total,
                    }
                    for threshold in SPLICEAI_THRESHOLDS:
                        high = sum(1 for row in items if float(row["spliceai"]) >= threshold)
                        key = f"spliceai_ge_{str(threshold).replace('.', '_')}"
                        record[f"{key}_count"] = high
                        record[f"{key}_percent"] = f"{high / total * 100:.2f}" if total else "0.00"
                    output.append(record)
    return output


def relative_position_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for gene in ["BRCA1", "BRCA2"]:
        for site_type in SITE_TYPES:
            for group in GROUP_ORDER:
                for rel in RELATIVE_POSITIONS:
                    items = [
                        row
                        for row in rows
                        if row["gene"] == gene
                        and row["group"] == group
                        and row.get("nearest_site_type") == site_type
                        and row.get("relative_boundary_position") == rel
                    ]
                    total = len(items)
                    high = sum(1 for row in items if float(row["spliceai"]) >= 0.20)
                    output.append(
                        {
                            "gene": gene,
                            "site_type": site_type,
                            "group": group,
                            "relative_boundary_position": rel,
                            "total_variants": total,
                            "spliceai_ge_0_20": high,
                            "percent": f"{high / total * 100:.2f}" if total else "0.00",
                        }
                    )
    return output


def exon_boundary_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    by_key = sorted({(row["gene"], row.get("nearest_boundary_exon"), row.get("nearest_site_type")) for row in rows if row.get("nearest_boundary_exon")})
    for gene, exon, site_type in by_key:
        near = [
            row
            for row in rows
            if row["gene"] == gene
            and row.get("nearest_boundary_exon") == exon
            and row.get("nearest_site_type") == site_type
            and isinstance(row.get("splice_boundary_distance"), int)
            and int(row["splice_boundary_distance"]) <= 10
        ]
        total = len(near)
        high = sum(1 for row in near if float(row["spliceai"]) >= 0.20)
        record = {
            "gene": gene,
            "cds_exon": exon,
            "site_type": site_type,
            "near_10bp_total": total,
            "spliceai_ge_0_20": high,
            "percent": f"{high / total * 100:.2f}" if total else "0.00",
        }
        for group in GROUP_ORDER:
            group_items = [row for row in near if row["group"] == group]
            group_high = sum(1 for row in group_items if float(row["spliceai"]) >= 0.20)
            record[f"{group}_total"] = len(group_items)
            record[f"{group}_high_spliceai"] = group_high
            record[f"{group}_percent"] = f"{group_high / len(group_items) * 100:.2f}" if group_items else "0.00"
        output.append(record)
    return output


def variant_type_boundary_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for gene in ["BRCA1", "BRCA2"]:
        for site_type in SITE_TYPES:
            for vtype in VARIANT_TYPES:
                for bin_name in DISTANCE_BINS:
                    items = [
                        row
                        for row in rows
                        if row["gene"] == gene
                        and row["variant_type"] == vtype
                        and row.get("nearest_site_type") == site_type
                        and row.get("splice_boundary_distance_bin") == bin_name
                    ]
                    total = len(items)
                    high = sum(1 for row in items if float(row["spliceai"]) >= 0.20)
                    output.append(
                        {
                            "gene": gene,
                            "site_type": site_type,
                            "variant_type": vtype,
                            "distance_bin": bin_name,
                            "total_variants": total,
                            "spliceai_ge_0_20": high,
                            "percent": f"{high / total * 100:.2f}" if total else "0.00",
                        }
                    )
    return output


def boundary_conflicts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for row in rows:
        distance = row.get("splice_boundary_distance")
        if not isinstance(distance, int):
            continue
        reasons = []
        if row["group"] == "benign" and float(row["spliceai"]) >= 0.20 and distance <= 2:
            reasons.append("benign_high_spliceai_near_boundary")
        if row["group"] == "vus" and float(row["spliceai"]) >= 0.80 and distance <= 2:
            reasons.append("vus_very_high_spliceai_near_boundary")
        if row["group"] == "pathogenic" and float(row["spliceai"]) < 0.10 and distance <= 2:
            reasons.append("pathogenic_low_spliceai_near_boundary")
        if row["group"] == "vus" and float(row["spliceai"]) >= 0.20 and distance > 100:
            reasons.append("vus_high_spliceai_far_from_boundary")
        if not reasons:
            continue
        output.append(
            {
                "boundary_conflict_reason": ";".join(reasons),
                "gene": row["gene"],
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "variant_type": row["variant_type"],
                "group": row["group"],
                "predicted_class": row["predicted_class"],
                "spliceai_score": row["spliceai_score"],
                "nearest_site_type": row.get("nearest_site_type", ""),
                "nearest_boundary_exon": row.get("nearest_boundary_exon", ""),
                "relative_boundary_position": row.get("relative_boundary_position", ""),
                "splice_boundary_distance": row.get("splice_boundary_distance", ""),
                "criteria": row.get("criteria", ""),
            }
        )
    return sorted(output, key=lambda row: (row["boundary_conflict_reason"], -float(row["spliceai_score"] or 0)))


def strongest_boundary_variants(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    selected = [
        row
        for row in rows
        if row["group"] in {"vus", "pathogenic"}
        and isinstance(row.get("splice_boundary_distance"), int)
        and int(row["splice_boundary_distance"]) <= 2
        and float(row["spliceai"]) >= 0.20
    ]
    selected = sorted(selected, key=lambda row: (float(row["spliceai"]), row["group"] == "pathogenic"), reverse=True)
    output = []
    for row in selected:
        output.append(
            {
                "gene": row["gene"],
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "variant_type": row["variant_type"],
                "group": row["group"],
                "predicted_class": row["predicted_class"],
                "spliceai_score": row["spliceai_score"],
                "nearest_site_type": row.get("nearest_site_type", ""),
                "nearest_boundary_exon": row.get("nearest_boundary_exon", ""),
                "relative_boundary_position": row.get("relative_boundary_position", ""),
                "splice_boundary_distance": row.get("splice_boundary_distance", ""),
                "criteria": row.get("criteria", ""),
            }
        )
    return output


def draw_site_type_rate_plot(summary_rows: list[dict[str, object]]) -> None:
    width = 1180
    height = 520
    left = 90
    top = 80
    panel_w = 460
    panel_h = 300
    panel_gap = 80
    bar_gap = 12
    body = ['<text x="28" y="34" class="title">SpliceAI >=0.20 rate by donor-like and acceptor-like boundary</text>']
    body.append('<text x="28" y="55" class="small">Rates are grouped by nearest internal CDS splice boundary and grouped class</text>')
    for pi, gene in enumerate(["BRCA1", "BRCA2"]):
        x0 = left + pi * (panel_w + panel_gap)
        body.append(f'<text x="{x0}" y="{top - 14}" class="label">{gene}</text>')
        body.append(f'<rect x="{x0}" y="{top}" width="{panel_w}" height="{panel_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
        for tick in [0, 20, 40, 60, 80, 100]:
            y = top + panel_h - tick / 100 * panel_h
            body.append(f'<line x1="{x0}" y1="{y}" x2="{x0 + panel_w}" y2="{y}" class="grid"/>')
            if pi == 0:
                body.append(f'<text x="42" y="{y + 4}" class="small">{tick}%</text>')
        slots = [(site, group) for site in SITE_TYPES for group in GROUP_ORDER]
        slot_w = (panel_w - 50) / len(slots)
        for si, (site_type, group) in enumerate(slots):
            match = next(
                row
                for row in summary_rows
                if row["gene"] == gene and row["site_type"] == site_type and row["group"] == group and row["distance_bin"] in {"0", "1-2"}
            )
            near_rows = [
                row
                for row in summary_rows
                if row["gene"] == gene
                and row["site_type"] == site_type
                and row["group"] == group
                and row["distance_bin"] in {"0", "1-2"}
            ]
            total = sum(int(row["total_variants"]) for row in near_rows)
            high = sum(int(row["spliceai_ge_0_20"]) for row in near_rows)
            value = high / total * 100 if total else 0
            x = x0 + 25 + si * slot_w
            bar_h = panel_h * value / 100
            y = top + panel_h - bar_h
            body.append(f'<rect x="{x}" y="{y}" width="{slot_w - bar_gap}" height="{bar_h}" fill="{GROUP_COLORS[group]}"/>')
            body.append(f'<text x="{x}" y="{top + panel_h + 18}" class="small">{site_type[:1]}-{group[:1]}</text>')
    legend_x = left
    legend_y = height - 80
    for i, group in enumerate(GROUP_ORDER):
        x = legend_x + i * 230
        body.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{GROUP_COLORS[group]}"/>')
        body.append(f'<text x="{x + 20}" y="{legend_y + 12}" class="small">{esc(GROUP_LABELS[group])}</text>')
    body.append(f'<text x="{left}" y="{height - 36}" class="small">x labels: a = acceptor-like, d = donor-like; b/v/p = benign/VUS/pathogenic</text>')
    write_svg(PLOT_DIR / "donor_acceptor_near_boundary_rate.svg", "\n".join(body), width, height)


def draw_distance_rate_heatmap(summary_rows: list[dict[str, object]]) -> None:
    width = 1120
    height = 620
    left = 210
    top = 70
    cell_w = 88
    cell_h = 27
    rows = []
    for gene in ["BRCA1", "BRCA2"]:
        for site_type in SITE_TYPES:
            for group in GROUP_ORDER:
                rows.append((gene, site_type, group))
    body = ['<text x="28" y="34" class="title">Boundary-distance SpliceAI rate heatmap</text>']
    for xi, bin_name in enumerate(DISTANCE_BINS):
        x = left + xi * cell_w
        body.append(f'<text x="{x + 4}" y="{top - 10}" class="small">{esc(bin_name)}</text>')
    for yi, (gene, site_type, group) in enumerate(rows):
        y = top + yi * cell_h
        label = f"{gene} {site_type} {group}"
        body.append(f'<text x="20" y="{y + 18}" class="small">{esc(label)}</text>')
        for xi, bin_name in enumerate(DISTANCE_BINS):
            row = next(
                item
                for item in summary_rows
                if item["gene"] == gene
                and item["site_type"] == site_type
                and item["group"] == group
                and item["distance_bin"] == bin_name
            )
            value = float(row["percent"])
            color = heat_color(value, 100.0, (124, 58, 237))
            x = left + xi * cell_w
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h - 2}" fill="{color}"/>')
            if value:
                body.append(f'<text x="{x + 5}" y="{y + 17}" class="small">{value:.1f}</text>')
    write_svg(PLOT_DIR / "boundary_distance_rate_heatmap.svg", "\n".join(body), width, height)


def draw_relative_position_plot(relative_rows: list[dict[str, object]]) -> None:
    width = 1120
    height = 500
    left = 82
    top = 70
    panel_w = 455
    panel_h = 290
    panel_gap = 70
    body = ['<text x="28" y="34" class="title">Exact relative position around splice boundaries</text>']
    body.append('<text x="28" y="55" class="small">Position 0 is the nearest internal CDS splice boundary; rates use SpliceAI >=0.20</text>')
    colors = {"acceptor_like": "#7c3aed", "donor_like": "#0f766e"}
    for pi, gene in enumerate(["BRCA1", "BRCA2"]):
        x0 = left + pi * (panel_w + panel_gap)
        body.append(f'<text x="{x0}" y="{top - 12}" class="label">{gene}</text>')
        body.append(f'<rect x="{x0}" y="{top}" width="{panel_w}" height="{panel_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
        for tick in [0, 20, 40, 60, 80, 100]:
            y = top + panel_h - tick / 100 * panel_h
            body.append(f'<line x1="{x0}" y1="{y}" x2="{x0 + panel_w}" y2="{y}" class="grid"/>')
            if pi == 0:
                body.append(f'<text x="38" y="{y + 4}" class="small">{tick}%</text>')
        for site_type in SITE_TYPES:
            points = []
            for rel in RELATIVE_POSITIONS:
                items = [
                    row
                    for row in relative_rows
                    if row["gene"] == gene
                    and row["site_type"] == site_type
                    and row["relative_boundary_position"] == rel
                ]
                total = sum(int(row["total_variants"]) for row in items)
                high = sum(int(row["spliceai_ge_0_20"]) for row in items)
                value = high / total * 100 if total else 0
                x = x0 + (rel - RELATIVE_POSITIONS[0]) / (RELATIVE_POSITIONS[-1] - RELATIVE_POSITIONS[0]) * panel_w
                y = top + panel_h - value / 100 * panel_h
                points.append(f"{x:.1f},{y:.1f}")
            body.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{colors[site_type]}" stroke-width="2"/>')
        x_mid = x0 + (-RELATIVE_POSITIONS[0]) / (RELATIVE_POSITIONS[-1] - RELATIVE_POSITIONS[0]) * panel_w
        body.append(f'<line x1="{x_mid}" y1="{top}" x2="{x_mid}" y2="{top + panel_h}" stroke="#111827" stroke-width="1" opacity="0.4"/>')
    legend_y = height - 78
    for i, site_type in enumerate(SITE_TYPES):
        x = left + i * 190
        body.append(f'<line x1="{x}" y1="{legend_y}" x2="{x + 26}" y2="{legend_y}" stroke="{colors[site_type]}" stroke-width="3"/>')
        body.append(f'<text x="{x + 34}" y="{legend_y + 4}" class="small">{site_type}</text>')
    body.append(f'<text x="{left + 390}" y="{height - 35}" class="label">Relative coding bases from nearest boundary</text>')
    write_svg(PLOT_DIR / "relative_position_spliceai_rate.svg", "\n".join(body), width, height)


def draw_threshold_sensitivity(threshold_rows: list[dict[str, object]]) -> None:
    width = 1100
    height = 450
    left = 85
    top = 70
    panel_w = 455
    panel_h = 260
    panel_gap = 72
    colors = {0.10: "#60a5fa", 0.20: "#7c3aed", 0.50: "#f59e0b", 0.80: "#dc2626"}
    body = ['<text x="28" y="34" class="title">SpliceAI threshold sensitivity near boundaries</text>']
    body.append('<text x="28" y="55" class="small">Near boundary means distance 0 or 1-2 bp, all groups combined</text>')
    for pi, gene in enumerate(["BRCA1", "BRCA2"]):
        x0 = left + pi * (panel_w + panel_gap)
        body.append(f'<text x="{x0}" y="{top - 12}" class="label">{gene}</text>')
        body.append(f'<rect x="{x0}" y="{top}" width="{panel_w}" height="{panel_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
        slots = SITE_TYPES
        slot_w = panel_w / len(slots)
        for si, site_type in enumerate(slots):
            x_slot = x0 + si * slot_w
            near_rows = [
                row
                for row in threshold_rows
                if row["gene"] == gene and row["site_type"] == site_type and row["distance_bin"] in {"0", "1-2"}
            ]
            for ti, threshold in enumerate(SPLICEAI_THRESHOLDS):
                key = f"spliceai_ge_{str(threshold).replace('.', '_')}"
                total = sum(int(row["total_variants"]) for row in near_rows)
                high = sum(int(row[f"{key}_count"]) for row in near_rows)
                value = high / total * 100 if total else 0
                bar_w = 28
                x = x_slot + 35 + ti * (bar_w + 8)
                h = panel_h * value / 100
                y = top + panel_h - h
                body.append(f'<rect x="{x}" y="{y}" width="{bar_w}" height="{h}" fill="{colors[threshold]}"/>')
            body.append(f'<text x="{x_slot + 40}" y="{top + panel_h + 20}" class="small">{site_type}</text>')
        for tick in [0, 20, 40, 60, 80, 100]:
            y = top + panel_h - tick / 100 * panel_h
            body.append(f'<line x1="{x0}" y1="{y}" x2="{x0 + panel_w}" y2="{y}" class="grid"/>')
            if pi == 0:
                body.append(f'<text x="36" y="{y + 4}" class="small">{tick}%</text>')
    legend_y = height - 70
    for i, threshold in enumerate(SPLICEAI_THRESHOLDS):
        x = left + i * 145
        body.append(f'<rect x="{x}" y="{legend_y - 10}" width="14" height="14" fill="{colors[threshold]}"/>')
        body.append(f'<text x="{x + 20}" y="{legend_y + 2}" class="small">>= {threshold:.2f}</text>')
    write_svg(PLOT_DIR / "threshold_sensitivity_near_boundary.svg", "\n".join(body), width, height)


def draw_exon_boundary_heatmap(exon_rows: list[dict[str, object]]) -> None:
    for gene in ["BRCA1", "BRCA2"]:
        rows = [row for row in exon_rows if row["gene"] == gene]
        exons = sorted({int(row["cds_exon"]) for row in rows})
        cell_w = 32
        cell_h = 30
        left = 125
        top = 75
        width = left + len(exons) * cell_w + 45
        height = top + len(SITE_TYPES) * cell_h + 85
        body = [f'<text x="28" y="34" class="title">{gene} exon boundary high SpliceAI rate</text>']
        body.append('<text x="28" y="55" class="small">Near-boundary variants within 10 bp, rate uses SpliceAI >=0.20</text>')
        for xi, exon in enumerate(exons):
            x = left + xi * cell_w
            body.append(f'<text x="{x + 4}" y="{top - 10}" class="small">E{exon}</text>')
        for yi, site_type in enumerate(SITE_TYPES):
            y = top + yi * cell_h
            body.append(f'<text x="22" y="{y + 20}" class="small">{site_type}</text>')
            for xi, exon in enumerate(exons):
                match = next((row for row in rows if int(row["cds_exon"]) == exon and row["site_type"] == site_type), None)
                value = float(match["percent"]) if match else 0
                color = heat_color(value, 100, (124, 58, 237))
                x = left + xi * cell_w
                body.append(f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h - 2}" fill="{color}"/>')
        write_svg(PLOT_DIR / f"{gene.lower()}_exon_boundary_heatmap.svg", "\n".join(body), width, height)


def draw_variant_type_heatmap(type_rows: list[dict[str, object]]) -> None:
    width = 1120
    height = 620
    left = 210
    top = 70
    cell_w = 88
    cell_h = 27
    rows = [(gene, site_type, vtype) for gene in ["BRCA1", "BRCA2"] for site_type in SITE_TYPES for vtype in VARIANT_TYPES]
    body = ['<text x="28" y="34" class="title">Boundary SpliceAI rate by variant type</text>']
    for xi, bin_name in enumerate(DISTANCE_BINS):
        x = left + xi * cell_w
        body.append(f'<text x="{x + 4}" y="{top - 10}" class="small">{esc(bin_name)}</text>')
    for yi, (gene, site_type, vtype) in enumerate(rows):
        y = top + yi * cell_h
        body.append(f'<text x="20" y="{y + 18}" class="small">{gene} {site_type} {vtype}</text>')
        for xi, bin_name in enumerate(DISTANCE_BINS):
            row = next(
                item
                for item in type_rows
                if item["gene"] == gene
                and item["site_type"] == site_type
                and item["variant_type"] == vtype
                and item["distance_bin"] == bin_name
            )
            value = float(row["percent"])
            color = heat_color(value, 100.0, (14, 116, 144))
            x = left + xi * cell_w
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h - 2}" fill="{color}"/>')
            if value:
                body.append(f'<text x="{x + 5}" y="{y + 17}" class="small">{value:.1f}</text>')
    write_svg(PLOT_DIR / "variant_type_boundary_rate_heatmap.svg", "\n".join(body), width, height)


def report_table(rows: list[dict[str, object]], columns: list[str], limit: int = 12) -> str:
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows[:limit]:
        body.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join([header, sep] + body)


def write_report(
    site_rows: list[dict[str, object]],
    distance_rows: list[dict[str, object]],
    relative_rows: list[dict[str, object]],
    conflicts: list[dict[str, object]],
    strongest: list[dict[str, object]],
) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    near_rows = []
    for row in distance_rows:
        if row["distance_bin"] not in {"0", "1-2"}:
            continue
        near_rows.append(row)
    top_near = sorted(near_rows, key=lambda row: float(row["percent"]), reverse=True)
    top_relative = sorted(
        [row for row in relative_rows if int(row["total_variants"]) >= 3],
        key=lambda row: float(row["percent"]),
        reverse=True,
    )
    text = f"""# Boundary Enrichment Analysis

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

This analysis refines the earlier boundary-distance view. Internal CDS exon starts are treated as acceptor-like boundaries and internal CDS exon ends are treated as donor-like boundaries in transcript coordinate space. The first CDS start and last CDS end are not treated as splice junctions.

Because the snapshot contains coding SNVs only, this is an exonic-side approximation. It does not include intronic variants.

## Outputs

- `tables/boundary_enrichment/variant_boundary_annotations.csv`
- `tables/boundary_enrichment/boundary_distance_group_summary.csv`
- `tables/boundary_enrichment/donor_acceptor_site_summary.csv`
- `tables/boundary_enrichment/relative_position_summary.csv`
- `tables/boundary_enrichment/threshold_sensitivity_summary.csv`
- `tables/boundary_enrichment/exon_boundary_summary.csv`
- `tables/boundary_enrichment/variant_type_boundary_summary.csv`
- `tables/boundary_enrichment/boundary_conflict_sanity_checks.csv`
- `tables/boundary_enrichment/strongest_boundary_vus_pathogenic_variants.csv`
- `plots/09_boundary_enrichment/donor_acceptor_near_boundary_rate.svg`
- `plots/09_boundary_enrichment/boundary_distance_rate_heatmap.svg`
- `plots/09_boundary_enrichment/relative_position_spliceai_rate.svg`
- `plots/09_boundary_enrichment/threshold_sensitivity_near_boundary.svg`
- `plots/09_boundary_enrichment/brca1_exon_boundary_heatmap.svg`
- `plots/09_boundary_enrichment/brca2_exon_boundary_heatmap.svg`
- `plots/09_boundary_enrichment/variant_type_boundary_rate_heatmap.svg`

## Donor-like versus Acceptor-like Summary

{report_table(site_rows, ["gene", "site_type", "total_variants", "spliceai_ge_0_20", "percent", "benign_percent", "vus_percent", "pathogenic_percent"])}

## Highest Near-boundary Rates

Near-boundary here means distance bin 0 or 1-2 bp.

{report_table(top_near, ["gene", "site_type", "group", "distance_bin", "total_variants", "spliceai_ge_0_20", "percent"], limit=18)}

## Highest Exact Relative-position Rates

Only rows with at least 3 variants are shown here.

{report_table(top_relative, ["gene", "site_type", "group", "relative_boundary_position", "total_variants", "spliceai_ge_0_20", "percent"], limit=18)}

## Boundary Conflict Sanity Checks

Conflict rows: {len(conflicts)}

{report_table(conflicts, ["boundary_conflict_reason", "gene", "c_notation", "p_notation", "group", "spliceai_score", "nearest_site_type", "relative_boundary_position"], limit=18)}

## Strongest Boundary VUS/Pathogenic Variants

Variants are VUS or pathogenic-group, within 2 bp of an internal CDS splice boundary, and SpliceAI >=0.20.

{report_table(strongest, ["gene", "c_notation", "p_notation", "group", "spliceai_score", "nearest_site_type", "relative_boundary_position", "criteria"], limit=18)}

## Reading Notes

If donor-like and acceptor-like rates are similar, the main signal is general proximity to a splice boundary. If one side is consistently higher, that suggests side-specific splice disruption patterns worth reviewing against RNA or curated rules.
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    boundaries = build_splice_boundaries()
    rows = load_rows(boundaries)
    distance_rows = summarize(rows)
    site_rows = site_summary(rows)
    threshold_rows = threshold_sensitivity(rows)
    relative_rows = relative_position_summary(rows)
    exon_rows = exon_boundary_summary(rows)
    type_rows = variant_type_boundary_summary(rows)
    conflicts = boundary_conflicts(rows)
    strongest = strongest_boundary_variants(rows)
    write_csv(
        TABLE_DIR / "variant_boundary_annotations.csv",
        rows,
        [
            "gene",
            "c_notation",
            "p_notation",
            "cds_pos",
            "variant_type",
            "spliceai_score",
            "group",
            "predicted_class",
            "nearest_site_type",
            "nearest_boundary_cds_pos",
            "nearest_boundary_exon",
            "relative_boundary_position",
            "splice_boundary_distance",
            "splice_boundary_distance_bin",
            "criteria",
        ],
    )
    write_csv(
        TABLE_DIR / "boundary_distance_group_summary.csv",
        distance_rows,
        ["gene", "site_type", "group", "group_label", "distance_bin", "total_variants", "spliceai_ge_0_20", "percent"],
    )
    write_csv(
        TABLE_DIR / "donor_acceptor_site_summary.csv",
        site_rows,
        ["gene", "site_type", "total_variants", "spliceai_ge_0_20", "percent"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["total", "high_spliceai", "percent"]],
    )
    threshold_fields = ["gene", "site_type", "group", "distance_bin", "total_variants"]
    for threshold in SPLICEAI_THRESHOLDS:
        key = f"spliceai_ge_{str(threshold).replace('.', '_')}"
        threshold_fields.extend([f"{key}_count", f"{key}_percent"])
    write_csv(TABLE_DIR / "threshold_sensitivity_summary.csv", threshold_rows, threshold_fields)
    write_csv(
        TABLE_DIR / "relative_position_summary.csv",
        relative_rows,
        ["gene", "site_type", "group", "relative_boundary_position", "total_variants", "spliceai_ge_0_20", "percent"],
    )
    write_csv(
        TABLE_DIR / "exon_boundary_summary.csv",
        exon_rows,
        ["gene", "cds_exon", "site_type", "near_10bp_total", "spliceai_ge_0_20", "percent"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["total", "high_spliceai", "percent"]],
    )
    write_csv(
        TABLE_DIR / "variant_type_boundary_summary.csv",
        type_rows,
        ["gene", "site_type", "variant_type", "distance_bin", "total_variants", "spliceai_ge_0_20", "percent"],
    )
    write_csv(
        TABLE_DIR / "boundary_conflict_sanity_checks.csv",
        conflicts,
        [
            "boundary_conflict_reason",
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "group",
            "predicted_class",
            "spliceai_score",
            "nearest_site_type",
            "nearest_boundary_exon",
            "relative_boundary_position",
            "splice_boundary_distance",
            "criteria",
        ],
    )
    write_csv(
        TABLE_DIR / "strongest_boundary_vus_pathogenic_variants.csv",
        strongest,
        [
            "gene",
            "c_notation",
            "p_notation",
            "variant_type",
            "group",
            "predicted_class",
            "spliceai_score",
            "nearest_site_type",
            "nearest_boundary_exon",
            "relative_boundary_position",
            "splice_boundary_distance",
            "criteria",
        ],
    )
    draw_site_type_rate_plot(distance_rows)
    draw_distance_rate_heatmap(distance_rows)
    draw_relative_position_plot(relative_rows)
    draw_threshold_sensitivity(threshold_rows)
    draw_exon_boundary_heatmap(exon_rows)
    draw_variant_type_heatmap(type_rows)
    write_report(site_rows, distance_rows, relative_rows, conflicts, strongest)
    print(f"Wrote boundary enrichment analysis to {REPORT}")


if __name__ == "__main__":
    main()
