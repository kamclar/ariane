"""Generate dependency-free SVG plots for the precomputed BRCA SNV snapshot."""

from __future__ import annotations

import csv
import html
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
TRANSCRIPT_MAPS = ROOT / "variant_space_scan" / "coordinate_mapping" / "data" / "transcript_maps.json"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
PLOT_DIR = OUT_DIR / "plots"
OVERVIEW_PLOT_DIR = PLOT_DIR / "01_overview"
POSITION_PLOT_DIR = PLOT_DIR / "02_position"
BOUNDARY_PLOT_DIR = PLOT_DIR / "03_boundary_spliceai"
CLUSTER_PLOT_DIR = PLOT_DIR / "04_clusters"
TABLE_DIR = OUT_DIR / "tables"
REPORT = OUT_DIR / "visualization_report.md"

CLASS_LABELS = {
    "1": "Benign",
    "2": "Likely Benign",
    "3": "VUS",
    "4": "Likely Pathogenic",
    "5": "Pathogenic",
}
CLASS_COLORS = {
    "1": "#2b8cbe",
    "2": "#7bccc4",
    "3": "#8a8f9c",
    "4": "#fdae61",
    "5": "#d7301f",
}
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
GROUP_ORDER = ["benign", "vus", "pathogenic"]
CLUSTER_COLORS = ["#2563eb", "#0891b2", "#65a30d", "#f59e0b", "#dc2626"]
CRITERIA_ORDER = [
    "BA1",
    "BS1_Strong",
    "BS1_Supporting",
    "BS3",
    "BP1",
    "BP4",
    "BP5",
    "BP7",
    "PM2_Supporting",
    "PM5_PTC",
    "PP3",
    "PP4",
    "PS1",
    "PS3",
    "PVS1",
]
CDS_LENGTHS = {"BRCA1": 5592, "BRCA2": 10257}
C_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def parse_float(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_cds_pos(c_notation: str) -> int | None:
    match = C_RE.match(c_notation or "")
    if not match:
        return None
    return int(match.group(1))


def parse_criteria(criteria: str) -> set[str]:
    if not criteria:
        return set()
    codes = set()
    for item in criteria.split(";"):
        code = item.split(":", 1)[0]
        if code:
            codes.add(code)
    return codes


def splice_bin(score: float | None) -> str:
    if score is None:
        return "missing"
    if score == 0:
        return "0"
    if score < 0.10:
        return "0-0.099"
    if score < 0.20:
        return "0.10-0.199"
    if score < 0.50:
        return "0.20-0.499"
    if score < 0.80:
        return "0.50-0.799"
    return "0.80-1.00"


def boundary_bin(distance: int | None) -> str:
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


def boundary_bin_compact(distance: int | None) -> str:
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
    if distance <= 250:
        return "101-250"
    return ">250"


def splice_bin_compact(score: float | None) -> str:
    if score is None:
        return "missing"
    if score < 0.05:
        return "0-0.049"
    if score < 0.10:
        return "0.05-0.099"
    if score < 0.20:
        return "0.10-0.199"
    if score < 0.50:
        return "0.20-0.499"
    if score < 0.80:
        return "0.50-0.799"
    return "0.80-1.00"


def class_group(predicted_class: str | object) -> str:
    cls = str(predicted_class)
    if cls in {"1", "2"}:
        return "benign"
    if cls in {"4", "5"}:
        return "pathogenic"
    return "vus"


def heat_color(value: float, max_value: float, base: tuple[int, int, int] = (37, 99, 235)) -> str:
    if max_value <= 0:
        intensity = 0.0
    else:
        intensity = math.sqrt(max(value, 0) / max_value)
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
.axis {{ stroke: #9aa4b2; stroke-width: 1; }}
.grid {{ stroke: #e2e8f0; stroke-width: 1; }}
.panel {{ fill: #f3f4f6; }}
.tile {{ stroke: #ffffff; stroke-width: 1; }}
.legend-title {{ font-size: 11px; font-weight: 700; }}
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def load_rows() -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["cds_pos"] = parse_cds_pos(row["c_notation"])
            row["spliceai"] = parse_float(row["spliceai_score"])
            row["criteria_codes"] = parse_criteria(row.get("criteria", ""))
            rows.append(row)
    return rows


def cds_boundaries_from_maps() -> dict[str, list[int]]:
    data = json.loads(TRANSCRIPT_MAPS.read_text(encoding="utf-8"))
    boundaries = {}
    for gene, info in data["assemblies"]["GRCh38"].items():
        cds_features = list(info["cds"])
        strand = info["strand"]
        if strand == "+":
            ordered = sorted(cds_features, key=lambda item: item["start"])
        else:
            ordered = sorted(cds_features, key=lambda item: item["start"], reverse=True)
        pos = 1
        gene_boundaries = []
        for feature in ordered:
            length = feature["end"] - feature["start"] + 1
            start = pos
            end = pos + length - 1
            gene_boundaries.extend([start, end])
            pos = end + 1
        boundaries[gene] = sorted(set(x for x in gene_boundaries if 1 <= x <= CDS_LENGTHS[gene]))
    return boundaries


def add_boundary_distances(rows: list[dict[str, object]], boundaries: dict[str, list[int]]) -> None:
    for row in rows:
        gene = str(row["gene"])
        pos = row.get("cds_pos")
        if not isinstance(pos, int):
            row["boundary_distance"] = None
            continue
        row["boundary_distance"] = min(abs(pos - boundary) for boundary in boundaries[gene])


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def bar_class_distribution(rows: list[dict[str, object]]) -> None:
    counts = Counter(str(row["predicted_class"]) for row in rows)
    max_count = max(counts.values())
    width = 760
    height = 380
    left = 150
    top = 60
    bar_h = 42
    gap = 18
    plot_w = 520
    body = [f'<text x="30" y="32" class="title">Predicted class distribution</text>']
    for i, cls in enumerate(sorted(CLASS_LABELS, key=int)):
        y = top + i * (bar_h + gap)
        count = counts[cls]
        w = round(plot_w * count / max_count)
        body.append(f'<text x="30" y="{y + 27}" class="label">{cls} {CLASS_LABELS[cls]}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="{CLASS_COLORS[cls]}"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 27}" class="label">{count}</text>')
    write_svg(OVERVIEW_PLOT_DIR / "class_distribution.svg", "\n".join(body), width, height)


def bar_group_distribution(rows: list[dict[str, object]]) -> None:
    counts = Counter(class_group(row["predicted_class"]) for row in rows)
    max_count = max(counts.values())
    width = 760
    height = 280
    left = 210
    top = 60
    bar_h = 42
    gap = 18
    plot_w = 470
    body = [f'<text x="30" y="32" class="title">Grouped predicted class distribution</text>']
    for i, group in enumerate(GROUP_ORDER):
        y = top + i * (bar_h + gap)
        count = counts[group]
        w = round(plot_w * count / max_count)
        body.append(f'<text x="30" y="{y + 27}" class="label">{GROUP_LABELS[group]}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="{GROUP_COLORS[group]}"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 27}" class="label">{count}</text>')
    write_svg(OVERVIEW_PLOT_DIR / "grouped_class_distribution.svg", "\n".join(body), width, height)


def heatmap(
    *,
    path: Path,
    title: str,
    x_labels: list[str],
    y_labels: list[str],
    values: dict[tuple[str, str], int],
    x_title: str = "",
) -> None:
    cell_w = 92
    cell_h = 30
    left = 150
    top = 96
    width = left + cell_w * len(x_labels) + 128
    height = top + cell_h * len(y_labels) + 70
    max_value = max(values.values()) if values else 0
    body = [f'<text x="30" y="32" class="title">{esc(title)}</text>']
    if x_title:
        body.append(f'<text x="{left}" y="58" class="small">{esc(x_title)}</text>')
    panel_w = cell_w * len(x_labels)
    panel_h = cell_h * len(y_labels)
    body.append(f'<rect x="{left - 2}" y="{top - 2}" width="{panel_w + 2}" height="{panel_h + 2}" class="panel"/>')
    for xi, x_label in enumerate(x_labels):
        x = left + xi * cell_w
        body.append(f'<text x="{x + 4}" y="{top - 16}" class="small">{esc(x_label)}</text>')
    for yi, y_label in enumerate(y_labels):
        y = top + yi * cell_h
        body.append(f'<text x="20" y="{y + 20}" class="small">{esc(y_label)}</text>')
        for xi, x_label in enumerate(x_labels):
            x = left + xi * cell_w
            value = values.get((y_label, x_label), 0)
            color = heat_color(value, max_value)
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{color}" class="tile"/>')
            if value:
                body.append(f'<text x="{x + 6}" y="{y + 19}" class="small">{value}</text>')
    legend_x = left + panel_w + 32
    legend_y = top + 4
    body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">count</text>')
    legend_steps = 6
    for step in range(legend_steps):
        frac = 1 - step / (legend_steps - 1)
        color = heat_color(frac * max_value, max_value)
        y = legend_y + step * 18
        body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{color}" class="tile"/>')
    body.append(f'<text x="{legend_x + 26}" y="{legend_y + 12}" class="small">{max_value:.0f}</text>')
    body.append(f'<text x="{legend_x + 26}" y="{legend_y + (legend_steps - 1) * 18 + 12}" class="small">0</text>')
    write_svg(path, "\n".join(body), width, height)


def percent_heatmap(
    *,
    path: Path,
    title: str,
    x_labels: list[str],
    y_labels: list[str],
    values: dict[tuple[str, str], float],
    x_title: str = "",
    note: str = "",
) -> None:
    cell_w = 112
    cell_h = 34
    left = 190
    title_words = title.split()
    title_lines: list[str] = []
    current = ""
    for word in title_words:
        candidate = f"{current} {word}".strip()
        if len(candidate) > 48 and current:
            title_lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        title_lines.append(current)
    title_lines = title_lines[:2]

    title_y = 28
    note_y = title_y + 22 * len(title_lines) + 2
    x_title_y = note_y + (22 if note else 0)
    top = x_title_y + 36
    width = left + cell_w * len(x_labels) + 130
    height = top + cell_h * len(y_labels) + 78
    max_value = max(values.values()) if values else 0.0
    body = []
    for idx, line in enumerate(title_lines):
        body.append(f'<text x="30" y="{title_y + idx * 22}" class="title">{esc(line)}</text>')
    if note:
        body.append(f'<text x="30" y="{note_y}" class="small">{esc(note)}</text>')
    if x_title:
        body.append(f'<text x="{left}" y="{x_title_y}" class="small">{esc(x_title)}</text>')
    panel_w = cell_w * len(x_labels)
    panel_h = cell_h * len(y_labels)
    body.append(f'<rect x="{left - 2}" y="{top - 2}" width="{panel_w + 2}" height="{panel_h + 2}" class="panel"/>')
    for xi, x_label in enumerate(x_labels):
        x = left + xi * cell_w
        body.append(f'<text x="{x + 6}" y="{top - 16}" class="small">{esc(x_label)}</text>')
    for yi, y_label in enumerate(y_labels):
        y = top + yi * cell_h
        body.append(f'<text x="20" y="{y + 22}" class="small">{esc(y_label)}</text>')
        for xi, x_label in enumerate(x_labels):
            x = left + xi * cell_w
            value = values.get((y_label, x_label), 0.0)
            color = heat_color(value, max_value)
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{color}" class="tile"/>')
            if value:
                body.append(f'<text x="{x + 6}" y="{y + 21}" class="small">{value:.2f}%</text>')
    legend_x = left + panel_w + 34
    legend_y = top + 4
    body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">percent</text>')
    legend_steps = 6
    for step in range(legend_steps):
        frac = 1 - step / (legend_steps - 1)
        color = heat_color(frac * max_value, max_value)
        y = legend_y + step * 18
        body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{color}" class="tile"/>')
    body.append(f'<text x="{legend_x + 26}" y="{legend_y + 12}" class="small">{max_value:.1f}%</text>')
    body.append(f'<text x="{legend_x + 26}" y="{legend_y + (legend_steps - 1) * 18 + 12}" class="small">0</text>')
    write_svg(path, "\n".join(body), width, height)


def plot_criteria_heatmap(rows: list[dict[str, object]]) -> None:
    values = Counter()
    for row in rows:
        cls = str(row["predicted_class"])
        for code in row["criteria_codes"]:
            values[(code, cls)] += 1
    heatmap(
        path=OVERVIEW_PLOT_DIR / "criteria_by_class_heatmap.svg",
        title="Criteria by predicted class",
        x_labels=sorted(CLASS_LABELS, key=int),
        y_labels=[code for code in CRITERIA_ORDER if any((code, cls) in values for cls in CLASS_LABELS)],
        values=values,
        x_title="Class",
    )


def plot_splice_bins_heatmap(rows: list[dict[str, object]]) -> None:
    bins = ["0", "0-0.099", "0.10-0.199", "0.20-0.499", "0.50-0.799", "0.80-1.00", "missing"]
    values = Counter((splice_bin(row["spliceai"]), str(row["predicted_class"])) for row in rows)
    heatmap(
        path=OVERVIEW_PLOT_DIR / "spliceai_bins_by_class_heatmap.svg",
        title="SpliceAI score bins by predicted class",
        x_labels=sorted(CLASS_LABELS, key=int),
        y_labels=bins,
        values=values,
        x_title="Class",
    )


def plot_grouped_splice_bins_heatmap(rows: list[dict[str, object]]) -> None:
    bins = ["0", "0-0.099", "0.10-0.199", "0.20-0.499", "0.50-0.799", "0.80-1.00", "missing"]
    values = Counter((splice_bin(row["spliceai"]), class_group(row["predicted_class"])) for row in rows)
    heatmap(
        path=OVERVIEW_PLOT_DIR / "spliceai_bins_by_group_heatmap.svg",
        title="SpliceAI score bins by grouped class",
        x_labels=GROUP_ORDER,
        y_labels=bins,
        values=values,
        x_title="Grouped class",
    )


def plot_gene_normalized_overview(rows: list[dict[str, object]]) -> None:
    genes = ["BRCA1", "BRCA2"]
    totals = Counter(str(row["gene"]) for row in rows)

    class_counts = Counter((str(row["predicted_class"]), str(row["gene"])) for row in rows)
    class_rates = {
        (cls, gene): (class_counts[(cls, gene)] / totals[gene] * 100 if totals[gene] else 0.0)
        for cls in sorted(CLASS_LABELS, key=int)
        for gene in genes
    }
    percent_heatmap(
        path=OVERVIEW_PLOT_DIR / "gene_normalized_class_distribution_heatmap.svg",
        title="Classes by gene",
        x_labels=genes,
        y_labels=[f"{cls} {CLASS_LABELS[cls]}" for cls in sorted(CLASS_LABELS, key=int)],
        values={(f"{cls} {CLASS_LABELS[cls]}", gene): class_rates[(cls, gene)] for cls in sorted(CLASS_LABELS, key=int) for gene in genes},
        x_title="Gene",
        note="Cells show percent of generated coding SNVs within each gene.",
    )

    group_counts = Counter((class_group(row["predicted_class"]), str(row["gene"])) for row in rows)
    group_rates = {
        (GROUP_LABELS[group], gene): (group_counts[(group, gene)] / totals[gene] * 100 if totals[gene] else 0.0)
        for group in GROUP_ORDER
        for gene in genes
    }
    percent_heatmap(
        path=OVERVIEW_PLOT_DIR / "gene_normalized_grouped_class_distribution_heatmap.svg",
        title="Grouped classes by gene",
        x_labels=genes,
        y_labels=[GROUP_LABELS[group] for group in GROUP_ORDER],
        values=group_rates,
        x_title="Gene",
        note="Benign combines class 1/2; pathogenic combines class 4/5.",
    )

    variant_types = sorted({str(row["variant_type"]) for row in rows})
    type_counts = Counter((str(row["variant_type"]), str(row["gene"])) for row in rows)
    type_rates = {
        (vtype, gene): (type_counts[(vtype, gene)] / totals[gene] * 100 if totals[gene] else 0.0)
        for vtype in variant_types
        for gene in genes
    }
    percent_heatmap(
        path=OVERVIEW_PLOT_DIR / "gene_normalized_variant_type_heatmap.svg",
        title="Variant types by gene",
        x_labels=genes,
        y_labels=variant_types,
        values=type_rates,
        x_title="Gene",
        note="Cells show percent of generated coding SNVs within each gene.",
    )

    criteria_counts = Counter()
    for row in rows:
        gene = str(row["gene"])
        for code in row["criteria_codes"]:
            criteria_counts[(code, gene)] += 1
    criteria_rates = {
        (code, gene): (criteria_counts[(code, gene)] / totals[gene] * 100 if totals[gene] else 0.0)
        for code in CRITERIA_ORDER
        for gene in genes
    }
    percent_heatmap(
        path=OVERVIEW_PLOT_DIR / "gene_normalized_criteria_rate_heatmap.svg",
        title="Criteria rates by gene",
        x_labels=genes,
        y_labels=[code for code in CRITERIA_ORDER if any(criteria_counts[(code, gene)] for gene in genes)],
        values=criteria_rates,
        x_title="Gene",
        note="Cells show percent of generated coding SNVs with each criterion.",
    )

    bins = ["0", "0-0.099", "0.10-0.199", "0.20-0.499", "0.50-0.799", "0.80-1.00", "missing"]
    splice_counts = Counter((splice_bin(row["spliceai"]), str(row["gene"])) for row in rows)
    splice_rates = {
        (bin_name, gene): (splice_counts[(bin_name, gene)] / totals[gene] * 100 if totals[gene] else 0.0)
        for bin_name in bins
        for gene in genes
    }
    percent_heatmap(
        path=OVERVIEW_PLOT_DIR / "gene_normalized_spliceai_bins_heatmap.svg",
        title="SpliceAI bins by gene",
        x_labels=genes,
        y_labels=bins,
        values=splice_rates,
        x_title="Gene",
        note="Cells show percent of generated coding SNVs within each gene.",
    )


def plot_boundary_heatmap(rows: list[dict[str, object]], *, gene: str | None = None, filename_prefix: str = "") -> None:
    if gene is not None:
        rows = [row for row in rows if row.get("gene") == gene]
    bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100", "missing"]
    values = Counter((boundary_bin(row["boundary_distance"]), str(row["predicted_class"])) for row in rows)
    heatmap(
        path=BOUNDARY_PLOT_DIR / f"{filename_prefix}boundary_distance_by_class_heatmap.svg",
        title=f"{gene + ': ' if gene else ''}Boundary distance by class",
        x_labels=sorted(CLASS_LABELS, key=int),
        y_labels=bins,
        values=values,
        x_title="Class",
    )
    high = Counter()
    total = Counter()
    for row in rows:
        b = boundary_bin(row["boundary_distance"])
        total[b] += 1
        if row["spliceai"] is not None and row["spliceai"] >= 0.20:
            high[b] += 1
    table_name = "spliceai_high_by_boundary_distance.csv" if not filename_prefix else f"{filename_prefix.rstrip('_')}_spliceai_high_by_boundary_distance.csv"
    write_csv(
        TABLE_DIR / table_name,
        ["boundary_distance_bin", "total_variants", "spliceai_ge_0_20", "percent"],
        [[b, total[b], high[b], f"{(high[b] / total[b] * 100):.2f}" if total[b] else "0.00"] for b in bins],
    )


def plot_grouped_boundary_heatmap(rows: list[dict[str, object]], *, gene: str | None = None, filename_prefix: str = "") -> None:
    if gene is not None:
        rows = [row for row in rows if row.get("gene") == gene]
    bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100", "missing"]
    values = Counter((boundary_bin(row["boundary_distance"]), class_group(row["predicted_class"])) for row in rows)
    heatmap(
        path=BOUNDARY_PLOT_DIR / f"{filename_prefix}boundary_distance_by_group_heatmap.svg",
        title=f"{gene + ': ' if gene else ''}Boundary distance by group",
        x_labels=GROUP_ORDER,
        y_labels=bins,
        values=values,
        x_title="Grouped class",
    )


def plot_position_tracks(rows: list[dict[str, object]], boundaries: dict[str, list[int]]) -> None:
    for gene in ["BRCA1", "BRCA2"]:
        gene_rows = [row for row in rows if row["gene"] == gene and isinstance(row.get("cds_pos"), int)]
        length = CDS_LENGTHS[gene]
        bins = 120
        class_counts = {cls: [0] * bins for cls in CLASS_LABELS}
        splice_max = [0.0] * bins
        splice_high = [0] * bins
        for row in gene_rows:
            pos = int(row["cds_pos"])
            idx = min(bins - 1, max(0, int((pos - 1) / length * bins)))
            cls = str(row["predicted_class"])
            class_counts[cls][idx] += 1
            score = row["spliceai"]
            if score is not None:
                splice_max[idx] = max(splice_max[idx], float(score))
                if score >= 0.20:
                    splice_high[idx] += 1
        write_position_csv(gene, class_counts, splice_max, splice_high, bins, length)
        draw_position_svg(gene, class_counts, splice_max, splice_high, boundaries[gene], bins, length)


def plot_grouped_position_tracks(rows: list[dict[str, object]], boundaries: dict[str, list[int]]) -> None:
    for gene in ["BRCA1", "BRCA2"]:
        gene_rows = [row for row in rows if row["gene"] == gene and isinstance(row.get("cds_pos"), int)]
        length = CDS_LENGTHS[gene]
        bins = 120
        group_counts = {group: [0] * bins for group in GROUP_ORDER}
        splice_max = [0.0] * bins
        splice_high = [0] * bins
        for row in gene_rows:
            pos = int(row["cds_pos"])
            idx = min(bins - 1, max(0, int((pos - 1) / length * bins)))
            group_counts[class_group(row["predicted_class"])][idx] += 1
            score = row["spliceai"]
            if score is not None:
                splice_max[idx] = max(splice_max[idx], float(score))
                if score >= 0.20:
                    splice_high[idx] += 1
        write_grouped_position_csv(gene, group_counts, splice_max, splice_high, bins, length)
        draw_grouped_position_svg(gene, group_counts, splice_max, splice_high, boundaries[gene], bins, length)


def write_position_csv(
    gene: str,
    class_counts: dict[str, list[int]],
    splice_max: list[float],
    splice_high: list[int],
    bins: int,
    length: int,
) -> None:
    rows = []
    for idx in range(bins):
        start = int(idx * length / bins) + 1
        end = int((idx + 1) * length / bins)
        row = [idx + 1, start, end, f"{splice_max[idx]:.3f}", splice_high[idx]]
        row.extend(class_counts[cls][idx] for cls in sorted(CLASS_LABELS, key=int))
        rows.append(row)
    write_csv(
        TABLE_DIR / f"{gene.lower()}_position_bins.csv",
        ["bin", "cds_start", "cds_end", "max_spliceai", "spliceai_ge_0_20"]
        + [f"class_{cls}" for cls in sorted(CLASS_LABELS, key=int)],
        rows,
    )


def draw_position_svg(
    gene: str,
    class_counts: dict[str, list[int]],
    splice_max: list[float],
    splice_high: list[int],
    boundaries: list[int],
    bins: int,
    length: int,
) -> None:
    cell_w = 7
    row_h = 26
    left = 75
    top = 80
    plot_w = cell_w * bins
    width = left + plot_w + 40
    height = top + row_h * 7 + 95
    body = [f'<text x="28" y="34" class="title">{gene} CDS position overview</text>']
    body.append(f'<text x="{left}" y="58" class="small">CDS positions 1 to {length}; vertical ticks mark CDS exon boundaries</text>')
    labels = [f"class {cls}" for cls in sorted(CLASS_LABELS, key=int)] + ["max SpliceAI", "SpAI >=0.20"]
    data_rows = [class_counts[cls] for cls in sorted(CLASS_LABELS, key=int)] + [
        [round(value * 100) for value in splice_max],
        splice_high,
    ]
    max_values = [max(row) if row else 0 for row in data_rows]
    bases = [(43, 140, 190), (123, 204, 196), (138, 143, 156), (253, 174, 97), (215, 48, 31), (88, 80, 141), (49, 130, 189)]
    for yi, label in enumerate(labels):
        y = top + yi * row_h
        body.append(f'<text x="10" y="{y + 17}" class="small">{esc(label)}</text>')
        for xi, value in enumerate(data_rows[yi]):
            x = left + xi * cell_w
            color = heat_color(value, max_values[yi], bases[yi])
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h - 3}" fill="{color}"/>')
    for boundary in boundaries:
        x = left + int((boundary - 1) / length * plot_w)
        body.append(f'<line x1="{x}" y1="{top - 6}" x2="{x}" y2="{top + row_h * len(labels)}" stroke="#1f2937" stroke-width="0.5" opacity="0.35"/>')
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        x = left + frac * plot_w
        pos = round(1 + frac * (length - 1))
        body.append(f'<line x1="{x}" y1="{top + row_h * len(labels) + 4}" x2="{x}" y2="{top + row_h * len(labels) + 10}" class="axis"/>')
        body.append(f'<text x="{x - 16}" y="{top + row_h * len(labels) + 26}" class="small">{pos}</text>')
    write_svg(POSITION_PLOT_DIR / f"{gene.lower()}_cds_position_overview.svg", "\n".join(body), width, height)


def write_grouped_position_csv(
    gene: str,
    group_counts: dict[str, list[int]],
    splice_max: list[float],
    splice_high: list[int],
    bins: int,
    length: int,
) -> None:
    rows = []
    for idx in range(bins):
        start = int(idx * length / bins) + 1
        end = int((idx + 1) * length / bins)
        row = [idx + 1, start, end, f"{splice_max[idx]:.3f}", splice_high[idx]]
        row.extend(group_counts[group][idx] for group in GROUP_ORDER)
        rows.append(row)
    write_csv(
        TABLE_DIR / f"{gene.lower()}_grouped_position_bins.csv",
        ["bin", "cds_start", "cds_end", "max_spliceai", "spliceai_ge_0_20"]
        + [f"group_{group}" for group in GROUP_ORDER],
        rows,
    )


def draw_grouped_position_svg(
    gene: str,
    group_counts: dict[str, list[int]],
    splice_max: list[float],
    splice_high: list[int],
    boundaries: list[int],
    bins: int,
    length: int,
) -> None:
    cell_w = 7
    row_h = 30
    left = 110
    top = 80
    plot_w = cell_w * bins
    width = left + plot_w + 40
    height = top + row_h * 5 + 95
    body = [f'<text x="28" y="34" class="title">{gene} grouped CDS position overview</text>']
    body.append(f'<text x="{left}" y="58" class="small">Benign combines classes 1-2; pathogenic combines classes 4-5</text>')
    labels = [GROUP_LABELS[group] for group in GROUP_ORDER] + ["max SpliceAI", "SpAI >=0.20"]
    data_rows = [group_counts[group] for group in GROUP_ORDER] + [
        [round(value * 100) for value in splice_max],
        splice_high,
    ]
    max_values = [max(row) if row else 0 for row in data_rows]
    bases = [(43, 140, 190), (138, 143, 156), (215, 48, 31), (88, 80, 141), (49, 130, 189)]
    for yi, label in enumerate(labels):
        y = top + yi * row_h
        body.append(f'<text x="10" y="{y + 19}" class="small">{esc(label)}</text>')
        for xi, value in enumerate(data_rows[yi]):
            x = left + xi * cell_w
            color = heat_color(value, max_values[yi], bases[yi])
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h - 4}" fill="{color}"/>')
    for boundary in boundaries:
        x = left + int((boundary - 1) / length * plot_w)
        body.append(f'<line x1="{x}" y1="{top - 6}" x2="{x}" y2="{top + row_h * len(labels)}" stroke="#1f2937" stroke-width="0.5" opacity="0.35"/>')
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        x = left + frac * plot_w
        pos = round(1 + frac * (length - 1))
        body.append(f'<line x1="{x}" y1="{top + row_h * len(labels) + 4}" x2="{x}" y2="{top + row_h * len(labels) + 10}" class="axis"/>')
        body.append(f'<text x="{x - 16}" y="{top + row_h * len(labels) + 26}" class="small">{pos}</text>')
    write_svg(POSITION_PLOT_DIR / f"{gene.lower()}_grouped_cds_position_overview.svg", "\n".join(body), width, height)


def euclidean(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def kmeans(features: list[list[float]], k: int = 5, iterations: int = 50) -> list[int]:
    if not features:
        return []
    seeds = [0, len(features) // 4, len(features) // 2, 3 * len(features) // 4, len(features) - 1]
    centers = [features[index] for index in seeds[:k]]
    assignments = [0] * len(features)
    for _ in range(iterations):
        changed = False
        for i, feature in enumerate(features):
            distances = [euclidean(feature, center) for center in centers]
            assigned = distances.index(min(distances))
            if assigned != assignments[i]:
                assignments[i] = assigned
                changed = True
        new_centers = []
        for cluster in range(k):
            members = [feature for feature, assigned in zip(features, assignments) if assigned == cluster]
            if not members:
                new_centers.append(centers[cluster])
                continue
            new_centers.append([sum(values) / len(members) for values in zip(*members)])
        centers = new_centers
        if not changed:
            break
    return assignments


def clustered_position_analysis(
    rows: list[dict[str, object]],
    boundaries: dict[str, list[int]],
    *,
    include_spliceai_features: bool,
    suffix: str,
    title_suffix: str,
) -> None:
    summary_rows = []
    segment_rows = []
    for gene in ["BRCA1", "BRCA2"]:
        length = CDS_LENGTHS[gene]
        bins = 120
        bin_rows = [[] for _ in range(bins)]
        for row in rows:
            if row["gene"] != gene or not isinstance(row.get("cds_pos"), int):
                continue
            idx = min(bins - 1, max(0, int((int(row["cds_pos"]) - 1) / length * bins)))
            bin_rows[idx].append(row)
        features = []
        table_rows = []
        for idx, items in enumerate(bin_rows):
            total = len(items) or 1
            benign = sum(1 for row in items if class_group(row["predicted_class"]) == "benign") / total
            vus = sum(1 for row in items if class_group(row["predicted_class"]) == "vus") / total
            pathogenic = sum(1 for row in items if class_group(row["predicted_class"]) == "pathogenic") / total
            high_splice = sum(1 for row in items if row["spliceai"] is not None and row["spliceai"] >= 0.20) / total
            mean_splice = sum(float(row["spliceai"] or 0.0) for row in items) / total
            near_boundary = sum(1 for row in items if isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 10) / total
            pvs1 = sum(1 for row in items if "PVS1" in row["criteria_codes"]) / total
            pp3 = sum(1 for row in items if "PP3" in row["criteria_codes"]) / total
            bs3 = sum(1 for row in items if "BS3" in row["criteria_codes"]) / total
            bp1 = sum(1 for row in items if "BP1" in row["criteria_codes"]) / total
            feature_row = [benign, vus, pathogenic, near_boundary, pvs1, pp3, bs3, bp1]
            if include_spliceai_features:
                feature_row = [benign, vus, pathogenic, high_splice, mean_splice, near_boundary, pvs1, pp3, bs3, bp1]
            features.append(feature_row)
        assignments = kmeans(features, k=5)
        for idx, cluster in enumerate(assignments):
            start = int(idx * length / bins) + 1
            end = int((idx + 1) * length / bins)
            items = bin_rows[idx]
            total = len(items) or 1
            benign = sum(1 for row in items if class_group(row["predicted_class"]) == "benign") / total
            vus = sum(1 for row in items if class_group(row["predicted_class"]) == "vus") / total
            pathogenic = sum(1 for row in items if class_group(row["predicted_class"]) == "pathogenic") / total
            high_splice = sum(1 for row in items if row["spliceai"] is not None and row["spliceai"] >= 0.20) / total
            mean_splice = sum(float(row["spliceai"] or 0.0) for row in items) / total
            near_boundary = sum(1 for row in items if isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 10) / total
            pvs1 = sum(1 for row in items if "PVS1" in row["criteria_codes"]) / total
            pp3 = sum(1 for row in items if "PP3" in row["criteria_codes"]) / total
            bs3 = sum(1 for row in items if "BS3" in row["criteria_codes"]) / total
            bp1 = sum(1 for row in items if "BP1" in row["criteria_codes"]) / total
            output_features = [benign, vus, pathogenic, high_splice, mean_splice, near_boundary, pvs1, pp3, bs3, bp1]
            table_rows.append([gene, idx + 1, start, end, cluster + 1] + [f"{value:.4f}" for value in output_features])
        write_csv(
            TABLE_DIR / f"{gene.lower()}_position_clusters_{suffix}.csv",
            [
                "gene",
                "bin",
                "cds_start",
                "cds_end",
                "cluster",
                "benign_fraction",
                "vus_fraction",
                "pathogenic_fraction",
                "spliceai_ge_0_20_fraction",
                "mean_spliceai",
                "near_boundary_fraction",
                "pvs1_fraction",
                "pp3_fraction",
                "bs3_fraction",
                "bp1_fraction",
            ],
            table_rows,
        )
        summary_rows.extend(cluster_summary_rows(gene, table_rows))
        segment_rows.extend(cluster_segment_rows(gene, assignments, bins, length))
        draw_cluster_svg(gene, assignments, boundaries[gene], bins, length, suffix, title_suffix)
    write_csv(
        TABLE_DIR / f"position_cluster_summary_{suffix}.csv",
        [
            "gene",
            "cluster",
            "n_bins",
            "cds_start_min",
            "cds_end_max",
            "mean_benign_fraction",
            "mean_vus_fraction",
            "mean_pathogenic_fraction",
            "mean_spliceai_ge_0_20_fraction",
            "mean_near_boundary_fraction",
            "mean_pvs1_fraction",
            "mean_pp3_fraction",
            "mean_bs3_fraction",
            "mean_bp1_fraction",
        ],
        summary_rows,
    )
    write_csv(
        TABLE_DIR / f"position_cluster_segments_{suffix}.csv",
        ["gene", "cluster", "start_bin", "end_bin", "cds_start", "cds_end", "n_bins"],
        segment_rows,
    )
    draw_cluster_profile_heatmap(summary_rows, suffix, title_suffix)


def compare_cluster_runs() -> list[list[object]]:
    rows = []
    for gene in ["brca1", "brca2"]:
        with_path = TABLE_DIR / f"{gene}_position_clusters_with_spliceai.csv"
        without_path = TABLE_DIR / f"{gene}_position_clusters_without_spliceai.csv"
        with_rows = list(csv.DictReader(with_path.open(newline="", encoding="utf-8")))
        without_rows = list(csv.DictReader(without_path.open(newline="", encoding="utf-8")))
        total = min(len(with_rows), len(without_rows))
        same = sum(1 for left, right in zip(with_rows, without_rows) if left["cluster"] == right["cluster"])
        changed = total - same
        rows.append([gene.upper(), total, same, changed, f"{same / total * 100:.2f}" if total else "0.00"])
    write_csv(
        TABLE_DIR / "cluster_spliceai_ablation_comparison.csv",
        ["gene", "n_bins", "same_cluster_bins", "changed_cluster_bins", "same_percent"],
        rows,
    )
    return rows


def pathogenic_driver_counts(rows: list[dict[str, object]]) -> dict[str, int]:
    pathogenic_rows = [row for row in rows if str(row["predicted_class"]) in {"4", "5"}]
    return {
        "class_4_5_total": len(pathogenic_rows),
        "with_pp3": sum(1 for row in pathogenic_rows if "PP3" in row["criteria_codes"]),
        "with_spliceai_ge_0_20": sum(1 for row in pathogenic_rows if row["spliceai"] is not None and row["spliceai"] >= 0.20),
        "with_pvs1": sum(1 for row in pathogenic_rows if "PVS1" in row["criteria_codes"]),
        "with_pm5_ptc": sum(1 for row in pathogenic_rows if "PM5_PTC" in row["criteria_codes"]),
        "nonsense": sum(1 for row in pathogenic_rows if row["variant_type"] == "nonsense"),
        "missense": sum(1 for row in pathogenic_rows if row["variant_type"] == "missense"),
        "synonymous": sum(1 for row in pathogenic_rows if row["variant_type"] == "synonymous"),
        "initiation_codon": sum(1 for row in pathogenic_rows if row["variant_type"] == "initiation_codon"),
    }


def clustered_position_analysis_no_splice_no_class(rows: list[dict[str, object]], boundaries: dict[str, list[int]]) -> None:
    suffix = "no_splice_no_class"
    title_suffix = "without class and SpliceAI-derived features"
    summary_rows = []
    segment_rows = []
    for gene in ["BRCA1", "BRCA2"]:
        length = CDS_LENGTHS[gene]
        bins = 120
        bin_rows = [[] for _ in range(bins)]
        for row in rows:
            if row["gene"] != gene or not isinstance(row.get("cds_pos"), int):
                continue
            idx = min(bins - 1, max(0, int((int(row["cds_pos"]) - 1) / length * bins)))
            bin_rows[idx].append(row)
        features = []
        table_rows = []
        for items in bin_rows:
            total = len(items) or 1
            missense = sum(1 for row in items if row["variant_type"] == "missense") / total
            nonsense = sum(1 for row in items if row["variant_type"] == "nonsense") / total
            synonymous = sum(1 for row in items if row["variant_type"] == "synonymous") / total
            initiation = sum(1 for row in items if row["variant_type"] == "initiation_codon") / total
            near_boundary = sum(1 for row in items if isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 10) / total
            pvs1 = sum(1 for row in items if "PVS1" in row["criteria_codes"]) / total
            pm5_ptc = sum(1 for row in items if "PM5_PTC" in row["criteria_codes"]) / total
            ps3 = sum(1 for row in items if "PS3" in row["criteria_codes"]) / total
            bs3 = sum(1 for row in items if "BS3" in row["criteria_codes"]) / total
            bp1 = sum(1 for row in items if "BP1" in row["criteria_codes"]) / total
            pm2 = sum(1 for row in items if "PM2_Supporting" in row["criteria_codes"]) / total
            features.append([missense, nonsense, synonymous, initiation, near_boundary, pvs1, pm5_ptc, ps3, bs3, bp1, pm2])
        assignments = kmeans(features, k=5)
        for idx, cluster in enumerate(assignments):
            start = int(idx * length / bins) + 1
            end = int((idx + 1) * length / bins)
            items = bin_rows[idx]
            total = len(items) or 1
            benign = sum(1 for row in items if class_group(row["predicted_class"]) == "benign") / total
            vus = sum(1 for row in items if class_group(row["predicted_class"]) == "vus") / total
            pathogenic = sum(1 for row in items if class_group(row["predicted_class"]) == "pathogenic") / total
            high_splice = sum(1 for row in items if row["spliceai"] is not None and row["spliceai"] >= 0.20) / total
            mean_splice = sum(float(row["spliceai"] or 0.0) for row in items) / total
            near_boundary = sum(1 for row in items if isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 10) / total
            pvs1 = sum(1 for row in items if "PVS1" in row["criteria_codes"]) / total
            pp3 = sum(1 for row in items if "PP3" in row["criteria_codes"]) / total
            bs3 = sum(1 for row in items if "BS3" in row["criteria_codes"]) / total
            bp1 = sum(1 for row in items if "BP1" in row["criteria_codes"]) / total
            output_features = [benign, vus, pathogenic, high_splice, mean_splice, near_boundary, pvs1, pp3, bs3, bp1]
            table_rows.append([gene, idx + 1, start, end, cluster + 1] + [f"{value:.4f}" for value in output_features])
        write_csv(
            TABLE_DIR / f"{gene.lower()}_position_clusters_{suffix}.csv",
            [
                "gene",
                "bin",
                "cds_start",
                "cds_end",
                "cluster",
                "benign_fraction",
                "vus_fraction",
                "pathogenic_fraction",
                "spliceai_ge_0_20_fraction",
                "mean_spliceai",
                "near_boundary_fraction",
                "pvs1_fraction",
                "pp3_fraction",
                "bs3_fraction",
                "bp1_fraction",
            ],
            table_rows,
        )
        summary_rows.extend(cluster_summary_rows(gene, table_rows))
        segment_rows.extend(cluster_segment_rows(gene, assignments, bins, length))
        draw_cluster_svg(gene, assignments, boundaries[gene], bins, length, suffix, title_suffix)
    write_csv(
        TABLE_DIR / f"position_cluster_summary_{suffix}.csv",
        [
            "gene",
            "cluster",
            "n_bins",
            "cds_start_min",
            "cds_end_max",
            "mean_benign_fraction",
            "mean_vus_fraction",
            "mean_pathogenic_fraction",
            "mean_spliceai_ge_0_20_fraction",
            "mean_near_boundary_fraction",
            "mean_pvs1_fraction",
            "mean_pp3_fraction",
            "mean_bs3_fraction",
            "mean_bp1_fraction",
        ],
        summary_rows,
    )
    write_csv(
        TABLE_DIR / f"position_cluster_segments_{suffix}.csv",
        ["gene", "cluster", "start_bin", "end_bin", "cds_start", "cds_end", "n_bins"],
        segment_rows,
    )
    draw_cluster_profile_heatmap(summary_rows, suffix, title_suffix)


def cluster_summary_rows(gene: str, table_rows: list[list[object]]) -> list[list[object]]:
    by_cluster: dict[int, list[list[object]]] = defaultdict(list)
    for row in table_rows:
        by_cluster[int(row[4])].append(row)
    summary = []
    for cluster, rows in sorted(by_cluster.items()):
        n = len(rows)
        start_min = min(int(row[2]) for row in rows)
        end_max = max(int(row[3]) for row in rows)
        means = []
        for col in [5, 6, 7, 8, 10, 11, 12, 13, 14]:
            means.append(sum(float(row[col]) for row in rows) / n)
        summary.append([gene, cluster, n, start_min, end_max] + [f"{value:.4f}" for value in means])
    return summary


def cluster_segment_rows(gene: str, assignments: list[int], bins: int, length: int) -> list[list[object]]:
    segments = []
    if not assignments:
        return segments
    start_idx = 0
    current = assignments[0]
    for idx, cluster in enumerate(assignments[1:], start=1):
        if cluster == current:
            continue
        segments.append(cluster_segment_row(gene, current + 1, start_idx, idx - 1, bins, length))
        start_idx = idx
        current = cluster
    segments.append(cluster_segment_row(gene, current + 1, start_idx, len(assignments) - 1, bins, length))
    return segments


def cluster_segment_row(gene: str, cluster: int, start_idx: int, end_idx: int, bins: int, length: int) -> list[object]:
    cds_start = int(start_idx * length / bins) + 1
    cds_end = int((end_idx + 1) * length / bins)
    return [gene, cluster, start_idx + 1, end_idx + 1, cds_start, cds_end, end_idx - start_idx + 1]


def draw_cluster_profile_heatmap(summary_rows: list[list[object]], suffix: str, title_suffix: str) -> None:
    features = [
        ("benign", 5),
        ("vus", 6),
        ("pathogenic", 7),
        ("SpAI >=0.20", 8),
        ("near boundary", 9),
        ("PVS1", 10),
        ("PP3", 11),
        ("BS3", 12),
        ("BP1", 13),
    ]
    x_labels = [f"{row[0]} C{row[1]}" for row in summary_rows]
    values = {}
    for row, x_label in zip(summary_rows, x_labels):
        for label, col in features:
            values[(label, x_label)] = round(float(row[col]), 3)
    heatmap(
        path=CLUSTER_PLOT_DIR / f"cluster_profile_heatmap_{suffix}.svg",
        title=f"Cluster mean feature profile {title_suffix}",
        x_labels=x_labels,
        y_labels=[label for label, _ in features],
        values=values,
        x_title="Gene and cluster",
    )
    for gene in ["BRCA1", "BRCA2"]:
        gene_rows = [row for row in summary_rows if row[0] == gene]
        draw_gene_cluster_profile_heatmap(gene_rows, suffix, title_suffix)


def draw_gene_cluster_profile_heatmap(summary_rows: list[list[object]], suffix: str, title_suffix: str) -> None:
    if not summary_rows:
        return
    gene = str(summary_rows[0][0])
    features = [
        ("benign", 5),
        ("vus", 6),
        ("pathogenic", 7),
        ("SpAI >=0.20", 8),
        ("near boundary", 9),
        ("PVS1", 10),
        ("PP3", 11),
        ("BS3", 12),
        ("BP1", 13),
    ]
    cell_w = 96
    cell_h = 30
    left = 150
    top = 116
    x_labels = [f"C{row[1]}" for row in summary_rows]
    y_labels = [label for label, _ in features]
    width = left + cell_w * len(x_labels) + 130
    height = top + cell_h * len(y_labels) + 75
    max_value = max(float(row[col]) for row in summary_rows for _, col in features)
    body = [f'<text x="30" y="32" class="title">{gene} cluster profile {esc(title_suffix)}</text>']
    body.append(f'<text x="{left}" y="58" class="small">Colored strip maps cluster IDs to the position cluster plot</text>')
    panel_w = cell_w * len(x_labels)
    panel_h = cell_h * len(y_labels)
    body.append(f'<rect x="{left - 2}" y="{top - 2}" width="{panel_w + 2}" height="{panel_h + 2}" class="panel"/>')
    for xi, (row, x_label) in enumerate(zip(summary_rows, x_labels)):
        x = left + xi * cell_w
        cluster = int(row[1])
        body.append(f'<text x="{x + 6}" y="78" class="small">{esc(x_label)}</text>')
        body.append(f'<rect x="{x}" y="86" width="{cell_w}" height="16" fill="{CLUSTER_COLORS[cluster - 1]}" class="tile"/>')
    for yi, y_label in enumerate(y_labels):
        y = top + yi * cell_h
        body.append(f'<text x="20" y="{y + 20}" class="small">{esc(y_label)}</text>')
        for xi, row in enumerate(summary_rows):
            x = left + xi * cell_w
            value = float(row[features[yi][1]])
            color = heat_color(value, max_value)
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{color}" class="tile"/>')
            if value:
                body.append(f'<text x="{x + 6}" y="{y + 19}" class="small">{value:.2f}</text>')
    legend_x = left + panel_w + 34
    legend_y = top + 4
    body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">mean</text>')
    for step in range(6):
        frac = 1 - step / 5
        color = heat_color(frac * max_value, max_value)
        y = legend_y + step * 18
        body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{color}" class="tile"/>')
    body.append(f'<text x="{legend_x + 26}" y="{legend_y + 12}" class="small">{max_value:.2f}</text>')
    body.append(f'<text x="{legend_x + 26}" y="{legend_y + 102}" class="small">0</text>')
    write_svg(CLUSTER_PLOT_DIR / f"{gene.lower()}_cluster_profile_heatmap_{suffix}.svg", "\n".join(body), width, height)


def draw_cluster_svg(
    gene: str,
    assignments: list[int],
    boundaries: list[int],
    bins: int,
    length: int,
    suffix: str,
    title_suffix: str,
) -> None:
    cell_w = 7
    left = 75
    top = 80
    plot_w = cell_w * bins
    width = left + plot_w + 40
    height = 205
    body = [f'<text x="28" y="34" class="title">{gene} k-means clusters along CDS {esc(title_suffix)}</text>']
    body.append(f'<text x="{left}" y="58" class="small">Compare with and without SpliceAI features; vertical ticks mark CDS exon boundaries</text>')
    for idx, cluster in enumerate(assignments):
        x = left + idx * cell_w
        body.append(f'<rect x="{x}" y="{top}" width="{cell_w}" height="42" fill="{CLUSTER_COLORS[cluster]}"/>')
    for boundary in boundaries:
        x = left + int((boundary - 1) / length * plot_w)
        body.append(f'<line x1="{x}" y1="{top - 8}" x2="{x}" y2="{top + 50}" stroke="#1f2937" stroke-width="0.5" opacity="0.35"/>')
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        x = left + frac * plot_w
        pos = round(1 + frac * (length - 1))
        body.append(f'<line x1="{x}" y1="{top + 48}" x2="{x}" y2="{top + 54}" class="axis"/>')
        body.append(f'<text x="{x - 16}" y="{top + 70}" class="small">{pos}</text>')
    legend_y = top + 100
    for cluster in range(5):
        x = left + cluster * 110
        body.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{CLUSTER_COLORS[cluster]}"/>')
        body.append(f'<text x="{x + 20}" y="{legend_y + 12}" class="small">cluster {cluster + 1}</text>')
    write_svg(CLUSTER_PLOT_DIR / f"{gene.lower()}_position_clusters_{suffix}.svg", "\n".join(body), width, height)


def plot_splice_boundary_scatter(rows: list[dict[str, object]]) -> None:
    points = [
        row
        for row in rows
        if isinstance(row.get("boundary_distance"), int) and row["spliceai"] is not None and row["spliceai"] >= 0.05
    ]
    width = 860
    height = 460
    left = 70
    top = 55
    plot_w = 720
    plot_h = 330
    max_dist = 250
    body = [f'<text x="28" y="32" class="title">SpliceAI versus distance to nearest CDS exon boundary</text>']
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    for tick in [0, 50, 100, 150, 200, 250]:
        x = left + tick / max_dist * plot_w
        body.append(f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + plot_h}" class="grid"/>')
        body.append(f'<text x="{x - 8}" y="{top + plot_h + 20}" class="small">{tick}</text>')
    for tick in [0, 0.2, 0.5, 0.8, 1.0]:
        y = top + plot_h - tick * plot_h
        body.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_w}" y2="{y}" class="grid"/>')
        body.append(f'<text x="32" y="{y + 4}" class="small">{tick:.1f}</text>')
    for row in points:
        dist = min(int(row["boundary_distance"]), max_dist)
        score = float(row["spliceai"])
        x = left + dist / max_dist * plot_w
        y = top + plot_h - score * plot_h
        cls = str(row["predicted_class"])
        body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.2" fill="{CLASS_COLORS[cls]}" opacity="0.55"/>')
    body.append(f'<text x="{left + 225}" y="{height - 25}" class="label">Distance to nearest CDS exon boundary, capped at 250 bp</text>')
    body.append(f'<text x="16" y="{top + 155}" class="label" transform="rotate(-90 16 {top + 155})">SpliceAI score</text>')
    legend_x = left + plot_w - 130
    for i, cls in enumerate(sorted(CLASS_LABELS, key=int)):
        y = top + 18 + i * 20
        body.append(f'<rect x="{legend_x}" y="{y - 10}" width="12" height="12" fill="{CLASS_COLORS[cls]}"/>')
        body.append(f'<text x="{legend_x + 18}" y="{y}" class="small">{cls} {CLASS_LABELS[cls]}</text>')
    write_svg(BOUNDARY_PLOT_DIR / "spliceai_vs_boundary_distance.svg", "\n".join(body), width, height)


def plot_grouped_splice_boundary_scatter(
    rows: list[dict[str, object]],
    *,
    gene: str | None = None,
    filename_prefix: str = "",
) -> list[list[object]]:
    points_by_group = {group: [] for group in GROUP_ORDER}
    for row in rows:
        if gene is not None and row.get("gene") != gene:
            continue
        if not isinstance(row.get("boundary_distance"), int) or row["spliceai"] is None:
            continue
        if float(row["spliceai"]) < 0.05:
            continue
        points_by_group[class_group(row["predicted_class"])].append(row)

    width = 1180
    height = 460
    left = 65
    top = 92
    panel_w = 330
    panel_h = 265
    panel_gap = 38
    max_dist = 250
    title_gene = f"{gene}: " if gene else ""
    body = [f'<text x="28" y="34" class="title">{title_gene}SpliceAI versus boundary distance</text>']
    body.append(f'<text x="28" y="56" class="small">Grouped by generated class. Only points with SpliceAI >= 0.05 are shown; distance is capped at 250 bp.</text>')
    for gi, group in enumerate(GROUP_ORDER):
        x0 = left + gi * (panel_w + panel_gap)
        y0 = top
        body.append(f'<text x="{x0}" y="{y0 - 18}" class="label">{esc(GROUP_LABELS[group])}</text>')
        body.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
        for tick in [0, 50, 100, 150, 200, 250]:
            x = x0 + tick / max_dist * panel_w
            body.append(f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y0 + panel_h}" class="grid"/>')
            body.append(f'<text x="{x - 8}" y="{y0 + panel_h + 18}" class="small">{tick}</text>')
        for tick in [0, 0.2, 0.5, 0.8, 1.0]:
            y = y0 + panel_h - tick * panel_h
            body.append(f'<line x1="{x0}" y1="{y}" x2="{x0 + panel_w}" y2="{y}" class="grid"/>')
            if gi == 0:
                body.append(f'<text x="28" y="{y + 4}" class="small">{tick:.1f}</text>')
        for row in points_by_group[group]:
            dist = min(int(row["boundary_distance"]), max_dist)
            score = float(row["spliceai"])
            x = x0 + dist / max_dist * panel_w
            y = y0 + panel_h - score * panel_h
            body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.0" fill="{GROUP_COLORS[group]}" opacity="0.52"/>')
        n_points = len(points_by_group[group])
        n_high = sum(1 for row in points_by_group[group] if float(row["spliceai"]) >= 0.20)
        body.append(f'<text x="{x0}" y="{y0 + panel_h + 40}" class="small">shown: {n_points}; SpAI >=0.20: {n_high}</text>')
    body.append(f'<text x="{left + 380}" y="{height - 20}" class="label">Distance to nearest CDS exon boundary, bp</text>')
    body.append(f'<text x="16" y="{top + 145}" class="label" transform="rotate(-90 16 {top + 145})">SpliceAI score</text>')
    write_svg(BOUNDARY_PLOT_DIR / f"{filename_prefix}spliceai_vs_boundary_distance_grouped_classes.svg", "\n".join(body), width, height)

    bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"]
    summary_rows = []
    for group in GROUP_ORDER:
        group_rows = [
            row
            for row in rows
            if (gene is None or row.get("gene") == gene)
            and class_group(row["predicted_class"]) == group
            and isinstance(row.get("boundary_distance"), int)
        ]
        for b in bins:
            in_bin = [row for row in group_rows if boundary_bin(row["boundary_distance"]) == b]
            high = [row for row in in_bin if row["spliceai"] is not None and float(row["spliceai"]) >= 0.20]
            total = len(in_bin)
            summary_rows.append(
                [
                    gene or "both",
                    group,
                    GROUP_LABELS[group],
                    b,
                    total,
                    len(high),
                    f"{(len(high) / total * 100):.2f}" if total else "0.00",
                ]
            )
    suffix = filename_prefix.rstrip("_")
    table_name = "spliceai_boundary_distance_by_grouped_class.csv" if not suffix else f"{suffix}_spliceai_boundary_distance_by_grouped_class.csv"
    write_csv(
        TABLE_DIR / table_name,
        ["gene", "group", "label", "boundary_distance_bin", "total_variants", "spliceai_ge_0_20", "percent"],
        summary_rows,
    )
    plot_grouped_boundary_rate_bars(summary_rows, gene=gene, filename_prefix=filename_prefix)
    plot_grouped_boundary_distribution_bars(summary_rows, gene=gene, filename_prefix=filename_prefix)
    return summary_rows


def plot_grouped_boundary_rate_bars(
    summary_rows: list[list[object]],
    *,
    gene: str | None = None,
    filename_prefix: str = "",
) -> None:
    bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"]
    values = {
        (str(row[1]), str(row[3])): float(row[6])
        for row in summary_rows
    }
    width = 1080
    height = 470
    left = 90
    top = 70
    plot_w = 870
    plot_h = 300
    group_gap = 10
    bin_w = plot_w / len(bins)
    bar_w = (bin_w - group_gap) / len(GROUP_ORDER)
    max_value = max(values.values()) if values else 1.0
    title_gene = f"{gene}: " if gene else ""
    body = [f'<text x="28" y="34" class="title">{title_gene}SpliceAI >=0.20 rate by distance</text>']
    body.append(f'<text x="28" y="54" class="small">This is normalized within each distance bin and group, so it is more informative than the scatter</text>')
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    for tick in [0, 20, 40, 60, 80, 100]:
        y = top + plot_h - (tick / 100) * plot_h
        body.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_w}" y2="{y}" class="grid"/>')
        body.append(f'<text x="42" y="{y + 4}" class="small">{tick}%</text>')
    for bi, bin_name in enumerate(bins):
        x_bin = left + bi * bin_w
        body.append(f'<text x="{x_bin + 4}" y="{top + plot_h + 20}" class="small">{esc(bin_name)}</text>')
        for gi, group in enumerate(GROUP_ORDER):
            value = values.get((group, bin_name), 0.0)
            bar_h = (value / max_value) * plot_h if max_value else 0
            x = x_bin + gi * bar_w + 4
            y = top + plot_h - bar_h
            body.append(f'<rect x="{x}" y="{y}" width="{bar_w - 2}" height="{bar_h}" fill="{GROUP_COLORS[group]}" opacity="0.85"/>')
    legend_x = left + plot_w - 245
    for i, group in enumerate(GROUP_ORDER):
        y = top + 18 + i * 20
        body.append(f'<rect x="{legend_x}" y="{y - 10}" width="12" height="12" fill="{GROUP_COLORS[group]}"/>')
        body.append(f'<text x="{legend_x + 18}" y="{y}" class="small">{esc(GROUP_LABELS[group])}</text>')
    body.append(f'<text x="{left + 315}" y="{height - 22}" class="label">Distance to nearest CDS exon boundary, bp</text>')
    write_svg(BOUNDARY_PLOT_DIR / f"{filename_prefix}spliceai_high_rate_by_boundary_grouped_classes.svg", "\n".join(body), width, height)


def plot_grouped_boundary_distribution_bars(
    summary_rows: list[list[object]],
    *,
    gene: str | None = None,
    filename_prefix: str = "",
) -> None:
    bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"]
    totals_by_group = {
        group: sum(int(row[4]) for row in summary_rows if row[1] == group)
        for group in GROUP_ORDER
    }
    values = {}
    for row in summary_rows:
        group = str(row[1])
        bin_name = str(row[3])
        total = totals_by_group[group]
        values[(group, bin_name)] = int(row[4]) / total * 100 if total else 0.0
    width = 1080
    height = 470
    left = 90
    top = 70
    plot_w = 870
    plot_h = 300
    group_gap = 10
    bin_w = plot_w / len(bins)
    bar_w = (bin_w - group_gap) / len(GROUP_ORDER)
    max_value = max(values.values()) if values else 1.0
    title_gene = f"{gene}: " if gene else ""
    body = [f'<text x="28" y="34" class="title">{title_gene}Boundary-distance distribution</text>']
    body.append(f'<text x="28" y="54" class="small">Percent of all variants in each group that fall into each distance bin</text>')
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    for tick in [0, 20, 40, 60, 80, 100]:
        y = top + plot_h - (tick / 100) * plot_h
        body.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_w}" y2="{y}" class="grid"/>')
        body.append(f'<text x="42" y="{y + 4}" class="small">{tick}%</text>')
    for bi, bin_name in enumerate(bins):
        x_bin = left + bi * bin_w
        body.append(f'<text x="{x_bin + 4}" y="{top + plot_h + 20}" class="small">{esc(bin_name)}</text>')
        for gi, group in enumerate(GROUP_ORDER):
            value = values.get((group, bin_name), 0.0)
            bar_h = (value / max_value) * plot_h if max_value else 0
            x = x_bin + gi * bar_w + 4
            y = top + plot_h - bar_h
            body.append(f'<rect x="{x}" y="{y}" width="{bar_w - 2}" height="{bar_h}" fill="{GROUP_COLORS[group]}" opacity="0.85"/>')
    legend_x = left + plot_w - 245
    for i, group in enumerate(GROUP_ORDER):
        y = top + 18 + i * 20
        body.append(f'<rect x="{legend_x}" y="{y - 10}" width="12" height="12" fill="{GROUP_COLORS[group]}"/>')
        body.append(f'<text x="{legend_x + 18}" y="{y}" class="small">{esc(GROUP_LABELS[group])}</text>')
    body.append(f'<text x="{left + 315}" y="{height - 22}" class="label">Distance to nearest CDS exon boundary, bp</text>')
    write_svg(BOUNDARY_PLOT_DIR / f"{filename_prefix}boundary_distance_distribution_grouped_classes.svg", "\n".join(body), width, height)


def plot_gene_boundary_comparison_heatmaps(summary_rows: list[list[object]]) -> None:
    bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"]
    x_labels = [f"{gene} {group}" for gene in ["BRCA1", "BRCA2"] for group in GROUP_ORDER]

    high_rate_values = {}
    distribution_values = {}
    totals_by_gene_group = Counter()
    for row in summary_rows:
        gene = str(row[0])
        group = str(row[1])
        totals_by_gene_group[(gene, group)] += int(row[4])
    for row in summary_rows:
        gene = str(row[0])
        group = str(row[1])
        bin_name = str(row[3])
        x_label = f"{gene} {group}"
        high_rate_values[(bin_name, x_label)] = float(row[6])
        denom = totals_by_gene_group[(gene, group)]
        distribution_values[(bin_name, x_label)] = int(row[4]) / denom * 100 if denom else 0.0

    percent_heatmap(
        path=BOUNDARY_PLOT_DIR / "gene_normalized_spliceai_high_rate_by_boundary_group_heatmap.svg",
        title="SpliceAI high-rate by gene and group",
        x_labels=x_labels,
        y_labels=bins,
        values=high_rate_values,
        x_title="Gene and grouped class",
        note="Cells show percent of variants with SpliceAI >=0.20 within each gene, group, and boundary-distance bin.",
    )
    percent_heatmap(
        path=BOUNDARY_PLOT_DIR / "gene_normalized_boundary_distribution_by_group_heatmap.svg",
        title="Boundary distribution by gene and group",
        x_labels=x_labels,
        y_labels=bins,
        values=distribution_values,
        x_title="Gene and grouped class",
        note="Cells show percent of each gene/group that falls into each boundary-distance bin.",
    )


def plot_boundary_spliceai_binned_heatmap(rows: list[dict[str, object]], gene: str) -> None:
    gene_rows = [
        row
        for row in rows
        if row.get("gene") == gene
        and isinstance(row.get("boundary_distance"), int)
        and row.get("spliceai") is not None
    ]
    boundary_bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", "101-250", ">250"]
    splice_bins = ["0-0.049", "0.05-0.099", "0.10-0.199", "0.20-0.499", "0.50-0.799", "0.80-1.00"]
    cell_w = 58
    cell_h = 28
    panel_gap = 42
    left = 95
    top = 112
    panel_w = cell_w * len(boundary_bins)
    panel_h = cell_h * len(splice_bins)
    width = left + len(GROUP_ORDER) * panel_w + (len(GROUP_ORDER) - 1) * panel_gap + 80
    height = top + panel_h + 112

    counts = Counter()
    totals = Counter()
    for row in gene_rows:
        group = class_group(row["predicted_class"])
        boundary = boundary_bin_compact(row["boundary_distance"])
        splice = splice_bin_compact(row["spliceai"])
        counts[(group, splice, boundary)] += 1
        totals[group] += 1
    percentages = {
        key: count / totals[key[0]] * 100 if totals[key[0]] else 0.0
        for key, count in counts.items()
    }
    max_percent = max(percentages.values()) if percentages else 0.0

    body = [
        f'<text x="28" y="32" class="title">{gene}: SpliceAI and boundary context</text>',
        '<text x="28" y="54" class="small">Each tile is percent of variants within the generated class group. This replaces the dense scatter view.</text>',
        f'<text x="{left}" y="82" class="small">Boundary distance bin, bp</text>',
    ]
    for gi, group in enumerate(GROUP_ORDER):
        x0 = left + gi * (panel_w + panel_gap)
        body.append(f'<text x="{x0}" y="{top - 18}" class="label">{esc(GROUP_LABELS[group])} (n={totals[group]})</text>')
        body.append(f'<rect x="{x0 - 2}" y="{top - 2}" width="{panel_w + 2}" height="{panel_h + 2}" class="panel"/>')
        for xi, boundary in enumerate(boundary_bins):
            x = x0 + xi * cell_w
            body.append(f'<text x="{x + 4}" y="{top - 4}" class="small">{esc(boundary)}</text>')
        for yi, splice in enumerate(splice_bins):
            y = top + yi * cell_h
            if gi == 0:
                body.append(f'<text x="18" y="{y + 18}" class="small">{esc(splice)}</text>')
            for xi, boundary in enumerate(boundary_bins):
                x = x0 + xi * cell_w
                value = percentages.get((group, splice, boundary), 0.0)
                color = heat_color(value, max_percent)
                body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{color}" class="tile"/>')
                if value >= 0.5:
                    body.append(f'<text x="{x + 4}" y="{y + 18}" class="small">{value:.1f}</text>')

    body.append(f'<text x="16" y="{top + 94}" class="label" transform="rotate(-90 16 {top + 94})">SpliceAI score bin</text>')
    legend_x = width - 62
    legend_y = top
    body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">percent</text>')
    for step in range(6):
        frac = 1 - step / 5
        color = heat_color(frac * max_percent, max_percent)
        y = legend_y + step * 18
        body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{color}" class="tile"/>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 12}" class="small">{max_percent:.1f}</text>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 102}" class="small">0</text>')
    body.append(f'<text x="{left + panel_w}" y="{height - 26}" class="label">Distance to nearest CDS exon boundary, bp</text>')

    write_svg(BOUNDARY_PLOT_DIR / f"{gene.lower()}_boundary_spliceai_binned_heatmap.svg", "\n".join(body), width, height)

    table_rows = []
    for group in GROUP_ORDER:
        for splice in splice_bins:
            for boundary in boundary_bins:
                count = counts[(group, splice, boundary)]
                table_rows.append(
                    [
                        gene,
                        group,
                        GROUP_LABELS[group],
                        splice,
                        boundary,
                        count,
                        totals[group],
                        f"{(count / totals[group] * 100):.2f}" if totals[group] else "0.00",
                    ]
                )
    write_csv(
        TABLE_DIR / f"{gene.lower()}_boundary_spliceai_binned_heatmap.csv",
        ["gene", "group", "label", "spliceai_bin", "boundary_distance_bin", "count", "group_total", "percent_within_group"],
        table_rows,
    )


def main() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    boundaries = cds_boundaries_from_maps()
    add_boundary_distances(rows, boundaries)
    bar_class_distribution(rows)
    bar_group_distribution(rows)
    plot_criteria_heatmap(rows)
    plot_splice_bins_heatmap(rows)
    plot_grouped_splice_bins_heatmap(rows)
    plot_gene_normalized_overview(rows)
    plot_boundary_heatmap(rows)
    plot_grouped_boundary_heatmap(rows)
    for gene in ["BRCA1", "BRCA2"]:
        prefix = f"{gene.lower()}_"
        plot_boundary_heatmap(rows, gene=gene, filename_prefix=prefix)
        plot_grouped_boundary_heatmap(rows, gene=gene, filename_prefix=prefix)
    plot_position_tracks(rows, boundaries)
    plot_grouped_position_tracks(rows, boundaries)
    clustered_position_analysis(
        rows,
        boundaries,
        include_spliceai_features=True,
        suffix="with_spliceai",
        title_suffix="with SpliceAI features",
    )
    clustered_position_analysis(
        rows,
        boundaries,
        include_spliceai_features=False,
        suffix="without_spliceai",
        title_suffix="without SpliceAI features",
    )
    clustered_position_analysis_no_splice_no_class(rows, boundaries)
    cluster_comparison = compare_cluster_runs()
    plot_splice_boundary_scatter(rows)
    plot_grouped_splice_boundary_scatter(rows)
    gene_boundary_summaries = []
    for gene in ["BRCA1", "BRCA2"]:
        gene_boundary_summaries.extend(
            plot_grouped_splice_boundary_scatter(rows, gene=gene, filename_prefix=f"{gene.lower()}_")
        )
        plot_boundary_spliceai_binned_heatmap(rows, gene)
    plot_gene_boundary_comparison_heatmaps(gene_boundary_summaries)
    write_report(rows, cluster_comparison)
    print(f"Wrote plots to {PLOT_DIR}")


def write_report(rows: list[dict[str, object]], cluster_comparison: list[list[object]]) -> None:
    bins = ["0", "1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"]
    total_by_bin = Counter()
    high_by_bin = Counter()
    for row in rows:
        b = boundary_bin(row.get("boundary_distance"))
        total_by_bin[b] += 1
        if row["spliceai"] is not None and row["spliceai"] >= 0.20:
            high_by_bin[b] += 1
    high_near = sum(
        1
        for row in rows
        if row["spliceai"] is not None
        and row["spliceai"] >= 0.20
        and isinstance(row.get("boundary_distance"), int)
        and row["boundary_distance"] <= 10
    )
    high_total = sum(1 for row in rows if row["spliceai"] is not None and row["spliceai"] >= 0.20)
    overall_rate = high_total / len(rows) if rows else 0
    near_total = sum(
        1
        for row in rows
        if isinstance(row.get("boundary_distance"), int) and row["boundary_distance"] <= 10
    )
    bin_lines = []
    for b in bins:
        total = total_by_bin[b]
        high = high_by_bin[b]
        rate = high / total if total else 0
        enrichment = rate / overall_rate if overall_rate else 0
        bin_lines.append(f"| {b} | {total} | {high} | {rate * 100:.2f}% | {enrichment:.1f}x |")
    bin_table = "\n".join(bin_lines)
    driver_counts = pathogenic_driver_counts(rows)
    comparison_lines = "\n".join(
        f"| {row[0]} | {row[1]} | {row[2]} | {row[3]} | {row[4]}% |"
        for row in cluster_comparison
    )
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# Visualization Report: Precomputed BRCA Module 1 SNV Snapshot

Generated: {generated}

These plots are generated without external plotting dependencies. They use the precomputed Module 1 coding SNV snapshot and the local GRCh38 RefSeq CDS map for NM_007294.4 and NM_000059.4.

## Splice Boundary Summary

- Variants with SpliceAI >= 0.20: {high_total}
- Variants with SpliceAI >= 0.20 and within 10 bp of a CDS exon boundary: {high_near}
- All variants within 10 bp of a CDS exon boundary: {near_total}

This uses distance to the nearest CDS exon boundary, so it is a coding-SNV approximation of splice-site proximity. It does not include intronic variants.

Overall SpliceAI >= 0.20 rate across the coding SNV snapshot: {overall_rate * 100:.2f}%.

| Distance to CDS exon boundary | Total variants | SpliceAI >= 0.20 | Percent | Enrichment vs overall |
|---|---:|---:|---:|---:|
{bin_table}

## Plots

- `plots/01_overview/class_distribution.svg`
- `plots/01_overview/grouped_class_distribution.svg`
- `plots/01_overview/criteria_by_class_heatmap.svg`
- `plots/01_overview/spliceai_bins_by_class_heatmap.svg`
- `plots/01_overview/spliceai_bins_by_group_heatmap.svg`
- `plots/02_position/brca1_cds_position_overview.svg`
- `plots/02_position/brca2_cds_position_overview.svg`
- `plots/02_position/brca1_grouped_cds_position_overview.svg`
- `plots/02_position/brca2_grouped_cds_position_overview.svg`
- `plots/03_boundary_spliceai/boundary_distance_by_class_heatmap.svg`
- `plots/03_boundary_spliceai/boundary_distance_by_group_heatmap.svg`
- `plots/03_boundary_spliceai/spliceai_vs_boundary_distance.svg`
- `plots/03_boundary_spliceai/spliceai_vs_boundary_distance_grouped_classes.svg`
- `plots/03_boundary_spliceai/spliceai_high_rate_by_boundary_grouped_classes.svg`
- `plots/03_boundary_spliceai/boundary_distance_distribution_grouped_classes.svg`
- `plots/04_clusters/brca1_position_clusters_with_spliceai.svg`
- `plots/04_clusters/brca2_position_clusters_with_spliceai.svg`
- `plots/04_clusters/brca1_position_clusters_without_spliceai.svg`
- `plots/04_clusters/brca2_position_clusters_without_spliceai.svg`
- `plots/04_clusters/brca1_position_clusters_no_splice_no_class.svg`
- `plots/04_clusters/brca2_position_clusters_no_splice_no_class.svg`
- `plots/04_clusters/cluster_profile_heatmap_with_spliceai.svg`
- `plots/04_clusters/cluster_profile_heatmap_without_spliceai.svg`
- `plots/04_clusters/cluster_profile_heatmap_no_splice_no_class.svg`
- `plots/04_clusters/brca1_cluster_profile_heatmap_with_spliceai.svg`
- `plots/04_clusters/brca2_cluster_profile_heatmap_with_spliceai.svg`
- `plots/04_clusters/brca1_cluster_profile_heatmap_without_spliceai.svg`
- `plots/04_clusters/brca2_cluster_profile_heatmap_without_spliceai.svg`
- `plots/04_clusters/brca1_cluster_profile_heatmap_no_splice_no_class.svg`
- `plots/04_clusters/brca2_cluster_profile_heatmap_no_splice_no_class.svg`

## Derived Tables

- `tables/spliceai_high_by_boundary_distance.csv`
- `tables/spliceai_boundary_distance_by_grouped_class.csv`
- `tables/brca1_position_bins.csv`
- `tables/brca2_position_bins.csv`
- `tables/brca1_grouped_position_bins.csv`
- `tables/brca2_grouped_position_bins.csv`
- `tables/brca1_position_clusters_with_spliceai.csv`
- `tables/brca2_position_clusters_with_spliceai.csv`
- `tables/position_cluster_summary_with_spliceai.csv`
- `tables/position_cluster_segments_with_spliceai.csv`
- `tables/brca1_position_clusters_without_spliceai.csv`
- `tables/brca2_position_clusters_without_spliceai.csv`
- `tables/position_cluster_summary_without_spliceai.csv`
- `tables/position_cluster_segments_without_spliceai.csv`
- `tables/brca1_position_clusters_no_splice_no_class.csv`
- `tables/brca2_position_clusters_no_splice_no_class.csv`
- `tables/position_cluster_summary_no_splice_no_class.csv`
- `tables/position_cluster_segments_no_splice_no_class.csv`
- `tables/cluster_spliceai_ablation_comparison.csv`

## Cluster Analysis

The cluster plots group 120 CDS position bins per gene using a simple k-means model with k=5. The `with_spliceai` run uses grouped class fractions, high SpliceAI fraction, mean SpliceAI, near-boundary fraction, and selected criteria fractions. The `without_spliceai` run removes high SpliceAI and mean SpliceAI from the clustering feature space, but still reports them afterward as annotations. This is exploratory pattern finding, not a validated classifier.

Clusters can recur in several non-contiguous CDS regions. Use the `position_cluster_segments_*` tables for contiguous stretches and the `cluster_profile_heatmap_*` plots for the mean feature profile of each cluster.

The stricter `no_splice_no_class` run removes predicted class fractions and SpliceAI-derived criteria from the clustering feature space. It clusters on variant type, near-boundary fraction, and non-SpliceAI criteria fractions, then reports class and SpliceAI afterward as annotations.

## Pathogenic Driver Check

For class 4/5 variants in this coding SNV snapshot:

| Signal | Count |
|---|---:|
| Total class 4/5 | {driver_counts["class_4_5_total"]} |
| PP3 present | {driver_counts["with_pp3"]} |
| SpliceAI >= 0.20 | {driver_counts["with_spliceai_ge_0_20"]} |
| PVS1 present | {driver_counts["with_pvs1"]} |
| PM5_PTC present | {driver_counts["with_pm5_ptc"]} |
| Nonsense variants | {driver_counts["nonsense"]} |
| Missense variants | {driver_counts["missense"]} |
| Synonymous variants | {driver_counts["synonymous"]} |
| Initiation codon variants | {driver_counts["initiation_codon"]} |

In this automated coding SNV landscape, most class 4/5 calls are PVS1/PM5_PTC driven, not PP3/SpliceAI driven. SpliceAI is still a strong local signal around splice boundaries and is important for specific variants, but it is not the dominant source of pathogenic classifications across the full coding SNV snapshot.

## SpliceAI Ablation Check

| Gene | Bins compared | Same cluster | Changed cluster | Same percent |
|---|---:|---:|---:|---:|
{comparison_lines}

In this run, removing SpliceAI from the clustering feature space did not change the bin assignments. The SpliceAI signal visible in cluster profiles should therefore be interpreted as an annotation correlated with the class/criteria/boundary structure, not as an independent driver of the k-means split.

## Reading Notes

The position overview plots use 120 bins along the coding sequence. Vertical lines are CDS exon boundaries. Darker cells mean more variants in a class or higher maximum SpliceAI signal in that bin.

The boundary-distance plot caps distance at 250 bp for readability. Points with larger distances are placed at the right edge.
"""
    REPORT.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
