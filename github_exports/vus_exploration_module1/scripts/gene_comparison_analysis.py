"""Domain and broad-region context comparison for the coding SNV snapshot."""

from __future__ import annotations

import csv
import math
import re
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
TABLE_DIR = OUT_DIR / "tables" / "gene_comparison"
PLOT_DIR = OUT_DIR / "plots" / "07_gene_comparison"
REPORT = OUT_DIR / "gene_comparison_report.md"

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
VARIANT_TYPES = ["initiation_codon", "missense", "nonsense", "synonymous"]
CRITERIA = ["PVS1", "PM5_PTC", "PP3", "PS3", "BS3", "BP1", "BP7", "PM2_Supporting"]
DOMAIN_FEATURES = {
    "BRCA1": [
        ("RING zinc-binding/E3 ligase", 2, 101),
        ("coiled-coil PALB2 interaction", 1391, 1424),
        ("BRCT phosphopeptide-binding", 1650, 1857),
    ],
    "BRCA2": [
        ("PALB2 binding N-terminal", 10, 40),
        ("BRC repeats / RAD51-binding", 1000, 2080),
        ("DNA-binding domain / helical-OB-DSS1", 2481, 3186),
        ("C-terminal RAD51/NLS context", 3187, 3418),
    ],
}
PROTEIN_LENGTHS = {gene: length // 3 for gene, length in CDS_LENGTHS.items()}


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


def aa_position(p_notation: str) -> int | None:
    match = re.search(r"[A-Z][a-z]{2}(\d+)", p_notation or "")
    if match:
        return int(match.group(1))
    return None


def domain_feature(gene: str, pos: int | None) -> tuple[str, str]:
    if pos is None:
        return "unknown protein position", "unknown"
    for label, start, end in DOMAIN_FEATURES.get(gene, []):
        if start <= pos <= end:
            return label, f"{start}-{end}"
    return "outside mapped domains", "outside"


def domain_aa_length(gene: str, domain: str) -> int:
    for label, start, end in DOMAIN_FEATURES.get(gene, []):
        if label == domain:
            return end - start + 1
    if domain == "outside mapped domains":
        mapped = sum(end - start + 1 for _, start, end in DOMAIN_FEATURES.get(gene, []))
        return max(PROTEIN_LENGTHS[gene] - mapped, 0)
    return 0


def per_100aa(count: int, aa_length: int) -> str:
    return f"{count / aa_length * 100:.2f}" if aa_length else "0.00"


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


def load_rows() -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["cds_pos"] = parse_cds_pos(row["c_notation"])
            row["spliceai"] = parse_float(row["spliceai_score"])
            row["criteria_codes"] = parse_criteria(row.get("criteria", ""))
            row["group"] = class_group(row["predicted_class"])
            row["aa_pos"] = aa_position(row.get("p_notation", ""))
            row["domain_context"], row["domain_range"] = domain_feature(row["gene"], row["aa_pos"])
            rows.append(row)
    return rows


def pct(part: int, total: int) -> str:
    return f"{part / total * 100:.2f}" if total else "0.00"


def gene_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for gene in ["BRCA1", "BRCA2"]:
        items = [row for row in rows if row["gene"] == gene]
        total = len(items)
        record = {"gene": gene, "n_variants": total, "cds_length": CDS_LENGTHS[gene]}
        for group in GROUP_ORDER:
            count = sum(1 for row in items if row["group"] == group)
            record[f"{group}_count"] = count
            record[f"{group}_percent"] = pct(count, total)
        for vtype in VARIANT_TYPES:
            count = sum(1 for row in items if row["variant_type"] == vtype)
            record[f"{vtype}_count"] = count
            record[f"{vtype}_percent"] = pct(count, total)
        high = sum(1 for row in items if row["spliceai"] >= 0.20)
        record["high_spliceai_count"] = high
        record["high_spliceai_percent"] = pct(high, total)
        for criterion in CRITERIA:
            count = sum(1 for row in items if criterion in row["criteria_codes"])
            record[f"{criterion.lower()}_count"] = count
            record[f"{criterion.lower()}_percent"] = pct(count, total)
        output.append(record)
    return output


def domain_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    domain_order = {
        gene: [feature[0] for feature in features] + ["outside mapped domains", "unknown protein position"]
        for gene, features in DOMAIN_FEATURES.items()
    }
    for gene in ["BRCA1", "BRCA2"]:
        for domain in domain_order[gene]:
            items = [row for row in rows if row["gene"] == gene and row["domain_context"] == domain]
            if not items:
                continue
            total = len(items)
            aa_length = domain_aa_length(gene, domain)
            record = {
                "gene": gene,
                "domain_context": domain,
                "domain_range": next((row["domain_range"] for row in items if row["domain_context"] == domain), ""),
                "aa_length": aa_length,
                "n_variants": total,
                "variants_per_100aa": per_100aa(total, aa_length),
            }
            for group in GROUP_ORDER:
                count = sum(1 for row in items if row["group"] == group)
                record[f"{group}_count"] = count
                record[f"{group}_percent"] = pct(count, total)
                record[f"{group}_per_100aa"] = per_100aa(count, aa_length)
            high = sum(1 for row in items if row["spliceai"] >= 0.20)
            record["high_spliceai_count"] = high
            record["high_spliceai_percent"] = pct(high, total)
            record["high_spliceai_per_100aa"] = per_100aa(high, aa_length)
            for criterion in ["PVS1", "PM5_PTC", "PP3", "PS3", "BS3", "BP1"]:
                count = sum(1 for row in items if criterion in row["criteria_codes"])
                record[f"{criterion.lower()}_count"] = count
                record[f"{criterion.lower()}_percent"] = pct(count, total)
                record[f"{criterion.lower()}_per_100aa"] = per_100aa(count, aa_length)
            output.append(record)
    return output


def domain_vus_examples(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    selected = []
    for row in rows:
        if row["group"] != "vus" or row["domain_context"] in {"outside mapped domains", "unknown protein position"}:
            continue
        score = 0
        reasons = []
        if row["spliceai"] >= 0.20:
            score += 30
            reasons.append("SpliceAI>=0.20")
        for criterion, points in [("PS3", 25), ("PP3", 15), ("PVS1", 20), ("BS3", 8)]:
            if criterion in row["criteria_codes"]:
                score += points
                reasons.append(criterion)
        if score == 0:
            continue
        selected.append(
            {
                "gene": row["gene"],
                "domain_context": row["domain_context"],
                "score": score,
                "reasons": ";".join(reasons),
                "c_notation": row["c_notation"],
                "p_notation": row["p_notation"],
                "spliceai_score": row["spliceai_score"],
                "criteria": row["criteria"],
            }
        )
    return sorted(selected, key=lambda row: (-int(row["score"]), row["gene"], row["domain_context"], row["c_notation"]))


def normalized_bins(rows: list[dict[str, object]], bins: int = 120) -> list[dict[str, object]]:
    output = []
    for gene in ["BRCA1", "BRCA2"]:
        length = CDS_LENGTHS[gene]
        by_bin = [[] for _ in range(bins)]
        for row in rows:
            if row["gene"] != gene or not isinstance(row.get("cds_pos"), int):
                continue
            idx = min(bins - 1, max(0, int((int(row["cds_pos"]) - 1) / length * bins)))
            by_bin[idx].append(row)
        for idx, items in enumerate(by_bin):
            total = len(items)
            record = {
                "gene": gene,
                "bin": idx + 1,
                "relative_start_percent": f"{idx / bins * 100:.2f}",
                "relative_end_percent": f"{(idx + 1) / bins * 100:.2f}",
                "cds_start": int(idx * length / bins) + 1,
                "cds_end": int((idx + 1) * length / bins),
                "n_variants": total,
            }
            for group in GROUP_ORDER:
                count = sum(1 for row in items if row["group"] == group)
                record[f"{group}_count"] = count
                record[f"{group}_percent"] = pct(count, total)
            high = sum(1 for row in items if row["spliceai"] >= 0.20)
            record["high_spliceai_count"] = high
            record["high_spliceai_percent"] = pct(high, total)
            for vtype in VARIANT_TYPES:
                count = sum(1 for row in items if row["variant_type"] == vtype)
                record[f"{vtype}_percent"] = pct(count, total)
            output.append(record)
    return output


def bin_differences(bin_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_gene = {
        gene: {int(row["bin"]): row for row in bin_rows if row["gene"] == gene}
        for gene in ["BRCA1", "BRCA2"]
    }
    output = []
    for idx in range(1, 121):
        left = by_gene["BRCA1"][idx]
        right = by_gene["BRCA2"][idx]
        record = {
            "bin": idx,
            "relative_start_percent": left["relative_start_percent"],
            "relative_end_percent": left["relative_end_percent"],
        }
        for metric in ["benign_percent", "vus_percent", "pathogenic_percent", "high_spliceai_percent", "missense_percent", "nonsense_percent", "synonymous_percent"]:
            record[f"brca1_{metric}"] = left[metric]
            record[f"brca2_{metric}"] = right[metric]
            record[f"diff_{metric}"] = f"{float(left[metric]) - float(right[metric]):.2f}"
        output.append(record)
    return output


def draw_gene_group_barplot(summary: list[dict[str, object]]) -> None:
    width = 820
    height = 360
    left = 120
    top = 68
    plot_w = 560
    bar_h = 38
    gap = 28
    body = ['<text x="28" y="34" class="title">BRCA1 vs BRCA2 grouped class composition</text>']
    for gi, gene_row in enumerate(summary):
        y = top + gi * (bar_h + gap)
        body.append(f'<text x="28" y="{y + 24}" class="label">{gene_row["gene"]}</text>')
        x = left
        for group in GROUP_ORDER:
            value = float(gene_row[f"{group}_percent"])
            w = plot_w * value / 100
            body.append(f'<rect x="{x}" y="{y}" width="{w}" height="{bar_h}" fill="{GROUP_COLORS[group]}"/>')
            if w > 35:
                body.append(f'<text x="{x + 4}" y="{y + 24}" class="small">{value:.1f}%</text>')
            x += w
    legend_y = top + 2 * (bar_h + gap) + 22
    for i, group in enumerate(GROUP_ORDER):
        x = left + i * 180
        body.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{GROUP_COLORS[group]}"/>')
        body.append(f'<text x="{x + 20}" y="{legend_y + 12}" class="small">{esc(GROUP_LABELS[group])}</text>')
    write_svg(PLOT_DIR / "gene_grouped_class_composition.svg", "\n".join(body), width, height)


def draw_gene_variant_type_barplot(summary: list[dict[str, object]]) -> None:
    colors = {
        "initiation_codon": "#f59e0b",
        "missense": "#2563eb",
        "nonsense": "#dc2626",
        "synonymous": "#16a34a",
    }
    width = 860
    height = 360
    left = 120
    top = 68
    plot_w = 580
    bar_h = 38
    gap = 28
    body = ['<text x="28" y="34" class="title">BRCA1 vs BRCA2 variant type composition</text>']
    for gi, gene_row in enumerate(summary):
        y = top + gi * (bar_h + gap)
        body.append(f'<text x="28" y="{y + 24}" class="label">{gene_row["gene"]}</text>')
        x = left
        for vtype in VARIANT_TYPES:
            value = float(gene_row[f"{vtype}_percent"])
            w = plot_w * value / 100
            body.append(f'<rect x="{x}" y="{y}" width="{w}" height="{bar_h}" fill="{colors[vtype]}"/>')
            if w > 38:
                body.append(f'<text x="{x + 4}" y="{y + 24}" class="small">{value:.1f}%</text>')
            x += w
    legend_y = top + 2 * (bar_h + gap) + 22
    for i, vtype in enumerate(VARIANT_TYPES):
        x = left + i * 155
        body.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{colors[vtype]}"/>')
        body.append(f'<text x="{x + 20}" y="{legend_y + 12}" class="small">{esc(vtype)}</text>')
    write_svg(PLOT_DIR / "gene_variant_type_composition.svg", "\n".join(body), width, height)


def draw_normalized_profile(bin_rows: list[dict[str, object]], metric: str, title: str, filename: str) -> None:
    width = 1020
    height = 390
    left = 70
    top = 62
    plot_w = 850
    plot_h = 250
    body = [f'<text x="28" y="34" class="title">{esc(title)}</text>']
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    for tick in [0, 20, 40, 60, 80, 100]:
        y = top + plot_h - tick / 100 * plot_h
        body.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_w}" y2="{y}" class="grid"/>')
        body.append(f'<text x="34" y="{y + 4}" class="small">{tick}%</text>')
    colors = {"BRCA1": "#2563eb", "BRCA2": "#dc2626"}
    for gene in ["BRCA1", "BRCA2"]:
        rows = sorted([row for row in bin_rows if row["gene"] == gene], key=lambda row: int(row["bin"]))
        points = []
        for row in rows:
            x = left + (int(row["bin"]) - 0.5) / 120 * plot_w
            y = top + plot_h - float(row[metric]) / 100 * plot_h
            points.append(f"{x:.1f},{y:.1f}")
        body.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{colors[gene]}" stroke-width="2"/>')
    for tick in [0, 25, 50, 75, 100]:
        x = left + tick / 100 * plot_w
        body.append(f'<line x1="{x}" y1="{top + plot_h + 4}" x2="{x}" y2="{top + plot_h + 10}" class="axis"/>')
        body.append(f'<text x="{x - 10}" y="{top + plot_h + 27}" class="small">{tick}%</text>')
    legend_x = left + plot_w - 135
    for i, gene in enumerate(["BRCA1", "BRCA2"]):
        y = top + 20 + i * 22
        body.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 22}" y2="{y}" stroke="{colors[gene]}" stroke-width="3"/>')
        body.append(f'<text x="{legend_x + 30}" y="{y + 4}" class="small">{gene}</text>')
    body.append(f'<text x="{left + 350}" y="{height - 22}" class="label">Relative CDS position</text>')
    write_svg(PLOT_DIR / filename, "\n".join(body), width, height)


def draw_difference_heatmap(diff_rows: list[dict[str, object]]) -> None:
    metrics = ["vus_percent", "pathogenic_percent", "high_spliceai_percent", "missense_percent", "nonsense_percent", "synonymous_percent"]
    labels = ["VUS", "Pathogenic", "SpliceAI >=0.20", "Missense", "Nonsense", "Synonymous"]
    width = 1030
    height = 330
    left = 135
    top = 70
    cell_w = 7
    row_h = 34
    body = ['<text x="28" y="34" class="title">BRCA1 minus BRCA2 normalized bin differences</text>']
    body.append(f'<text x="{left}" y="55" class="small">Blue means higher in BRCA1, red means higher in BRCA2</text>')
    for yi, (metric, label) in enumerate(zip(metrics, labels)):
        y = top + yi * row_h
        body.append(f'<text x="22" y="{y + 20}" class="small">{esc(label)}</text>')
        for xi, row in enumerate(diff_rows):
            value = float(row[f"diff_{metric}"])
            magnitude = min(abs(value), 80.0)
            base = (37, 99, 235) if value >= 0 else (220, 38, 38)
            color = heat_color(magnitude, 80.0, base)
            x = left + xi * cell_w
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h - 4}" fill="{color}"/>')
    write_svg(PLOT_DIR / "normalized_bin_difference_heatmap.svg", "\n".join(body), width, height)


def draw_domain_signal_heatmap(domain_rows: list[dict[str, object]]) -> None:
    metrics = [
        ("benign 1/2", "benign_percent"),
        ("VUS 3", "vus_percent"),
        ("pathogenic 4/5", "pathogenic_percent"),
        ("SpAI >=0.20", "high_spliceai_percent"),
        ("PVS1", "pvs1_percent"),
        ("PP3", "pp3_percent"),
        ("PS3", "ps3_percent"),
        ("BS3", "bs3_percent"),
        ("BP1", "bp1_percent"),
    ]
    rows = [row for row in domain_rows if row["domain_context"] != "unknown protein position"]
    cell_w = 88
    cell_h = 30
    left = 255
    top = 112
    width = left + cell_w * len(metrics) + 145
    height = top + cell_h * len(rows) + 95
    max_value = max(float(row[col]) for row in rows for _, col in metrics) if rows else 1.0
    body = ['<text x="28" y="34" class="title">Domain and region context, Module 1 signals</text>']
    body.append('<text x="28" y="56" class="small">Rows are mapped BRCA1/2 domains or broad functional regions. Cells show percent within region.</text>')
    body.append(f'<rect x="{left - 2}" y="{top - 2}" width="{cell_w * len(metrics) + 2}" height="{cell_h * len(rows) + 2}" class="panel"/>')
    for xi, (label, _) in enumerate(metrics):
        x = left + xi * cell_w
        body.append(f'<text x="{x + 4}" y="{top - 14}" class="small">{esc(label)}</text>')
    for yi, row in enumerate(rows):
        y = top + yi * cell_h
        label = f"{row['gene']} {row['domain_context']}"
        body.append(f'<text x="20" y="{y + 20}" class="small">{esc(label[:38])}</text>')
        for xi, (_, col) in enumerate(metrics):
            x = left + xi * cell_w
            value = float(row[col])
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{heat_color(value, max_value, (37, 99, 235))}" class="tile"/>')
            if value >= 10:
                body.append(f'<text x="{x + 5}" y="{y + 19}" class="small">{value:.0f}</text>')
    legend_x = left + cell_w * len(metrics) + 32
    legend_y = top
    body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">percent</text>')
    for step in range(6):
        frac = 1 - step / 5
        y = legend_y + step * 18
        body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{heat_color(frac * max_value, max_value, (37, 99, 235))}" class="tile"/>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 12}" class="small">{max_value:.0f}</text>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 102}" class="small">0</text>')
    write_svg(PLOT_DIR / "domain_region_signal_heatmap.svg", "\n".join(body), width, height)


def draw_domain_normalized_density_heatmap(domain_rows: list[dict[str, object]]) -> None:
    metrics = [
        ("benign 1/2", "benign_per_100aa"),
        ("VUS 3", "vus_per_100aa"),
        ("pathogenic 4/5", "pathogenic_per_100aa"),
        ("SpAI >=0.20", "high_spliceai_per_100aa"),
        ("PVS1", "pvs1_per_100aa"),
        ("PP3", "pp3_per_100aa"),
        ("PS3", "ps3_per_100aa"),
        ("BS3", "bs3_per_100aa"),
    ]
    rows = [row for row in domain_rows if row["domain_context"] != "unknown protein position"]
    cell_w = 92
    cell_h = 30
    left = 270
    top = 116
    width = left + cell_w * len(metrics) + 145
    height = top + cell_h * len(rows) + 95
    max_value = max(float(row[col]) for row in rows for _, col in metrics) if rows else 1.0
    body = ['<text x="28" y="34" class="title">BRCA1 domains vs BRCA2 domains, normalized density</text>']
    body.append('<text x="28" y="56" class="small">Cells show generated Module 1 signal counts per 100 amino acids in each mapped domain or broad region.</text>')
    body.append('<text x="28" y="74" class="small">This complements within-domain percentages and reduces the effect of different domain lengths.</text>')
    for xi, (label, _) in enumerate(metrics):
        x = left + xi * cell_w
        body.append(f'<text x="{x + 4}" y="{top - 14}" class="small">{esc(label)}</text>')
    for yi, row in enumerate(rows):
        y = top + yi * cell_h
        label = f"{row['gene']} {row['domain_context']}"
        body.append(f'<text x="20" y="{y + 20}" class="small">{esc(label[:40])}</text>')
        for xi, (_, col) in enumerate(metrics):
            x = left + xi * cell_w
            value = float(row[col])
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{heat_color(value, max_value, (21, 128, 61))}" class="tile"/>')
            if value >= 5:
                body.append(f'<text x="{x + 5}" y="{y + 19}" class="small">{value:.0f}</text>')
    legend_x = left + cell_w * len(metrics) + 32
    legend_y = top
    body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">per 100 aa</text>')
    for step in range(6):
        frac = 1 - step / 5
        y = legend_y + step * 18
        body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{heat_color(frac * max_value, max_value, (21, 128, 61))}" class="tile"/>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 12}" class="small">{max_value:.0f}</text>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 102}" class="small">0</text>')
    write_svg(PLOT_DIR / "domain_vs_domain_normalized_density_heatmap.svg", "\n".join(body), width, height)


def draw_domain_gene_split_percent_heatmap(domain_rows: list[dict[str, object]]) -> None:
    metrics = [
        ("benign 1/2", "benign_percent"),
        ("VUS 3", "vus_percent"),
        ("pathogenic 4/5", "pathogenic_percent"),
        ("SpAI >=0.20", "high_spliceai_percent"),
        ("PP3", "pp3_percent"),
        ("PS3", "ps3_percent"),
        ("BS3", "bs3_percent"),
        ("BP1", "bp1_percent"),
    ]
    rows = [row for row in domain_rows if row["domain_context"] != "unknown protein position"]
    rows.sort(key=lambda row: (row["gene"], row["domain_context"] == "outside mapped domains", row["domain_context"]))
    cell_w = 90
    cell_h = 30
    left = 280
    top = 118
    width = left + cell_w * len(metrics) + 135
    height = top + cell_h * len(rows) + 95
    max_value = max(float(row[col]) for row in rows for _, col in metrics) if rows else 1.0
    body = ['<text x="28" y="34" class="title">BRCA1 domains vs BRCA2 domains, within-domain percent</text>']
    body.append('<text x="28" y="56" class="small">Cells show percent of generated SNVs within each domain or region. This is normalized within that domain.</text>')
    body.append('<text x="28" y="74" class="small">Rows are descriptive contexts, not homologous regions and not independent ACMG/ENIGMA evidence.</text>')
    current_gene = None
    for xi, (label, _) in enumerate(metrics):
        x = left + xi * cell_w
        body.append(f'<text x="{x + 4}" y="{top - 14}" class="small">{esc(label)}</text>')
    for yi, row in enumerate(rows):
        y = top + yi * cell_h
        if row["gene"] != current_gene:
            current_gene = row["gene"]
            body.append(f'<text x="20" y="{y - 4}" class="legend-title">{current_gene}</text>')
        label = f"{row['domain_context']} ({row['domain_range']})"
        body.append(f'<text x="58" y="{y + 20}" class="small">{esc(label[:38])}</text>')
        for xi, (_, col) in enumerate(metrics):
            x = left + xi * cell_w
            value = float(row[col])
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{heat_color(value, max_value, (126, 87, 194))}" class="tile"/>')
            if value >= 10:
                body.append(f'<text x="{x + 5}" y="{y + 19}" class="small">{value:.0f}</text>')
    legend_x = left + cell_w * len(metrics) + 32
    legend_y = top
    body.append(f'<text x="{legend_x}" y="{legend_y - 10}" class="legend-title">percent</text>')
    for step in range(6):
        frac = 1 - step / 5
        y = legend_y + step * 18
        body.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{heat_color(frac * max_value, max_value, (126, 87, 194))}" class="tile"/>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 12}" class="small">{max_value:.0f}</text>')
    body.append(f'<text x="{legend_x + 24}" y="{legend_y + 102}" class="small">0</text>')
    write_svg(PLOT_DIR / "domain_vs_domain_within_domain_percent_heatmap.svg", "\n".join(body), width, height)


def draw_domain_grouped_class_barplot(domain_rows: list[dict[str, object]]) -> None:
    rows = [row for row in domain_rows if row["domain_context"] != "unknown protein position"]
    width = 1100
    left = 300
    top = 84
    bar_h = 24
    gap = 12
    plot_w = 620
    height = top + len(rows) * (bar_h + gap) + 82
    body = ['<text x="28" y="34" class="title">Generated class mix by mapped domain or region</text>']
    body.append('<text x="28" y="56" class="small">Stacked percent within region. Domains are not ACMG/ENIGMA evidence by themselves.</text>')
    for idx, row in enumerate(rows):
        y = top + idx * (bar_h + gap)
        label = f"{row['gene']} {row['domain_context']}"
        body.append(f'<text x="24" y="{y + 17}" class="small">{esc(label[:42])}</text>')
        x = left
        for group in GROUP_ORDER:
            value = float(row[f"{group}_percent"])
            w = plot_w * value / 100
            body.append(f'<rect x="{x}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{GROUP_COLORS[group]}"/>')
            x += w
    legend_y = height - 42
    for idx, group in enumerate(GROUP_ORDER):
        x = left + idx * 185
        body.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{GROUP_COLORS[group]}"/>')
        body.append(f'<text x="{x + 20}" y="{legend_y + 12}" class="small">{esc(GROUP_LABELS[group])}</text>')
    write_svg(PLOT_DIR / "domain_region_grouped_class_mix.svg", "\n".join(body), width, height)


def top_differences(diff_rows: list[dict[str, object]], metric: str, limit: int = 15) -> list[dict[str, object]]:
    ranked = sorted(diff_rows, key=lambda row: abs(float(row[f"diff_{metric}"])), reverse=True)
    return ranked[:limit]


def report_table(rows: list[dict[str, object]], columns: list[str], limit: int = 12) -> str:
    if not rows:
        return "| none |\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows[:limit]:
        body.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join([header, sep] + body)


def write_report(summary: list[dict[str, object]], domain_rows: list[dict[str, object]], vus_examples: list[dict[str, object]], diff_rows: list[dict[str, object]]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    brca1, brca2 = summary
    text = f"""# Domain And Region Context Comparison

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

This analysis replaces a direct BRCA1-versus-BRCA2 whole-gene comparison with a
domain and broad-region view. Direct whole-gene comparison is limited because
BRCA1 and BRCA2 are different proteins with different lengths, domain
architecture, exon structure, and structured evidence coverage.

The domain view is still descriptive. Domain membership is not an ACMG/ENIGMA
criterion by itself. It helps explain where generated Module 1 signals are
concentrated and where benign counterexamples exist.

## Outputs

- `tables/gene_comparison/gene_summary.csv`
- `tables/gene_comparison/domain_region_summary.csv`
- `tables/gene_comparison/domain_region_vus_examples.csv`
- `tables/gene_comparison/normalized_position_bins.csv`
- `tables/gene_comparison/brca1_vs_brca2_bin_differences.csv`
- `plots/07_gene_comparison/domain_region_signal_heatmap.svg`
- `plots/07_gene_comparison/domain_region_grouped_class_mix.svg`
- `plots/07_gene_comparison/domain_vs_domain_within_domain_percent_heatmap.svg`
- `plots/07_gene_comparison/domain_vs_domain_normalized_density_heatmap.svg`

## Whole-Gene Summary

| Metric | BRCA1 | BRCA2 |
| --- | ---: | ---: |
| Variants | {brca1["n_variants"]} | {brca2["n_variants"]} |
| Benign percent | {brca1["benign_percent"]}% | {brca2["benign_percent"]}% |
| VUS percent | {brca1["vus_percent"]}% | {brca2["vus_percent"]}% |
| Pathogenic percent | {brca1["pathogenic_percent"]}% | {brca2["pathogenic_percent"]}% |
| High SpliceAI percent | {brca1["high_spliceai_percent"]}% | {brca2["high_spliceai_percent"]}% |
| Missense percent | {brca1["missense_percent"]}% | {brca2["missense_percent"]}% |
| Nonsense percent | {brca1["nonsense_percent"]}% | {brca2["nonsense_percent"]}% |
| Synonymous percent | {brca1["synonymous_percent"]}% | {brca2["synonymous_percent"]}% |

## Domain And Region Summary

{report_table(domain_rows, ["gene", "domain_context", "domain_range", "aa_length", "n_variants", "variants_per_100aa", "benign_percent", "vus_percent", "pathogenic_percent", "high_spliceai_percent", "ps3_percent", "bs3_percent"], limit=20)}

## BRCA1 Domains Versus BRCA2 Domains

The main comparison is domain-to-domain and region-to-region, not whole-gene
BRCA1 versus whole-gene BRCA2. Two normalizations are useful:

1. Within-domain percent asks: among generated SNVs in this domain, what
   fraction is benign, VUS, pathogenic, high-SpliceAI, PS3, BS3, and so on?
2. Density per 100 amino acids asks: how many generated signals occur per
   protein length unit in that domain?

The first view is best for class mix. The second view is useful when comparing
small regions such as BRCA1 RING or BRCA2 PALB2-binding context with larger
regions such as BRCA2 BRC repeats or the DNA-binding domain.

## VUS Examples Inside Mapped Domains

These are examples with a nonzero heuristic signal inside mapped domains. They
are not clinical recommendations and do not imply pathogenicity.

{report_table(vus_examples, ["gene", "domain_context", "score", "reasons", "c_notation", "p_notation", "spliceai_score"], limit=20)}

## Legacy Relative-CDS Bin Comparison

The older normalized-bin comparison is retained as a descriptive quality check,
but it should not be interpreted as homologous domain comparison. Positive
difference means BRCA1 higher; negative difference means BRCA2 higher at the
same relative CDS percentage.

{report_table(top_differences(diff_rows, "vus_percent"), ["bin", "relative_start_percent", "relative_end_percent", "brca1_vus_percent", "brca2_vus_percent", "diff_vus_percent"], limit=8)}

## Reading Notes

The most useful interpretation is within-region rather than gene-versus-gene.
For example, BRCA1 BRCT and BRCA2 DNA-binding/BRC regions are not equivalent
domains, but both can be inspected as functionally important contexts with
different Module 1 signal patterns. Domain context should always be combined
with exact variant consequence, accepted ACMG/ENIGMA evidence, functional/RNA
metadata, and benign counterexamples.
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    summary = gene_summary(rows)
    domains = domain_summary(rows)
    vus_examples = domain_vus_examples(rows)
    bins = normalized_bins(rows)
    diffs = bin_differences(bins)
    write_csv(
        TABLE_DIR / "gene_summary.csv",
        summary,
        ["gene", "n_variants", "cds_length"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent"]]
        + [f"{vtype}_{suffix}" for vtype in VARIANT_TYPES for suffix in ["count", "percent"]]
        + ["high_spliceai_count", "high_spliceai_percent"]
        + [f"{criterion.lower()}_{suffix}" for criterion in CRITERIA for suffix in ["count", "percent"]],
    )
    write_csv(
        TABLE_DIR / "domain_region_summary.csv",
        domains,
        ["gene", "domain_context", "domain_range", "aa_length", "n_variants", "variants_per_100aa"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent", "per_100aa"]]
        + ["high_spliceai_count", "high_spliceai_percent", "high_spliceai_per_100aa"]
        + [f"{criterion.lower()}_{suffix}" for criterion in ["PVS1", "PM5_PTC", "PP3", "PS3", "BS3", "BP1"] for suffix in ["count", "percent", "per_100aa"]],
    )
    write_csv(
        TABLE_DIR / "domain_region_vus_examples.csv",
        vus_examples,
        ["gene", "domain_context", "score", "reasons", "c_notation", "p_notation", "spliceai_score", "criteria"],
    )
    write_csv(
        TABLE_DIR / "normalized_position_bins.csv",
        bins,
        ["gene", "bin", "relative_start_percent", "relative_end_percent", "cds_start", "cds_end", "n_variants"]
        + [f"{group}_{suffix}" for group in GROUP_ORDER for suffix in ["count", "percent"]]
        + ["high_spliceai_count", "high_spliceai_percent"]
        + [f"{vtype}_percent" for vtype in VARIANT_TYPES],
    )
    write_csv(
        TABLE_DIR / "brca1_vs_brca2_bin_differences.csv",
        diffs,
        ["bin", "relative_start_percent", "relative_end_percent"]
        + [
            item
            for metric in ["benign_percent", "vus_percent", "pathogenic_percent", "high_spliceai_percent", "missense_percent", "nonsense_percent", "synonymous_percent"]
            for item in [f"brca1_{metric}", f"brca2_{metric}", f"diff_{metric}"]
        ],
    )
    draw_gene_group_barplot(summary)
    draw_gene_variant_type_barplot(summary)
    draw_domain_signal_heatmap(domains)
    draw_domain_gene_split_percent_heatmap(domains)
    draw_domain_normalized_density_heatmap(domains)
    draw_domain_grouped_class_barplot(domains)
    draw_normalized_profile(bins, "vus_percent", "Normalized VUS profile", "normalized_vus_profile.svg")
    draw_normalized_profile(bins, "pathogenic_percent", "Normalized pathogenic profile", "normalized_pathogenic_profile.svg")
    draw_normalized_profile(bins, "high_spliceai_percent", "Normalized high SpliceAI profile", "normalized_spliceai_profile.svg")
    draw_difference_heatmap(diffs)
    write_report(summary, domains, vus_examples, diffs)
    print(f"Wrote domain/region comparison to {REPORT}")


if __name__ == "__main__":
    main()
