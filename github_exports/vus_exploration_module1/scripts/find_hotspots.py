"""Find exploratory hotspot and coldspot regions in the BRCA coding SNV snapshot."""

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
TABLE_DIR = OUT_DIR / "tables"
PLOT_DIR = OUT_DIR / "plots" / "05_hotspots"
REPORT = OUT_DIR / "hotspot_coldspot_report.md"

CDS_LENGTHS = {"BRCA1": 5592, "BRCA2": 10257}
C_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")
SIGNALS = [
    "pathogenic",
    "vus",
    "high_spliceai",
    "pvs1",
    "pm5_ptc",
    "bs3",
    "bp1",
    "pp3",
]
SIGNAL_LABELS = {
    "pathogenic": "Pathogenic/Likely Pathogenic",
    "vus": "VUS",
    "high_spliceai": "SpliceAI >= 0.20",
    "pvs1": "PVS1",
    "pm5_ptc": "PM5_PTC",
    "bs3": "BS3",
    "bp1": "BP1",
    "pp3": "PP3",
}
SIGNAL_COLORS = {
    "pathogenic": (215, 48, 31),
    "vus": (138, 143, 156),
    "high_spliceai": (88, 80, 141),
    "pvs1": (234, 88, 12),
    "pm5_ptc": (245, 158, 11),
    "bs3": (43, 140, 190),
    "bp1": (22, 163, 74),
    "pp3": (124, 58, 237),
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


def group(predicted_class: str) -> str:
    if predicted_class in {"1", "2"}:
        return "benign"
    if predicted_class in {"4", "5"}:
        return "pathogenic"
    return "vus"


def load_rows() -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["cds_pos"] = parse_cds_pos(row["c_notation"])
            row["spliceai"] = parse_float(row["spliceai_score"])
            row["criteria_codes"] = parse_criteria(row.get("criteria", ""))
            row["group"] = group(row["predicted_class"])
            rows.append(row)
    return rows


def signal_present(row: dict[str, object], signal: str) -> bool:
    codes = row["criteria_codes"]
    if signal == "pathogenic":
        return row["group"] == "pathogenic"
    if signal == "vus":
        return row["group"] == "vus"
    if signal == "high_spliceai":
        return float(row["spliceai"]) >= 0.20
    if signal == "pvs1":
        return "PVS1" in codes
    if signal == "pm5_ptc":
        return "PM5_PTC" in codes
    if signal == "bs3":
        return "BS3" in codes
    if signal == "bp1":
        return "BP1" in codes
    if signal == "pp3":
        return "PP3" in codes
    raise ValueError(f"Unknown signal: {signal}")


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def poisson_tail(k: int, expected: float) -> float:
    if expected <= 0:
        return 0.0 if k > 0 else 1.0
    term = math.exp(-expected)
    cumulative = term
    for i in range(1, k):
        term *= expected / i
        cumulative += term
    return max(0.0, min(1.0, 1.0 - cumulative))


def build_bin_rows(rows: list[dict[str, object]], bins: int = 120) -> tuple[list[dict[str, object]], dict[str, int]]:
    totals = {signal: sum(1 for row in rows if signal_present(row, signal)) for signal in SIGNALS}
    total_variants = len(rows)
    bin_rows = []
    for gene in ["BRCA1", "BRCA2"]:
        length = CDS_LENGTHS[gene]
        gene_rows = [row for row in rows if row["gene"] == gene and isinstance(row.get("cds_pos"), int)]
        by_bin = [[] for _ in range(bins)]
        for row in gene_rows:
            idx = min(bins - 1, max(0, int((int(row["cds_pos"]) - 1) / length * bins)))
            by_bin[idx].append(row)
        for idx, items in enumerate(by_bin):
            start = int(idx * length / bins) + 1
            end = int((idx + 1) * length / bins)
            record = {
                "gene": gene,
                "bin": idx + 1,
                "cds_start": start,
                "cds_end": end,
                "n_variants": len(items),
            }
            for signal in SIGNALS:
                count = sum(1 for row in items if signal_present(row, signal))
                rate = count / len(items) if items else 0.0
                background = totals[signal] / total_variants if total_variants else 0.0
                expected = len(items) * background
                enrichment = rate / background if background else 0.0
                record[f"{signal}_count"] = count
                record[f"{signal}_rate"] = rate
                record[f"{signal}_expected"] = expected
                record[f"{signal}_enrichment"] = enrichment
                record[f"{signal}_poisson_tail"] = poisson_tail(count, expected)
            bin_rows.append(record)
    return bin_rows, totals


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def esc(value: object) -> str:
    text = str(value)
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def heat_color(value: float, max_value: float, base: tuple[int, int, int]) -> str:
    if max_value <= 0:
        intensity = 0.0
    else:
        intensity = min(1.0, math.sqrt(max(value, 0.0) / max_value))
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
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def fieldnames() -> list[str]:
    fields = ["gene", "bin", "cds_start", "cds_end", "n_variants"]
    for signal in SIGNALS:
        fields.extend(
            [
                f"{signal}_count",
                f"{signal}_rate",
                f"{signal}_expected",
                f"{signal}_enrichment",
                f"{signal}_poisson_tail",
            ]
        )
    return fields


def hotspot_rows(bin_rows: list[dict[str, object]], signal: str, min_count: int = 5) -> list[dict[str, object]]:
    selected = [
        row
        for row in bin_rows
        if int(row[f"{signal}_count"]) >= min_count and float(row[f"{signal}_enrichment"]) >= 2.0
    ]
    return sorted(
        selected,
        key=lambda row: (float(row[f"{signal}_enrichment"]), int(row[f"{signal}_count"])),
        reverse=True,
    )


def coldspot_rows(bin_rows: list[dict[str, object]], signal: str, min_expected: float = 5.0) -> list[dict[str, object]]:
    selected = [
        row
        for row in bin_rows
        if float(row[f"{signal}_expected"]) >= min_expected and int(row[f"{signal}_count"]) == 0
    ]
    return sorted(selected, key=lambda row: float(row[f"{signal}_expected"]), reverse=True)


def compact_rows(rows: list[dict[str, object]], signal: str, limit: int = 25) -> list[dict[str, object]]:
    compact = []
    for row in rows[:limit]:
        compact.append(
            {
                "gene": row["gene"],
                "bin": row["bin"],
                "cds_start": row["cds_start"],
                "cds_end": row["cds_end"],
                "n_variants": row["n_variants"],
                "signal": signal,
                "count": row[f"{signal}_count"],
                "rate": f"{float(row[f'{signal}_rate']) * 100:.2f}",
                "expected": f"{float(row[f'{signal}_expected']):.2f}",
                "enrichment": f"{float(row[f'{signal}_enrichment']):.2f}",
                "poisson_tail": f"{float(row[f'{signal}_poisson_tail']):.4g}",
            }
        )
    return compact


def write_signal_tables(bin_rows: list[dict[str, object]]) -> dict[str, dict[str, int]]:
    summary = {}
    compact_fields = [
        "gene",
        "bin",
        "cds_start",
        "cds_end",
        "n_variants",
        "signal",
        "count",
        "rate",
        "expected",
        "enrichment",
        "poisson_tail",
    ]
    for signal in SIGNALS:
        hot = hotspot_rows(bin_rows, signal)
        cold = coldspot_rows(bin_rows, signal)
        write_csv(TABLE_DIR / f"hotspots_{signal}.csv", compact_rows(hot, signal, limit=50), compact_fields)
        write_csv(TABLE_DIR / f"coldspots_{signal}.csv", compact_rows(cold, signal, limit=50), compact_fields)
        summary[signal] = {"hotspots": len(hot), "coldspots": len(cold)}
    return summary


def draw_signal_heatmap(bin_rows: list[dict[str, object]], signal: str) -> None:
    genes = ["BRCA1", "BRCA2"]
    bins = 120
    cell_w = 7
    row_h = 34
    left = 72
    top = 78
    plot_w = bins * cell_w
    width = left + plot_w + 55
    height = top + len(genes) * row_h + 92
    values = [
        float(row[f"{signal}_enrichment"])
        for row in bin_rows
        if int(row[f"{signal}_count"]) > 0
    ]
    max_value = min(max(values) if values else 1.0, 16.0)
    body = [f'<text x="28" y="34" class="title">{esc(SIGNAL_LABELS[signal])} hotspot enrichment</text>']
    body.append(f'<text x="{left}" y="58" class="small">Each gene is split into 120 CDS bins; darker means stronger enrichment, capped at 16x</text>')
    base = SIGNAL_COLORS[signal]
    for yi, gene in enumerate(genes):
        y = top + yi * row_h
        body.append(f'<text x="24" y="{y + 21}" class="label">{gene}</text>')
        gene_rows = sorted([row for row in bin_rows if row["gene"] == gene], key=lambda row: int(row["bin"]))
        for xi, row in enumerate(gene_rows):
            value = min(float(row[f"{signal}_enrichment"]), 16.0)
            count = int(row[f"{signal}_count"])
            color = "#f8fafc" if count == 0 else heat_color(value, max_value, base)
            x = left + xi * cell_w
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h - 5}" fill="{color}"/>')
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        x = left + frac * plot_w
        body.append(f'<line x1="{x}" y1="{top + len(genes) * row_h + 4}" x2="{x}" y2="{top + len(genes) * row_h + 10}" class="axis"/>')
        body.append(f'<text x="{x - 10}" y="{top + len(genes) * row_h + 27}" class="small">{round(frac * 100)}%</text>')
    write_svg(PLOT_DIR / f"hotspot_heatmap_{signal}.svg", "\n".join(body), width, height)


def draw_top_hotspot_bars(bin_rows: list[dict[str, object]], signal: str, limit: int = 12) -> None:
    rows = hotspot_rows(bin_rows, signal)[:limit]
    width = 980
    left = 220
    top = 66
    bar_h = 24
    gap = 10
    plot_w = 650
    height = top + len(rows) * (bar_h + gap) + 60
    body = [f'<text x="28" y="34" class="title">Top {esc(SIGNAL_LABELS[signal])} hotspots</text>']
    if not rows:
        body.append('<text x="28" y="80" class="label">No hotspots passed the current thresholds.</text>')
        write_svg(PLOT_DIR / f"top_hotspots_{signal}.svg", "\n".join(body), width, 140)
        return
    max_enrichment = max(float(row[f"{signal}_enrichment"]) for row in rows)
    base = SIGNAL_COLORS[signal]
    for i, row in enumerate(rows):
        y = top + i * (bar_h + gap)
        label = f"{row['gene']} bin {row['bin']} c.{row['cds_start']}-{row['cds_end']}"
        enrichment = float(row[f"{signal}_enrichment"])
        count = int(row[f"{signal}_count"])
        rate = float(row[f"{signal}_rate"]) * 100
        w = round(plot_w * enrichment / max_enrichment) if max_enrichment else 0
        body.append(f'<text x="26" y="{y + 17}" class="small">{esc(label)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="{bar_h}" fill="{heat_color(enrichment, max_enrichment, base)}"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 17}" class="small">{enrichment:.2f}x, n={count}, {rate:.1f}%</text>')
    write_svg(PLOT_DIR / f"top_hotspots_{signal}.svg", "\n".join(body), width, height)


def draw_combined_signal_heatmap(bin_rows: list[dict[str, object]]) -> None:
    genes = ["BRCA1", "BRCA2"]
    bins = 120
    cell_w = 7
    row_h = 25
    left = 150
    top = 78
    rows = [(gene, signal) for gene in genes for signal in SIGNALS]
    width = left + bins * cell_w + 55
    height = top + len(rows) * row_h + 82
    body = ['<text x="28" y="34" class="title">Hotspot enrichment overview</text>']
    body.append(f'<text x="{left}" y="58" class="small">Enrichment capped at 16x; empty bins are white</text>')
    for yi, (gene, signal) in enumerate(rows):
        y = top + yi * row_h
        body.append(f'<text x="18" y="{y + 16}" class="small">{gene} {esc(signal)}</text>')
        gene_rows = sorted([row for row in bin_rows if row["gene"] == gene], key=lambda row: int(row["bin"]))
        base = SIGNAL_COLORS[signal]
        for xi, row in enumerate(gene_rows):
            count = int(row[f"{signal}_count"])
            value = min(float(row[f"{signal}_enrichment"]), 16.0)
            color = "#f8fafc" if count == 0 else heat_color(value, 16.0, base)
            x = left + xi * cell_w
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h - 3}" fill="{color}"/>')
    write_svg(PLOT_DIR / "hotspot_heatmap_all_signals.svg", "\n".join(body), width, height)


def draw_hotspot_plots(bin_rows: list[dict[str, object]]) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    draw_combined_signal_heatmap(bin_rows)
    for signal in SIGNALS:
        draw_signal_heatmap(bin_rows, signal)
        draw_top_hotspot_bars(bin_rows, signal)


def report_table(rows: list[dict[str, object]]) -> str:
    if not rows:
        return "| none | | | | | | | | |\n"
    lines = []
    for row in rows[:10]:
        lines.append(
            f"| {row['gene']} | {row['bin']} | {row['cds_start']}-{row['cds_end']} | "
            f"{row['n_variants']} | {row['count']} | {row['rate']}% | "
            f"{row['expected']} | {row['enrichment']}x | {row['poisson_tail']} |"
        )
    return "\n".join(lines)


def write_report(bin_rows: list[dict[str, object]], totals: dict[str, int], summary: dict[str, dict[str, int]]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    sections = []
    for signal in ["pathogenic", "vus", "high_spliceai", "pvs1", "pm5_ptc", "bs3", "bp1", "pp3"]:
        hot = compact_rows(hotspot_rows(bin_rows, signal), signal, limit=10)
        cold = compact_rows(coldspot_rows(bin_rows, signal), signal, limit=10)
        sections.append(
            f"""## {signal}

Total signal count: {totals[signal]}

Hotspot bins found: {summary[signal]["hotspots"]}
Coldspot bins found: {summary[signal]["coldspots"]}

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{report_table(hot)}

Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{report_table(cold)}
"""
        )
    text = f"""# Hotspot and Coldspot Report

Generated: {generated}

Input: `{INPUT_CSV.relative_to(ROOT)}`

This is an exploratory bin-level scan over BRCA1 and BRCA2 coding SNVs. Each gene is split into 120 CDS bins. Hotspots are bins with at least 5 signal variants and at least 2x enrichment over the full coding-SNV background. Coldspots are bins with expected signal count at least 5 but observed count 0.

Signals:

- pathogenic: predicted classes 4 and 5
- vus: predicted class 3
- high_spliceai: SpliceAI >= 0.20
- pvs1, pm5_ptc, bs3, bp1, pp3: criterion present

Interpretation caveat: this is an automated Module 1 landscape, not final expert classification. Bin-level enrichment is descriptive and should be followed by variant-level review.

## Output Files

- `tables/hotspot_coldspot_bins_all_signals.csv`
- `tables/hotspots_<signal>.csv`
- `tables/coldspots_<signal>.csv`
- `plots/05_hotspots/hotspot_heatmap_all_signals.svg`
- `plots/05_hotspots/hotspot_heatmap_<signal>.svg`
- `plots/05_hotspots/top_hotspots_<signal>.svg`

{chr(10).join(sections)}
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    bin_rows, totals = build_bin_rows(rows)
    fields = fieldnames()
    formatted = []
    for row in bin_rows:
        out = dict(row)
        for signal in SIGNALS:
            for suffix in ["rate", "expected", "enrichment", "poisson_tail"]:
                key = f"{signal}_{suffix}"
                out[key] = f"{float(out[key]):.6g}"
        formatted.append(out)
    write_csv(TABLE_DIR / "hotspot_coldspot_bins_all_signals.csv", formatted, fields)
    summary = write_signal_tables(bin_rows)
    draw_hotspot_plots(bin_rows)
    write_report(bin_rows, totals, summary)
    print(f"Wrote hotspot report to {REPORT}")


if __name__ == "__main__":
    main()
