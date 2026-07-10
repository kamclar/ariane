"""VUS proximity to pathogenic versus benign regions."""

from __future__ import annotations

import csv
import math
import re
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
PATHOGENIC_HOTSPOTS = ANALYSIS_DIR / "tables" / "hotspots_pathogenic.csv"
TABLE_DIR = ANALYSIS_DIR / "tables" / "vus_pathogenic_regions"
PLOT_DIR = ANALYSIS_DIR / "plots" / "18_vus_pathogenic_regions"
REPORT = ANALYSIS_DIR / "vus_pathogenic_region_report.md"

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


def group(predicted_class: str) -> str:
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
            row["cds_pos_int"] = parse_cds_pos(row["c_notation"])
            row["group"] = group(row["predicted_class"])
            row["spliceai_float"] = parse_float(row.get("spliceai_score"))
            row["criteria_codes"] = parse_criteria(row.get("criteria", ""))
            if row["cds_pos_int"] is not None:
                rows.append(row)
    return rows


def load_pathogenic_hotspots() -> list[dict]:
    if not PATHOGENIC_HOTSPOTS.exists():
        return []
    rows = []
    with PATHOGENIC_HOTSPOTS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "gene": row["gene"],
                    "cds_start": int(row["cds_start"]),
                    "cds_end": int(row["cds_end"]),
                    "bin": row["bin"],
                    "enrichment": row["enrichment"],
                    "count": row["count"],
                }
            )
    return rows


def nearest_distance(pos: int, positions: list[int]) -> int | None:
    if not positions:
        return None
    return min(abs(pos - other) for other in positions)


def count_within(pos: int, positions: list[int], window: int) -> int:
    start = pos - window
    end = pos + window
    return sum(1 for other in positions if start <= other <= end)


def hotspot_hit(gene: str, pos: int, hotspots: list[dict]) -> dict | None:
    for row in hotspots:
        if row["gene"] == gene and row["cds_start"] <= pos <= row["cds_end"]:
            return row
    return None


def broad_pathogenic_context_label(row: dict, hotspot: dict | None) -> str:
    labels = []
    if hotspot:
        labels.append(f"pathogenic_hotspot_bin_{hotspot['bin']}")
    if int(row["pathogenic_count_20bp"]) >= 5:
        labels.append("pathogenic_dense_20bp")
    if int(row["pathogenic_count_50bp"]) >= 10:
        labels.append("pathogenic_dense_50bp")
    if int(row["pathogenic_count_100bp"]) >= 20:
        labels.append("pathogenic_dense_100bp")
    return ";".join(labels) if labels else "not_pathogenic_region_by_current_thresholds"


def strong_pathogenic_region_label(row: dict, hotspot: dict | None) -> str:
    labels = []
    path20 = int(row["pathogenic_count_20bp"])
    benign20 = int(row["benign_count_20bp"])
    path50 = int(row["pathogenic_count_50bp"])
    benign50 = int(row["benign_count_50bp"])
    if hotspot:
        labels.append(f"pathogenic_hotspot_bin_{hotspot['bin']}")
    if path20 >= 10 and path20 > benign20:
        labels.append("strong_pathogenic_dense_20bp")
    if path50 >= 20 and path50 >= math.ceil(1.5 * max(benign50, 1)):
        labels.append("strong_pathogenic_dense_50bp")
    return ";".join(labels) if labels else "not_strong_pathogenic_region"


def closeness_label(path_dist: int | None, benign_dist: int | None) -> str:
    if path_dist is None and benign_dist is None:
        return "no_pathogenic_or_benign_neighbor"
    if path_dist is None:
        return "closer_to_benign"
    if benign_dist is None:
        return "closer_to_pathogenic"
    if path_dist < benign_dist:
        return "closer_to_pathogenic"
    if benign_dist < path_dist:
        return "closer_to_benign"
    return "equal_distance"


def proximity_score(row: dict) -> int:
    score = 0
    score += int(row["pathogenic_count_20bp"]) * 5
    score += int(row["pathogenic_count_50bp"]) * 2
    score += int(row["pathogenic_count_100bp"])
    score -= int(row["benign_count_20bp"]) * 5
    score -= int(row["benign_count_50bp"]) * 2
    score -= int(row["benign_count_100bp"])
    if row["closer_to_group"] == "closer_to_pathogenic":
        score += 20
    if row["closer_to_group"] == "closer_to_benign":
        score -= 20
    if row["strong_pathogenic_region"] != "not_strong_pathogenic_region":
        score += 30
    if row["spliceai_float"] >= 0.20:
        score += 10
    if "PP3" in row["criteria_codes"]:
        score += 5
    if "PS3" in row["criteria_codes"] or "PS1" in row["criteria_codes"]:
        score += 8
    return score


def analyze() -> tuple[list[dict], list[dict], list[dict]]:
    rows = load_rows()
    hotspots = load_pathogenic_hotspots()
    positions: dict[tuple[str, str], list[int]] = {}
    for gene in ["BRCA1", "BRCA2"]:
        for grp in ["benign", "pathogenic", "vus"]:
            positions[(gene, grp)] = sorted(
                int(row["cds_pos_int"])
                for row in rows
                if row["gene"] == gene and row["group"] == grp
            )

    annotated = []
    for row in rows:
        if row["group"] != "vus":
            continue
        gene = row["gene"]
        pos = int(row["cds_pos_int"])
        path_positions = positions[(gene, "pathogenic")]
        benign_positions = positions[(gene, "benign")]
        path_dist = nearest_distance(pos, path_positions)
        benign_dist = nearest_distance(pos, benign_positions)
        hit = hotspot_hit(gene, pos, hotspots)
        out = {
            "gene": gene,
            "c_notation": row["c_notation"],
            "p_notation": row["p_notation"],
            "cds_pos": pos,
            "variant_type": row["variant_type"],
            "spliceai_score": row["spliceai_score"],
            "criteria": row["criteria"],
            "total_points": row["total_points"],
            "predicted_class": row["predicted_class"],
            "nearest_pathogenic_distance": "" if path_dist is None else path_dist,
            "nearest_benign_distance": "" if benign_dist is None else benign_dist,
            "pathogenic_count_20bp": count_within(pos, path_positions, 20),
            "benign_count_20bp": count_within(pos, benign_positions, 20),
            "pathogenic_count_50bp": count_within(pos, path_positions, 50),
            "benign_count_50bp": count_within(pos, benign_positions, 50),
            "pathogenic_count_100bp": count_within(pos, path_positions, 100),
            "benign_count_100bp": count_within(pos, benign_positions, 100),
            "pathogenic_hotspot_bin": "" if not hit else hit["bin"],
            "pathogenic_hotspot_enrichment": "" if not hit else hit["enrichment"],
            "pathogenic_hotspot_count": "" if not hit else hit["count"],
            "closer_to_group": closeness_label(path_dist, benign_dist),
            "spliceai_float": row["spliceai_float"],
            "criteria_codes": row["criteria_codes"],
        }
        out["broad_pathogenic_context"] = broad_pathogenic_context_label(out, hit)
        out["strong_pathogenic_region"] = strong_pathogenic_region_label(out, hit)
        out["pathogenic_proximity_score"] = proximity_score(out)
        annotated.append(out)

    annotated.sort(
        key=lambda item: (
            item["broad_pathogenic_context"] != "not_pathogenic_region_by_current_thresholds",
            item["strong_pathogenic_region"] != "not_strong_pathogenic_region",
            item["closer_to_group"] == "closer_to_pathogenic",
            int(item["pathogenic_proximity_score"]),
            int(item["pathogenic_count_20bp"]),
            float(item["spliceai_float"]),
        ),
        reverse=True,
    )
    summaries = build_summaries(annotated)
    return annotated, summaries, top_candidates(annotated)


def build_summaries(rows: list[dict]) -> list[dict]:
    counts = Counter(row["closer_to_group"] for row in rows)
    region_counts = Counter(
        "broad_pathogenic_context"
        if row["broad_pathogenic_context"] != "not_pathogenic_region_by_current_thresholds"
        else "no_broad_pathogenic_context"
        for row in rows
    )
    strong_counts = Counter(
        "strong_pathogenic_region"
        if row["strong_pathogenic_region"] != "not_strong_pathogenic_region"
        else "not_strong_pathogenic_region"
        for row in rows
    )
    output = []
    for label, count in counts.most_common():
        output.append({"summary_type": "closer_to_group", "label": label, "count": count})
    for label, count in region_counts.most_common():
        output.append({"summary_type": "broad_pathogenic_context", "label": label, "count": count})
    for label, count in strong_counts.most_common():
        output.append({"summary_type": "strong_pathogenic_region", "label": label, "count": count})
    output.append(
        {
            "summary_type": "intersection",
            "label": "strong_pathogenic_region_and_closer_to_pathogenic",
            "count": sum(
                1
                for row in rows
                if row["strong_pathogenic_region"] != "not_strong_pathogenic_region"
                and row["closer_to_group"] == "closer_to_pathogenic"
            ),
        }
    )
    return output


def top_candidates(rows: list[dict]) -> list[dict]:
    return [
        row
        for row in rows
        if row["strong_pathogenic_region"] != "not_strong_pathogenic_region"
        and row["closer_to_group"] in {"closer_to_pathogenic", "equal_distance"}
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


def draw_summary_plot(summary: list[dict]) -> None:
    rows = summary
    width = 980
    height = 70 + len(rows) * 34
    left = 360
    top = 48
    plot_w = 520
    max_count = max((int(row["count"]) for row in rows), default=1)
    body = ['<text x="28" y="30" class="title">VUS pathogenic-region and proximity summary</text>']
    for i, row in enumerate(rows):
        y = top + i * 34
        count = int(row["count"])
        w = plot_w * count / max_count
        color = "#b91c1c" if "pathogenic" in row["label"] else "#64748b"
        label = f"{row['summary_type']}: {row['label']}"
        body.append(f'<text x="28" y="{y + 18}" class="small">{esc(label)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w}" height="22" fill="{color}"/>')
        body.append(f'<text x="{left + w + 8}" y="{y + 17}" class="small">{count}</text>')
    write_svg(PLOT_DIR / "vus_pathogenic_region_summary.svg", "\n".join(body), width, height)


def draw_scatter(rows: list[dict], window: int) -> None:
    width = 1040
    height = 740
    left = 122
    top = 92
    plot_w = 650
    plot_h = 460
    path_key = f"pathogenic_count_{window}bp"
    benign_key = f"benign_count_{window}bp"
    max_count = max(
        max(int(row[path_key]), int(row[benign_key]))
        for row in rows
    )
    max_count = max(1, max_count)
    tick_step = max(1, math.ceil(max_count / 5))
    ticks = list(range(0, max_count + 1, tick_step))
    if ticks[-1] != max_count:
        ticks.append(max_count)
    body = [
        f'<text x="28" y="30" class="title">VUS local neighborhood, +/-{window} bp window</text>',
        f'<text x="28" y="54" class="small">Each dot is one VUS. Count nearby generated variants in a coding window centered on that VUS.</text>',
        f'<text x="28" y="72" class="small">Right = more class 1/2 neighbors. Up = more class 4/5 neighbors. Numbers on axes are counts, not genomic positions.</text>',
    ]
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    body.append(f'<polygon points="{left},{top} {left + plot_w},{top} {left},{top + plot_h}" fill="#fee2e2" opacity="0.50"/>')
    body.append(f'<polygon points="{left},{top + plot_h} {left + plot_w},{top} {left + plot_w},{top + plot_h}" fill="#dcfce7" opacity="0.55"/>')
    for tick in ticks:
        x = left + plot_w * tick / max_count
        y = top + plot_h - plot_h * tick / max_count
        body.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" stroke="#e2e8f0"/>')
        body.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e2e8f0"/>')
        body.append(f'<line x1="{x:.1f}" y1="{top + plot_h}" x2="{x:.1f}" y2="{top + plot_h + 6}" stroke="#64748b"/>')
        body.append(f'<line x1="{left - 6}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" stroke="#64748b"/>')
        body.append(f'<text x="{x - 4:.1f}" y="{top + plot_h + 22}" class="small">{tick}</text>')
        body.append(f'<text x="{left - 28}" y="{y + 4:.1f}" class="small">{tick}</text>')
    body.append(f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top}" stroke="#94a3b8" stroke-dasharray="5 5"/>')
    body.append(f'<rect x="{left + 16}" y="{top + 14}" width="162" height="22" rx="5" fill="#ffffff" opacity="0.82"/>')
    body.append(f'<text x="{left + 24}" y="{top + 30}" class="small" fill="#7f1d1d">pathogenic-enriched</text>')
    body.append(f'<rect x="{left + plot_w - 164}" y="{top + plot_h - 78}" width="148" height="22" rx="5" fill="#ffffff" opacity="0.82"/>')
    body.append(f'<text x="{left + plot_w - 156}" y="{top + plot_h - 62}" class="small" fill="#14532d">benign-enriched</text>')
    for row in rows:
        path_count = int(row[path_key])
        benign_count = int(row[benign_key])
        x = left + plot_w * int(row[benign_key]) / max_count
        y = top + plot_h - plot_h * int(row[path_key]) / max_count
        if path_count > benign_count:
            color = "#dc2626"
        elif benign_count > path_count:
            color = "#16a34a"
        else:
            color = "#64748b"
        radius = 4 if row["closer_to_group"] == "closer_to_pathogenic" else 2.5
        stroke = "#111827" if row["strong_pathogenic_region"] != "not_strong_pathogenic_region" else "none"
        stroke_w = 0.7 if stroke != "none" else 0
        body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius}" fill="{color}" stroke="{stroke}" stroke-width="{stroke_w}" opacity="0.58"/>')
    if window == 20:
        cluster = [
            row
            for row in rows
            if row["gene"] == "BRCA1" and 5044 <= int(row["cds_pos"]) <= 5143
        ]
        if cluster:
            xs = [left + plot_w * int(row[benign_key]) / max_count for row in cluster]
            ys = [top + plot_h - plot_h * int(row[path_key]) / max_count for row in cluster]
            min_x, max_x = min(xs), max(xs)
            min_y, max_y = min(ys), max(ys)
            pad = 14
            body.append(
                f'<rect x="{min_x - pad:.1f}" y="{min_y - pad:.1f}" '
                f'width="{max_x - min_x + 2 * pad:.1f}" height="{max_y - min_y + 2 * pad:.1f}" '
                'rx="12" fill="none" stroke="#7c3aed" stroke-width="2" stroke-dasharray="6 4"/>'
            )
            label_x = min(max_x + 22, left + plot_w - 240)
            label_y = max(min_y - 34, top + 58)
            body.append(
                f'<rect x="{label_x:.1f}" y="{label_y:.1f}" width="238" height="46" '
                'rx="7" fill="#ffffff" opacity="0.90" stroke="#7c3aed"/>'
            )
            body.append(f'<text x="{label_x + 9:.1f}" y="{label_y + 17:.1f}" class="small" fill="#4c1d95">BRCA1 BRCT1 mixed cluster</text>')
            body.append(f'<text x="{label_x + 9:.1f}" y="{label_y + 34:.1f}" class="small" fill="#4c1d95">c.5044-c.5143 / p.1682-p.1715</text>')
            body.append(
                f'<line x1="{label_x:.1f}" y1="{label_y + 40:.1f}" '
                f'x2="{(min_x + max_x) / 2:.1f}" y2="{min_y - 4:.1f}" '
                'stroke="#7c3aed" stroke-width="1.5"/>'
            )
    body.append(f'<text x="{left + 150}" y="{top + plot_h + 54}" class="label">Number of class 1/2 variants within +/-{window} bp around this VUS</text>')
    y_label_x = 34
    y_label_y = top + plot_h / 2 + 190
    body.append(f'<text x="{y_label_x}" y="{y_label_y:.1f}" class="label" transform="rotate(-90 {y_label_x} {y_label_y:.1f})">Number of class 4/5 variants within +/-{window} bp around this VUS</text>')
    legend_x = left + plot_w + 28
    legend_y = top + 18
    body.append(f'<text x="{legend_x}" y="{legend_y}" class="label">How to read it</text>')
    body.append(f'<circle cx="{legend_x + 8}" cy="{legend_y + 26}" r="4" fill="#dc2626" opacity="0.7"/>')
    body.append(f'<text x="{legend_x + 22}" y="{legend_y + 30}" class="small">VUS has more 4/5 neighbors</text>')
    body.append(f'<circle cx="{legend_x + 8}" cy="{legend_y + 50}" r="4" fill="#16a34a" opacity="0.7"/>')
    body.append(f'<text x="{legend_x + 22}" y="{legend_y + 54}" class="small">VUS has more 1/2 neighbors</text>')
    body.append(f'<circle cx="{legend_x + 8}" cy="{legend_y + 74}" r="4" fill="#64748b" opacity="0.7"/>')
    body.append(f'<text x="{legend_x + 22}" y="{legend_y + 78}" class="small">equal counts</text>')
    body.append(f'<circle cx="{legend_x + 8}" cy="{legend_y + 98}" r="4" fill="#ffffff" stroke="#111827" stroke-width="1"/>')
    body.append(f'<text x="{legend_x + 22}" y="{legend_y + 102}" class="small">black outline: strong context flag</text>')
    window_y = legend_y + 150
    window_x = legend_x + 18
    body.append(f'<text x="{legend_x}" y="{window_y - 18}" class="label">Window definition</text>')
    body.append(f'<line x1="{window_x}" y1="{window_y}" x2="{window_x + 180}" y2="{window_y}" stroke="#334155" stroke-width="2"/>')
    body.append(f'<line x1="{window_x}" y1="{window_y - 8}" x2="{window_x}" y2="{window_y + 8}" stroke="#334155"/>')
    body.append(f'<line x1="{window_x + 90}" y1="{window_y - 12}" x2="{window_x + 90}" y2="{window_y + 12}" stroke="#334155" stroke-width="2"/>')
    body.append(f'<line x1="{window_x + 180}" y1="{window_y - 8}" x2="{window_x + 180}" y2="{window_y + 8}" stroke="#334155"/>')
    body.append(f'<text x="{window_x - 12}" y="{window_y + 26}" class="small">-{window} bp</text>')
    body.append(f'<text x="{window_x + 78}" y="{window_y + 26}" class="small">VUS</text>')
    body.append(f'<text x="{window_x + 160}" y="{window_y + 26}" class="small">+{window} bp</text>')
    write_svg(PLOT_DIR / f"vus_pathogenic_vs_benign_neighborhood_{window}bp.svg", "\n".join(body), width, height)


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


def write_report(rows: list[dict], summary: list[dict], top: list[dict]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# VUS in Pathogenic Regions

Generated: {generated}

Input: `{INPUT.relative_to(ROOT)}`

## Purpose

This analysis asks which VUS lie in regions dominated by generated
pathogenic/likely pathogenic variants and which VUS are locally closer to the
pathogenic group than to the benign/likely benign group.

This is a prioritization and interpretation layer only. It does not change any
classification.

## Definitions

- pathogenic group: predicted classes 4 and 5
- benign group: predicted classes 1 and 2
- VUS group: predicted class 3
- in the local-neighborhood scatter plots, each dot is one VUS
- the VUS is the center of the window; `+/-20 bp` means 20 coding bases
  upstream and 20 coding bases downstream of that VUS position
- x-axis: number of benign-group variants within that window around the VUS
- y-axis: number of pathogenic-group variants within that window around the VUS
- broad pathogenic context: VUS is in a pathogenic hotspot bin, or has at least 5
  pathogenic-group variants within +/-20 coding bases, or at least 10 within
  +/-50 coding bases, or at least 20 within +/-100 coding bases
- strong pathogenic region: VUS is in a pathogenic hotspot bin, or the local
  pathogenic group clearly outnumbers the benign group by stricter thresholds
- closer to pathogenic: nearest class 4/5 variant is closer than nearest class
  1/2 variant in coding position

## Summary

{markdown_table(summary, ["summary_type", "label", "count"], 20)}

## Top VUS In Pathogenic Regions And Closer To Pathogenic

{markdown_table(top, ["pathogenic_proximity_score", "gene", "c_notation", "p_notation", "variant_type", "spliceai_score", "strong_pathogenic_region", "broad_pathogenic_context", "closer_to_group", "nearest_pathogenic_distance", "nearest_benign_distance", "pathogenic_count_20bp", "benign_count_20bp", "criteria"], 30)}

## Interpretation

These variants are useful manual-review candidates because their local generated
classification landscape is more pathogenic than benign. The signal can come
from dense pathogenic neighborhoods, pathogenic hotspot bins, or shorter
distance to class 4/5 variants than to class 1/2 variants.

This does not prove pathogenicity. It can reflect real biological constraint,
rule behavior in the automated snapshot, local splice effects, or clustering of
criteria such as PS3, PP3, PVS1, or PM5_PTC.

## Outputs

- `tables/vus_pathogenic_regions/vus_pathogenic_region_annotated.csv`
- `tables/vus_pathogenic_regions/vus_pathogenic_region_top_candidates.csv`
- `tables/vus_pathogenic_regions/vus_pathogenic_region_summary.csv`
- `plots/18_vus_pathogenic_regions/vus_pathogenic_region_summary.svg`
- `plots/18_vus_pathogenic_regions/vus_pathogenic_vs_benign_neighborhood_20bp.svg`
- `plots/18_vus_pathogenic_regions/vus_pathogenic_vs_benign_neighborhood_50bp.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    annotated, summary, top = analyze()
    fields = [
        "pathogenic_proximity_score",
        "gene",
        "c_notation",
        "p_notation",
        "cds_pos",
        "variant_type",
        "spliceai_score",
        "strong_pathogenic_region",
        "broad_pathogenic_context",
        "closer_to_group",
        "nearest_pathogenic_distance",
        "nearest_benign_distance",
        "pathogenic_count_20bp",
        "benign_count_20bp",
        "pathogenic_count_50bp",
        "benign_count_50bp",
        "pathogenic_count_100bp",
        "benign_count_100bp",
        "pathogenic_hotspot_bin",
        "pathogenic_hotspot_enrichment",
        "pathogenic_hotspot_count",
        "criteria",
        "total_points",
        "predicted_class",
    ]
    write_csv(TABLE_DIR / "vus_pathogenic_region_annotated.csv", annotated, fields)
    write_csv(TABLE_DIR / "vus_pathogenic_region_top_candidates.csv", top, fields)
    write_csv(TABLE_DIR / "vus_pathogenic_region_summary.csv", summary, ["summary_type", "label", "count"])
    draw_summary_plot(summary)
    draw_scatter(annotated, 20)
    draw_scatter(annotated, 50)
    write_report(annotated, summary, top)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
