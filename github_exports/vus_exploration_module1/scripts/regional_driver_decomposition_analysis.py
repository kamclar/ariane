"""Regional driver decomposition for the precomputed Module 1 map."""

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
TABLE_DIR = ANALYSIS_DIR / "tables" / "regional_driver_decomposition"
PLOT_DIR = ANALYSIS_DIR / "plots" / "25_regional_driver_decomposition"
REPORT = ANALYSIS_DIR / "regional_driver_decomposition_report.md"

CDS_LENGTHS = {"BRCA1": 5592, "BRCA2": 10257}
BIN_SIZE = 100
C_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")
P_RE = re.compile(r"[A-Z][a-z]{2}(\d+)")

GROUPS = ["benign", "vus", "pathogenic"]
GROUP_COLORS = {"benign": "#16a34a", "vus": "#64748b", "pathogenic": "#b91c1c"}
DRIVERS = [
    "truncation_driver",
    "splice_driver",
    "pathogenic_functional_driver",
    "benign_functional_driver",
    "benign_rule_driver",
    "frequency_driver",
    "low_information",
]
DRIVER_LABELS = {
    "truncation_driver": "PTC/PVS1/PM5",
    "splice_driver": "PP3/SpliceAI",
    "pathogenic_functional_driver": "PS3/PS1/PP4",
    "benign_functional_driver": "BS3",
    "benign_rule_driver": "BP/BA/BS1",
    "frequency_driver": "PM2/BS1/BA1",
    "low_information": "PM2-only/none",
}
DRIVER_COLORS = {
    "truncation_driver": "#b91c1c",
    "splice_driver": "#f97316",
    "pathogenic_functional_driver": "#7c3aed",
    "benign_functional_driver": "#16a34a",
    "benign_rule_driver": "#0f766e",
    "frequency_driver": "#2563eb",
    "low_information": "#94a3b8",
}


STRUCTURE_FEATURES = {
    "BRCA1": [
        ("RING zinc-binding/E3 ligase region", 2, 101),
        ("coiled-coil PALB2 interaction region", 1391, 1424),
        ("BRCT phosphopeptide-binding region", 1650, 1857),
    ],
    "BRCA2": [
        ("PALB2 binding N-terminal region", 10, 40),
        ("BRC repeats / RAD51-binding region", 1000, 2080),
        ("DNA-binding domain / helical-OB-DSS1 region", 2481, 3186),
        ("C-terminal RAD51-binding/nuclear localization context", 3187, 3418),
    ],
}


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


def feature_hits(gene: str, aa_pos: int | None) -> list[str]:
    if aa_pos is None:
        return []
    return [
        name
        for name, start, end in STRUCTURE_FEATURES.get(gene, [])
        if start <= aa_pos <= end
    ]


def load_rows() -> list[dict]:
    boundaries = cds_boundaries()
    rows = []
    with INPUT.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pos = parse_cds_pos(row["c_notation"])
            if pos is None:
                continue
            criteria = parse_criteria(row.get("criteria", ""))
            score = parse_float(row.get("spliceai_score"))
            aa_pos = aa_position(row.get("p_notation", ""))
            row["cds_pos"] = pos
            row["aa_pos"] = aa_pos
            row["group"] = group(row["predicted_class"])
            row["criteria_codes"] = criteria
            row["spliceai_float"] = score
            row["boundary_distance"] = nearest_boundary_distance(row["gene"], pos, boundaries)
            row["feature_hits"] = feature_hits(row["gene"], aa_pos)
            row["truncation_driver"] = row["variant_type"] == "nonsense" or bool(criteria & {"PVS1", "PM5_PTC"})
            row["splice_driver"] = "PP3" in criteria or score >= 0.20
            row["pathogenic_functional_driver"] = bool(criteria & {"PS3", "PS1", "PP4"})
            row["benign_functional_driver"] = "BS3" in criteria
            row["benign_rule_driver"] = bool(criteria & {"BP1", "BP4", "BP5", "BP7", "BA1", "BS1_Strong", "BS1_Supporting"})
            row["frequency_driver"] = bool(criteria & {"PM2_Supporting", "BA1", "BS1_Strong", "BS1_Supporting"})
            row["low_information"] = criteria <= {"PM2_Supporting"} or not criteria
            rows.append(row)
    return rows


def bin_id(pos: int) -> int:
    return (pos - 1) // BIN_SIZE + 1


def bin_range(gene: str, idx: int) -> tuple[int, int]:
    start = (idx - 1) * BIN_SIZE + 1
    end = min(idx * BIN_SIZE, CDS_LENGTHS[gene])
    return start, end


def fraction(num: int, den: int) -> float:
    return num / den if den else 0.0


def pathogenic_fraction(rows: list[dict]) -> float:
    return fraction(sum(1 for row in rows if row["group"] == "pathogenic"), len(rows))


def classify_region(row: dict) -> str:
    path = int(row["pathogenic_count"])
    if path == 0:
        if int(row["benign_functional_count"]) or int(row["benign_rule_count"]):
            return "benign_evidence_dominant"
        return "no_pathogenic_signal"
    driver_fracs = {
        "truncation_driven": float(row["truncation_pathogenic_fraction"]),
        "splice_driven": float(row["splice_pathogenic_fraction"]),
        "functional_pathogenic_driven": float(row["functional_pathogenic_fraction"]),
    }
    best, best_value = max(driver_fracs.items(), key=lambda item: item[1])
    if best_value >= 0.50:
        return best
    if float(row["residual_pathogenic_fraction"]) >= 0.10 and path >= 5:
        return "residual_after_major_drivers"
    if int(row["mixed_position_count"]) or int(row["mixed_codon_count"]):
        return "mixed_position_or_codon"
    return "mixed_or_low_specificity"


def summarize_bins(rows: list[dict]) -> list[dict]:
    by_bin: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in rows:
        by_bin[(row["gene"], bin_id(int(row["cds_pos"])))].append(row)

    position_groups: dict[tuple[str, int], set[str]] = defaultdict(set)
    codon_groups: dict[tuple[str, int], set[str]] = defaultdict(set)
    for row in rows:
        position_groups[(row["gene"], int(row["cds_pos"]))].add(row["group"])
        if row["aa_pos"]:
            codon_groups[(row["gene"], int(row["aa_pos"]))].add(row["group"])

    output = []
    for gene in ["BRCA1", "BRCA2"]:
        max_bin = math.ceil(CDS_LENGTHS[gene] / BIN_SIZE)
        for idx in range(1, max_bin + 1):
            bucket = by_bin.get((gene, idx), [])
            start, end = bin_range(gene, idx)
            groups = Counter(row["group"] for row in bucket)
            pathogenic = [row for row in bucket if row["group"] == "pathogenic"]
            pathogenic_no_trunc = [row for row in bucket if row["group"] == "pathogenic" and not row["truncation_driver"]]
            no_trunc = [row for row in bucket if not row["truncation_driver"]]
            no_splice = [row for row in bucket if not row["splice_driver"]]
            no_functional = [row for row in bucket if not row["pathogenic_functional_driver"] and not row["benign_functional_driver"]]
            residual = [
                row for row in bucket
                if not row["truncation_driver"]
                and not row["splice_driver"]
                and not row["pathogenic_functional_driver"]
                and not row["benign_functional_driver"]
            ]
            feature_counts = Counter(feature for row in bucket for feature in row["feature_hits"])
            mixed_positions = sum(
                1
                for pos in range(start, end + 1)
                if len(position_groups.get((gene, pos), set())) > 1
            )
            aa_start = (start - 1) // 3 + 1
            aa_end = (end - 1) // 3 + 1
            mixed_codons = sum(
                1
                for aa in range(aa_start, aa_end + 1)
                if len(codon_groups.get((gene, aa), set())) > 1
            )
            row = {
                "gene": gene,
                "bin": idx,
                "cds_start": start,
                "cds_end": end,
                "aa_start": aa_start,
                "aa_end": aa_end,
                "total_count": len(bucket),
                "benign_count": groups["benign"],
                "vus_count": groups["vus"],
                "pathogenic_count": groups["pathogenic"],
                "pathogenic_fraction": f"{pathogenic_fraction(bucket):.4f}",
                "pathogenic_fraction_no_truncation": f"{pathogenic_fraction(no_trunc):.4f}",
                "pathogenic_fraction_no_splice": f"{pathogenic_fraction(no_splice):.4f}",
                "pathogenic_fraction_no_functional": f"{pathogenic_fraction(no_functional):.4f}",
                "residual_pathogenic_fraction": f"{pathogenic_fraction(residual):.4f}",
                "truncation_pathogenic_count": sum(1 for item in pathogenic if item["truncation_driver"]),
                "splice_pathogenic_count": sum(1 for item in pathogenic if item["splice_driver"]),
                "functional_pathogenic_count": sum(1 for item in pathogenic if item["pathogenic_functional_driver"]),
                "benign_functional_count": sum(1 for item in bucket if item["benign_functional_driver"]),
                "benign_rule_count": sum(1 for item in bucket if item["benign_rule_driver"]),
                "high_spliceai_count": sum(1 for item in bucket if item["spliceai_float"] >= 0.20),
                "pp3_count": sum(1 for item in bucket if "PP3" in item["criteria_codes"]),
                "pvs1_pm5_count": sum(1 for item in bucket if item["truncation_driver"]),
                "ps3_ps1_pp4_count": sum(1 for item in bucket if item["pathogenic_functional_driver"]),
                "bs3_count": sum(1 for item in bucket if item["benign_functional_driver"]),
                "bp_benign_rule_count": sum(1 for item in bucket if item["benign_rule_driver"]),
                "pm2_only_or_none_count": sum(1 for item in bucket if item["low_information"]),
                "near_boundary_count": sum(1 for item in bucket if item["boundary_distance"] <= 10),
                "mixed_position_count": mixed_positions,
                "mixed_codon_count": mixed_codons,
                "top_structure_feature": feature_counts.most_common(1)[0][0] if feature_counts else "none",
            }
            row["truncation_pathogenic_fraction"] = f"{fraction(int(row['truncation_pathogenic_count']), len(pathogenic)):.4f}"
            row["splice_pathogenic_fraction"] = f"{fraction(int(row['splice_pathogenic_count']), len(pathogenic)):.4f}"
            row["functional_pathogenic_fraction"] = f"{fraction(int(row['functional_pathogenic_count']), len(pathogenic)):.4f}"
            row["region_category"] = classify_region(row)
            output.append(row)
    return output


def top_regions(rows: list[dict]) -> list[dict]:
    return sorted(
        rows,
        key=lambda row: (
            float(row["pathogenic_fraction"]),
            int(row["pathogenic_count"]),
            int(row["total_count"]),
        ),
        reverse=True,
    )[:40]


def driver_summary(rows: list[dict]) -> list[dict]:
    counts = Counter(row["region_category"] for row in rows)
    return [{"region_category": key, "count": value} for key, value in counts.most_common()]


def ablation_summary(rows: list[dict]) -> list[dict]:
    output = []
    for gene in ["BRCA1", "BRCA2", "both"]:
        selected = rows if gene == "both" else [row for row in rows if row["gene"] == gene]
        for metric in [
            "pathogenic_fraction",
            "pathogenic_fraction_no_truncation",
            "pathogenic_fraction_no_splice",
            "pathogenic_fraction_no_functional",
            "residual_pathogenic_fraction",
        ]:
            values = [float(row[metric]) for row in selected if int(row["total_count"]) > 0]
            output.append(
                {
                    "gene": gene,
                    "metric": metric,
                    "mean": f"{(sum(values) / len(values)):.4f}" if values else "0.0000",
                    "max": f"{max(values):.4f}" if values else "0.0000",
                    "bins_ge_0_10": sum(1 for value in values if value >= 0.10),
                    "bins_ge_0_25": sum(1 for value in values if value >= 0.25),
                }
            )
    return output


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


def heat_color(value: float, max_value: float, base: tuple[int, int, int]) -> str:
    intensity = 0.0 if max_value <= 0 else min(1.0, math.sqrt(max(value, 0.0) / max_value))
    bg = (248, 250, 252)
    rgb = tuple(round(bg[i] * (1 - intensity) + base[i] * intensity) for i in range(3))
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def driver_track_svg(rows: list[dict], gene: str) -> None:
    gene_rows = [row for row in rows if row["gene"] == gene]
    metrics = [
        ("pathogenic_fraction", "pathogenic fraction", (185, 28, 28)),
        ("truncation_pathogenic_fraction", "PTC/PVS1/PM5 among path", (185, 28, 28)),
        ("splice_pathogenic_fraction", "PP3/SpAI among path", (249, 115, 22)),
        ("functional_pathogenic_fraction", "PS3/PS1/PP4 among path", (124, 58, 237)),
        ("bs3_count", "BS3 count", (22, 163, 74)),
        ("mixed_codon_count", "mixed codons", (37, 99, 235)),
        ("residual_pathogenic_fraction", "residual path fraction", (15, 118, 110)),
    ]
    cell_w = 8
    row_h = 28
    left = 185
    top = 70
    plot_w = cell_w * len(gene_rows)
    width = left + plot_w + 40
    height = top + row_h * len(metrics) + 90
    body = [f'<text x="28" y="32" class="title">{gene} regional driver decomposition, {BIN_SIZE} bp bins</text>']
    body.append(f'<text x="{left}" y="52" class="small">Darker cells indicate stronger regional signal in the precomputed Module 1 map.</text>')
    for mi, (metric, label, base) in enumerate(metrics):
        y = top + mi * row_h
        values = [float(row[metric]) if "fraction" in metric else float(row[metric]) for row in gene_rows]
        max_value = max(values) if values else 1.0
        body.append(f'<text x="18" y="{y + 18}" class="small">{esc(label)}</text>')
        for xi, value in enumerate(values):
            x = left + xi * cell_w
            body.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h - 3}" fill="{heat_color(value, max_value, base)}"/>')
    for feature, start, end in STRUCTURE_FEATURES.get(gene, []):
        cds_start = (start - 1) * 3 + 1
        cds_end = end * 3
        x = left + max(0, (cds_start - 1) / BIN_SIZE) * cell_w
        w = max(3, (cds_end - cds_start + 1) / BIN_SIZE * cell_w)
        y = top + row_h * len(metrics) + 16
        body.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="10" fill="#ddd6fe" stroke="#7c3aed"/>')
        body.append(f'<text x="{x:.1f}" y="{y + 24}" class="small">{esc(feature.split("/")[0])}</text>')
    write_svg(PLOT_DIR / f"{gene.lower()}_regional_driver_tracks.svg", "\n".join(body), width, height)


def category_bar_svg(summary: list[dict]) -> None:
    width = 960
    height = 80 + len(summary) * 34
    left = 300
    top = 48
    plot_w = 520
    max_count = max((int(row["count"]) for row in summary), default=1)
    body = ['<text x="28" y="30" class="title">Regional signal categories</text>']
    for i, row in enumerate(summary):
        y = top + i * 34
        count = int(row["count"])
        w = plot_w * count / max_count
        body.append(f'<text x="28" y="{y + 17}" class="small">{esc(row["region_category"])}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w:.1f}" height="22" fill="#7c3aed" opacity="0.82"/>')
        body.append(f'<text x="{left + w + 8:.1f}" y="{y + 17}" class="small">{count}</text>')
    write_svg(PLOT_DIR / "regional_signal_categories.svg", "\n".join(body), width, height)


def ablation_svg(rows: list[dict]) -> None:
    metrics = [
        ("pathogenic_fraction", "all"),
        ("pathogenic_fraction_no_truncation", "no PTC"),
        ("pathogenic_fraction_no_splice", "no PP3/SpAI"),
        ("pathogenic_fraction_no_functional", "no PS3/BS3"),
        ("residual_pathogenic_fraction", "residual"),
    ]
    genes = ["BRCA1", "BRCA2"]
    width = 900
    height = 410
    left = 80
    top = 70
    panel_w = 340
    panel_h = 250
    gap = 70
    body = ['<text x="28" y="32" class="title">Mean regional pathogenic fraction after driver exclusions</text>']
    for gi, gene in enumerate(genes):
        selected = [row for row in rows if row["gene"] == gene and int(row["total_count"]) > 0]
        x0 = left + gi * (panel_w + gap)
        body.append(f'<text x="{x0}" y="{top - 14}" class="label">{gene}</text>')
        body.append(f'<rect x="{x0}" y="{top}" width="{panel_w}" height="{panel_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
        means = []
        for metric, _ in metrics:
            values = [float(row[metric]) for row in selected]
            means.append(sum(values) / len(values) if values else 0.0)
        max_value = max(means + [0.01])
        bar_w = panel_w / len(metrics) - 12
        for i, ((_, label), value) in enumerate(zip(metrics, means)):
            h = panel_h * value / max_value
            x = x0 + 8 + i * (bar_w + 12)
            y = top + panel_h - h
            body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="#b91c1c" opacity="0.78"/>')
            body.append(f'<text x="{x:.1f}" y="{top + panel_h + 18}" class="small">{esc(label)}</text>')
            body.append(f'<text x="{x:.1f}" y="{y - 5:.1f}" class="tiny">{value:.2f}</text>')
    write_svg(PLOT_DIR / "pathogenic_fraction_driver_ablation.svg", "\n".join(body), width, height)


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


def write_report(rows: list[dict], categories: list[dict], ablation: list[dict], top: list[dict]) -> None:
    text = f"""# Regional Driver Decomposition

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Purpose

This analysis asks what drives local pathogenic enrichment in the precomputed
ARIANE Module 1 coding SNV map. It is not a classifier and it does not add new
ACMG/ENIGMA evidence. It decomposes regional signal into rule/mechanism
components so that local context can be interpreted more safely.

Core interpretation:

> Local pathogenic enrichment is a context signal, not an ACMG/ENIGMA evidence
> criterion.

## Method

The CDS of BRCA1 and BRCA2 was split into {BIN_SIZE} bp bins. For each bin, the
analysis counted benign/VUS/pathogenic generated classes and several driver
signals:

- truncation driver: nonsense, `PVS1`, or `PM5_PTC`
- splice driver: `PP3` or SpliceAI >= 0.20
- pathogenic functional driver: `PS3`, `PS1`, or `PP4`
- benign functional driver: `BS3`
- benign rule driver: `BP1`, `BP4`, `BP5`, `BP7`, `BA1`, or `BS1`
- frequency/absence driver: `PM2`, `BS1`, or `BA1`
- low-information variants: `PM2_Supporting` only or no criteria

The analysis also recalculated regional pathogenic fraction after excluding
major drivers:

- all variants
- excluding truncation drivers
- excluding splice drivers
- excluding functional evidence drivers
- excluding truncation, splice, and functional evidence together

## Regional Categories

{markdown_table(categories, ["region_category", "count"], 20)}

## Driver Ablation Summary

{markdown_table(ablation, ["gene", "metric", "mean", "max", "bins_ge_0_10", "bins_ge_0_25"], 30)}

## Top Pathogenic-Enriched Bins

{markdown_table(top, ["gene", "bin", "cds_start", "cds_end", "aa_start", "aa_end", "pathogenic_fraction", "pathogenic_count", "total_count", "region_category", "top_structure_feature", "truncation_pathogenic_fraction", "splice_pathogenic_fraction", "functional_pathogenic_fraction", "residual_pathogenic_fraction", "mixed_codon_count"], 30)}

## Interpretation

Regional pathogenic enrichment is not one biological phenomenon. In this
precomputed Module 1 map it can arise from several different sources:

- truncating possibilities and `PVS1`/`PM5_PTC`
- splice prediction and `PP3`
- functional evidence such as `PS3` or benign evidence such as `BS3`
- benign rule structure such as `BP1`, `BP4`, or `BP7`
- mixed codons and positions where different substitutions activate different
  criteria

The practical use is triage: a region with high pathogenic enrichment should
tell the curator what to inspect first, not how to classify the VUS.

## Outputs

- `tables/regional_driver_decomposition/regional_driver_bins.csv`
- `tables/regional_driver_decomposition/regional_signal_categories.csv`
- `tables/regional_driver_decomposition/regional_driver_ablation_summary.csv`
- `tables/regional_driver_decomposition/top_pathogenic_enriched_bins.csv`
- `plots/25_regional_driver_decomposition/brca1_regional_driver_tracks.svg`
- `plots/25_regional_driver_decomposition/brca2_regional_driver_tracks.svg`
- `plots/25_regional_driver_decomposition/regional_signal_categories.svg`
- `plots/25_regional_driver_decomposition/pathogenic_fraction_driver_ablation.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = load_rows()
    regional = summarize_bins(rows)
    categories = driver_summary(regional)
    ablation = ablation_summary(regional)
    top = top_regions(regional)
    fields = [
        "gene",
        "bin",
        "cds_start",
        "cds_end",
        "aa_start",
        "aa_end",
        "total_count",
        "benign_count",
        "vus_count",
        "pathogenic_count",
        "pathogenic_fraction",
        "pathogenic_fraction_no_truncation",
        "pathogenic_fraction_no_splice",
        "pathogenic_fraction_no_functional",
        "residual_pathogenic_fraction",
        "truncation_pathogenic_count",
        "splice_pathogenic_count",
        "functional_pathogenic_count",
        "benign_functional_count",
        "benign_rule_count",
        "high_spliceai_count",
        "pp3_count",
        "pvs1_pm5_count",
        "ps3_ps1_pp4_count",
        "bs3_count",
        "bp_benign_rule_count",
        "pm2_only_or_none_count",
        "near_boundary_count",
        "mixed_position_count",
        "mixed_codon_count",
        "top_structure_feature",
        "truncation_pathogenic_fraction",
        "splice_pathogenic_fraction",
        "functional_pathogenic_fraction",
        "region_category",
    ]
    write_csv(TABLE_DIR / "regional_driver_bins.csv", regional, fields)
    write_csv(TABLE_DIR / "regional_signal_categories.csv", categories, ["region_category", "count"])
    write_csv(TABLE_DIR / "regional_driver_ablation_summary.csv", ablation, ["gene", "metric", "mean", "max", "bins_ge_0_10", "bins_ge_0_25"])
    write_csv(TABLE_DIR / "top_pathogenic_enriched_bins.csv", top, fields)
    driver_track_svg(regional, "BRCA1")
    driver_track_svg(regional, "BRCA2")
    category_bar_svg(categories)
    ablation_svg(regional)
    write_report(regional, categories, ablation, top)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
