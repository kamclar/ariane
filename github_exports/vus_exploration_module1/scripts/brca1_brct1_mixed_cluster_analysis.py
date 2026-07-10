"""Focused analysis of the BRCA1 BRCT1 mixed VUS neighborhood."""

from __future__ import annotations

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
VUS_NEIGHBORHOODS = ANALYSIS_DIR / "tables" / "vus_pathogenic_regions" / "vus_pathogenic_region_annotated.csv"
BOUNDARY = ANALYSIS_DIR / "tables" / "boundary_enrichment" / "variant_boundary_annotations.csv"
TABLE_DIR = ANALYSIS_DIR / "tables" / "brca1_brct1_mixed_cluster"
PLOT_DIR = ANALYSIS_DIR / "plots" / "23_brca1_brct1_mixed_cluster"
REPORT = ANALYSIS_DIR / "brca1_brct1_mixed_cluster_report.md"

GENE = "BRCA1"
CDS_START = 5044
CDS_END = 5143
AA_START = 1682
AA_END = 1715


def parse_float(value: str | None) -> float:
    if not value:
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def parse_criteria(value: str) -> set[str]:
    return {item.split(":", 1)[0] for item in value.split(";") if item}


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


def in_cluster(row: dict) -> bool:
    if row.get("gene") != GENE:
        return False
    try:
        pos = int(row.get("cds_pos") or row.get("c_pos") or row.get("cds_pos_int") or "")
    except ValueError:
        return False
    return CDS_START <= pos <= CDS_END


def load_full_cluster() -> list[dict]:
    rows = []
    with INPUT.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            c_notation = row["c_notation"]
            try:
                cds_pos = int(c_notation.split(".", 1)[1].split(">", 1)[0][:-1])
            except (IndexError, ValueError):
                continue
            if row["gene"] == GENE and CDS_START <= cds_pos <= CDS_END:
                row["cds_pos"] = cds_pos
                row["aa_region"] = f"p.{AA_START}-p.{AA_END}"
                row["class_group"] = class_group(row["predicted_class"])
                row["spliceai_float"] = parse_float(row.get("spliceai_score"))
                row["criteria_codes"] = "+".join(sorted(parse_criteria(row.get("criteria", "")))) or "none"
                row["domain_context"] = "BRCA1 BRCT phosphopeptide-binding region; UniProt BRCT 1 subdomain"
                rows.append(row)
    return rows


def load_vus_cluster() -> list[dict]:
    rows = []
    with VUS_NEIGHBORHOODS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if in_cluster(row):
                rows.append(row)
    return rows


def load_boundary_cluster() -> list[dict]:
    if not BOUNDARY.exists():
        return []
    rows = []
    with BOUNDARY.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("gene") != GENE:
                continue
            try:
                cds_pos = int(row.get("cds_pos", ""))
            except ValueError:
                continue
            if CDS_START <= cds_pos <= CDS_END:
                rows.append(row)
    return rows


def count_rows(rows: list[dict], key: str) -> list[dict]:
    return [
        {"field": key, "value": value, "count": count}
        for value, count in Counter(row.get(key, "") for row in rows).most_common()
    ]


def criteria_counts(rows: list[dict]) -> list[dict]:
    counts = Counter()
    for row in rows:
        for code in parse_criteria(row.get("criteria", "")):
            counts[code] += 1
    return [{"criterion": key, "count": value} for key, value in counts.most_common()]


def splice_summary(rows: list[dict]) -> list[dict]:
    bins = Counter()
    for row in rows:
        score = parse_float(row.get("spliceai_score"))
        if score >= 0.50:
            bins[">=0.50"] += 1
        elif score >= 0.20:
            bins["0.20-0.49"] += 1
        elif score >= 0.10:
            bins["0.10-0.19"] += 1
        else:
            bins["<0.10"] += 1
    return [{"spliceai_bin": key, "count": bins[key]} for key in [">=0.50", "0.20-0.49", "0.10-0.19", "<0.10"]]


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


def class_profile_svg(rows: list[dict]) -> None:
    by_pos: dict[int, Counter] = {}
    for row in rows:
        pos = int(row["cds_pos"])
        by_pos.setdefault(pos, Counter())[row["class_group"]] += 1
    width = 1180
    height = 390
    left = 78
    top = 64
    plot_w = 980
    plot_h = 220
    max_count = max((sum(counter.values()) for counter in by_pos.values()), default=1)
    colors = {"benign": "#16a34a", "vus": "#64748b", "pathogenic": "#b91c1c"}
    body = [
        '<text x="28" y="30" class="title">BRCA1 BRCT1 mixed cluster: generated classes by CDS position</text>',
        '<text x="28" y="50" class="small">Region c.5044-c.5143 / p.1682-p.1715. Stacked bars show all possible generated SNVs at each CDS position.</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>',
    ]
    for tick in [CDS_START, 5074, CDS_END]:
        x = left + plot_w * (tick - CDS_START) / (CDS_END - CDS_START)
        body.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" stroke="#cbd5e1"/>')
        body.append(f'<text x="{x - 18:.1f}" y="{top + plot_h + 20}" class="small">c.{tick}</text>')
    boundary_x = left + plot_w * (5074 - CDS_START) / (CDS_END - CDS_START)
    body.append(f'<line x1="{boundary_x:.1f}" y1="{top - 8}" x2="{boundary_x:.1f}" y2="{top + plot_h}" stroke="#7c3aed" stroke-width="2" stroke-dasharray="5 4"/>')
    body.append(f'<text x="{boundary_x + 8:.1f}" y="{top - 14}" class="small" fill="#4c1d95">nearest donor-like boundary c.5074</text>')
    bar_w = max(2.0, plot_w / (CDS_END - CDS_START + 1) - 1)
    for pos in sorted(by_pos):
        x = left + plot_w * (pos - CDS_START) / (CDS_END - CDS_START + 1)
        y_bottom = top + plot_h
        for grp in ["benign", "vus", "pathogenic"]:
            count = by_pos[pos][grp]
            if not count:
                continue
            h = plot_h * count / max_count
            y_bottom -= h
            body.append(f'<rect x="{x:.1f}" y="{y_bottom:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{colors[grp]}" opacity="0.82"/>')
    legend_x = left + plot_w + 22
    legend_y = top + 20
    for i, grp in enumerate(["benign", "vus", "pathogenic"]):
        y = legend_y + i * 24
        body.append(f'<rect x="{legend_x}" y="{y - 11}" width="13" height="13" fill="{colors[grp]}"/>')
        body.append(f'<text x="{legend_x + 20}" y="{y}" class="small">{grp}</text>')
    body.append(f'<text x="{left}" y="{top + plot_h + 48}" class="label">CDS position</text>')
    body.append(f'<text x="24" y="{top + 20}" class="label">SNV count</text>')
    write_svg(PLOT_DIR / "brca1_brct1_mixed_cluster_class_profile.svg", "\n".join(body), width, height)


def markdown_table(rows: list[dict], columns: list[str], limit: int = 30) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def build_report(full_rows: list[dict], vus_rows: list[dict], boundary_rows: list[dict]) -> None:
    full_by_group = count_rows(full_rows, "class_group")
    full_by_class = count_rows(full_rows, "predicted_class")
    full_by_type = count_rows(full_rows, "variant_type")
    full_criteria = criteria_counts(full_rows)
    vus_criteria = criteria_counts(vus_rows)
    full_splice = splice_summary(full_rows)
    vus_splice = splice_summary(vus_rows)
    boundary_bins = count_rows(boundary_rows, "splice_boundary_distance_bin") if boundary_rows else []
    donor_acceptor = count_rows(boundary_rows, "nearest_site_type") if boundary_rows else []

    top_vus = sorted(
        vus_rows,
        key=lambda row: (
            int(row.get("benign_count_20bp", 0)),
            int(row.get("pathogenic_count_20bp", 0)),
            parse_float(row.get("spliceai_score")),
        ),
        reverse=True,
    )

    text = f"""# BRCA1 BRCT1 Mixed Local Neighborhood

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Why This Region Was Flagged

In the VUS local-neighborhood scatter plot, a visible green cluster appears in
the benign-enriched sector. Those dots are VUS with many neighboring generated
variants of both directions, but with more class 1/2 than class 4/5 neighbors
within +/-20 coding bp.

The cluster maps mainly to:

- gene: `BRCA1`
- coding range: `c.{CDS_START}-c.{CDS_END}`
- protein range: approximately `p.{AA_START}-p.{AA_END}`
- local feature context: BRCA1 BRCT phosphopeptide-binding region, inside the
  UniProt BRCT 1 subdomain

This is a prioritization/interpretation finding only. It does not reclassify
any VUS.

## What The Region Represents

This range is inside the BRCA1 C-terminal BRCT region. The generated landscape
is mixed because different substitutions in the same short region can trigger
different automated evidence patterns:

- many missense substitutions receive benign-direction evidence such as `BS3`
  or no strong pathogenic evidence
- some substitutions have `PP3` because of local SpliceAI signal
- some substitutions have `PS3`
- nearby nonsense substitutions receive truncation-related pathogenic evidence
  such as `PVS1` and `PM5_PTC`

So the region is not uniformly benign or pathogenic. It is a mixed BRCT1
neighborhood where the exact nucleotide/protein consequence matters.

## All Generated SNVs In The Region

Total generated SNVs in region: {len(full_rows)}

Grouped classes:

{markdown_table(full_by_group, ["field", "value", "count"], 10)}

Exact classes:

{markdown_table(full_by_class, ["field", "value", "count"], 10)}

Variant types:

{markdown_table(full_by_type, ["field", "value", "count"], 10)}

SpliceAI bins:

{markdown_table(full_splice, ["spliceai_bin", "count"], 10)}

Most frequent criteria:

{markdown_table(full_criteria, ["criterion", "count"], 20)}

## VUS In The Region

VUS in region: {len(vus_rows)}

VUS SpliceAI bins:

{markdown_table(vus_splice, ["spliceai_bin", "count"], 10)}

VUS criteria:

{markdown_table(vus_criteria, ["criterion", "count"], 20)}

Top VUS by local benign/pathogenic neighbor counts:

{markdown_table(top_vus, ["gene", "c_notation", "p_notation", "variant_type", "spliceai_score", "benign_count_20bp", "pathogenic_count_20bp", "criteria"], 25)}

## Boundary Context

Boundary annotation rows in this region: {len(boundary_rows)}

Nearest boundary bin summary:

{markdown_table(boundary_bins, ["field", "value", "count"], 10)}

Donor/acceptor-like context:

{markdown_table(donor_acceptor, ["field", "value", "count"], 10)}

## Interpretation

This cluster is interesting because it sits in an important functional BRCT
region but is locally benign-enriched in the generated landscape. The green
position in the scatter plot should therefore be interpreted as:

- not "these VUS are benign"
- but "this BRCT1 subregion has many generated benign/likely benign neighbors,
  despite also containing pathogenic-direction evidence nearby"

The practical manual-review question is whether each VUS is driven by:

- protein functional evidence (`PS3`/`BS3`)
- splice prediction (`PP3`/SpliceAI)
- nearby truncating consequences
- or a mixed rule pattern that needs curated functional assay metadata

## Outputs

- `tables/brca1_brct1_mixed_cluster/all_region_variants.csv`
- `tables/brca1_brct1_mixed_cluster/vus_region_variants.csv`
- `tables/brca1_brct1_mixed_cluster/boundary_region_variants.csv`
- `plots/23_brca1_brct1_mixed_cluster/brca1_brct1_mixed_cluster_class_profile.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    full_rows = load_full_cluster()
    vus_rows = load_vus_cluster()
    boundary_rows = load_boundary_cluster()
    write_csv(
        TABLE_DIR / "all_region_variants.csv",
        full_rows,
        [
            "gene",
            "c_notation",
            "p_notation",
            "cds_pos",
            "aa_region",
            "variant_type",
            "spliceai_score",
            "criteria",
            "total_points",
            "predicted_class",
            "predicted_label",
            "class_group",
            "criteria_codes",
            "domain_context",
        ],
    )
    write_csv(
        TABLE_DIR / "vus_region_variants.csv",
        vus_rows,
        [
            "gene",
            "c_notation",
            "p_notation",
            "cds_pos",
            "variant_type",
            "spliceai_score",
            "pathogenic_count_20bp",
            "benign_count_20bp",
            "pathogenic_count_50bp",
            "benign_count_50bp",
            "criteria",
            "total_points",
            "predicted_class",
        ],
    )
    if boundary_rows:
        write_csv(TABLE_DIR / "boundary_region_variants.csv", boundary_rows, list(boundary_rows[0]))
    class_profile_svg(full_rows)
    build_report(full_rows, vus_rows, boundary_rows)
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
