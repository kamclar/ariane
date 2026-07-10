from __future__ import annotations

import csv
import html
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
RESIDUES = ROOT / "backend" / "data" / "clinically_important_residues.json"
TABLE_DIR = ANALYSIS_DIR / "tables" / "pathogenic_mechanisms"
PLOT_DIR = ANALYSIS_DIR / "plots" / "16_pathogenic_mechanisms"
REPORT = ANALYSIS_DIR / "pathogenic_mechanism_report.md"


PATHOGENIC_CLASSES = {"4", "5"}
TRUNCATING_TYPES = {"nonsense", "frameshift", "splice_site"}


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_float(value: str) -> Optional[float]:
    try:
        if value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def aa_position(p_notation: str) -> Optional[int]:
    match = re.search(r"[A-Z][a-z]{2}(\d+)", p_notation or "")
    if match:
        return int(match.group(1))
    return None


def parse_criteria(value: str) -> Dict[str, dict]:
    criteria: Dict[str, dict] = {}
    if not value:
        return criteria
    for item in value.split(";"):
        parts = item.split(":")
        if len(parts) < 3:
            continue
        code = parts[0].strip()
        if code:
            criteria[code] = {
                "strength": parts[1].strip(),
                "points": parse_int(parts[2].strip()),
            }
    return criteria


def criteria_combo(criteria: Dict[str, dict]) -> str:
    if not criteria:
        return "none"
    return "+".join(sorted(criteria))


def load_rows() -> List[dict]:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_domain_data() -> tuple[dict, dict]:
    payload = json.loads(RESIDUES.read_text(encoding="utf-8"))
    domains = {}
    pathogenic_positions = {}
    for gene, gene_domains in payload.get("domains", {}).items():
        domains[gene] = {}
        pathogenic_positions[gene] = {}
        for name, data in gene_domains.items():
            start, end = data["range"]
            domains[gene][name] = (int(start), int(end))
            for residue in data.get("pathogenic_residues", []):
                pathogenic_positions[gene].setdefault(int(residue["pos"]), []).append(
                    {
                        "domain": name,
                        "aa": residue.get("aa"),
                        "c": residue.get("c"),
                    }
                )
    return domains, pathogenic_positions


def domain_for(gene: str, pos: Optional[int], domains: dict) -> str:
    if pos is None:
        return "unknown_position"
    for name, (start, end) in domains.get(gene, {}).items():
        if start <= pos <= end:
            return name
    return "outside_enigma_functional_domain"


def pathogenic_residue_hit(gene: str, pos: Optional[int], pathogenic_positions: dict) -> str:
    if pos is None:
        return ""
    hits = pathogenic_positions.get(gene, {}).get(pos, [])
    if not hits:
        return ""
    return ";".join(f"{hit['domain']}:{hit['aa']}:{hit['c']}" for hit in hits)


def mechanism_for(row: dict, criteria: Dict[str, dict], domain: str, residue_hit: str) -> str:
    variant_type = row["variant_type"]
    spliceai = parse_float(row.get("spliceai_score", ""))
    has_pvs1 = "PVS1" in criteria
    has_pm5_ptc = "PM5_PTC" in criteria
    has_ps3 = "PS3" in criteria
    has_pp4 = "PP4" in criteria
    has_ps1 = "PS1" in criteria
    has_pp3 = "PP3" in criteria

    if has_pvs1 and has_pm5_ptc:
        return "truncating_pvs1_pm5_ptc"
    if has_pvs1:
        return "truncating_pvs1_only"
    if has_pm5_ptc:
        return "truncating_pm5_ptc_only"
    if has_ps3 and has_pp4:
        return "functional_multifactorial_pathogenic"
    if has_ps3 and has_pp3:
        return "functional_plus_splice_or_computational"
    if has_ps3:
        return "functional_ps3_only_or_with_pm2"
    if has_pp4:
        return "multifactorial_pp4"
    if has_ps1:
        return "same_amino_acid_ps1"
    if has_pp3 and spliceai is not None and spliceai >= 0.20:
        return "spliceai_pp3_high_signal"
    if residue_hit:
        return "known_pathogenic_residue_context"
    if domain not in ("outside_enigma_functional_domain", "unknown_position"):
        return "enigma_functional_domain_context"
    if variant_type in TRUNCATING_TYPES:
        return "truncating_without_pvs1"
    return "other_pathogenic_mechanism"


def aberrant_protein_label(row: dict, criteria: Dict[str, dict]) -> str:
    variant_type = row["variant_type"]
    if "PVS1" in criteria or "PM5_PTC" in criteria:
        if variant_type == "nonsense":
            return "premature_stop_or_truncating_protein"
        if variant_type == "frameshift":
            return "frameshift_truncating_protein"
        if variant_type == "splice_site":
            return "splice_related_loss_of_function"
        return "loss_of_function_or_truncating_rule"
    if variant_type == "missense":
        return "missense_protein_function_or_domain"
    if variant_type == "synonymous":
        return "synonymous_splice_or_functional_evidence"
    if variant_type == "initiation_codon":
        return "translation_initiation_effect"
    return "other_or_unclear"


def write_csv(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def bar_svg(path: Path, rows: List[dict], label_key: str, value_key: str, title: str, max_items: int = 20) -> None:
    rows = rows[:max_items]
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 1050
    left = 390
    top = 55
    row_h = 30
    chart_w = 560
    height = top + row_h * len(rows) + 45
    max_value = max((int(row[value_key]) for row in rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:13px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".small{font-size:11px;fill:#555}",
        ".bar{fill:#b24a3a}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = int(row[value_key])
        bar_w = 0 if max_value == 0 else int(chart_w * value / max_value)
        parts.append(f'<text x="24" y="{y + 18}">{html.escape(str(row[label_key]))}</text>')
        parts.append(f'<rect x="{left}" y="{y + 5}" width="{bar_w}" height="18" class="bar"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 19}" class="small">{value}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def heatmap_svg(
    path: Path,
    rows: List[dict],
    row_key: str,
    col_key: str,
    value_key: str,
    title: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row_labels = sorted({str(row[row_key]) for row in rows})
    col_labels = sorted({str(row[col_key]) for row in rows})
    matrix = {(str(row[row_key]), str(row[col_key])): int(row[value_key]) for row in rows}
    max_value = max(matrix.values(), default=1)
    cell_w = 125
    cell_h = 34
    left = 250
    top = 85
    width = left + cell_w * len(col_labels) + 80
    height = top + cell_h * len(row_labels) + 80

    def color(value: int) -> str:
        if value <= 0:
            return "#f7f7f7"
        frac = value / max_value
        # light to dark red
        r = 255 - int(75 * frac)
        g = 238 - int(155 * frac)
        b = 230 - int(170 * frac)
        return f"#{r:02x}{g:02x}{b:02x}"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:12px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".small{font-size:10px;fill:#555}",
        ".cell{stroke:#fff;stroke-width:1}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
        f'<text x="24" y="52" class="small">Counts of generated ARIANE class 4/5 variants, not independently validated biology.</text>',
    ]
    for j, col in enumerate(col_labels):
        x = left + j * cell_w + cell_w / 2
        parts.append(
            f'<text x="{x}" y="{top - 12}" text-anchor="middle">{html.escape(col)}</text>'
        )
    for i, row_label in enumerate(row_labels):
        y = top + i * cell_h
        parts.append(f'<text x="24" y="{y + 22}">{html.escape(row_label)}</text>')
        for j, col in enumerate(col_labels):
            x = left + j * cell_w
            value = matrix.get((row_label, col), 0)
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" class="cell" fill="{color(value)}"/>'
            )
            if value:
                parts.append(
                    f'<text x="{x + cell_w / 2}" y="{y + 22}" text-anchor="middle">{value}</text>'
                )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def markdown_table(rows: List[dict], columns: List[str], limit: int = 25) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(html.escape(str(row.get(col, ""))) for col in columns) + " |")
    return "\n".join(lines)


def analyze() -> tuple[list[dict], dict[str, list[dict]]]:
    domains, pathogenic_positions = load_domain_data()
    rows = [row for row in load_rows() if row["predicted_class"] in PATHOGENIC_CLASSES]
    annotated = []
    counters = {
        "mechanism": Counter(),
        "aberrant": Counter(),
        "variant_type": Counter(),
        "domain": Counter(),
        "gene_domain": Counter(),
        "criteria_combo": Counter(),
        "residue_hits": Counter(),
        "missense_domain": Counter(),
        "class_mechanism": Counter(),
        "gene_mechanism": Counter(),
        "missense_domain_gene": Counter(),
    }

    for row in rows:
        criteria = parse_criteria(row["criteria"])
        pos = aa_position(row["p_notation"])
        domain = domain_for(row["gene"], pos, domains)
        residue_hit = pathogenic_residue_hit(row["gene"], pos, pathogenic_positions)
        mechanism = mechanism_for(row, criteria, domain, residue_hit)
        aberrant = aberrant_protein_label(row, criteria)
        combo = criteria_combo(criteria)
        annotated_row = {
            "gene": row["gene"],
            "c_notation": row["c_notation"],
            "p_notation": row["p_notation"],
            "aa_position": "" if pos is None else pos,
            "variant_type": row["variant_type"],
            "predicted_class": row["predicted_class"],
            "predicted_label": row["predicted_label"],
            "total_points": row["total_points"],
            "mechanism_bucket": mechanism,
            "aberrant_protein_bucket": aberrant,
            "enigma_domain": domain,
            "known_pathogenic_residue_context": residue_hit,
            "spliceai_score": row["spliceai_score"],
            "criteria_combo": combo,
            "criteria": row["criteria"],
        }
        annotated.append(annotated_row)
        counters["mechanism"][mechanism] += 1
        counters["aberrant"][aberrant] += 1
        counters["variant_type"][row["variant_type"]] += 1
        counters["domain"][domain] += 1
        counters["gene_domain"][(row["gene"], domain)] += 1
        counters["criteria_combo"][combo] += 1
        counters["class_mechanism"][(row["predicted_class"], mechanism)] += 1
        counters["gene_mechanism"][(row["gene"], mechanism)] += 1
        if row["variant_type"] == "missense":
            counters["missense_domain"][(row["gene"], domain, bool(residue_hit))] += 1
            counters["missense_domain_gene"][(row["gene"], domain)] += 1
        if residue_hit:
            counters["residue_hits"][(row["gene"], row["p_notation"], residue_hit)] += 1

    tables = {
        "mechanism_summary": [
            {"mechanism_bucket": k, "count": v}
            for k, v in counters["mechanism"].most_common()
        ],
        "aberrant_protein_summary": [
            {"aberrant_protein_bucket": k, "count": v}
            for k, v in counters["aberrant"].most_common()
        ],
        "variant_type_summary": [
            {"variant_type": k, "count": v}
            for k, v in counters["variant_type"].most_common()
        ],
        "domain_summary": [
            {"enigma_domain": k, "count": v}
            for k, v in counters["domain"].most_common()
        ],
        "gene_domain_summary": [
            {"gene": k[0], "enigma_domain": k[1], "count": v}
            for k, v in counters["gene_domain"].most_common()
        ],
        "criteria_combo_summary": [
            {"criteria_combo": k, "count": v}
            for k, v in counters["criteria_combo"].most_common()
        ],
        "known_pathogenic_residue_context": [
            {"gene": k[0], "p_notation": k[1], "known_context": k[2], "count": v}
            for k, v in counters["residue_hits"].most_common()
        ],
        "missense_domain_summary": [
            {
                "gene": k[0],
                "enigma_domain": k[1],
                "known_pathogenic_residue_position": str(k[2]),
                "count": v,
            }
            for k, v in counters["missense_domain"].most_common()
        ],
        "class_mechanism_summary": [
            {"predicted_class": k[0], "mechanism_bucket": k[1], "count": v}
            for k, v in counters["class_mechanism"].most_common()
        ],
        "gene_mechanism_summary": [
            {"gene": k[0], "mechanism_bucket": k[1], "count": v}
            for k, v in counters["gene_mechanism"].most_common()
        ],
        "missense_domain_gene_summary": [
            {"gene": k[0], "enigma_domain": k[1], "count": v}
            for k, v in counters["missense_domain_gene"].most_common()
        ],
    }
    return annotated, tables


def build_report(annotated: list[dict], tables: dict[str, list[dict]]) -> None:
    total = len(annotated)
    truncating = sum(1 for row in annotated if row["aberrant_protein_bucket"] in {
        "premature_stop_or_truncating_protein",
        "frameshift_truncating_protein",
        "splice_related_loss_of_function",
        "loss_of_function_or_truncating_rule",
    })
    domain_context = sum(1 for row in annotated if row["enigma_domain"] not in ("outside_enigma_functional_domain", "unknown_position"))
    residue_context = sum(1 for row in annotated if row["known_pathogenic_residue_context"])
    missense = sum(1 for row in annotated if row["variant_type"] == "missense")

    text = f"""# Pathogenic Mechanism Exploration

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Why We Did This

This analysis looks only at generated class 4 and class 5 variants and asks what
kind of biological or rule-based mechanism appears to drive pathogenicity in the
current ARIANE Module 1 coding SNV snapshot.

It separates:

- truncating or loss-of-function-like outcomes, mostly `PVS1` and `PM5_PTC`
- functional evidence, mostly `PS3` and `PP4`
- same amino-acid evidence, `PS1`
- computational splice or protein evidence, mostly `PP3`
- location in ENIGMA-defined functional domains
- overlap with ENIGMA clinically important pathogenic residue positions

This is descriptive. It does not prove molecular mechanism for each variant.

When using this analysis in continuous text, keep two layers separate. One layer
is the analysis of ARIANE's generated classification and its criteria. The other
layer is the biological-looking pattern in the generated data landscape. The
mechanism buckets below describe what the current automated classification
appears to be driven by; they should not be presented as independently proven
biological mechanisms unless validated against external curated evidence.

## Important Limitation

ARIANE currently has ENIGMA-defined clinically important domains and selected
pathogenic residue positions. It does not yet have a full protein-function map
with all active sites, interaction interfaces, phosphorylation sites, DNA
contact residues, PALB2/RAD51/BARD1 mechanistic annotations, or structural
features. Therefore, "active site" in a strict biochemical sense is not
available from the current local data.

## Dataset

- Class 4/5 variants analyzed: {total}
- Missense class 4/5 variants: {missense}
- Truncating or loss-of-function-like bucket: {truncating}
- Variants inside an ENIGMA functional domain: {domain_context}
- Variants overlapping a known ENIGMA pathogenic residue position: {residue_context}

## Mechanism Buckets

{markdown_table(tables["mechanism_summary"], ["mechanism_bucket", "count"], 20)}

## Aberrant Protein Buckets

{markdown_table(tables["aberrant_protein_summary"], ["aberrant_protein_bucket", "count"], 20)}

## Variant Types

{markdown_table(tables["variant_type_summary"], ["variant_type", "count"], 20)}

## ENIGMA Domain Context

{markdown_table(tables["gene_domain_summary"], ["gene", "enigma_domain", "count"], 30)}

## Missense Domain Context

This table is the most relevant one for questions about functional domains or
active-site-like protein effects. It excludes nonsense variants and asks where
class 4/5 missense variants fall relative to the ENIGMA domain/residue map.

{markdown_table(tables["missense_domain_summary"], ["gene", "enigma_domain", "known_pathogenic_residue_position", "count"], 30)}

## Known Pathogenic Residue Context

These are positions where the generated class 4/5 variant falls at a residue
position listed in the local ENIGMA clinically important residue table. This is
position context, not proof that every alternate amino acid at that position has
the same evidence.

{markdown_table(tables["known_pathogenic_residue_context"], ["gene", "p_notation", "known_context", "count"], 30)}

## Most Common Criteria Combinations

{markdown_table(tables["criteria_combo_summary"], ["criteria_combo", "count"], 25)}

## Mechanism By Class

{markdown_table(tables["class_mechanism_summary"], ["predicted_class", "mechanism_bucket", "count"], 30)}

## Interpretation

- In this coding SNV snapshot, class 4/5 pathogenicity is dominated by
  truncating mechanisms: nonsense variants with `PVS1` and often `PM5_PTC`.
- Missense class 4/5 variants are much fewer and are mostly driven by curated
  functional or multifactorial evidence (`PS3`, `PP4`, `PS1`) and by location in
  ENIGMA functional domains or known pathogenic residue positions.
- SpliceAI/PP3 can contribute to pathogenic combinations, but by itself it does
  not explain most class 4/5 variants.
- We can describe "functional domain context" with current local data. We
  cannot yet describe all protein active sites or detailed structural mechanism
  without adding a curated protein feature source.

## Outputs

- `tables/pathogenic_mechanisms/pathogenic_mechanism_annotated_variants.csv`
- `tables/pathogenic_mechanisms/mechanism_summary.csv`
- `tables/pathogenic_mechanisms/aberrant_protein_summary.csv`
- `tables/pathogenic_mechanisms/variant_type_summary.csv`
- `tables/pathogenic_mechanisms/domain_summary.csv`
- `tables/pathogenic_mechanisms/gene_domain_summary.csv`
- `tables/pathogenic_mechanisms/criteria_combo_summary.csv`
- `tables/pathogenic_mechanisms/known_pathogenic_residue_context.csv`
- `tables/pathogenic_mechanisms/missense_domain_summary.csv`
- `tables/pathogenic_mechanisms/class_mechanism_summary.csv`
- `tables/pathogenic_mechanisms/gene_mechanism_summary.csv`
- `tables/pathogenic_mechanisms/missense_domain_gene_summary.csv`
- `plots/16_pathogenic_mechanisms/mechanism_summary.svg`
- `plots/16_pathogenic_mechanisms/aberrant_protein_summary.svg`
- `plots/16_pathogenic_mechanisms/gene_domain_summary.svg`
- `plots/16_pathogenic_mechanisms/gene_by_mechanism_heatmap.svg`
- `plots/16_pathogenic_mechanisms/gene_by_domain_heatmap.svg`
- `plots/16_pathogenic_mechanisms/class_by_mechanism_heatmap.svg`
- `plots/16_pathogenic_mechanisms/missense_domain_heatmap.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    annotated, tables = analyze()

    write_csv(
        TABLE_DIR / "pathogenic_mechanism_annotated_variants.csv",
        annotated,
        [
            "gene",
            "c_notation",
            "p_notation",
            "aa_position",
            "variant_type",
            "predicted_class",
            "predicted_label",
            "total_points",
            "mechanism_bucket",
            "aberrant_protein_bucket",
            "enigma_domain",
            "known_pathogenic_residue_context",
            "spliceai_score",
            "criteria_combo",
            "criteria",
        ],
    )
    write_csv(TABLE_DIR / "mechanism_summary.csv", tables["mechanism_summary"], ["mechanism_bucket", "count"])
    write_csv(TABLE_DIR / "aberrant_protein_summary.csv", tables["aberrant_protein_summary"], ["aberrant_protein_bucket", "count"])
    write_csv(TABLE_DIR / "variant_type_summary.csv", tables["variant_type_summary"], ["variant_type", "count"])
    write_csv(TABLE_DIR / "domain_summary.csv", tables["domain_summary"], ["enigma_domain", "count"])
    write_csv(TABLE_DIR / "gene_domain_summary.csv", tables["gene_domain_summary"], ["gene", "enigma_domain", "count"])
    write_csv(TABLE_DIR / "criteria_combo_summary.csv", tables["criteria_combo_summary"], ["criteria_combo", "count"])
    write_csv(
        TABLE_DIR / "missense_domain_summary.csv",
        tables["missense_domain_summary"],
        ["gene", "enigma_domain", "known_pathogenic_residue_position", "count"],
    )
    write_csv(
        TABLE_DIR / "class_mechanism_summary.csv",
        tables["class_mechanism_summary"],
        ["predicted_class", "mechanism_bucket", "count"],
    )
    write_csv(
        TABLE_DIR / "gene_mechanism_summary.csv",
        tables["gene_mechanism_summary"],
        ["gene", "mechanism_bucket", "count"],
    )
    write_csv(
        TABLE_DIR / "missense_domain_gene_summary.csv",
        tables["missense_domain_gene_summary"],
        ["gene", "enigma_domain", "count"],
    )
    write_csv(
        TABLE_DIR / "known_pathogenic_residue_context.csv",
        tables["known_pathogenic_residue_context"],
        ["gene", "p_notation", "known_context", "count"],
    )

    bar_svg(PLOT_DIR / "mechanism_summary.svg", tables["mechanism_summary"], "mechanism_bucket", "count", "Class 4/5 Mechanism Buckets")
    bar_svg(PLOT_DIR / "aberrant_protein_summary.svg", tables["aberrant_protein_summary"], "aberrant_protein_bucket", "count", "Class 4/5 Aberrant Protein Buckets")
    gene_domain_plot = [
        {"label": f"{row['gene']} / {row['enigma_domain']}", "count": row["count"]}
        for row in tables["gene_domain_summary"]
    ]
    bar_svg(PLOT_DIR / "gene_domain_summary.svg", gene_domain_plot, "label", "count", "Class 4/5 ENIGMA Domain Context")
    heatmap_svg(
        PLOT_DIR / "gene_by_mechanism_heatmap.svg",
        tables["gene_mechanism_summary"],
        "mechanism_bucket",
        "gene",
        "count",
        "Class 4/5 Mechanism By Gene",
    )
    heatmap_svg(
        PLOT_DIR / "gene_by_domain_heatmap.svg",
        tables["gene_domain_summary"],
        "enigma_domain",
        "gene",
        "count",
        "Class 4/5 Domain Context By Gene",
    )
    heatmap_svg(
        PLOT_DIR / "class_by_mechanism_heatmap.svg",
        tables["class_mechanism_summary"],
        "mechanism_bucket",
        "predicted_class",
        "count",
        "Class 4 Versus Class 5 Mechanism",
    )
    heatmap_svg(
        PLOT_DIR / "missense_domain_heatmap.svg",
        tables["missense_domain_gene_summary"],
        "enigma_domain",
        "gene",
        "count",
        "Missense Class 4/5 Domain Context",
    )
    build_report(annotated, tables)


if __name__ == "__main__":
    main()
