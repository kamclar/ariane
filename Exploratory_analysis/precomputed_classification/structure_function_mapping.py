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
TABLE_DIR = ANALYSIS_DIR / "tables" / "structure_function_mapping"
PLOT_DIR = ANALYSIS_DIR / "plots" / "17_structure_function_mapping"
REPORT = ANALYSIS_DIR / "structure_function_mapping_report.md"
ACTIVE_SITE_ANNOTATIONS = TABLE_DIR / "curated_active_site_interface_annotations.csv"


STRUCTURE_FEATURES = {
    "BRCA1": [
        {
            "feature": "RING zinc-binding/E3 ligase region",
            "start": 2,
            "end": 101,
            "function": "BRCA1-BARD1 RING complex, zinc-coordinating fold, ubiquitin ligase activity and protein interactions",
            "structure_status": "experimental structures available for N-terminal RING/BARD1 region",
        },
        {
            "feature": "coiled-coil PALB2 interaction region",
            "start": 1391,
            "end": 1424,
            "function": "PALB2 interaction and homologous recombination complex assembly context",
            "structure_status": "local ENIGMA domain; detailed interface mapping not in local data",
        },
        {
            "feature": "BRCT phosphopeptide-binding region",
            "start": 1650,
            "end": 1857,
            "function": "phosphopeptide binding, DNA damage response recruitment and protein interactions",
            "structure_status": "many experimental BRCT structures available",
        },
    ],
    "BRCA2": [
        {
            "feature": "PALB2 binding N-terminal region",
            "start": 10,
            "end": 40,
            "function": "PALB2 binding context and BRCA1-PALB2-BRCA2 homologous recombination complex",
            "structure_status": "local ENIGMA domain; detailed interface mapping not in local data",
        },
        {
            "feature": "BRC repeats / RAD51-binding region",
            "start": 1000,
            "end": 2080,
            "function": "RAD51 binding and regulation; included as broad literature feature, not an ENIGMA BP1 domain",
            "structure_status": "broad feature range only; exact repeat/interface mapping should be curated before clinical use",
        },
        {
            "feature": "DNA-binding domain / helical-OB-DSS1 region",
            "start": 2481,
            "end": 3186,
            "function": "DNA binding, DSS1 interaction, OB/helical/tower domain context",
            "structure_status": "experimental structures available for C-terminal DNA-binding domain fragments",
        },
        {
            "feature": "C-terminal RAD51-binding/nuclear localization context",
            "start": 3187,
            "end": 3418,
            "function": "C-terminal RAD51-binding and nuclear localization context",
            "structure_status": "broad feature context only in this analysis",
        },
    ],
}

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


def load_pathogenic_positions() -> dict:
    payload = json.loads(RESIDUES.read_text(encoding="utf-8"))
    positions = {}
    for gene, domains in payload.get("domains", {}).items():
        positions[gene] = {}
        for domain_name, data in domains.items():
            for residue in data.get("pathogenic_residues", []):
                positions[gene].setdefault(int(residue["pos"]), []).append(
                    f"{domain_name}:{residue.get('aa')}:{residue.get('c')}"
                )
    return positions


def load_curated_active_site_annotations() -> dict[tuple[str, int], list[dict]]:
    if not ACTIVE_SITE_ANNOTATIONS.exists():
        return {}
    with ACTIVE_SITE_ANNOTATIONS.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    annotations: dict[tuple[str, int], list[dict]] = {}
    for row in rows:
        pos = parse_int(row.get("protein_position", ""), default=-1)
        if pos <= 0:
            continue
        annotations.setdefault((row["gene"], pos), []).append(row)
    return annotations


def load_curated_interface_regions() -> list[dict]:
    rows: list[dict] = []
    for path in sorted(TABLE_DIR.glob("uniprot_*_interface_regions.csv")):
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                row["protein_start"] = parse_int(row.get("protein_start", ""), default=-1)
                row["protein_end"] = parse_int(row.get("protein_end", ""), default=-1)
                if row["protein_start"] > 0 and row["protein_end"] >= row["protein_start"]:
                    rows.append(row)
    return rows


def feature_hits(gene: str, pos: Optional[int]) -> list[dict]:
    if pos is None:
        return []
    return [
        feature
        for feature in STRUCTURE_FEATURES.get(gene, [])
        if int(feature["start"]) <= pos <= int(feature["end"])
    ]


def structure_domain_context(features: list[dict], known_context: str) -> str:
    if features:
        return ";".join(feature["feature"] for feature in features)
    if known_context:
        return "known_ENIGMA_pathogenic_residue_without_local_domain_annotation"
    return "outside_mapped_structure_function_features"


def curated_active_site_status(
    gene: str,
    pos: Optional[int],
    annotations: dict[tuple[str, int], list[dict]],
    regions: list[dict],
) -> str:
    if pos is None:
        return "unknown_position"
    exact = annotations.get((gene, pos), [])
    region_hits = [
        row
        for row in regions
        if row["gene"] == gene and int(row["protein_start"]) <= pos <= int(row["protein_end"])
    ]
    labels = []
    if exact:
        labels.extend(
            f"{row['annotation_type']}:{row['annotation_label']} [{row['curation_status']}]"
            for row in exact
        )
    if region_hits:
        labels.extend(
            f"region_level:{row['annotation_type']}:{row['annotation_label']} "
            f"({row['protein_start']}-{row['protein_end']}) [{row['curation_status']}]"
            for row in region_hits
        )
    if labels:
        return ";".join(labels)
    return "not_curated_in_current_dataset"


def is_nonaberrant_no_splice_pathogenic(row: dict, criteria: Dict[str, dict]) -> bool:
    if row["predicted_class"] not in {"4", "5"}:
        return False
    if row["variant_type"] in {"nonsense", "frameshift", "splice_site"}:
        return False
    if "PVS1" in criteria or "PM5_PTC" in criteria:
        return False
    spliceai = parse_float(row.get("spliceai_score", ""))
    if spliceai is not None and spliceai >= 0.20:
        return False
    if "PP3" in criteria:
        return False
    return True


def write_csv(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def bar_svg(path: Path, rows: List[dict], label_key: str, value_key: str, title: str, max_items: int = 20) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = rows[:max_items]
    width = 1100
    left = 440
    top = 55
    row_h = 32
    chart_w = 560
    height = top + row_h * len(rows) + 55
    max_value = max((int(row[value_key]) for row in rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:13px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".small{font-size:11px;fill:#555}",
        ".bar{fill:#6f58a8}",
        "</style>",
        f'<text x="24" y="30" class="title">{html.escape(title)}</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = int(row[value_key])
        bar_w = 0 if max_value == 0 else int(chart_w * value / max_value)
        parts.append(f'<text x="24" y="{y + 20}">{html.escape(str(row[label_key]))}</text>')
        parts.append(f'<rect x="{left}" y="{y + 6}" width="{bar_w}" height="19" class="bar"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 21}" class="small">{value}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def lollipop_svg(path: Path, rows: List[dict], gene: str) -> None:
    gene_rows = [row for row in rows if row["gene"] == gene and row["aa_position"]]
    if not gene_rows:
        return
    length = 1863 if gene == "BRCA1" else 3418
    width = 1200
    height = 240
    left = 80
    right = 40
    axis_y = 120
    scale = (width - left - right) / length
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,sans-serif;font-size:12px;fill:#222}",
        ".title{font-size:20px;font-weight:bold}",
        ".axis{stroke:#333;stroke-width:3}",
        ".feature{fill:#e8e2f4;stroke:#8b76bd}",
        ".v4{fill:#f0a03a;stroke:#222}",
        ".v5{fill:#b23a48;stroke:#222}",
        ".small{font-size:10px;fill:#555}",
        "</style>",
        f'<text x="24" y="30" class="title">{gene} non-truncating no-splice class 4/5 variants</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width - right}" y2="{axis_y}" class="axis"/>',
    ]
    for feature in STRUCTURE_FEATURES.get(gene, []):
        x = left + int(feature["start"]) * scale
        w = max(3, (int(feature["end"]) - int(feature["start"]) + 1) * scale)
        parts.append(f'<rect x="{x}" y="{axis_y - 28}" width="{w}" height="18" class="feature"/>')
        parts.append(f'<text x="{x}" y="{axis_y - 35}" class="small">{html.escape(feature["feature"].split("/")[0])}</text>')
    for row in gene_rows:
        pos = int(row["aa_position"])
        x = left + pos * scale
        y2 = axis_y - 50 if row["predicted_class"] == "5" else axis_y + 50
        klass = "v5" if row["predicted_class"] == "5" else "v4"
        parts.append(f'<line x1="{x}" y1="{axis_y}" x2="{x}" y2="{y2}" stroke="#777"/>')
        parts.append(f'<circle cx="{x}" cy="{y2}" r="5" class="{klass}"/>')
    parts.append(f'<text x="{left}" y="{axis_y + 42}" class="small">0</text>')
    parts.append(f'<text x="{width - right - 35}" y="{axis_y + 42}" class="small">{length} aa</text>')
    parts.append('<text x="24" y="215" class="small">Upper dots: class 5. Lower dots: class 4. Purple boxes: broad structural/function features.</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def markdown_table(rows: List[dict], columns: List[str], limit: int = 30) -> str:
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
    pathogenic_positions = load_pathogenic_positions()
    active_site_annotations = load_curated_active_site_annotations()
    interface_regions = load_curated_interface_regions()
    annotated = []
    counters = {
        "feature": Counter(),
        "gene_feature": Counter(),
        "criteria_combo": Counter(),
        "variant_type": Counter(),
        "known_residue": Counter(),
        "curated_active_site_status": Counter(),
    }
    for row in load_rows():
        criteria = parse_criteria(row["criteria"])
        if not is_nonaberrant_no_splice_pathogenic(row, criteria):
            continue
        pos = aa_position(row["p_notation"])
        hits = feature_hits(row["gene"], pos)
        feature_names = [hit["feature"] for hit in hits] or ["outside_mapped_structure_function_features"]
        known_context = ";".join(pathogenic_positions.get(row["gene"], {}).get(pos or -1, []))
        domain_context = structure_domain_context(hits, known_context)
        active_status = curated_active_site_status(row["gene"], pos, active_site_annotations, interface_regions)
        annotated_row = {
            "gene": row["gene"],
            "c_notation": row["c_notation"],
            "p_notation": row["p_notation"],
            "aa_position": "" if pos is None else pos,
            "variant_type": row["variant_type"],
            "predicted_class": row["predicted_class"],
            "predicted_label": row["predicted_label"],
            "total_points": row["total_points"],
            "spliceai_score": row["spliceai_score"],
            "structure_features": ";".join(feature_names),
            "structure_domain_context": domain_context,
            "curated_active_site_status": active_status,
            "known_pathogenic_residue_context": known_context,
            "criteria_combo": criteria_combo(criteria),
            "criteria": row["criteria"],
        }
        annotated.append(annotated_row)
        counters["criteria_combo"][annotated_row["criteria_combo"]] += 1
        counters["variant_type"][row["variant_type"]] += 1
        if known_context:
            counters["known_residue"]["known_pathogenic_residue_position"] += 1
        else:
            counters["known_residue"]["no_known_pathogenic_residue_position"] += 1
        counters["curated_active_site_status"][active_status] += 1
        for feature in feature_names:
            counters["feature"][feature] += 1
            counters["gene_feature"][(row["gene"], feature)] += 1

    tables = {
        "feature_summary": [{"feature": k, "count": v} for k, v in counters["feature"].most_common()],
        "gene_feature_summary": [
            {"gene": k[0], "feature": k[1], "count": v}
            for k, v in counters["gene_feature"].most_common()
        ],
        "criteria_combo_summary": [
            {"criteria_combo": k, "count": v}
            for k, v in counters["criteria_combo"].most_common()
        ],
        "variant_type_summary": [
            {"variant_type": k, "count": v}
            for k, v in counters["variant_type"].most_common()
        ],
        "known_residue_summary": [
            {"known_pathogenic_residue_context": k, "count": v}
            for k, v in counters["known_residue"].most_common()
        ],
        "curated_active_site_status_summary": [
            {"curated_active_site_status": k, "count": v}
            for k, v in counters["curated_active_site_status"].most_common()
        ],
    }
    return annotated, tables


def build_report(annotated: list[dict], tables: dict[str, list[dict]]) -> None:
    total = len(annotated)
    missense = sum(1 for row in annotated if row["variant_type"] == "missense")
    known = sum(1 for row in annotated if row["known_pathogenic_residue_context"])
    text = f"""# Structure Function Mapping

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Input: `{INPUT.relative_to(ROOT)}`

## Why We Did This

This analysis focuses on generated class 4/5 variants that are pathogenic
without an obvious aberrant-protein or splice driver in the current snapshot.
It excludes:

- nonsense/frameshift/splice-site variants
- variants with `PVS1` or `PM5_PTC`
- variants with SpliceAI >= 0.20
- variants with `PP3`

The goal is to ask: among the remaining pathogenic variants, do they map to
protein regions where a missense or initiation-codon change could plausibly
affect protein function?

## Dataset

- Non-truncating, no-strong-splice class 4/5 variants: {total}
- Missense variants in this subset: {missense}
- Variants at a known ENIGMA pathogenic residue position: {known}

## Important Interpretation Boundary

This is not full 3D structural modeling. It maps variants to broad
structure/function features using local ENIGMA domains plus a small curated
feature scaffold from public protein annotations. A real 3D analysis would need
explicit PDB or AlphaFold residue coordinates, structure confidence, partner
interfaces, ligand/DNA contacts, and probably separate models for BRCA1-BARD1
and BRCA2-DSS1/RAD51/PALB2 complexes.

## Feature Summary

{markdown_table(tables["feature_summary"], ["feature", "count"], 20)}

## Gene Feature Summary

{markdown_table(tables["gene_feature_summary"], ["gene", "feature", "count"], 30)}

## Criteria Driving This Subset

{markdown_table(tables["criteria_combo_summary"], ["criteria_combo", "count"], 20)}

## Variant Types

{markdown_table(tables["variant_type_summary"], ["variant_type", "count"], 10)}

## Known Pathogenic Residue Context

{markdown_table(tables["known_residue_summary"], ["known_pathogenic_residue_context", "count"], 10)}

## Curated Active Site Or Interface Status

The analysis now reads residue-level active-site/interface annotations from
`tables/structure_function_mapping/curated_active_site_interface_annotations.csv`.
It also reads UniProt region-level interface annotations from
`tables/structure_function_mapping/uniprot_brca1_interface_regions.csv` and
`tables/structure_function_mapping/uniprot_brca2_interface_regions.csv`.
This is deliberately separate from broad domain context. A variant can sit in a
BRCT or DBD domain without being annotated as an exact binding/interface residue.
Rows prefixed with `region_level` should be interpreted as regional context, not
as proof that the exact residue directly contacts a partner.

{markdown_table(tables["curated_active_site_status_summary"], ["curated_active_site_status", "count"], 20)}

## Top Variants

{markdown_table(annotated, ["gene", "c_notation", "p_notation", "variant_type", "predicted_class", "structure_features", "structure_domain_context", "curated_active_site_status", "known_pathogenic_residue_context", "criteria_combo"], 40)}

## Interpretation

- The class 4/5 variants without truncation and without strong splice/PP3 signal
  form a small subset.
- Most of this subset is driven by curated functional/multifactorial evidence,
  especially `PS3` and `PP4`, rather than by rules inferred from structure.
- The strongest structure/function concentration is in BRCA1 RING and BRCT
  regions, with some BRCA2 DNA-binding-domain variants.
- For truly structural interpretation, the next step is to add residue-level
  coordinates from selected experimental PDB structures and/or AlphaFold models
  and annotate partner interfaces.

## Active Site Annotation Curation

The first residue-level table was generated from downloaded UniProt BRCA1/BRCA2
snapshots plus a source-tracked seed for the BRCA1 RING zinc-coordinating core.
UniProt provides strong domain and interaction context, but not a complete
ready-made residue-level interface table for all BRCA1/2 functional surfaces.
Some BRCT and BRCA2 positions are now extracted from UniProt `Mutagenesis` or
`Natural variant` records with interaction effects. Remaining BRCT
phosphopeptide-binding surface details and BRCA2 DNA/DSS1/RAD51/PALB2 contact
residues still need explicit PDB/literature curation before we should treat them
as exact active-site/interface annotations.

Minimum useful fields for that table:

- `gene`
- `protein_position`
- `annotation_type`
- `annotation_label`
- `evidence_source`
- `source_id`
- `structure_id`
- `notes`

Candidate sources to curate from:

- UniProt feature annotations for BRCA1 and BRCA2
- selected experimental PDB structures
- AlphaFold models as spatial context, with confidence filtering
- InterPro/Pfam domain annotations
- literature for BRCA1-BARD1 RING, BRCA1 BRCT phosphopeptide binding, BRCA2
  DNA/DSS1/RAD51/PALB2 interfaces

## Suggested Next 3D Step

Start with the small subset in
`tables/structure_function_mapping/nontruncating_no_splice_pathogenic_variants.csv`.
For each variant, map the amino-acid position to a structure source:

- BRCA1 RING: experimental BRCA1-BARD1 N-terminal structures
- BRCA1 BRCT: experimental BRCT phosphopeptide-binding structures
- BRCA2 DBD: experimental DNA-binding/DSS1 region structures

Then annotate distance to zinc-coordinating residues, phosphopeptide-binding
surface, DNA/DSS1/RAD51/PALB2 interface, or other curated functional surfaces.

## Sources

- UniProt BRCA1 P38398: https://rest.uniprot.org/uniprotkb/P38398.txt
- UniProt BRCA2 P51587: https://rest.uniprot.org/uniprotkb/P51587.txt
- AlphaFold BRCA1 P38398: https://alphafold.ebi.ac.uk/entry/P38398
- AlphaFold BRCA2 P51587: https://alphafold.ebi.ac.uk/entry/P51587

## Outputs

- `tables/structure_function_mapping/nontruncating_no_splice_pathogenic_variants.csv`
- `tables/structure_function_mapping/feature_summary.csv`
- `tables/structure_function_mapping/gene_feature_summary.csv`
- `tables/structure_function_mapping/criteria_combo_summary.csv`
- `tables/structure_function_mapping/variant_type_summary.csv`
- `tables/structure_function_mapping/known_residue_summary.csv`
- `tables/structure_function_mapping/curated_active_site_status_summary.csv`
- `tables/structure_function_mapping/curated_active_site_interface_annotations.csv`
- `tables/structure_function_mapping/uniprot_brca1_interface_regions.csv`
- `tables/structure_function_mapping/uniprot_brca2_interface_regions.csv`
- `plots/17_structure_function_mapping/feature_summary.svg`
- `plots/17_structure_function_mapping/gene_feature_summary.svg`
- `plots/17_structure_function_mapping/curated_active_site_status_summary.svg`
- `plots/17_structure_function_mapping/brca1_lollipop.svg`
- `plots/17_structure_function_mapping/brca2_lollipop.svg`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    annotated, tables = analyze()
    fields = [
        "gene",
        "c_notation",
        "p_notation",
        "aa_position",
        "variant_type",
        "predicted_class",
        "predicted_label",
        "total_points",
        "spliceai_score",
        "structure_features",
        "structure_domain_context",
        "curated_active_site_status",
        "known_pathogenic_residue_context",
        "criteria_combo",
        "criteria",
    ]
    write_csv(TABLE_DIR / "nontruncating_no_splice_pathogenic_variants.csv", annotated, fields)
    write_csv(TABLE_DIR / "feature_summary.csv", tables["feature_summary"], ["feature", "count"])
    write_csv(TABLE_DIR / "gene_feature_summary.csv", tables["gene_feature_summary"], ["gene", "feature", "count"])
    write_csv(TABLE_DIR / "criteria_combo_summary.csv", tables["criteria_combo_summary"], ["criteria_combo", "count"])
    write_csv(TABLE_DIR / "variant_type_summary.csv", tables["variant_type_summary"], ["variant_type", "count"])
    write_csv(TABLE_DIR / "known_residue_summary.csv", tables["known_residue_summary"], ["known_pathogenic_residue_context", "count"])
    write_csv(TABLE_DIR / "curated_active_site_status_summary.csv", tables["curated_active_site_status_summary"], ["curated_active_site_status", "count"])
    bar_svg(PLOT_DIR / "feature_summary.svg", tables["feature_summary"], "feature", "count", "Non-truncating No-splice Pathogenic Variants By Protein Feature")
    bar_svg(
        PLOT_DIR / "curated_active_site_status_summary.svg",
        tables["curated_active_site_status_summary"],
        "curated_active_site_status",
        "count",
        "Curated Active Site Or Interface Status",
    )
    gene_feature_plot = [
        {"label": f"{row['gene']} / {row['feature']}", "count": row["count"]}
        for row in tables["gene_feature_summary"]
    ]
    bar_svg(PLOT_DIR / "gene_feature_summary.svg", gene_feature_plot, "label", "count", "Protein Feature By Gene")
    lollipop_svg(PLOT_DIR / "brca1_lollipop.svg", annotated, "BRCA1")
    lollipop_svg(PLOT_DIR / "brca2_lollipop.svg", annotated, "BRCA2")
    build_report(annotated, tables)


if __name__ == "__main__":
    main()
