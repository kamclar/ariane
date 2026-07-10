"""Prepare a clean GitHub export for the ARIANE Module 1 VUS exploration."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "Exploratory_analysis" / "precomputed_classification"
OUTPUT = ROOT / "github_exports" / "vus_exploration_module1"
SNAPSHOT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
SNAPSHOT_SUMMARY = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.summary.json"
MANIFEST = ROOT / "variant_space_scan" / "outputs" / "brca_snv_manifest.csv"
COORDS = ROOT / "variant_space_scan" / "outputs" / "brca_snv_manifest.with_coordinates.local_full.csv"
SPLICEAI = ROOT / "variant_space_scan" / "outputs" / "brca_snv_manifest.with_spliceai.ref_block_033_047_final.csv"
EXPORT_PLOT_ALLOWLIST = {
    "07_gene_comparison": {
        "domain_vs_domain_within_domain_percent_heatmap.svg",
        "domain_vs_domain_normalized_density_heatmap.svg",
        "domain_region_signal_heatmap.svg",
        "domain_region_grouped_class_mix.svg",
    },
}


KEY_REPORTS = [
    "README.md",
    "vus_prioritization_module1_exploratory_report.md",
    "vus_prioritization_module1_exploratory_report_updated_20260623_v22.docx",
    "findlay_sge_vus_tier_report.md",
    "vus_evidence_action_plan_report.md",
    "public_classification_snapshot_report.md",
    "priority_queue_synthesis_report.md",
    "cross_gene_normalized_report.md",
    "gene_comparison_report.md",
    "vus_bottleneck_report.md",
    "vus_manuscript_critique_response_report.md",
    "bs3_functional_source_audit_report.md",
]


KEY_TABLES = [
    "class_distribution.csv",
    "gene_by_class.csv",
    "criteria_counts.csv",
    "criteria_by_class.csv",
    "variant_type_by_class.csv",
    "spliceai_bins_by_class.csv",
    "gnomad_coarse_status.csv",
    "vus_prioritization/vus_priority_annotated.csv",
    "vus_prioritization/vus_priority_top_by_category.csv",
    "vus_bottleneck/vus_bottleneck_summary.csv",
    "vus_bottleneck/vus_bottleneck_combinations.csv",
    "vus_evidence_action_plan/vus_evidence_action_plan_by_action_group.csv",
    "vus_evidence_action_plan/vus_evidence_action_plan_by_bottleneck.csv",
    "vus_evidence_action_plan/vus_evidence_action_plan_top_candidates.csv",
    "findlay_sge_vus_tier/findlay_sge_by_priority_tier.csv",
    "findlay_sge_vus_tier/findlay_sge_lof_gap_by_tier.csv",
    "findlay_sge_vus_tier/findlay_sge_tier3_lof_action_overlap.csv",
    "public_classification_snapshot/public_classification_snapshot_by_category.csv",
    "public_classification_snapshot/public_classification_snapshot_by_discordance.csv",
    "public_classification_snapshot/public_classification_snapshot_variants.csv",
    "priority_queue_synthesis/priority_queue_counts.csv",
    "priority_queue_synthesis/priority_queue_overlap.csv",
    "priority_queue_synthesis/priority_queue_by_action_group.csv",
    "priority_queue_synthesis/priority_queue_top_multiqueue_examples.csv",
    "priority_queue_synthesis/priority_queue_synthesis_variants.csv",
    "cross_gene_normalized/gene_normalized_class_distribution.csv",
    "cross_gene_normalized/gene_normalized_grouped_class_distribution.csv",
    "cross_gene_normalized/gene_normalized_criteria_counts.csv",
    "cross_gene_normalized/gene_normalized_vus_priority_tiers.csv",
    "cross_gene_normalized/gene_normalized_evidence_action_groups.csv",
    "cross_gene_normalized/gene_normalized_priority_queues.csv",
    "gene_comparison/domain_region_summary.csv",
    "gene_comparison/domain_region_vus_examples.csv",
    "brca1_boundary_spliceai_binned_heatmap.csv",
    "brca2_boundary_spliceai_binned_heatmap.csv",
    "vus_manuscript_critique_response/neighborhood_increment_summary.csv",
    "vus_manuscript_critique_response/tier_transition_without_neighborhood.csv",
    "bs3_functional_source_audit/brca1_bs3_findlay_match_summary.csv",
]


def read_count(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def copy_file(src: Path, dst: Path) -> None:
    if not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def reset_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for child in path.iterdir():
        try:
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
        except PermissionError:
            # Windows can keep a viewed directory locked. Keep it and overwrite
            # files during copy; stale files are acceptable for the export
            # workflow and safer than failing the refresh.
            continue


def copy_tree_filtered(src: Path, dst: Path, suffixes: set[str]) -> None:
    if not src.exists():
        return
    for path in src.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in suffixes:
            continue
        rel = path.relative_to(src)
        allowed_names = EXPORT_PLOT_ALLOWLIST.get(rel.parts[0])
        if allowed_names is not None and path.name not in allowed_names:
            continue
        copy_file(path, dst / rel)


def gzip_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open("rb") as source, gzip.open(dst, "wb", compresslevel=9) as target:
        shutil.copyfileobj(source, target)


def snapshot_stats(path: Path) -> dict[str, object]:
    genes: Counter[str] = Counter()
    classes: Counter[str] = Counter()
    variant_types: Counter[str] = Counter()
    criteria: Counter[str] = Counter()
    high_splice = 0
    rows = 0
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows += 1
            genes[row["gene"]] += 1
            classes[row["predicted_class"]] += 1
            variant_types[row["variant_type"]] += 1
            try:
                if float(row.get("spliceai_score") or 0) >= 0.20:
                    high_splice += 1
            except ValueError:
                pass
            for item in (row.get("criteria") or "").split(";"):
                if item:
                    criteria[item.split(":")[0]] += 1
    return {
        "rows": rows,
        "genes": dict(genes),
        "classes": dict(classes),
        "variant_types": dict(variant_types),
        "vus_count": classes.get("3", 0),
        "high_spliceai_ge_0_20": high_splice,
        "criteria_counts": dict(criteria),
    }


def csv_to_markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def count_rows_from_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_data_generation_report(stats: dict[str, object]) -> None:
    class_rows = [
        {"class": key, "count": value}
        for key, value in sorted(stats["classes"].items(), key=lambda item: int(item[0]))
    ]
    gene_rows = [
        {"gene": key, "count": value}
        for key, value in sorted(stats["genes"].items())
    ]
    type_rows = [
        {"variant_type": key, "count": value}
        for key, value in sorted(stats["variant_types"].items())
    ]
    criteria_rows = [
        {"criterion": key, "count": value}
        for key, value in sorted(stats["criteria_counts"].items(), key=lambda item: (-item[1], item[0]))
    ]

    action_rows = count_rows_from_csv(ANALYSIS / "tables" / "vus_evidence_action_plan" / "vus_evidence_action_plan_by_action_group.csv")
    public_rows = count_rows_from_csv(ANALYSIS / "tables" / "public_classification_snapshot" / "public_classification_snapshot_by_category.csv")
    findlay_rows = count_rows_from_csv(ANALYSIS / "tables" / "findlay_sge_vus_tier" / "findlay_sge_by_priority_tier.csv")
    queue_rows = count_rows_from_csv(ANALYSIS / "tables" / "priority_queue_synthesis" / "priority_queue_counts.csv")
    normalized_group_rows = count_rows_from_csv(ANALYSIS / "tables" / "cross_gene_normalized" / "gene_normalized_grouped_class_distribution.csv")

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# Data Generation Report

Generated: {generated}

## Scope

This repository contains an exploratory analysis of a precomputed ARIANE
Module 1 dataset covering coding SNVs in BRCA1 and BRCA2.

Module 1 applies only automatable ACMG/AMP and ENIGMA VCEP criteria. The
resulting classes therefore describe the behavior of the current computational
framework and must not be interpreted as complete clinical classifications or
biological truth.

Even classifications based on the full set of applicable criteria remain
evidence-based interpretations rather than direct observations of biological
effects. This framework may later be extended to include broader criteria
coverage and additional evidence sources.

## Generation Pipeline

1. Generate the BRCA1/2 coding SNV manifest.
2. Resolve transcript variants to genomic coordinates.
3. Attach reference-transcript SpliceAI scores from the local precomputed cache.
4. Run ARIANE Module 1 automated ENIGMA/ACMG rules over the coding SNV manifest.
5. Export the full classification snapshot.
6. Run exploratory scripts over the snapshot to produce VUS tiers, bottlenecks,
   public assertion context, functional SGE checks, and figures.

## Main Input And Output Files

| File | Rows | Notes |
| --- | ---: | --- |
| `variant_space_scan/outputs/brca_snv_manifest.csv` | {read_count(MANIFEST) if MANIFEST.exists() else "missing"} | initial coding SNV manifest |
| `variant_space_scan/outputs/brca_snv_manifest.with_coordinates.local_full.csv` | {read_count(COORDS) if COORDS.exists() else "missing"} | manifest with local coordinate mapping |
| `variant_space_scan/outputs/brca_snv_manifest.with_spliceai.ref_block_033_047_final.csv` | {read_count(SPLICEAI) if SPLICEAI.exists() else "missing"} | manifest with reference-transcript SpliceAI scores |
| `variant_space_scan/outputs/brca_module1_full_snv_classification.csv` | {stats["rows"]} | final Module 1 classification snapshot |

The GitHub export stores the final classification snapshot as
`data/brca_module1_full_snv_classification.csv.gz`.

## Dataset Summary

Total variants: `{stats["rows"]}`

### By Gene

{csv_to_markdown_table(gene_rows, ["gene", "count"])}

### By Generated Class

{csv_to_markdown_table(class_rows, ["class", "count"])}

### By Variant Type

{csv_to_markdown_table(type_rows, ["variant_type", "count"])}

### VUS and SpliceAI Summary

Class 3 VUS: `{stats["vus_count"]}`

Variants with reference-transcript SpliceAI score >= 0.20:
`{stats["high_spliceai_ge_0_20"]}`

## Automatically Applied Criteria

The following table reports how often each criterion was assigned by ARIANE
Module 1. These counts reflect the implemented automated rules and the
availability of their required input data. They do not represent complete
variant assessments because non-automatable evidence and expert review are not
included. Multiple criteria may be assigned to the same variant.

{csv_to_markdown_table(criteria_rows, ["criterion", "count"])}

## Key VUS Exploration Outputs

### Evidence Action Groups

{csv_to_markdown_table(action_rows, ["action_group", "count", "high_priority_count", "high_spliceai_count", "median_points"])}

### Findlay SGE BRCA1 VUS Check

{csv_to_markdown_table(findlay_rows, ["priority_tier", "count", "median_function_score", "lof_count", "int_count", "func_count"])}

### Public Classification Snapshot

{csv_to_markdown_table(public_rows, ["public_snapshot_category", "count"])}

### Priority Queue Synthesis

{csv_to_markdown_table(queue_rows, ["queue_label", "count"])}

### Gene-Normalized Grouped Classes

{csv_to_markdown_table(normalized_group_rows, ["gene", "label", "count", "percent_within_gene", "per_1000"])}

## Included Export Contents

- `data/`: compressed final classification snapshot and summary JSON
- `docs/`: main exploratory report and key supporting reports
- `scripts/`: analysis scripts used to generate tables and figures
- `tables/`: selected summary tables and review queues
- `plots/`: SVG plots generated by the exploratory analysis

## Excluded From Export

Large or redundant intermediate tables are not copied into this export. The
full working directory contains many development-stage outputs, repeated Word
exports, and intermediate validation files. The export keeps the final snapshot,
summary tables, scripts, and documentation needed to understand or rebuild the
analysis.
"""
    (OUTPUT / "DATA_GENERATION_REPORT.md").write_text(text, encoding="utf-8")


def write_export_readme() -> None:
    text = """# ARIANE Module 1 BRCA1/2 VUS Exploration

This repository contains an exploratory analysis of a precomputed ARIANE
Module 1 dataset covering coding SNVs in BRCA1 and BRCA2.

Module 1 applies only automatable ACMG/AMP and ENIGMA VCEP criteria. The
resulting classes therefore describe the behavior of the current computational
framework and must not be interpreted as complete clinical classifications or
biological truth.

Even classifications based on the full set of applicable criteria remain
evidence-based interpretations rather than direct observations of biological
effects. This framework may later be extended to include broader criteria
coverage and additional evidence sources.

The main purpose is to show how a precomputed Module 1 map can be converted
into an auditable manual-review queue for VUS. The analysis separates
variant-level signals, neighborhood-driven context, evidence-conflict cases,
public assertion context, and functional SGE checks.

Start here:

- `DATA_GENERATION_REPORT.md`
- `docs/vus_prioritization_module1_exploratory_report.md`
- `docs/vus_prioritization_module1_exploratory_report_updated_20260623_v22.docx`

## How This GitHub Export Is Maintained

This directory is intended to be the GitHub-ready copy of the analysis. The
working analysis stays in the main ARIANE workspace, while this export contains
only the files that should be published or shared:

- compressed final Module 1 snapshot
- selected tables
- SVG figures
- analysis scripts
- documentation and reports

To refresh the GitHub-ready directory after new analysis results are generated,
run from the ARIANE workspace root:

```powershell
python scripts\\prepare_vus_github_export.py
```

The script recreates the whole `github_exports/vus_exploration_module1`
directory. This keeps the public copy clean and avoids manually copying files.
A zip archive can still be created for sharing, but it is not the primary
working format.

Important boundary:

This is a review-triage and method-development export. It is not a clinical
classification release and should not be used as standalone medical evidence.
"""
    (OUTPUT / "README.md").write_text(text, encoding="utf-8")


def write_exports_readme() -> None:
    text = """# GitHub Exports

This directory contains clean, shareable exports generated from the ARIANE
working analysis.

## Main Export

- `vus_exploration_module1/`: GitHub-ready directory for the BRCA1/2 Module 1
  VUS exploration.
- `vus_exploration_module1.zip`: optional archive for sending or storing a
  snapshot. The directory is the preferred source for GitHub.

## Refresh Workflow

From the ARIANE workspace root:

```powershell
python scripts\\prepare_vus_github_export.py
```

This rebuilds `github_exports/vus_exploration_module1` from the current
analysis outputs. Commit or upload that directory to GitHub.

Only create the zip when a separate archive is needed:

```powershell
Compress-Archive -Path github_exports\\vus_exploration_module1\\* -DestinationPath github_exports\\vus_exploration_module1.zip -Force
```
"""
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    (OUTPUT.parent / "README.md").write_text(text, encoding="utf-8")


def write_external_sources_note() -> None:
    text = """# External Sources

This export does not redistribute third-party PDFs or spreadsheets from journal
supplements. The working analysis used local copies of the following sources:

- Findlay et al. 2018 BRCA1 saturation genome editing Supplementary Table 1
- Findlay et al. 2018 Supplementary Information
- Dace/Findlay 2023 interim update as cited in local ENIGMA Table 9 source notes
- ClinVar VCV XML fetched through NCBI eutils
- ClinGen Evidence Repository API

Derived outputs included in this export:

- `tables/findlay_sge_vus_tier/findlay_sge_by_priority_tier.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_lof_gap_by_tier.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_tier3_lof_action_overlap.csv`
- `tables/bs3_functional_source_audit/brca1_bs3_findlay_match_summary.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_variants.csv`
- `tables/priority_queue_synthesis/priority_queue_counts.csv`
- `tables/priority_queue_synthesis/priority_queue_overlap.csv`
- `tables/priority_queue_synthesis/priority_queue_synthesis_variants.csv`
- `tables/cross_gene_normalized/gene_normalized_grouped_class_distribution.csv`
- `tables/cross_gene_normalized/gene_normalized_priority_queues.csv`
- `tables/cross_gene_normalized/gene_normalized_evidence_action_groups.csv`
- `docs/findlay_sge_vus_tier_report.md`
- `docs/public_classification_snapshot_report.md`
- `docs/priority_queue_synthesis_report.md`
- `docs/cross_gene_normalized_report.md`

If full reproducibility from raw external sources is needed, download the
source files from their original providers and place them into the paths
documented by the analysis scripts.
"""
    (OUTPUT / "EXTERNAL_SOURCES.md").write_text(text, encoding="utf-8")


def write_gitignore() -> None:
    text = """__pycache__/
*.pyc
.DS_Store
Thumbs.db
*.tmp
"""
    (OUTPUT / ".gitignore").write_text(text, encoding="utf-8")


def main() -> None:
    reset_output_dir(OUTPUT)

    stats = snapshot_stats(SNAPSHOT)

    gzip_file(SNAPSHOT, OUTPUT / "data" / "brca_module1_full_snv_classification.csv.gz")
    copy_file(SNAPSHOT_SUMMARY, OUTPUT / "data" / SNAPSHOT_SUMMARY.name)

    for report in KEY_REPORTS:
        copy_file(ANALYSIS / report, OUTPUT / "docs" / report)

    for script in ANALYSIS.glob("*.py"):
        if script.name == "parse_xlsx_to_csv.py":
            continue
        copy_file(script, OUTPUT / "scripts" / script.name)

    for rel in KEY_TABLES:
        copy_file(ANALYSIS / "tables" / rel, OUTPUT / "tables" / rel)

    copy_tree_filtered(ANALYSIS / "plots", OUTPUT / "plots", {".svg"})

    write_data_generation_report(stats)
    write_export_readme()
    write_exports_readme()
    write_external_sources_note()
    write_gitignore()

    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(OUTPUT.relative_to(ROOT)),
        "snapshot_rows": stats["rows"],
        "vus_count": stats["vus_count"],
        "included_reports": KEY_REPORTS,
        "included_tables": KEY_TABLES,
    }
    (OUTPUT / "export_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
