# ARIANE Module 1 BRCA1/2 VUS Exploration

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
python scripts\prepare_vus_github_export.py
```

The script recreates the whole `github_exports/vus_exploration_module1`
directory. This keeps the public copy clean and avoids manually copying files.
A zip archive can still be created for sharing, but it is not the primary
working format.

Important boundary:

This is a review-triage and method-development export. It is not a clinical
classification release and should not be used as standalone medical evidence.
