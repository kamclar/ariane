# Data Generation Report

Generated: 2026-06-24 10:57

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
| `variant_space_scan/outputs/brca_snv_manifest.csv` | 47547 | initial coding SNV manifest |
| `variant_space_scan/outputs/brca_snv_manifest.with_coordinates.local_full.csv` | 47547 | manifest with local coordinate mapping |
| `variant_space_scan/outputs/brca_snv_manifest.with_spliceai.ref_block_033_047_final.csv` | 47547 | manifest with reference-transcript SpliceAI scores |
| `variant_space_scan/outputs/brca_module1_full_snv_classification.csv` | 47547 | final Module 1 classification snapshot |

The GitHub export stores the final classification snapshot as
`data/brca_module1_full_snv_classification.csv.gz`.

## Dataset Summary

Total variants: `47547`

### By Gene

| gene | count |
| --- | --- |
| BRCA1 | 16776 |
| BRCA2 | 30771 |

### By Generated Class

| class | count |
| --- | --- |
| 1 | 802 |
| 2 | 36150 |
| 3 | 8154 |
| 4 | 141 |
| 5 | 2300 |

### By Variant Type

| variant_type | count |
| --- | --- |
| initiation_codon | 18 |
| missense | 35241 |
| nonsense | 2397 |
| synonymous | 9891 |

### VUS and SpliceAI Summary

Class 3 VUS: `8154`

Variants with reference-transcript SpliceAI score >= 0.20:
`1205`

## Automatically Applied Criteria

The following table reports how often each criterion was assigned by ARIANE
Module 1. These counts reflect the implemented automated rules and the
availability of their required input data. They do not represent complete
variant assessments because non-automatable evidence and expert review are not
included. Multiple criteria may be assigned to the same variant.

| criterion | count |
| --- | --- |
| PM2_Supporting | 43166 |
| BP1 | 34855 |
| BS3 | 2528 |
| PVS1 | 2327 |
| PM5_PTC | 2241 |
| BP4 | 1404 |
| BP7 | 1404 |
| PP3 | 1060 |
| PS3 | 699 |
| BP5 | 514 |
| BS1_Supporting | 141 |
| BS1_Strong | 101 |
| PP4 | 58 |
| BA1 | 48 |
| PS1 | 7 |

## Key VUS Exploration Outputs

### Evidence Action Groups

| action_group | count | high_priority_count | high_spliceai_count | median_points |
| --- | --- | --- | --- | --- |
| background_absence_only | 5048 | 0 | 10 | 1.0 |
| conflicting_computational_splice_context | 1276 | 0 | 0 | -1.0 |
| computational_signal_needs_noncomputational_support | 842 | 812 | 842 | 2.0 |
| near_pathogenic_threshold | 446 | 38 | 5 | 5.0 |
| module1_no_signal | 331 | 0 | 1 | 0 |
| near_benign_threshold | 191 | 2 | 0 | -4 |
| manual_conflict_adjudication | 20 | 7 | 8 | -0.5 |

### Findlay SGE BRCA1 VUS Check

| priority_tier | count | median_function_score | lof_count | int_count | func_count |
| --- | --- | --- | --- | --- | --- |
| tier1_urgent | 9 | -1.951 | 6 | 3 | 0 |
| tier2_high | 46 | -1.647 | 30 | 15 | 1 |
| tier3_medium | 365 | -1.891 | 303 | 15 | 47 |
| tier4_low | 253 | -0.828 | 14 | 139 | 100 |

### Public Classification Snapshot

| public_snapshot_category | count |
| --- | --- |
| multi-submitter assertion | 28 |
| single-submitter assertion | 28 |
| conflicting public assertions | 11 |
| no public assertion | 7 |
| panel-level public assertion | 4 |
| ClinGen/ENIGMA assertion | 2 |

### Priority Queue Synthesis

| queue_label | count |
| --- | --- |
| neighborhood context | 6649 |
| low information | 5379 |
| splice/RNA | 3291 |
| functional | 836 |
| near threshold | 637 |
| public assertion | 72 |
| evidence conflict | 30 |

### Gene-Normalized Grouped Classes

| gene | label | count | percent_within_gene | per_1000 |
| --- | --- | --- | --- | --- |
| BRCA1 | Benign/Likely Benign | 14274 | 85.09 | 850.86 |
| BRCA2 | Benign/Likely Benign | 22678 | 73.70 | 736.99 |
| BRCA1 | VUS | 1625 | 9.69 | 96.86 |
| BRCA2 | VUS | 6529 | 21.22 | 212.18 |
| BRCA1 | Likely Pathogenic/Pathogenic | 877 | 5.23 | 52.28 |
| BRCA2 | Likely Pathogenic/Pathogenic | 1564 | 5.08 | 50.83 |

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
