# Exploratory Analysis: Precomputed BRCA Module 1 SNV Classification

Generated: 2026-06-23 13:10

Input snapshot: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

This is an exploratory analysis of the automated Module 1 overlay. It is not an expert final classification and it does not add non-automated ACMG/ENIGMA criteria such as PS4, PM3, PP1, BS2, BS4, or curated RNA conclusions.

## Dataset

- Total variants: 47547
- Coordinate status: {'ok': 47547}
- Variants with SpliceAI score: 47547

## Class Distribution

| Class | Label | Count | Percent |
|---|---:|---:|---:|
| 1 | Benign | 802 | 1.69% |
| 2 | Likely Benign | 36150 | 76.03% |
| 3 | VUS | 8154 | 17.15% |
| 4 | Likely Pathogenic | 141 | 0.30% |
| 5 | Pathogenic | 2300 | 4.84% |

## Gene Coverage

| Gene | Count | Percent |
|---|---:|---:|
| BRCA1 | 16776 | 35.28% |
| BRCA2 | 30771 | 64.72% |

## Gene-Normalized Class Distribution

Because BRCA1 and BRCA2 contribute different numbers of coding SNVs to the
snapshot, cross-gene comparisons should use within-gene percentages or rates
per 1000 generated variants, not raw combined counts.

| Gene | Class | Label | Count | Percent within gene |
|---|---:|---|---:|---:|
| BRCA1 | 1 | Benign | 519 | 3.09% |
| BRCA1 | 2 | Likely Benign | 13755 | 81.99% |
| BRCA1 | 3 | VUS | 1625 | 9.69% |
| BRCA1 | 4 | Likely Pathogenic | 64 | 0.38% |
| BRCA1 | 5 | Pathogenic | 813 | 4.85% |
| BRCA2 | 1 | Benign | 283 | 0.92% |
| BRCA2 | 2 | Likely Benign | 22395 | 72.78% |
| BRCA2 | 3 | VUS | 6529 | 21.22% |
| BRCA2 | 4 | Likely Pathogenic | 77 | 0.25% |
| BRCA2 | 5 | Pathogenic | 1487 | 4.83% |

### Grouped Classes By Gene

| Gene | Grouped class | Count | Percent within gene |
|---|---|---:|---:|
| BRCA1 | Benign/Likely Benign | 14274 | 85.09% |
| BRCA1 | VUS | 1625 | 9.69% |
| BRCA1 | Likely Pathogenic/Pathogenic | 877 | 5.23% |
| BRCA2 | Benign/Likely Benign | 22678 | 73.70% |
| BRCA2 | VUS | 6529 | 21.22% |
| BRCA2 | Likely Pathogenic/Pathogenic | 1564 | 5.08% |

## Variant Types

| Variant type | Count | Percent |
|---|---:|---:|
| initiation_codon | 18 | 0.04% |
| missense | 35241 | 74.12% |
| nonsense | 2397 | 5.04% |
| synonymous | 9891 | 20.80% |

## Most Frequent Criteria

| Criterion | Count | Percent of all variants |
|---|---:|---:|
| PM2_Supporting | 43166 | 90.79% |
| BP1 | 34855 | 73.31% |
| BS3 | 2528 | 5.32% |
| PVS1 | 2327 | 4.89% |
| PM5_PTC | 2241 | 4.71% |
| BP4 | 1404 | 2.95% |
| BP7 | 1404 | 2.95% |
| PP3 | 1060 | 2.23% |
| PS3 | 699 | 1.47% |
| BP5 | 514 | 1.08% |
| BS1_Supporting | 141 | 0.30% |
| BS1_Strong | 101 | 0.21% |
| PP4 | 58 | 0.12% |
| BA1 | 48 | 0.10% |
| PS1 | 7 | 0.01% |

## SpliceAI Signal

Variants with reference-transcript SpliceAI score >= 0.20: 1205

| Class | Label | Count |
|---|---:|---:|
| 1 | Benign | 2 |
| 2 | Likely Benign | 126 |
| 3 | VUS | 866 |
| 4 | Likely Pathogenic | 75 |
| 5 | Pathogenic | 136 |

The full list is in `tables/high_spliceai_variants_ge_0_20.csv`.

## VUS Attention List

VUS variants flagged for manual review: 1479

The list is in `tables/vus_attention_list.csv`. Reasons include high SpliceAI, PVS1 still ending as VUS, functional evidence still ending as VUS, or PM5_PTC still ending as VUS.

## gnomAD Status

| gnomAD status | Count | Percent |
|---|---:|---:|
| absent_v2_only | 43641 | 91.78% |
| found | 3906 | 8.22% |

Interpret this cautiously: the current overlay uses the available local frequency lookup state. It is useful for spotting patterns, but it should not be treated as a complete population-frequency audit.

## Output Tables

- `tables/class_distribution.csv`
- `tables/gene_by_class.csv`
- `tables/gene_normalized_class_distribution.csv`
- `tables/gene_normalized_grouped_class_distribution.csv`
- `tables/gene_normalized_variant_type_distribution.csv`
- `tables/gene_normalized_criteria_counts.csv`
- `tables/gene_normalized_spliceai_bins.csv`
- `tables/gene_normalized_gnomad_coarse_status.csv`
- `tables/variant_type_by_class.csv`
- `tables/criteria_counts.csv`
- `tables/criteria_by_class.csv`
- `tables/spliceai_bins_by_class.csv`
- `tables/gnomad_coarse_status.csv`
- `tables/high_spliceai_variants_ge_0_20.csv`
- `tables/vus_attention_list.csv`
- `tables/likely_pathogenic_pathogenic_variants.csv`

## Next Analytical Questions

1. Compare VUS with high SpliceAI against exon/intron boundary context and known ENIGMA splice rules.
2. Review PVS1-only or PVS1-dominated VUS variants to see whether Table 3 combinations are behaving as expected.
3. Separate biological patterns from current evidence availability, especially for PM2 and sparse functional or curated evidence.
4. Add plots after the tabular checks are stable, for example class by CDS position, criteria heatmaps, and SpliceAI score density.
