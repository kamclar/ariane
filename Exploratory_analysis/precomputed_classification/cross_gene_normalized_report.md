# Cross-Gene Normalized Analysis

Generated: 2026-06-23T13:10:13

## Purpose

Several exploratory analyses initially used BRCA1 and BRCA2 together. That is
useful for a first overview, but raw combined counts can be misleading because
BRCA2 contributes more generated coding SNVs than BRCA1 in the current
snapshot.

This analysis adds per-gene normalized tables and plots. For full-snapshot
outputs, rates are normalized to all generated coding SNVs in the gene. For VUS
worklists, rates are normalized to VUS in the gene.

## Denominators

| Gene | All generated coding SNVs | VUS |
| --- | ---: | ---: |
| BRCA1 | 16776 | 1625 |
| BRCA2 | 30771 | 6529 |

## Outputs

- `tables/cross_gene_normalized/gene_normalized_class_distribution.csv`
- `tables/cross_gene_normalized/gene_normalized_grouped_class_distribution.csv`
- `tables/cross_gene_normalized/gene_normalized_criteria_counts.csv`
- `tables/cross_gene_normalized/gene_normalized_vus_priority_tiers.csv`
- `tables/cross_gene_normalized/gene_normalized_evidence_action_groups.csv`
- `tables/cross_gene_normalized/gene_normalized_priority_queues.csv`
- `plots/31_cross_gene_normalized/gene_normalized_grouped_class_per_1000.svg`
- `plots/31_cross_gene_normalized/gene_normalized_priority_queues_per_1000_vus.svg`
- `plots/31_cross_gene_normalized/gene_normalized_evidence_action_groups_per_1000_vus.svg`

## Interpretation Boundary

These normalized rates compare how the generated ARIANE Module 1 map behaves in
BRCA1 versus BRCA2. They do not prove biological differences by themselves.
Some differences reflect gene length, exon structure, available structured
functional data, curated lookup coverage, and how current Module 1 criteria are
implemented.

Findlay SGE-based analyses remain BRCA1-only in the current dataset and should
not be interpreted as a cross-gene comparison.
