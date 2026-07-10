# External Sources

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
