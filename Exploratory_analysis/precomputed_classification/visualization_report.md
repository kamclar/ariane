# Visualization Report: Precomputed BRCA Module 1 SNV Snapshot

Generated: 2026-06-23 14:16

These plots are generated without external plotting dependencies. They use the precomputed Module 1 coding SNV snapshot and the local GRCh38 RefSeq CDS map for NM_007294.4 and NM_000059.4.

## Splice Boundary Summary

- Variants with SpliceAI >= 0.20: 1205
- Variants with SpliceAI >= 0.20 and within 10 bp of a CDS exon boundary: 381
- All variants within 10 bp of a CDS exon boundary: 3168

This uses distance to the nearest CDS exon boundary, so it is a coding-SNV approximation of splice-site proximity. It does not include intronic variants.

Overall SpliceAI >= 0.20 rate across the coding SNV snapshot: 2.53%.

| Distance to CDS exon boundary | Total variants | SpliceAI >= 0.20 | Percent | Enrichment vs overall |
|---|---:|---:|---:|---:|
| 0 | 288 | 134 | 46.53% | 18.4x |
| 1-2 | 576 | 118 | 20.49% | 8.1x |
| 3-5 | 864 | 66 | 7.64% | 3.0x |
| 6-10 | 1440 | 63 | 4.38% | 1.7x |
| 11-20 | 2874 | 184 | 6.40% | 2.5x |
| 21-50 | 6630 | 338 | 5.10% | 2.0x |
| 51-100 | 5316 | 178 | 3.35% | 1.3x |
| >100 | 29559 | 124 | 0.42% | 0.2x |

## Plots

- `plots/01_overview/class_distribution.svg`
- `plots/01_overview/grouped_class_distribution.svg`
- `plots/01_overview/criteria_by_class_heatmap.svg`
- `plots/01_overview/spliceai_bins_by_class_heatmap.svg`
- `plots/01_overview/spliceai_bins_by_group_heatmap.svg`
- `plots/02_position/brca1_cds_position_overview.svg`
- `plots/02_position/brca2_cds_position_overview.svg`
- `plots/02_position/brca1_grouped_cds_position_overview.svg`
- `plots/02_position/brca2_grouped_cds_position_overview.svg`
- `plots/03_boundary_spliceai/boundary_distance_by_class_heatmap.svg`
- `plots/03_boundary_spliceai/boundary_distance_by_group_heatmap.svg`
- `plots/03_boundary_spliceai/spliceai_vs_boundary_distance.svg`
- `plots/03_boundary_spliceai/spliceai_vs_boundary_distance_grouped_classes.svg`
- `plots/03_boundary_spliceai/spliceai_high_rate_by_boundary_grouped_classes.svg`
- `plots/03_boundary_spliceai/boundary_distance_distribution_grouped_classes.svg`
- `plots/04_clusters/brca1_position_clusters_with_spliceai.svg`
- `plots/04_clusters/brca2_position_clusters_with_spliceai.svg`
- `plots/04_clusters/brca1_position_clusters_without_spliceai.svg`
- `plots/04_clusters/brca2_position_clusters_without_spliceai.svg`
- `plots/04_clusters/brca1_position_clusters_no_splice_no_class.svg`
- `plots/04_clusters/brca2_position_clusters_no_splice_no_class.svg`
- `plots/04_clusters/cluster_profile_heatmap_with_spliceai.svg`
- `plots/04_clusters/cluster_profile_heatmap_without_spliceai.svg`
- `plots/04_clusters/cluster_profile_heatmap_no_splice_no_class.svg`
- `plots/04_clusters/brca1_cluster_profile_heatmap_with_spliceai.svg`
- `plots/04_clusters/brca2_cluster_profile_heatmap_with_spliceai.svg`
- `plots/04_clusters/brca1_cluster_profile_heatmap_without_spliceai.svg`
- `plots/04_clusters/brca2_cluster_profile_heatmap_without_spliceai.svg`
- `plots/04_clusters/brca1_cluster_profile_heatmap_no_splice_no_class.svg`
- `plots/04_clusters/brca2_cluster_profile_heatmap_no_splice_no_class.svg`

## Derived Tables

- `tables/spliceai_high_by_boundary_distance.csv`
- `tables/spliceai_boundary_distance_by_grouped_class.csv`
- `tables/brca1_position_bins.csv`
- `tables/brca2_position_bins.csv`
- `tables/brca1_grouped_position_bins.csv`
- `tables/brca2_grouped_position_bins.csv`
- `tables/brca1_position_clusters_with_spliceai.csv`
- `tables/brca2_position_clusters_with_spliceai.csv`
- `tables/position_cluster_summary_with_spliceai.csv`
- `tables/position_cluster_segments_with_spliceai.csv`
- `tables/brca1_position_clusters_without_spliceai.csv`
- `tables/brca2_position_clusters_without_spliceai.csv`
- `tables/position_cluster_summary_without_spliceai.csv`
- `tables/position_cluster_segments_without_spliceai.csv`
- `tables/brca1_position_clusters_no_splice_no_class.csv`
- `tables/brca2_position_clusters_no_splice_no_class.csv`
- `tables/position_cluster_summary_no_splice_no_class.csv`
- `tables/position_cluster_segments_no_splice_no_class.csv`
- `tables/cluster_spliceai_ablation_comparison.csv`

## Cluster Analysis

The cluster plots group 120 CDS position bins per gene using a simple k-means model with k=5. The `with_spliceai` run uses grouped class fractions, high SpliceAI fraction, mean SpliceAI, near-boundary fraction, and selected criteria fractions. The `without_spliceai` run removes high SpliceAI and mean SpliceAI from the clustering feature space, but still reports them afterward as annotations. This is exploratory pattern finding, not a validated classifier.

Clusters can recur in several non-contiguous CDS regions. Use the `position_cluster_segments_*` tables for contiguous stretches and the `cluster_profile_heatmap_*` plots for the mean feature profile of each cluster.

The stricter `no_splice_no_class` run removes predicted class fractions and SpliceAI-derived criteria from the clustering feature space. It clusters on variant type, near-boundary fraction, and non-SpliceAI criteria fractions, then reports class and SpliceAI afterward as annotations.

## Pathogenic Driver Check

For class 4/5 variants in this coding SNV snapshot:

| Signal | Count |
|---|---:|
| Total class 4/5 | 2441 |
| PP3 present | 78 |
| SpliceAI >= 0.20 | 211 |
| PVS1 present | 2326 |
| PM5_PTC present | 2241 |
| Nonsense variants | 2326 |
| Missense variants | 103 |
| Synonymous variants | 9 |
| Initiation codon variants | 3 |

In this automated coding SNV landscape, most class 4/5 calls are PVS1/PM5_PTC driven, not PP3/SpliceAI driven. SpliceAI is still a strong local signal around splice boundaries and is important for specific variants, but it is not the dominant source of pathogenic classifications across the full coding SNV snapshot.

## SpliceAI Ablation Check

| Gene | Bins compared | Same cluster | Changed cluster | Same percent |
|---|---:|---:|---:|---:|
| BRCA1 | 120 | 120 | 0 | 100.00% |
| BRCA2 | 120 | 120 | 0 | 100.00% |

In this run, removing SpliceAI from the clustering feature space did not change the bin assignments. The SpliceAI signal visible in cluster profiles should therefore be interpreted as an annotation correlated with the class/criteria/boundary structure, not as an independent driver of the k-means split.

## Reading Notes

The position overview plots use 120 bins along the coding sequence. Vertical lines are CDS exon boundaries. Darker cells mean more variants in a class or higher maximum SpliceAI signal in that bin.

The boundary-distance plot caps distance at 250 bp for readability. Points with larger distances are placed at the right edge.
