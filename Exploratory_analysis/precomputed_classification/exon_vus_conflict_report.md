# Exon, VUS Priority, and Conflict Analysis

Generated: 2026-06-23 14:40

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

This report is based on the automated Module 1 coding SNV snapshot. It is a
descriptive exon-level audit of generated class mix, selected Module 1 signals,
and sanity-check patterns. It is not final expert classification and it is not
a validated manual-review prioritization criterion.

## Outputs

- `tables/exon_vus_conflict/exon_level_summary.csv`
- `tables/exon_vus_conflict/vus_priority_list.csv`
- `tables/exon_vus_conflict/conflict_sanity_checks.csv`
- `plots/06_exon_vus_conflict/brca1_exon_group_heatmap.svg`
- `plots/06_exon_vus_conflict/brca2_exon_group_heatmap.svg`
- `plots/06_exon_vus_conflict/brca1_exon_signal_heatmap.svg`
- `plots/06_exon_vus_conflict/brca2_exon_signal_heatmap.svg`
- `plots/06_exon_vus_conflict/heuristic_vus_score_examples.svg`
- `plots/06_exon_vus_conflict/conflict_sanity_check_patterns.svg`

## Summary

- Exon summary rows: 48
- VUS variants with heuristic score: 8154
- VUS variants with priority score >= 50: 859
- Conflict sanity-check rows: 243

## Conflict Types

- benign_with_pathogenic_evidence: 129
- benign_with_high_spliceai: 128
- vus_high_spliceai_far_from_boundary: 113
- benign_high_spliceai_at_boundary: 15
- vus_with_pvs1: 1

## Top Exons by VUS Percent

| gene | cds_exon | cds_start | cds_end | n_variants | vus_count | vus_percent |
| --- | --- | --- | --- | --- | --- | --- |
| BRCA2 | 20 | 8633 | 8754 | 366 | 334 | 91.26 |
| BRCA2 | 24 | 9257 | 9501 | 735 | 670 | 91.16 |
| BRCA2 | 23 | 9118 | 9256 | 417 | 379 | 90.89 |
| BRCA2 | 15 | 7618 | 7805 | 564 | 511 | 90.60 |
| BRCA2 | 18 | 8332 | 8487 | 468 | 424 | 90.60 |
| BRCA2 | 16 | 7806 | 7976 | 513 | 464 | 90.45 |
| BRCA2 | 17 | 7977 | 8331 | 1065 | 961 | 90.23 |
| BRCA2 | 22 | 8954 | 9117 | 492 | 440 | 89.43 |
| BRCA2 | 21 | 8755 | 8953 | 597 | 533 | 89.28 |
| BRCA2 | 14 | 7436 | 7617 | 546 | 482 | 88.28 |

## Top Exons by High SpliceAI Percent

| gene | cds_exon | cds_start | cds_end | n_variants | high_spliceai_count | high_spliceai_percent |
| --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | 4 | 213 | 301 | 267 | 84 | 31.46 |
| BRCA1 | 3 | 135 | 212 | 234 | 59 | 25.21 |
| BRCA1 | 8 | 594 | 670 | 231 | 51 | 22.08 |
| BRCA2 | 22 | 8954 | 9117 | 492 | 105 | 21.34 |
| BRCA1 | 15 | 4987 | 5074 | 264 | 49 | 18.56 |
| BRCA1 | 12 | 4358 | 4484 | 381 | 66 | 17.32 |
| BRCA2 | 23 | 9118 | 9256 | 417 | 60 | 14.39 |
| BRCA2 | 6 | 517 | 631 | 345 | 34 | 9.86 |
| BRCA2 | 4 | 426 | 475 | 150 | 14 | 9.33 |
| BRCA2 | 15 | 7618 | 7805 | 564 | 52 | 9.22 |

## Highest Heuristic VUS Score Examples

| priority_score | priority_reasons | gene | c_notation | p_notation | cds_exon | boundary_distance | spliceai_score |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 90 | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | 3 | 22 | 0.650 |
| 90 | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | 15 | 30 | 0.570 |
| 90 | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | 16 | 29 | 0.400 |
| 90 | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5453A&gt;G | p.(Asp1818Gly) | 21 | 14 | 0.350 |
| 90 | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | 16 | 29 | 0.220 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.1909G&gt;T | p.(Gly637Cys) | 9 | 0 | 0.970 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.4185G&gt;T | p.(Gln1395His) | 10 | 0 | 0.960 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.4185G&gt;C | p.(Gln1395His) | 10 | 0 | 0.950 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.1909G&gt;C | p.(Gly637Arg) | 9 | 0 | 0.950 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.6937G&gt;T | p.(Gly2313Cys) | 11 | 0 | 0.950 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.7435G&gt;T | p.(Asp2479Tyr) | 13 | 0 | 0.950 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.475G&gt;T | p.(Val159Leu) | 4 | 0 | 0.940 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.477G&gt;A | p.(Val159=) | 5 | 1 | 0.940 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.6937G&gt;A | p.(Gly2313Ser) | 11 | 0 | 0.940 |
| 85 | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.6937G&gt;C | p.(Gly2313Arg) | 11 | 0 | 0.940 |

## Top Conflict Sanity Checks

| conflict_reasons | gene | c_notation | p_notation | cds_exon | boundary_distance | spliceai_score | predicted_class |
| --- | --- | --- | --- | --- | --- | --- | --- |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.192T&gt;C | p.(Cys64=) | 3 | 20 | 0.870 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.193A&gt;C | p.(Lys65Gln) | 3 | 19 | 0.870 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.195G&gt;C | p.(Lys65Asn) | 3 | 17 | 0.870 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.195G&gt;A | p.(Lys65=) | 3 | 17 | 0.860 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.195G&gt;T | p.(Lys65Asn) | 3 | 17 | 0.860 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.194A&gt;G | p.(Lys65Arg) | 3 | 18 | 0.840 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.194A&gt;T | p.(Lys65Met) | 3 | 18 | 0.800 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.162G&gt;T | p.(Gln54His) | 3 | 27 | 0.700 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.187T&gt;G | p.(Leu63Val) | 3 | 25 | 0.700 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.5270A&gt;T | p.(Asp1757Val) | 18 | 7 | 0.690 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.194A&gt;C | p.(Lys65Thr) | 3 | 18 | 0.670 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.256C&gt;G | p.(Leu86Val) | 4 | 43 | 0.670 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.84G&gt;A | p.(Leu28=) | 2 | 3 | 0.650 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.165G&gt;T | p.(Lys55Asn) | 3 | 30 | 0.550 | 2 |
| benign_with_high_spliceai;benign_with_pathogenic_evidence | BRCA1 | c.231G&gt;T | p.(Thr77=) | 4 | 18 | 0.540 | 2 |
