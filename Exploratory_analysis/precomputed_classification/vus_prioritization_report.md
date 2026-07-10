# VUS Prioritization Report

Generated: 2026-06-19 21:03

Input: `Exploratory_analysis\precomputed_classification\tables\exon_vus_conflict\vus_priority_list.csv`

This is a focused manual-review layer over the automated Module 1 VUS set. It does not change classifications.

## What This Is For

This report helps prioritize VUS variants for manual curator review. It is not a classifier and it does not upgrade or downgrade any variant. The goal is to turn thousands of VUS variants into a smaller review queue by highlighting variants with stronger biological or rule-based reasons to look again.

The most useful output is `vus_priority_annotated.csv`, where every VUS gets a priority score, a tier, a review category, and explicit reasons. The `tier1_urgent` and `tier2_high` variants are the practical first-pass review list.

`Pathogenic density` means that many variants close to the VUS in coding-sequence position are already classified as likely pathogenic/pathogenic by the automated Module 1 snapshot. In this analysis, `near_pathogenic_dense` is assigned when at least 5 pathogenic-group variants occur within +/-20 coding bases of the VUS. This does not mean the VUS is pathogenic by itself. It means the local region has many pathogenic signals and deserves earlier curator attention.

Typical high-priority reasons include:

- strong or moderate SpliceAI signal
- proximity to a CDS exon boundary
- PP3, PS3, PS1, or PVS1 evidence still ending as VUS
- location in a dense pathogenic neighborhood
- high SpliceAI in a region otherwise rich in benign evidence, which is useful as a sanity check

## How This Was Derived

The starting point was the full precomputed BRCA1/BRCA2 coding SNV Module 1 snapshot. First, only variants classified as class 3 VUS by the automated snapshot were selected. Then each VUS was annotated with its CDS position, CDS exon, distance to the nearest CDS exon boundary, SpliceAI score, applied ACMG/ENIGMA Module 1 criteria, and local neighborhood context.

The priority score is a heuristic review score. Points are added for signals that make a VUS more useful to inspect manually:

- SpliceAI >= 0.20 adds strong splice-priority weight
- SpliceAI >= 0.10 adds moderate splice-priority weight
- distance <= 2 bp or <= 10 bp from a CDS exon boundary adds splice-site proximity weight
- PVS1, PS3, PS1, PP3, PP4, or BS3 add rule/evidence-based review weight
- `near_pathogenic_dense` adds weight when the local +/-20 coding-base neighborhood contains at least 5 pathogenic-group variants
- `high_spliceai_in_benign_dense_region` adds weight as a sanity-check flag when high SpliceAI appears in a locally benign-rich region

The tiers are then assigned from the total score: `tier1_urgent` for score >= 80, `tier2_high` for score >= 50, `tier3_medium` for score >= 25, and `tier4_low` below 25. Review categories are assigned from the dominant reason pattern, for example splice-boundary, splice-far-from-boundary, functional-plus-splice-prediction, or pathogenic-neighborhood.

This is intentionally conservative as a triage tool. It does not introduce new ACMG/ENIGMA criteria and it does not claim that a high-priority VUS should be reclassified. It only identifies variants that are more informative to open first.

## Outputs

- `tables/vus_prioritization/vus_priority_annotated.csv`
- `tables/vus_prioritization/vus_priority_top_by_category.csv`
- `tables/vus_prioritization/vus_priority_reason_counts.csv`
- `tables/vus_prioritization/vus_priority_category_counts.csv`
- `tables/vus_prioritization/vus_priority_exon_summary.csv`
- `plots/08_vus_prioritization/vus_priority_score_histogram.svg`
- `plots/08_vus_prioritization/vus_priority_reason_counts.svg`
- `plots/08_vus_prioritization/vus_priority_exon_barplot.svg`

## Tier Counts

| priority_tier | count |
| --- | --- |
| tier1_urgent | 99 |
| tier2_high | 760 |
| tier3_medium | 733 |
| tier4_low | 6562 |

## Review Categories

| review_category | count |
| --- | --- |
| pathogenic_neighborhood | 5740 |
| other_vus | 1805 |
| splice_far_from_boundary | 269 |
| splice_boundary | 185 |
| splice_other | 92 |
| functional_or_same_aa | 56 |
| functional_plus_splice_prediction | 6 |
| vus_with_pvs1 | 1 |

## Most Common Reasons

| reason | count |
| --- | --- |
| near_pathogenic_dense | 6080 |
| SpliceAI&gt;=0.10 | 1017 |
| SpliceAI&gt;=0.20 | 866 |
| PP3 | 854 |
| within_10bp_boundary | 846 |
| high_spliceai_in_benign_dense_region | 534 |
| PS3 | 448 |
| within_2bp_boundary | 426 |
| BS3 | 172 |
| PP4 | 6 |
| PS1 | 5 |
| PVS1 | 1 |

## Top Exons for High-priority VUS

| gene | cds_exon | vus_count | tier1_urgent | tier2_high | max_priority_score | mean_priority_score | high_spliceai_percent |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | 12 | 102 | 12 | 43 | 85 | 49.13 | 53.92 |
| BRCA1 | 11 | 233 | 10 | 2 | 85 | 16.41 | 5.15 |
| BRCA1 | 6 | 16 | 8 | 8 | 85 | 75.50 | 100.00 |
| BRCA1 | 9 | 222 | 7 | 74 | 85 | 38.58 | 36.49 |
| BRCA2 | 13 | 26 | 7 | 12 | 85 | 57.19 | 73.08 |
| BRCA2 | 4 | 15 | 7 | 5 | 85 | 66.13 | 80.00 |
| BRCA1 | 10 | 42 | 7 | 1 | 85 | 30.76 | 19.05 |
| BRCA1 | 13 | 66 | 6 | 18 | 85 | 44.24 | 36.36 |
| BRCA2 | 11 | 37 | 6 | 6 | 85 | 44.11 | 32.43 |
| BRCA2 | 9 | 28 | 5 | 14 | 85 | 54.68 | 67.86 |
| BRCA1 | 15 | 82 | 3 | 9 | 90 | 35.79 | 13.41 |
| BRCA2 | 10 | 36 | 2 | 16 | 85 | 44.72 | 50.00 |
| BRCA1 | 4 | 76 | 2 | 15 | 85 | 39.37 | 21.05 |
| BRCA2 | 2 | 151 | 2 | 9 | 85 | 19.09 | 7.28 |
| BRCA2 | 3 | 23 | 2 | 8 | 85 | 50.87 | 43.48 |

## Top VUS Variants

| priority_score | priority_tier | review_category | priority_reasons | gene | c_notation | p_notation | cds_exon | boundary_distance | spliceai_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | 3 | 22 | 0.650 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | 15 | 30 | 0.570 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | 16 | 29 | 0.400 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5453A&gt;G | p.(Asp1818Gly) | 21 | 14 | 0.350 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | 16 | 29 | 0.220 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.1909G&gt;T | p.(Gly637Cys) | 9 | 0 | 0.970 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.4185G&gt;T | p.(Gln1395His) | 10 | 0 | 0.960 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.4185G&gt;C | p.(Gln1395His) | 10 | 0 | 0.950 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.1909G&gt;C | p.(Gly637Arg) | 9 | 0 | 0.950 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.6937G&gt;T | p.(Gly2313Cys) | 11 | 0 | 0.950 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.7435G&gt;T | p.(Asp2479Tyr) | 13 | 0 | 0.950 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.475G&gt;T | p.(Val159Leu) | 4 | 0 | 0.940 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.477G&gt;A | p.(Val159=) | 5 | 1 | 0.940 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.6937G&gt;A | p.(Gly2313Ser) | 11 | 0 | 0.940 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.6937G&gt;C | p.(Gly2313Arg) | 11 | 0 | 0.940 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA2 | c.7435G&gt;C | p.(Asp2479His) | 13 | 0 | 0.940 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.443A&gt;C | p.(Gln148Pro) | 6 | 1 | 0.920 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.443A&gt;G | p.(Gln148Arg) | 6 | 1 | 0.920 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.443A&gt;T | p.(Gln148Leu) | 6 | 1 | 0.920 |
| 85 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.444G&gt;A | p.(Gln148=) | 6 | 2 | 0.920 |
