# Criteria Co-occurrence Analysis

Generated: 2026-06-20 08:18

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

This report describes which automated Module 1 criteria occur together and how their combinations relate to benign, VUS, and pathogenic grouped outcomes. It is descriptive and does not add new ACMG/ENIGMA criteria.

## Outputs

- `tables/criteria_cooccurrence/criteria_counts_by_group.csv`
- `tables/criteria_cooccurrence/criteria_pair_cooccurrence.csv`
- `tables/criteria_cooccurrence/criteria_combo_outcomes.csv`
- `tables/criteria_cooccurrence/vus_holding_combinations.csv`
- `tables/criteria_cooccurrence/pathogenic_dominant_combinations.csv`
- `tables/criteria_cooccurrence/criteria_sanity_checks.csv`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_all.svg`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_benign.svg`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_vus.svg`
- `plots/10_criteria_cooccurrence/criteria_pair_heatmap_pathogenic.svg`
- `plots/10_criteria_cooccurrence/top_vus_holding_combinations.svg`
- `plots/10_criteria_cooccurrence/top_pathogenic_combinations.svg`

## Most Frequent Criteria

| criterion | count | percent | benign_percent | vus_percent | pathogenic_percent |
| --- | --- | --- | --- | --- | --- |
| BA1 | 48 | 0.10 | 0.13 | 0.00 | 0.00 |
| BS1_Strong | 101 | 0.21 | 0.27 | 0.02 | 0.00 |
| BS1_Supporting | 141 | 0.30 | 0.35 | 0.13 | 0.04 |
| BS3 | 2528 | 5.32 | 6.38 | 2.11 | 0.00 |
| BP1 | 34855 | 73.31 | 94.33 | 0.00 | 0.00 |
| BP4 | 1404 | 2.95 | 0.35 | 15.65 | 0.00 |
| BP5 | 514 | 1.08 | 1.32 | 0.32 | 0.00 |
| BP7 | 1404 | 2.95 | 0.35 | 15.65 | 0.00 |
| PM2_Supporting | 43166 | 90.79 | 90.13 | 92.53 | 94.84 |
| PM5_PTC | 2241 | 4.71 | 0.00 | 0.00 | 91.81 |
| PP3 | 1060 | 2.23 | 0.35 | 10.47 | 3.20 |
| PP4 | 58 | 0.12 | 0.00 | 0.07 | 2.13 |
| PS1 | 7 | 0.01 | 0.00 | 0.06 | 0.08 |
| PS3 | 699 | 1.47 | 0.00 | 5.49 | 10.24 |
| PVS1 | 2327 | 4.89 | 0.00 | 0.01 | 95.29 |

## Top Criterion Pairs

| criterion_a | criterion_b | count | benign_count | vus_count | pathogenic_count |
| --- | --- | --- | --- | --- | --- |
| BP1 | PM2_Supporting | 31463 | 31463 | 0 | 0 |
| PM5_PTC | PVS1 | 2241 | 0 | 0 | 2241 |
| PM2_Supporting | PVS1 | 2213 | 0 | 1 | 2212 |
| PM2_Supporting | PM5_PTC | 2129 | 0 | 0 | 2129 |
| BS3 | PM2_Supporting | 2114 | 2112 | 2 | 0 |
| BP4 | BP7 | 1404 | 128 | 1276 | 0 |
| BP4 | PM2_Supporting | 1276 | 0 | 1276 | 0 |
| BP7 | PM2_Supporting | 1276 | 0 | 1276 | 0 |
| PM2_Supporting | PP3 | 984 | 119 | 795 | 70 |
| PM2_Supporting | PS3 | 645 | 0 | 416 | 229 |
| BS3 | BP1 | 452 | 452 | 0 | 0 |
| BP1 | BP5 | 415 | 415 | 0 | 0 |
| PS3 | PVS1 | 147 | 0 | 0 | 147 |
| BS3 | BP5 | 145 | 145 | 0 | 0 |
| BP5 | PM2_Supporting | 144 | 130 | 14 | 0 |

## VUS-holding Combinations

Combinations with at least 5 variants and at least 80% VUS.

| criteria_combo | count | vus_count | vus_percent_of_combo | benign_count | pathogenic_count |
| --- | --- | --- | --- | --- | --- |
| PM2_Supporting | 5046 | 5046 | 100.00 | 0 | 0 |
| BP4+BP7+PM2_Supporting | 1276 | 1276 | 100.00 | 0 | 0 |
| PM2_Supporting+PP3 | 791 | 791 | 100.00 | 0 | 0 |
| PM2_Supporting+PS3 | 407 | 407 | 100.00 | 0 | 0 |
| none | 331 | 331 | 100.00 | 0 | 0 |
| BS3 | 170 | 170 | 100.00 | 0 | 0 |
| PP3 | 51 | 51 | 100.00 | 0 | 0 |
| PS3 | 27 | 27 | 100.00 | 0 | 0 |
| BS1_Supporting | 8 | 8 | 100.00 | 0 | 0 |

## Pathogenic-dominant Combinations

Combinations with at least 5 variants and at least 80% pathogenic grouped outcome.

| criteria_combo | count | pathogenic_count | pathogenic_percent_of_combo | vus_count | benign_count |
| --- | --- | --- | --- | --- | --- |
| PM2_Supporting+PM5_PTC+PVS1 | 2010 | 2010 | 100.00 | 0 | 0 |
| PM2_Supporting+PM5_PTC+PS3+PVS1 | 119 | 119 | 100.00 | 0 | 0 |
| PM5_PTC+PVS1 | 100 | 100 | 100.00 | 0 | 0 |
| PM2_Supporting+PVS1 | 68 | 68 | 100.00 | 0 | 0 |
| PM2_Supporting+PP3+PS3 | 61 | 61 | 100.00 | 0 | 0 |
| PM2_Supporting+PS3+PVS1 | 15 | 15 | 100.00 | 0 | 0 |
| PM5_PTC+PS3+PVS1 | 11 | 11 | 100.00 | 0 | 0 |
| PM2_Supporting+PP3+PP4+PS3 | 5 | 5 | 100.00 | 0 | 0 |
| PM2_Supporting+PP4+PS3 | 30 | 27 | 90.00 | 3 | 0 |
| PP4+PS3 | 5 | 4 | 80.00 | 1 | 0 |

## Sanity-check Patterns

- PS3_but_VUS: 448
- pathogenic_and_benign_evidence: 140
- PP3_with_BS3: 128
- PVS1_but_VUS: 1

Top sanity-check rows:

| sanity_reason | gene | c_notation | p_notation | criteria_combo | predicted_class | spliceai_score |
| --- | --- | --- | --- | --- | --- | --- |
| PS3_but_VUS | BRCA1 | c.1A&gt;C | p.(Met1Leu) | PM2_Supporting+PS3 | 3 | 0.02 |
| PS3_but_VUS | BRCA1 | c.1A&gt;T | p.(Met1Leu) | PM2_Supporting+PS3 | 3 | 0.01 |
| PS3_but_VUS | BRCA1 | c.2T&gt;A | p.(Met1Lys) | PM2_Supporting+PS3 | 3 | 0.01 |
| PS3_but_VUS | BRCA1 | c.2T&gt;C | p.(Met1Thr) | PM2_Supporting+PS3 | 3 | 0.02 |
| PS3_but_VUS | BRCA1 | c.2T&gt;G | p.(Met1Arg) | PM2_Supporting+PS3 | 3 | 0.15 |
| PS3_but_VUS | BRCA1 | c.3G&gt;A | p.(Met1Ile) | PM2_Supporting+PS3 | 3 | 0.23 |
| PS3_but_VUS | BRCA1 | c.3G&gt;C | p.(Met1Ile) | PM2_Supporting+PS3 | 3 | 0.01 |
| PS3_but_VUS | BRCA1 | c.4G&gt;C | p.(Asp2His) | PM2_Supporting+PS3 | 3 | 0.12 |
| pathogenic_and_benign_evidence;PP3_with_BS3 | BRCA1 | c.26A&gt;C | p.(Glu9Ala) | BS3+PM2_Supporting+PP3 | 2 | 0.2 |
| PS3_but_VUS | BRCA1 | c.32T&gt;A | p.(Val11Glu) | PM2_Supporting+PS3 | 3 | 0.05 |
| PS3_but_VUS | BRCA1 | c.32T&gt;G | p.(Val11Gly) | PM2_Supporting+PS3 | 3 | 0.04 |
| PS3_but_VUS | BRCA1 | c.35A&gt;C | p.(Gln12Pro) | PM2_Supporting+PS3 | 3 | 0.02 |
| PS3_but_VUS | BRCA1 | c.40G&gt;T | p.(Val14Phe) | PM2_Supporting+PS3 | 3 | 0.02 |
| PS3_but_VUS | BRCA1 | c.41T&gt;G | p.(Val14Gly) | PM2_Supporting+PS3 | 3 | 0.03 |
| PS3_but_VUS | BRCA1 | c.44T&gt;A | p.(Ile15Asn) | PM2_Supporting+PS3 | 3 | 0.01 |
| PS3_but_VUS | BRCA1 | c.44T&gt;C | p.(Ile15Thr) | PM2_Supporting+PS3 | 3 | 0.0 |
| PS3_but_VUS | BRCA1 | c.44T&gt;G | p.(Ile15Ser) | PM2_Supporting+PS3 | 3 | 0.01 |
| PS3_but_VUS | BRCA1 | c.49G&gt;C | p.(Ala17Pro) | PM2_Supporting+PS3 | 3 | 0.01 |
| pathogenic_and_benign_evidence;PP3_with_BS3 | BRCA1 | c.50C&gt;G | p.(Ala17Gly) | BS3+PM2_Supporting+PP3 | 2 | 0.38 |
| PS3_but_VUS | BRCA1 | c.53T&gt;A | p.(Met18Lys) | PM2_Supporting+PS3 | 3 | 0.03 |

## Reading Notes

`PM2_Supporting` is common because much of this synthetic coding SNV landscape is absent from the available local population-frequency lookup. Pair counts involving `PM2_Supporting` should therefore be read as background context, not as independent biological structure. For the sanity-check table, `PM2_Supporting` is not counted as conflicting pathogenic evidence by itself; this avoids flagging expected combinations such as `BP1+PM2_Supporting` as false conflicts.
