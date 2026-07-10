# Findlay SGE Continuous Score Vs VUS Tier Analysis

Generated: 2026-06-22 14:32

## Purpose

This analysis compares ARIANE Module 1 VUS priority tiers with the continuous
Findlay et al. 2018 BRCA1 saturation genome editing `function.score.mean`.
This is intentionally different from comparing against the derived `BS3` label:
the continuous score is used as an external functional axis where available.

Scope is limited to BRCA1 regions covered by Findlay 2018, mainly RING and
BRCT. Lack of comparable dense BRCA2 functional data is itself a limitation.

## Matched VUS Count

Matched BRCA1 VUS with Findlay SGE rows: `673`

## By Priority Tier

| priority_tier | count | median_function_score | mean_function_score | lof_count | int_count | func_count |
| --- | --- | --- | --- | --- | --- | --- |
| tier1_urgent | 9 | -1.951 | -1.893 | 6 | 3 | 0 |
| tier2_high | 46 | -1.647 | -1.595 | 30 | 15 | 1 |
| tier3_medium | 365 | -1.891 | -1.796 | 303 | 15 | 47 |
| tier4_low | 253 | -0.828 | -0.721 | 14 | 139 | 100 |

## LOF Distribution And Priority Gap

| priority_tier | lof_count | with_ps3 | with_bs3 | with_splice_priority_reason | with_boundary_priority_reason | with_pathogenic_neighborhood_reason | top_criteria_combo | top_priority_reasons |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tier1_urgent | 6 | 6 | 0 | 6 | 1 | 5 | PP3:Supporting:1;PS3:Strong:4 | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region |
| tier2_high | 30 | 29 | 0 | 5 | 24 | 29 | PM2_Supporting:Supporting:1;PS3:Strong:4 | within_2bp_boundary;PS3;near_pathogenic_dense |
| tier3_medium | 303 | 302 | 0 | 37 | 54 | 277 | PM2_Supporting:Supporting:1;PS3:Strong:4 | PS3;near_pathogenic_dense |
| tier4_low | 14 | 0 | 0 | 1 | 3 | 13 | PM2_Supporting:Supporting:1 | near_pathogenic_dense |

The Findlay SGE axis does not support a monotonic interpretation of the current
priority tiers. Tier4 is clearly more function-preserved or benign-leaning, but
tier3 contains more LOF-class variants than tier1 and tier2 combined. In this
snapshot, most tier3 LOF variants already have `PS3`; they remain tier3 because
`PS3` alone, or `PS3+PM2_Supporting` with weak context, does not reach the
tier1/tier2 review threshold unless splice, boundary, PS1, PP4, or stronger
neighborhood signals are also present.

This is a useful critique of the current triage heuristic: it prioritizes
splice/boundary/context combinations, and may under-prioritize functionally
damaging missense VUS that have PS3 but lack additional review signals.

## Tier3 LOF Overlap With Evidence Action Plan

| action_group | bottleneck_category | count | top_criteria_combo | top_priority_reasons |
| --- | --- | --- | --- | --- |
| near_pathogenic_threshold | strong_pathogenic_combo_one_step_short | 302 | PM2_Supporting:Supporting:1;PS3:Strong:4 | PS3;near_pathogenic_dense |
| background_absence_only | PM2_plus_weak_context | 1 | PM2_Supporting:Supporting:1;PP4:Strong:4 | PP4;near_pathogenic_dense |

This connects the Findlay SGE axis to the action-plan analysis. The tier3 LOF
group is not a separate story: almost all of it overlaps with the
`near_pathogenic_threshold` / `strong_pathogenic_combo_one_step_short` queue.
This means review priority and distance-to-reclassification are related but not
identical. A variant can be close to likely pathogenic under the evidence model
while remaining only tier3 in the heuristic priority score if its signal is
mainly `PS3+PM2_Supporting`.

## By Region

| findlay_region | count | median_function_score | mean_function_score | lof_count | int_count | func_count |
| --- | --- | --- | --- | --- | --- | --- |
| BRCT_region | 461 | -1.412 | -1.329 | 244 | 105 | 112 |
| N_terminal_RING_region | 211 | -1.442 | -1.494 | 109 | 67 | 35 |
| outside_Findlay_Ring_BRCT_core | 1 | -0.610 | -0.610 | 0 | 0 | 1 |

## By Neighborhood-Driven Status

| neighborhood_driven | count | median_function_score | mean_function_score | lof_count | int_count | func_count |
| --- | --- | --- | --- | --- | --- | --- |
| no | 50 | -1.492 | -1.278 | 28 | 9 | 13 |
| yes | 623 | -1.409 | -1.388 | 325 | 163 | 135 |

## LOF-Class VUS Candidates

| priority_score | priority_tier | review_category | priority_reasons | gene | c_notation | p_notation | findlay_region | findlay_function_score_mean | findlay_mean_rna_score | spliceai_score | criteria |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5453A&gt;G | p.(Asp1818Gly) | BRCT_region | -3.163 | -4.430 | 0.350 | BP5:Moderate:-2;PM2_Supporting:Supporting:1;PP3:Supporting:1;PS3:Strong:4 |
| 83 | tier1_urgent | splice_boundary | SpliceAI&gt;=0.20;within_2bp_boundary;PS3;high_spliceai_in_benign_dense_region | BRCA1 | c.3G&gt;A | p.(Met1Ile) | N_terminal_RING_region | -2.502 | -0.981 | 0.230 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | BRCT_region | -2.229 | -3.081 | 0.400 | PP3:Supporting:1;PS3:Strong:4 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | BRCT_region | -2.008 | -0.548 | 0.220 | PP3:Supporting:1;PS3:Strong:4 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | N_terminal_RING_region | -1.951 | -5.421 | 0.650 | PP3:Supporting:1;PS3:Strong:4 |
| 90 | tier1_urgent | functional_plus_splice_prediction | SpliceAI&gt;=0.20;PS3;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | BRCT_region | -1.797 | -3.923 | 0.570 | BP5:Strong:-4;PM2_Supporting:Supporting:1;PP3:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5191G&gt;A | p.(Glu1731Lys) | BRCT_region | -2.770 |  | 0.010 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5279T&gt;C | p.(Ile1760Thr) | BRCT_region | -2.605 | 0.669 | 0.020 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | functional_or_same_aa | SpliceAI&gt;=0.10;within_2bp_boundary;PS3 | BRCA1 | c.2T&gt;G | p.(Met1Arg) | N_terminal_RING_region | -2.482 | -0.643 | 0.150 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 55 | tier2_high | pathogenic_neighborhood | PS3;PS1;near_pathogenic_dense | BRCA1 | c.5509T&gt;A | p.(Trp1837Arg) | BRCT_region | -2.461 | -0.852 | 0.000 | PM2_Supporting:Supporting:1;PS1:Strong:4;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.134A&gt;T | p.(Lys45Ile) | N_terminal_RING_region | -2.436 | -0.468 | 0.030 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.79T&gt;G | p.(Cys27Gly) | N_terminal_RING_region | -2.374 | -0.563 | 0.020 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.132C&gt;G | p.(Cys44Trp) | N_terminal_RING_region | -2.307 | -0.454 | 0.020 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | SpliceAI&gt;=0.10;within_10bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.218T&gt;G | p.(Leu73Arg) | N_terminal_RING_region | -2.289 | 0.278 | 0.160 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5279T&gt;G | p.(Ile1760Ser) | BRCT_region | -2.271 | 0.114 | 0.030 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 63 | tier2_high | pathogenic_neighborhood | within_10bp_boundary;PS3;PS1;near_pathogenic_dense | BRCA1 | c.131G&gt;C | p.(Cys44Ser) | N_terminal_RING_region | -2.268 | -0.255 | 0.000 | PM2_Supporting:Supporting:1;PS1:Strong:4;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5279T&gt;A | p.(Ile1760Asn) | BRCT_region | -2.182 | -0.823 | 0.020 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.79T&gt;A | p.(Cys27Ser) | N_terminal_RING_region | -1.980 | -0.758 | 0.040 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.79T&gt;C | p.(Cys27Arg) | N_terminal_RING_region | -1.959 | -0.881 | 0.030 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.81T&gt;G | p.(Cys27Trp) | N_terminal_RING_region | -1.923 | -0.683 | 0.030 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5152T&gt;G | p.(Trp1718Gly) | BRCT_region | -1.875 | -0.521 | 0.090 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | SpliceAI&gt;=0.10;within_10bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5145C&gt;G | p.(Ser1715Arg) | BRCT_region | -1.871 | -0.633 | 0.160 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5072C&gt;A | p.(Thr1691Lys) | BRCT_region | -1.855 | -0.798 | 0.050 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 52 | tier2_high | pathogenic_neighborhood | PS3;PP4;near_pathogenic_dense | BRCA1 | c.181T&gt;A | p.(Cys61Ser) | N_terminal_RING_region | -1.841 | -0.899 | 0.080 | PM2_Supporting:Supporting:1;PP4:Strong:4;PS3:Strong:4 |
| 55 | tier2_high | pathogenic_neighborhood | PS3;PS1;near_pathogenic_dense | BRCA1 | c.5212G&gt;C | p.(Gly1738Arg) | BRCT_region | -1.798 | -0.894 | 0.000 | PM2_Supporting:Supporting:1;PS1:Strong:4;PS3:Strong:4 |
| 63 | tier2_high | pathogenic_neighborhood | within_10bp_boundary;PS3;PS1;near_pathogenic_dense | BRCA1 | c.5145C&gt;A | p.(Ser1715Arg) | BRCT_region | -1.770 | -0.708 | 0.030 | PM2_Supporting:Supporting:1;PS1:Strong:4;PS3:Strong:4 |
| 52 | tier2_high | pathogenic_neighborhood | PS3;PP4;near_pathogenic_dense | BRCA1 | c.5207T&gt;G | p.(Val1736Gly) | BRCT_region | -1.733 | -0.659 | 0.000 | PM2_Supporting:Supporting:1;PP4:Strong:4;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5152T&gt;A | p.(Trp1718Arg) | BRCT_region | -1.678 | -0.521 | 0.020 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.5152T&gt;C | p.(Trp1718Arg) | BRCT_region | -1.650 | -0.399 | 0.040 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| 57 | tier2_high | pathogenic_neighborhood | within_2bp_boundary;PS3;near_pathogenic_dense | BRCA1 | c.83T&gt;C | p.(Leu28Pro) | N_terminal_RING_region | -1.644 | -0.323 | 0.000 | PM2_Supporting:Supporting:1;PS3:Strong:4 |

## Interpretation Boundary

Findlay SGE is not used as a clinical reference set here. It is a functional
measurement axis for a restricted BRCA1 subset. It is also not fully
independent of the generated Module 1 map, because BRCA1 RING/BRCT `PS3` and
`BS3` evidence in the local ENIGMA Table 9 lookup is partly derived from the
same Findlay SGE source. The least circular subset is therefore VUS with a
Findlay score but without applied `PS3` or `BS3`.

The main supported interpretation is coarse: tier4 is function-preserved or
benign-leaning compared with the other tiers. The data do not validate the fine
ordering tier1 over tier2 over tier3, and small counts in tier1 and tier2 make
their medians unstable.

## Outputs

- `tables/findlay_sge_vus_tier/brca1_vus_with_findlay_sge_scores.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_by_priority_tier.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_by_region.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_by_neighborhood_driven.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_lof_gap_by_tier.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_tier3_lof_action_overlap.csv`
- `tables/findlay_sge_vus_tier/findlay_sge_lof_vus_candidates.csv`
- `plots/28_findlay_sge_vus_tier/findlay_sge_score_by_priority_tier.svg`
