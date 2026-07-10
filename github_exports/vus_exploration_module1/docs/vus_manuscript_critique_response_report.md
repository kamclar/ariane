# VUS Manuscript Critique Response Analysis

Generated: 2026-06-22 09:26

Input: `Exploratory_analysis\precomputed_classification\tables\vus_prioritization\vus_priority_annotated.csv`

## Purpose

This analysis responds to manuscript critique about possible circularity and
redundancy of the local-neighborhood signal. It asks what happens if the two
local-neighborhood priority components are removed from the VUS priority score:

- `near_pathogenic_dense`: 12 points
- `high_spliceai_in_benign_dense_region`: 8 points

This does not validate the priority score externally. It only quantifies how
much the current tiering depends on local-neighborhood context.

## Neighborhood Contribution

| metric | count | percent_of_vus |
| --- | --- | --- |
| total_vus | 8154 | 100.00 |
| vus_with_any_neighborhood_reason | 6249 | 76.64 |
| vus_with_only_neighborhood_reasons | 5018 | 61.54 |
| tier1_or_tier2_original | 859 | 10.53 |
| tier1_or_tier2_without_neighborhood_points | 286 | 3.51 |
| lost_tier1_or_tier2_without_neighborhood_points | 573 | 7.03 |

## Tier Transitions

| original_tier | tier_without_neighborhood | count |
| --- | --- | --- |
| tier1_urgent | tier2_high | 99 |
| tier2_high | tier2_high | 187 |
| tier2_high | tier3_medium | 573 |
| tier3_medium | tier3_medium | 529 |
| tier3_medium | tier4_low | 204 |
| tier4_low | tier4_low | 6562 |

## Examples Affected By Neighborhood Removal

| priority_score | score_without_neighborhood | priority_tier | tier_without_neighborhood | neighborhood_points | priority_reasons | gene | c_notation | p_notation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.2213T&gt;G | p.(Val738Gly) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.224A&gt;T | p.(Glu75Val) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.2254T&gt;G | p.(Leu752Val) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.2342A&gt;T | p.(Glu781Val) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.238A&gt;C | p.(Ser80Arg) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.251A&gt;T | p.(Glu84Val) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.263A&gt;T | p.(Lys88Ile) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.2733A&gt;T | p.(Gly911=) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.283C&gt;G | p.(Leu95Val) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.28G&gt;A | p.(Glu10Lys) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.290C&gt;T | p.(Thr97Ile) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.2980T&gt;G | p.(Cys994Gly) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3378A&gt;T | p.(Pro1126=) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3380A&gt;T | p.(Tyr1127Phe) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3384G&gt;C | p.(Leu1128=) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3384G&gt;T | p.(Leu1128=) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3399A&gt;C | p.(Leu1133Phe) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3399A&gt;G | p.(Leu1133=) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3417T&gt;A | p.(Ser1139Arg) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3418A&gt;C | p.(Ser1140Arg) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3552A&gt;T | p.(Gly1184=) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3658G&gt;T | p.(Asp1220Tyr) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3693T&gt;A | p.(Phe1231Leu) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.379A&gt;G | p.(Ser127Gly) |
| 65 | 45 | tier2_high | tier3_medium | 20 | SpliceAI&gt;=0.20;PP3;near_pathogenic_dense;high_spliceai_in_benign_dense_region | BRCA1 | c.3875C&gt;A | p.(Ser1292Tyr) |

## Bottleneck Table Note

The manuscript bottleneck table intentionally shows the largest bottleneck
groups, not all mutually exclusive categories. The complete bottleneck table is
available in `tables/vus_bottleneck/vus_bottleneck_summary.csv`.

| metric | count |
| --- | --- |
| all_bottleneck_categories_total | 8154 |
| categories_shown_in_manuscript_table | 8111 |
| other_categories_not_shown_in_short_table | 43 |

## Interpretation

- The neighborhood signal should be described as a secondary triage and
  explanation signal, not as an independent classifier.
- Any VUS that remains tier1/tier2 after removing neighborhood points is
  prioritized mainly by variant-level criteria such as SpliceAI, boundary
  proximity, PS3, PP3, PS1, PVS1, PP4, or BS3.
- VUS that drop out of tier1/tier2 after removing neighborhood points should be
  labelled as neighborhood-driven candidates and should be interpreted more
  cautiously.
- External validation against ENIGMA, ClinVar high-review variants, EQA cases,
  or local expert review remains required before the score can be presented as
  a validated prioritization method.

## Outputs

- `tables/vus_manuscript_critique_response/neighborhood_increment_summary.csv`
- `tables/vus_manuscript_critique_response/tier_transition_without_neighborhood.csv`
- `tables/vus_manuscript_critique_response/neighborhood_affected_examples.csv`
- `tables/vus_manuscript_critique_response/bottleneck_short_table_note.csv`
- `plots/27_vus_manuscript_critique_response/neighborhood_increment_summary.svg`
- `plots/27_vus_manuscript_critique_response/tier_transition_without_neighborhood.svg`
