# VUS Bottleneck Deep Dive

Generated: 2026-06-20 09:47

Inputs:

- `tables/vus_bottleneck/vus_bottleneck_variant_examples.csv`
- `tables/vus_prioritization/vus_priority_annotated.csv`

## Why We Did This

The VUS bottleneck analysis showed three groups worth looking at more closely. The goal here is not to reclassify variants. The goal is to understand what kind of evidence is missing and which subsets are worth curator time.

We focused on:

1. `strong_pathogenic_combo_one_step_short`: VUS variants with strong pathogenic evidence that miss the exact ENIGMA combination threshold.
2. `PM2_only`: the largest VUS block, where absence from population data is the only automated evidence.
3. `BP4+BP7+PM2_Supporting`: benign splice prediction mixed with PM2 background, still ending as VUS.

## Overall Deep Buckets

| deep_bucket | parent_bottleneck | count | high_priority_count | near_pathogenic_dense_count | high_spliceai_count | top_combo |
| --- | --- | --- | --- | --- | --- | --- |
| PM2_only_pathogenic_neighborhood | PM2_only | 3692 | 0 | 3692 | 0 | PM2_Supporting |
| PM2_only_low_information_background | PM2_only | 1296 | 0 | 0 | 0 | PM2_Supporting |
| BP4_BP7_PM2_near_pathogenic_density | mixed_splice_benign_and_pathogenic | 930 | 0 | 930 | 0 | BP4+BP7+PM2_Supporting |
| PS3_PM2_needs_one_more_pathogenic_support | strong_pathogenic_combo_one_step_short | 407 | 25 | 356 | 1 | PM2_Supporting+PS3 |
| BP4_BP7_PM2_low_splice_benign_prediction | mixed_splice_benign_and_pathogenic | 340 | 0 | 0 | 0 | BP4+BP7+PM2_Supporting |
| PM2_only_splice_or_boundary_context | PM2_only | 58 | 0 | 0 | 10 | PM2_Supporting |
| PS3_only_needs_independent_support | strong_pathogenic_combo_one_step_short | 27 | 1 | 21 | 0 | PS3 |
| complex_strong_pathogenic_combo_review | strong_pathogenic_combo_one_step_short | 8 | 8 | 8 | 0 | PM2_Supporting+PS1+PS3 |
| BP4_BP7_PM2_boundary_splice_benign | mixed_splice_benign_and_pathogenic | 6 | 0 | 0 | 0 | BP4+BP7+PM2_Supporting |
| PS3_PP3_splice_or_functional_followup | strong_pathogenic_combo_one_step_short | 4 | 4 | 3 | 4 | PP3+PS3 |

## 1. One-Step-Short Pathogenic VUS

Rows reviewed: `446`

These are the most actionable pathogenic-side VUS. They often already have `PS3` and sometimes `PP3`, but they still do not satisfy the exact ENIGMA combination. This means the application is behaving conservatively: strong evidence alone, or strong plus one weak signal, is not enough.

What can be concluded:

- these are not automatic upgrades
- they are a good manual review queue
- the most useful missing evidence types are PP1, PM3, PS4, PP4, curated PS1, or another accepted supporting pathogenic item
- variants with both `PS3` and `PP3` are especially interesting because functional and computational/splice evidence point in the same direction

Top examples:

| deep_bucket | gene | c_notation | p_notation | criteria_combo | total_points | priority_score | spliceai_score | near_pathogenic_20bp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PS3_PP3_splice_or_functional_followup | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | PP3+PS3 | 5 | 90 | 0.65 | 18 |
| PS3_PP3_splice_or_functional_followup | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | PP3+PS3 | 5 | 90 | 0.4 | 15 |
| PS3_PP3_splice_or_functional_followup | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | PP3+PS3 | 5 | 90 | 0.22 | 15 |
| PS3_PP3_splice_or_functional_followup | BRCA2 | c.9218A&gt;G | p.(Asp3073Gly) | PP3+PS3 | 5 | 70 | 0.22 | 3 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.3G&gt;A | p.(Met1Ile) | PM2_Supporting+PS3 | 5 | 83 | 0.23 | 4 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.5077G&gt;T | p.(Ala1693Ser) | PM2_Supporting+PS3 | 5 | 69 | 0.1 | 20 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.5145C&gt;G | p.(Ser1715Arg) | PM2_Supporting+PS3 | 5 | 57 | 0.16 | 16 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.218T&gt;G | p.(Leu73Arg) | PM2_Supporting+PS3 | 5 | 57 | 0.16 | 12 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.2T&gt;G | p.(Met1Arg) | PM2_Supporting+PS3 | 5 | 57 | 0.15 | 4 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.5152T&gt;G | p.(Trp1718Gly) | PM2_Supporting+PS3 | 5 | 57 | 0.09 | 12 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.5153G&gt;T | p.(Trp1718Leu) | PM2_Supporting+PS3 | 5 | 57 | 0.07 | 13 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.5072C&gt;A | p.(Thr1691Lys) | PM2_Supporting+PS3 | 5 | 57 | 0.05 | 24 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.5152T&gt;C | p.(Trp1718Arg) | PM2_Supporting+PS3 | 5 | 57 | 0.04 | 12 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.79T&gt;A | p.(Cys27Ser) | PM2_Supporting+PS3 | 5 | 57 | 0.04 | 12 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA2 | c.7978T&gt;G | p.(Tyr2660Asp) | PM2_Supporting+PS3 | 5 | 57 | 0.03 | 14 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.79T&gt;C | p.(Cys27Arg) | PM2_Supporting+PS3 | 5 | 57 | 0.03 | 12 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.81T&gt;G | p.(Cys27Trp) | PM2_Supporting+PS3 | 5 | 57 | 0.03 | 12 |
| PS3_PM2_needs_one_more_pathogenic_support | BRCA1 | c.134A&gt;T | p.(Lys45Ile) | PM2_Supporting+PS3 | 5 | 57 | 0.03 | 10 |

## 2. PM2-Only VUS

Rows reviewed: `5046`

This is the largest group, but also the least informative by itself. `PM2_Supporting` means the variant is absent or very rare in the available population data. That does not tell us whether the variant is damaging.

What can be concluded:

- most PM2-only variants should stay low priority until external evidence appears
- pathogenic density alone is too broad in this full-SNV synthetic landscape and should not be treated as an action trigger
- a smaller subset becomes interesting if PM2-only also has splice/boundary context or another independent signal
- PM2-only is a data-gap group, not a biological mechanism

Top PM2-only examples selected by local context. These are not proposed upgrades; they are examples of how broad the PM2-only background remains even after adding local context:

| deep_bucket | gene | c_notation | p_notation | variant_type | priority_score | priority_reasons | cds_exon | near_pathogenic_20bp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.9120T&gt;G | p.(Val3040=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 23 | 12 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.318A&gt;C | p.(Gly106=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 3 | 8 |
| PM2_only_pathogenic_neighborhood | BRCA1 | c.548G&gt;A | p.(Gly183Glu) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 7 | 5 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.8631A&gt;T | p.(Glu2877Asp) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 19 | 7 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.8752G&gt;A | p.(Glu2918Lys) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 20 | 7 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.9256G&gt;A | p.(Gly3086Arg) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 23 | 5 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.9120T&gt;C | p.(Val3040=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 23 | 12 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.7436A&gt;T | p.(Asp2479Val) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 14 | 9 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.7804A&gt;G | p.(Arg2602Gly) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 15 | 7 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.8752G&gt;C | p.(Glu2918Gln) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 20 | 7 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.6843A&gt;C | p.(Gly2281=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 11 | 9 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.6843A&gt;G | p.(Gly2281=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 11 | 9 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.6843A&gt;T | p.(Gly2281=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 11 | 9 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.318A&gt;G | p.(Gly106=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 3 | 8 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.6844G&gt;A | p.(Glu2282Lys) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 11 | 8 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.6844G&gt;C | p.(Glu2282Gln) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 11 | 8 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.319A&gt;C | p.(Arg107=) | synonymous | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 3 | 7 |
| PM2_only_pathogenic_neighborhood | BRCA2 | c.319A&gt;G | p.(Arg107Gly) | missense | 44 | SpliceAI&gt;=0.10;within_2bp_boundary;near_pathogenic_dense | 3 | 7 |

## 3. BP4/BP7 Plus PM2

Rows reviewed: `1276`

This group has benign splice prediction (`BP4+BP7`) but also PM2 background, so the final point total usually remains VUS. This is not a contradiction in the biological sense; it is a weak benign signal mixed with weak pathogenic background.

What can be concluded:

- these are benign-leaning VUS, not likely benign by current automated rules
- one additional benign supporting item could often move them toward likely benign in the contradictory-evidence point path
- RNA evidence would be the cleanest way to strengthen the benign splice branch
- cases in pathogenic-dense regions deserve manual caution, because benign splice prediction does not exclude non-splice pathogenic mechanisms

Top BP4/BP7 plus PM2 examples:

| deep_bucket | gene | c_notation | p_notation | variant_type | priority_score | boundary_distance | near_pathogenic_20bp | spliceai_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BP4_BP7_PM2_near_pathogenic_density | BRCA1 | c.222A&gt;G | p.(Gln74=) | synonymous | 32 | 9 | 13 | 0.1 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.7971A&gt;G | p.(Lys2657=) | synonymous | 32 | 5 | 13 | 0.1 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8955T&gt;C | p.(Val2985=) | synonymous | 32 | 1 | 8 | 0.07 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8955T&gt;G | p.(Val2985=) | synonymous | 32 | 1 | 8 | 0.06 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.7804A&gt;C | p.(Arg2602=) | synonymous | 32 | 1 | 7 | 0.06 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8955T&gt;A | p.(Val2985=) | synonymous | 32 | 1 | 8 | 0.04 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.7618C&gt;T | p.(Leu2540=) | synonymous | 32 | 0 | 7 | 0.04 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.7806G&gt;A | p.(Arg2602=) | synonymous | 32 | 0 | 7 | 0.03 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.9120T&gt;A | p.(Val3040=) | synonymous | 32 | 2 | 12 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8757T&gt;A | p.(Gly2919=) | synonymous | 32 | 2 | 10 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8757T&gt;C | p.(Gly2919=) | synonymous | 32 | 2 | 10 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8757T&gt;G | p.(Gly2919=) | synonymous | 32 | 2 | 10 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.7620G&gt;C | p.(Leu2540=) | synonymous | 32 | 2 | 8 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.7803T&gt;C | p.(Tyr2601=) | synonymous | 32 | 2 | 7 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.69T&gt;C | p.(Asp23=) | synonymous | 32 | 1 | 5 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.70T&gt;C | p.(Leu24=) | synonymous | 32 | 2 | 5 | 0.02 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8952A&gt;C | p.(Ser2984=) | synonymous | 32 | 1 | 10 | 0.01 |
| BP4_BP7_PM2_near_pathogenic_density | BRCA2 | c.8952A&gt;G | p.(Ser2984=) | synonymous | 32 | 1 | 10 | 0.01 |

## Exon-Level Concentrations

| bottleneck_category | gene | cds_exon | count | high_priority_count | near_pathogenic_dense_count | max_spliceai | top_deep_bucket |
| --- | --- | --- | --- | --- | --- | --- | --- |
| strong_pathogenic_combo_one_step_short | BRCA1 | 16 | 60 | 8 | 60 | 0.40 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 1 | 34 | 5 | 24 | 0.23 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 2 | 32 | 5 | 32 | 0.04 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 17 | 8 | 4 | 8 | 0.07 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 19 | 17 | 3 | 17 | 0.03 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 18 | 42 | 2 | 42 | 0.07 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA2 | 16 | 28 | 2 | 25 | 0.17 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 3 | 18 | 2 | 18 | 0.65 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 22 | 33 | 1 | 33 | 0.01 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA2 | 17 | 25 | 1 | 16 | 0.03 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 4 | 24 | 1 | 24 | 0.16 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 15 | 23 | 1 | 23 | 0.19 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA1 | 20 | 19 | 1 | 19 | 0.01 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA2 | 23 | 8 | 1 | 2 | 0.22 | PS3_PM2_needs_one_more_pathogenic_support |
| strong_pathogenic_combo_one_step_short | BRCA2 | 12 | 1 | 1 | 1 | 0.18 | PS3_only_needs_independent_support |
| PM2_only | BRCA2 | 17 | 640 | 0 | 334 | 0.14 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 24 | 483 | 0 | 310 | 0.16 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 21 | 391 | 0 | 387 | 0.17 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 15 | 352 | 0 | 320 | 0.19 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 14 | 313 | 0 | 177 | 0.19 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 16 | 293 | 0 | 226 | 0.18 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 19 | 282 | 0 | 282 | 0.18 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 18 | 275 | 0 | 210 | 0.18 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 22 | 255 | 0 | 255 | 0.19 | PM2_only_pathogenic_neighborhood |
| PM2_only | BRCA2 | 23 | 232 | 0 | 102 | 0.19 | PM2_only_low_information_background |

## Practical Interpretation

For curation, the best next queue is not the largest group. It is the group where one credible additional evidence item can change the interpretation:

1. start with `PS3_PP3_splice_or_functional_followup` and `PS3_PM2_needs_one_more_pathogenic_support`
2. then review BP4/BP7 plus PM2 variants only when there is an additional reason beyond pathogenic density
3. keep most PM2-only variants as low-priority background unless another evidence source appears

For the application, this suggests that the UI could eventually expose a VUS explanation like: "VUS because evidence is one supporting pathogenic criterion short", "VUS because PM2 is the only evidence", or "benign-leaning VUS, needs additional benign confirmation."

## Outputs

- `tables/vus_bottleneck_deep_dive/deep_bucket_summary.csv`
- `tables/vus_bottleneck_deep_dive/deep_dive_action_queue.csv`
- `tables/vus_bottleneck_deep_dive/deep_dive_exon_summary.csv`
- `tables/vus_bottleneck_deep_dive/one_step_short_variants.csv`
- `tables/vus_bottleneck_deep_dive/pm2_only_variants.csv`
- `tables/vus_bottleneck_deep_dive/bp4_bp7_pm2_variants.csv`
- `plots/13_vus_bottleneck_deep_dive/deep_bucket_counts.svg`
- `plots/13_vus_bottleneck_deep_dive/action_queue_top_buckets.svg`
