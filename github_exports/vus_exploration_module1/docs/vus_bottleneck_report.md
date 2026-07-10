# VUS Decision Bottleneck Analysis

Generated: 2026-06-23 14:35

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Purpose

This analysis asks why automated Module 1 leaves variants as VUS and what kind
of evidence type would theoretically be needed next. It is not another
classifier and it does not change any variant class. It should be read as an
evidence-pattern audit, not as a validated list of variants that must be
manually reviewed.

The main idea is to group VUS variants by decision mechanism:

- not enough pathogenic evidence to meet ENIGMA combinations
- not enough benign evidence to meet benign combinations
- contradictory benign and pathogenic evidence
- no automated evidence beyond the input variant itself

This is more useful than a simple BRCA1 versus BRCA2 comparison because it points to what is missing.

## Overall VUS Count

Total VUS rows analyzed: `8154`

## Bottleneck Summary

| bottleneck_category | count | percent_of_vus | brca1_count | brca2_count | high_spliceai_count | top_combo |
| --- | --- | --- | --- | --- | --- | --- |
| PM2_only | 5046 | 61.88 | 725 | 4321 | 10 | PM2_Supporting |
| mixed_splice_benign_and_pathogenic | 1276 | 15.65 | 90 | 1186 | 0 | BP4+BP7+PM2_Supporting |
| computational_pathogenic_evidence_not_enough | 842 | 10.33 | 286 | 556 | 842 | PM2_Supporting+PP3 |
| strong_pathogenic_combo_one_step_short | 446 | 5.47 | 342 | 104 | 5 | PM2_Supporting+PS3 |
| no_automated_evidence | 331 | 4.06 | 51 | 280 | 1 | none |
| benign_functional_evidence_not_enough | 170 | 2.08 | 123 | 47 | 0 | BS3 |
| benign_evidence_not_enough | 21 | 0.26 | 3 | 18 | 0 | BP5 |
| other_conflicting_evidence | 19 | 0.23 | 4 | 15 | 8 | BP5+PM2_Supporting |
| PM2_plus_weak_context | 2 | 0.02 | 1 | 1 | 0 | PM2_Supporting+PP4 |
| conflicting_PVS1_BS3 | 1 | 0.01 | 0 | 1 | 0 | BS3+PM2_Supporting+PVS1 |

## Most Common Blocking Combinations

| bottleneck_category | criteria_combo | count | brca1_count | brca2_count | median_points | high_spliceai_count |
| --- | --- | --- | --- | --- | --- | --- |
| PM2_only | PM2_Supporting | 5046 | 725 | 4321 | 1.0 | 10 |
| mixed_splice_benign_and_pathogenic | BP4+BP7+PM2_Supporting | 1276 | 90 | 1186 | -1.0 | 0 |
| computational_pathogenic_evidence_not_enough | PM2_Supporting+PP3 | 791 | 273 | 518 | 2 | 791 |
| strong_pathogenic_combo_one_step_short | PM2_Supporting+PS3 | 407 | 321 | 86 | 5 | 1 |
| no_automated_evidence | none | 331 | 51 | 280 | 0 | 1 |
| benign_functional_evidence_not_enough | BS3 | 170 | 123 | 47 | -4.0 | 0 |
| computational_pathogenic_evidence_not_enough | PP3 | 51 | 13 | 38 | 1 | 51 |
| strong_pathogenic_combo_one_step_short | PS3 | 27 | 11 | 16 | 4 | 0 |
| benign_evidence_not_enough | BP5 | 11 | 2 | 9 | -2 | 0 |
| other_conflicting_evidence | BP5+PM2_Supporting | 10 | 0 | 10 | -1.0 | 0 |
| benign_evidence_not_enough | BS1_Supporting | 8 | 0 | 8 | -1.0 | 0 |
| strong_pathogenic_combo_one_step_short | PM2_Supporting+PS1+PS3 | 4 | 4 | 0 | 9.0 | 0 |
| strong_pathogenic_combo_one_step_short | PP3+PS3 | 4 | 3 | 1 | 5.0 | 4 |
| strong_pathogenic_combo_one_step_short | PM2_Supporting+PP4+PS3 | 3 | 2 | 1 | 9 | 0 |
| other_conflicting_evidence | BS1_Supporting+PP3 | 3 | 0 | 3 | 0 | 3 |
| benign_evidence_not_enough | BS1_Strong | 2 | 1 | 1 | -4.0 | 0 |
| other_conflicting_evidence | BP5+PM2_Supporting+PP3 | 2 | 1 | 1 | 0.0 | 2 |
| other_conflicting_evidence | BP5+PM2_Supporting+PP3+PS3 | 2 | 2 | 0 | 3.0 | 2 |
| PM2_plus_weak_context | PM2_Supporting+PP4 | 2 | 1 | 1 | 5.0 | 0 |
| strong_pathogenic_combo_one_step_short | PP4+PS3 | 1 | 1 | 0 | 8 | 0 |
| other_conflicting_evidence | BS3+PM2_Supporting+PS1 | 1 | 1 | 0 | 1 | 0 |
| other_conflicting_evidence | BP5+PP3 | 1 | 0 | 1 | -1 | 1 |
| conflicting_PVS1_BS3 | BS3+PM2_Supporting+PVS1 | 1 | 0 | 1 | 5 | 0 |

## Example Variants From Selected Bottleneck Groups

| bottleneck_category | gene | c_notation | p_notation | variant_type | criteria_combo | total_points | spliceai_score |
| --- | --- | --- | --- | --- | --- | --- | --- |
| conflicting_PVS1_BS3 | BRCA2 | c.9925G&gt;T | p.(Glu3309Ter) | nonsense | BS3+PM2_Supporting+PVS1 | 5 | 0.01 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | missense | PP3+PS3 | 5 | 0.65 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | missense | PP3+PS3 | 5 | 0.4 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.3G&gt;A | p.(Met1Ile) | initiation_codon | PM2_Supporting+PS3 | 5 | 0.23 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | missense | PP3+PS3 | 5 | 0.22 |
| strong_pathogenic_combo_one_step_short | BRCA2 | c.9218A&gt;G | p.(Asp3073Gly) | missense | PP3+PS3 | 5 | 0.22 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.190T&gt;C | p.(Cys64Arg) | missense | PM2_Supporting+PS3 | 5 | 0.19 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5056C&gt;A | p.(His1686Asn) | missense | PM2_Supporting+PS3 | 5 | 0.19 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5096G&gt;T | p.(Arg1699Leu) | missense | PM2_Supporting+PS3 | 5 | 0.19 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.53T&gt;G | p.(Met18Arg) | missense | PM2_Supporting+PS3 | 5 | 0.19 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5057A&gt;C | p.(His1686Pro) | missense | PM2_Supporting+PS3 | 5 | 0.18 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5434C&gt;G | p.(Pro1812Ala) | missense | PS3 | 4 | 0.18 |
| strong_pathogenic_combo_one_step_short | BRCA2 | c.7007G&gt;C | p.(Arg2336Pro) | missense | PS3 | 4 | 0.18 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5060T&gt;A | p.(Val1687Asp) | missense | PM2_Supporting+PS3 | 5 | 0.17 |
| strong_pathogenic_combo_one_step_short | BRCA2 | c.7522G&gt;C | p.(Gly2508Arg) | missense | PM2_Supporting+PS3 | 5 | 0.17 |
| strong_pathogenic_combo_one_step_short | BRCA2 | c.7857G&gt;C | p.(Trp2619Cys) | missense | PM2_Supporting+PS3 | 5 | 0.17 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.173C&gt;A | p.(Pro58His) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.218T&gt;G | p.(Leu73Arg) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.245T&gt;G | p.(Leu82Arg) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5009G&gt;A | p.(Arg1670Lys) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5056C&gt;G | p.(His1686Asp) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5059G&gt;T | p.(Val1687Phe) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5062G&gt;T | p.(Val1688Phe) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5116G&gt;C | p.(Gly1706Arg) | missense | PM2_Supporting+PS3 | 5 | 0.16 |
| strong_pathogenic_combo_one_step_short | BRCA1 | c.5145C&gt;G | p.(Ser1715Arg) | missense | PM2_Supporting+PS3 | 5 | 0.16 |

## Interpretation

The largest bottleneck is expected: variants with only `PM2_Supporting` remain VUS because absence from population data is not independent pathogenic evidence. These variants mostly need outside evidence, not more arithmetic.

The second major bottleneck is benign splice prediction plus PM2 background, commonly `BP4+BP7+PM2_Supporting`. This is a VUS because the benign evidence is not strong enough and PM2 pulls in the opposite direction. RNA/functional confirmation or stronger benign evidence would be needed to resolve these.

The most actionable pathogenic-side group is `strong_pathogenic_combo_one_step_short`, usually `PM2_Supporting+PS3`. These are close to likely pathogenic: one additional supporting pathogenic criterion, or one moderate pathogenic criterion, could move many cases. This group also catches rarer high-point combinations that still miss the exact ENIGMA Table 3 combination pattern.

The conflict groups are smaller but more important for quality control. `PP3+BS3` and `PVS1+BS3` should be reviewed as evidence-adjudication cases rather than treated as automatic upgrades or downgrades.

## Outputs

- `tables/vus_bottleneck/vus_bottleneck_summary.csv`
- `tables/vus_bottleneck/vus_bottleneck_combinations.csv`
- `tables/vus_bottleneck/vus_bottleneck_variant_examples.csv`
- `plots/12_vus_bottleneck/vus_bottleneck_categories.svg`
- `plots/12_vus_bottleneck/vus_bottleneck_combinations.svg`
