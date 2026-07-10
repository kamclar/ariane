# Criteria Sanity Audit

Generated: 2026-06-20 08:26

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

This report expands the criteria co-occurrence sanity checks. It asks whether flagged combinations are likely to be real rule conflicts, expected borderline outcomes, or useful candidates for manual review. It does not add or change ACMG/ENIGMA criteria.

## What Was Checked

The audit focuses on automatically detectable patterns:

- pathogenic and benign evidence in the same variant, excluding `PM2_Supporting` as a conflict driver
- `PP3` together with `BS3`
- `PVS1` variants that remain VUS
- `PS3` variants that remain VUS
- any `BS3` variant that would still classify as pathogenic

`PM2_Supporting` is intentionally treated as background context for conflict detection because it is very common in this synthetic SNV landscape. Counting it as a conflict driver would falsely flag many expected combinations such as `BP1+PM2_Supporting`.

## Main Result

Most sanity hits are not software errors. The largest group is `PS3_but_VUS`: variants with strong pathogenic evidence from the automated rule set, but without enough additional evidence to cross the ENIGMA Table 3 combination threshold. This is an expected Module 1 limitation and a useful triage group.

The most review-worthy groups are:

- the single `PVS1_but_VUS` row, because PVS1 normally contributes very strong evidence but can still remain VUS if the Table 3 combination is not met
- `PP3_with_BS3`, because a computational pathogenic signal and functional benign evidence point in opposite directions
- `pathogenic_and_benign_evidence`, after excluding PM2 background, because these are genuine mixed-direction evidence cases

## Summary By Pattern

| sanity_reason | count | brca1_count | brca2_count | benign_count | vus_count | pathogenic_count | top_variant_type | median_spliceai | max_spliceai |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PS3_but_VUS | 448 | 344 | 104 | 0 | 448 | 0 | missense | 0.01 | 0.65 |
| pathogenic_and_benign_evidence | 140 | 128 | 12 | 129 | 10 | 1 | missense | 0.28 | 0.99 |
| PP3_with_BS3 | 128 | 123 | 5 | 128 | 0 | 0 | missense | 0.28 | 0.99 |
| PVS1_but_VUS | 1 | 0 | 1 | 0 | 1 | 0 | nonsense | 0.01 | 0.01 |

## Interpretation Buckets

| bucket | count |
| --- | --- |
| expected_vus_single_strong_evidence | 434 |
| benign_functional_overrides_pp3_signal | 128 |
| vus_with_ps3_and_extra_context | 14 |
| mixed_direction_evidence | 9 |
| manual_review_single_pvs1_case | 1 |

## Top Detailed Pattern Splits

| sanity_reason | gene | variant_type | criteria_combo | group | count |
| --- | --- | --- | --- | --- | --- |
| PS3_but_VUS | BRCA1 | missense | PM2_Supporting+PS3 | vus | 312 |
| PS3_but_VUS | BRCA2 | missense | PM2_Supporting+PS3 | vus | 86 |
| pathogenic_and_benign_evidence | BRCA1 | missense | BS3+PM2_Supporting+PP3 | benign | 83 |
| PP3_with_BS3 | BRCA1 | missense | BS3+PM2_Supporting+PP3 | benign | 83 |
| pathogenic_and_benign_evidence | BRCA1 | synonymous | BS3+PM2_Supporting+PP3 | benign | 35 |
| PP3_with_BS3 | BRCA1 | synonymous | BS3+PM2_Supporting+PP3 | benign | 35 |
| PS3_but_VUS | BRCA2 | missense | PS3 | vus | 16 |
| PS3_but_VUS | BRCA1 | missense | PS3 | vus | 11 |
| PS3_but_VUS | BRCA1 | initiation_codon | PM2_Supporting+PS3 | vus | 7 |
| PS3_but_VUS | BRCA1 | missense | PM2_Supporting+PS1+PS3 | vus | 4 |
| pathogenic_and_benign_evidence | BRCA1 | missense | BS3+PP3 | benign | 4 |
| PP3_with_BS3 | BRCA1 | missense | BS3+PP3 | benign | 4 |
| PS3_but_VUS | BRCA1 | missense | PP3+PS3 | vus | 3 |
| pathogenic_and_benign_evidence | BRCA2 | missense | BS1_Supporting+PP3 | vus | 3 |
| PS3_but_VUS | BRCA1 | synonymous | PM2_Supporting+PS3 | vus | 2 |
| PS3_but_VUS | BRCA1 | missense | PM2_Supporting+PP4+PS3 | vus | 2 |
| pathogenic_and_benign_evidence | BRCA1 | missense | BP5+PM2_Supporting+PP3+PS3 | vus | 2 |
| PS3_but_VUS | BRCA1 | missense | BP5+PM2_Supporting+PP3+PS3 | vus | 2 |
| pathogenic_and_benign_evidence | BRCA2 | missense | BP5+BS3+PP3 | benign | 2 |
| PP3_with_BS3 | BRCA2 | missense | BP5+BS3+PP3 | benign | 2 |
| pathogenic_and_benign_evidence | BRCA1 | synonymous | BS3+PP3 | benign | 1 |
| PP3_with_BS3 | BRCA1 | synonymous | BS3+PP3 | benign | 1 |
| pathogenic_and_benign_evidence | BRCA1 | missense | BP5+PM2_Supporting+PP3 | vus | 1 |
| pathogenic_and_benign_evidence | BRCA1 | nonsense | BS1_Supporting+PM5_PTC+PVS1 | pathogenic | 1 |
| PS3_but_VUS | BRCA1 | missense | PP4+PS3 | vus | 1 |

## PS3 But VUS Examples

These are mostly expected VUS outcomes: `PS3` is strong evidence, but alone or with only `PM2_Supporting` it does not satisfy the pathogenic combination rule.

| gene | c_notation | p_notation | variant_type | criteria_combo | total_points | predicted_class | spliceai_score | interpretation_bucket |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | c.190T&gt;G | p.(Cys64Gly) | missense | PP3+PS3 | 5 | 3 | 0.65 | vus_with_ps3_and_extra_context |
| BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | missense | BP5+PM2_Supporting+PP3+PS3 | 2 | 3 | 0.57 | vus_with_ps3_and_extra_context |
| BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | missense | PP3+PS3 | 5 | 3 | 0.4 | vus_with_ps3_and_extra_context |
| BRCA1 | c.5453A&gt;G | p.(Asp1818Gly) | missense | BP5+PM2_Supporting+PP3+PS3 | 4 | 3 | 0.35 | vus_with_ps3_and_extra_context |
| BRCA1 | c.3G&gt;A | p.(Met1Ile) | initiation_codon | PM2_Supporting+PS3 | 5 | 3 | 0.23 | expected_vus_single_strong_evidence |
| BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | missense | PP3+PS3 | 5 | 3 | 0.22 | vus_with_ps3_and_extra_context |
| BRCA2 | c.9218A&gt;G | p.(Asp3073Gly) | missense | PP3+PS3 | 5 | 3 | 0.22 | vus_with_ps3_and_extra_context |
| BRCA1 | c.53T&gt;G | p.(Met18Arg) | missense | PM2_Supporting+PS3 | 5 | 3 | 0.19 | expected_vus_single_strong_evidence |
| BRCA1 | c.190T&gt;C | p.(Cys64Arg) | missense | PM2_Supporting+PS3 | 5 | 3 | 0.19 | expected_vus_single_strong_evidence |
| BRCA1 | c.5056C&gt;A | p.(His1686Asn) | missense | PM2_Supporting+PS3 | 5 | 3 | 0.19 | expected_vus_single_strong_evidence |
| BRCA1 | c.5096G&gt;T | p.(Arg1699Leu) | missense | PM2_Supporting+PS3 | 5 | 3 | 0.19 | expected_vus_single_strong_evidence |
| BRCA1 | c.5057A&gt;C | p.(His1686Pro) | missense | PM2_Supporting+PS3 | 5 | 3 | 0.18 | expected_vus_single_strong_evidence |

## PP3 With BS3 Examples

These are useful manual-review candidates because the splice/computational signal (`PP3`) and functional benign evidence (`BS3`) point in opposite directions. In the current point system the benign evidence usually dominates.

| gene | c_notation | p_notation | variant_type | criteria_combo | total_points | predicted_class | spliceai_score | interpretation_bucket |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | c.5196T&gt;G | p.(His1732Gln) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.99 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.5470A&gt;G | p.(Ile1824Val) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.93 | benign_functional_overrides_pp3_signal |
| BRCA2 | c.6842G&gt;T | p.(Gly2281Val) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.92 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.192T&gt;C | p.(Cys64=) | synonymous | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.87 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.193A&gt;C | p.(Lys65Gln) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.87 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.195G&gt;C | p.(Lys65Asn) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.87 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.195G&gt;A | p.(Lys65=) | synonymous | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.86 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.195G&gt;T | p.(Lys65Asn) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.86 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.194A&gt;G | p.(Lys65Arg) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.84 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.194A&gt;T | p.(Lys65Met) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.8 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.5073A&gt;G | p.(Thr1691=) | synonymous | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.71 | benign_functional_overrides_pp3_signal |
| BRCA1 | c.5196T&gt;A | p.(His1732Gln) | missense | BS3+PM2_Supporting+PP3 | -2 | 2 | 0.71 | benign_functional_overrides_pp3_signal |

## PVS1 But VUS

This is the single highest-priority sanity row to inspect manually.

| gene | c_notation | p_notation | variant_type | criteria_combo | total_points | predicted_class | spliceai_score | interpretation_bucket |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA2 | c.9925G&gt;T | p.(Glu3309Ter) | nonsense | BS3+PM2_Supporting+PVS1 | 5 | 3 | 0.01 | manual_review_single_pvs1_case |

## Practical Takeaways

1. `PS3_but_VUS` should be treated as a triage label, not as an error label. It identifies variants where Module 1 found meaningful evidence but the rule combination is still incomplete.
2. `PP3_with_BS3` is the best conflict set for curator review, especially where SpliceAI is high.
3. `PVS1_but_VUS` should be manually checked because it tests whether the PVS1 rule strength and final combination logic are behaving as intended.
4. No `BS3_but_pathogenic` pattern was observed, which is reassuring for the current automated evidence combination.

## Follow-up Case Audit

The single `PVS1_but_VUS` row was reviewed separately in `pvs1_but_vus_case_audit.md`. The short conclusion is that `BRCA2 c.9925G>T p.(Glu3309Ter)` looks like an expected edge case rather than an implementation bug: `PVS1 Very Strong` is applied because the stop is exactly at the BRCA2 p.3309 critical boundary, but `BS3 Strong` comes from Table 9 functional evidence and creates contradictory evidence. With `PM2_Supporting`, the point total is 5, which remains class 3 VUS in the contradictory-evidence path.

## Outputs

- `tables/criteria_sanity_audit/sanity_reason_summary.csv`
- `tables/criteria_sanity_audit/sanity_reason_gene_type_combo.csv`
- `tables/criteria_sanity_audit/sanity_spliceai_bins.csv`
- `tables/criteria_sanity_audit/sanity_variant_level_audit.csv`
- `tables/criteria_sanity_audit/sanity_examples_PS3_but_VUS.csv`
- `tables/criteria_sanity_audit/sanity_examples_PP3_with_BS3.csv`
- `tables/criteria_sanity_audit/sanity_examples_PVS1_but_VUS.csv`
- `tables/criteria_sanity_audit/pvs1_but_vus_case_audit.csv`
- `plots/11_criteria_sanity_audit/sanity_reason_counts.svg`
- `plots/11_criteria_sanity_audit/sanity_interpretation_buckets.svg`
