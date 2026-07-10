# Frequency Baseline Audit

Generated: 2026-06-20 14:11

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Why We Did This

Budapest 2026 notes indicate that frequency calibration is moving toward
gnomAD v4.1. Before changing data sources, this audit documents how much the
current precomputed snapshot depends on the existing local frequency evidence.

This is a baseline analysis, not a gnomAD v4.1 recalculation. The local cache
currently represents a BRCA regional gnomAD v2.1.1 exome source with real
coverage, while the snapshot often marks variants as `absent_v2_only`.

## Dataset

- Total variants: 47547
- Variants with `PM2_Supporting`: 43166
- Variants with `BS1_Strong` or `BS1_Supporting`: 242
- Variants with `BA1`: 48
- gnomAD status `absent_v2_only`: 43641
- gnomAD status `found`: 3906

## gnomAD Status By Class

| gnomad_status_group | class | count |
| --- | --- | --- |
| absent_v2_only | 1 Benign | 328 |
| absent_v2_only | 2 Likely Benign | 33428 |
| absent_v2_only | 3 VUS | 7545 |
| absent_v2_only | 4 Likely Pathogenic | 136 |
| absent_v2_only | 5 Pathogenic | 2204 |
| found | 1 Benign | 474 |
| found | 2 Likely Benign | 2722 |
| found | 3 VUS | 609 |
| found | 4 Likely Pathogenic | 5 |
| found | 5 Pathogenic | 96 |

## Frequency Criteria By Class

| frequency_criterion | class | count |
| --- | --- | --- |
| BA1 | 1 Benign | 48 |
| BS1_Strong | 1 Benign | 97 |
| BS1_Strong | 2 Likely Benign | 2 |
| BS1_Strong | 3 VUS | 2 |
| BS1_Supporting | 1 Benign | 29 |
| BS1_Supporting | 2 Likely Benign | 100 |
| BS1_Supporting | 3 VUS | 11 |
| BS1_Supporting | 5 Pathogenic | 1 |
| PM2_Supporting | 1 Benign | 328 |
| PM2_Supporting | 2 Likely Benign | 32978 |
| PM2_Supporting | 3 VUS | 7545 |
| PM2_Supporting | 4 Likely Pathogenic | 136 |
| PM2_Supporting | 5 Pathogenic | 2179 |

## Frequency Criterion Combos By Group

| frequency_combo | group | count |
| --- | --- | --- |
| PM2_Supporting | benign_group | 33306 |
| PM2_Supporting | vus | 7545 |
| PM2_Supporting | pathogenic_group | 2315 |
| BS1_Supporting | benign_group | 129 |
| BS1_Strong | benign_group | 99 |
| BA1 | benign_group | 48 |
| BS1_Supporting | vus | 11 |
| BS1_Strong | vus | 2 |
| BS1_Supporting | pathogenic_group | 1 |

## PM2 Loss Simulation

This simulation removes `PM2_Supporting` from each variant that has it and
re-runs the same ENIGMA combination logic. It approximates what could happen if
future gnomAD v4.1 calibration shows that PM2 should not be applied for some
variants.

| transition | count |
| --- | --- |
| 2 Likely Benign -&gt; 2 Likely Benign | 31270 |
| 3 VUS -&gt; 3 VUS | 6269 |
| 5 Pathogenic -&gt; 5 Pathogenic | 2126 |
| 2 Likely Benign -&gt; 3 VUS | 1708 |
| 3 VUS -&gt; 2 Likely Benign | 1276 |
| 1 Benign -&gt; 1 Benign | 328 |
| 4 Likely Pathogenic -&gt; 3 VUS | 134 |
| 5 Pathogenic -&gt; 4 Likely Pathogenic | 53 |
| 4 Likely Pathogenic -&gt; 4 Likely Pathogenic | 2 |

Variants whose class changes when PM2 is removed: 3171

Examples:

| gene | c_notation | p_notation | original_class | new_class_without_pm2 | original_points | new_points_without_pm2 | original_combo |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | c.100C&gt;G | p.(Pro34Ala) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.102T&gt;A | p.(Pro34=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.102T&gt;C | p.(Pro34=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.102T&gt;G | p.(Pro34=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.103G&gt;A | p.(Val35Ile) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.105C&gt;A | p.(Val35=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.105C&gt;G | p.(Val35=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.105C&gt;T | p.(Val35=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.106T&gt;A | p.(Ser36Thr) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.106T&gt;G | p.(Ser36Ala) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.109A&gt;C | p.(Thr37Pro) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.109A&gt;G | p.(Thr37Ala) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.109A&gt;T | p.(Thr37Ser) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.10T&gt;A | p.(Ser4Thr) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.10T&gt;C | p.(Ser4Pro) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.10T&gt;G | p.(Ser4Ala) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.111A&gt;C | p.(Thr37=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.111A&gt;G | p.(Thr37=) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.112A&gt;C | p.(Lys38Gln) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.112A&gt;G | p.(Lys38Glu) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.113A&gt;C | p.(Lys38Thr) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.113A&gt;G | p.(Lys38Arg) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.113A&gt;T | p.(Lys38Met) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.114G&gt;C | p.(Lys38Asn) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |
| BRCA1 | c.114G&gt;T | p.(Lys38Asn) | 2 | 3 | -3 | -4 | BS3+PM2_Supporting |

## Interpretation

- The current snapshot is highly dependent on PM2 because most synthetic coding
  SNVs are absent from the local population cache.
- PM2-only VUS variants remain the largest low-information group.
- Loss of PM2 is most important when PM2 completes a pathogenic combination or
  participates in contradictory evidence.
- The `VUS -> Likely Benign` transitions after removing PM2 should be read in
  terms of what PM2 means. PM2 is not direct proof that a variant is damaging;
  it means the variant is absent or very rare in population data and therefore
  lacks population-based benign reassurance. In combinations such as
  `BP4+BP7+PM2_Supporting`, PM2 acts as a cautionary counterweight to weak
  benign computational/splice-neutral evidence. If PM2 is not established, only
  the benign criteria remain, and the same combination logic can move the
  variant to Likely Benign.
- A real gnomAD v4.1 analysis should not simply replace this simulation. It
  should recompute BA1, BS1, and PM2 from v4.1/v4.1.1 frequency and allele
  number data, including exome/genome discordance flags and coverage/mappability
  context.

## Outputs

- `tables/frequency_baseline/gnomad_status_by_class.csv`
- `tables/frequency_baseline/frequency_criteria_by_class.csv`
- `tables/frequency_baseline/frequency_combo_by_group.csv`
- `tables/frequency_baseline/frequency_examples.csv`
- `tables/frequency_baseline/pm2_loss_transitions.csv`
- `tables/frequency_baseline/pm2_loss_changed_variants.csv`
- `plots/15_frequency_baseline/frequency_combo_by_group.svg`
- `plots/15_frequency_baseline/pm2_loss_transitions.svg`
