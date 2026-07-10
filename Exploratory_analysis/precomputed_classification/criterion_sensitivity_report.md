# Criterion Sensitivity Analysis

Generated: 2026-06-20 13:37

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Why We Did This

This analysis asks which automated criteria actually hold the generated
classification in place. It is not a new classifier and it does not change any
variant class.

For each variant, the analysis removes one applied criterion at a time and
re-runs the same ENIGMA combination logic on the remaining criteria. If the
class changes, that criterion is classification-sensitive for that variant.

For VUS variants, it also asks a separate hypothetical question: what would
happen if one additional accepted supporting or strong evidence item appeared?
This is a triage model only. The hypothetical evidence is not applied to the
real classification.

## Dataset

- Total variants: 47547
- VUS variants: 8154
- Single-criterion removals that changed class: 45461

## Criteria Most Often Changing Class When Removed

| criterion | applications | class_changed_when_removed | percent_changed |
| --- | --- | --- | --- |
| BP1 | 34855 | 34738 | 99.66 |
| PM2_Supporting | 43166 | 3171 | 7.35 |
| PVS1 | 2327 | 2327 | 100.00 |
| BS3 | 2528 | 2267 | 89.68 |
| PM5_PTC | 2241 | 2111 | 94.20 |
| BP5 | 514 | 246 | 47.86 |
| BP4 | 1404 | 119 | 8.48 |
| BP7 | 1404 | 119 | 8.48 |
| PS3 | 699 | 117 | 16.74 |
| PP3 | 1060 | 74 | 6.98 |
| BS1_Strong | 101 | 58 | 57.43 |
| PP4 | 58 | 52 | 89.66 |
| BA1 | 48 | 48 | 100.00 |
| BS1_Supporting | 141 | 11 | 7.80 |
| PS1 | 7 | 3 | 42.86 |

## Most Common Removal Transitions

| removed_criterion | transition | count |
| --- | --- | --- |
| BP1 | 2 Likely Benign -&gt; 3 VUS | 34140 |
| PM2_Supporting | 2 Likely Benign -&gt; 2 Likely Benign | 31270 |
| PM2_Supporting | 3 VUS -&gt; 3 VUS | 6269 |
| PVS1 | 5 Pathogenic -&gt; 3 VUS | 2257 |
| PM2_Supporting | 5 Pathogenic -&gt; 5 Pathogenic | 2126 |
| PM5_PTC | 5 Pathogenic -&gt; 4 Likely Pathogenic | 2011 |
| BS3 | 2 Likely Benign -&gt; 3 VUS | 1860 |
| PM2_Supporting | 2 Likely Benign -&gt; 3 VUS | 1708 |
| BP4 | 3 VUS -&gt; 3 VUS | 1276 |
| BP7 | 3 VUS -&gt; 3 VUS | 1276 |
| PM2_Supporting | 3 VUS -&gt; 2 Likely Benign | 1276 |
| PP3 | 3 VUS -&gt; 3 VUS | 854 |
| PS3 | 3 VUS -&gt; 3 VUS | 447 |
| BP1 | 1 Benign -&gt; 2 Likely Benign | 395 |
| BS3 | 1 Benign -&gt; 2 Likely Benign | 391 |
| PM2_Supporting | 1 Benign -&gt; 1 Benign | 328 |
| BP1 | 1 Benign -&gt; 3 VUS | 203 |
| BP5 | 1 Benign -&gt; 2 Likely Benign | 197 |
| BS3 | 3 VUS -&gt; 3 VUS | 171 |
| PS3 | 5 Pathogenic -&gt; 5 Pathogenic | 135 |

## VUS Hypothetical Evidence Impact

| scenario | transition | count |
| --- | --- | --- |
| pathogenic_supporting | 3 VUS -&gt; 3 VUS | 7558 |
| pathogenic_supporting | 3 VUS -&gt; 4 Likely Pathogenic | 414 |
| pathogenic_supporting | 3 VUS -&gt; 2 Likely Benign | 172 |
| pathogenic_supporting | 3 VUS -&gt; 5 Pathogenic | 7 |
| pathogenic_supporting | 3 VUS -&gt; 1 Benign | 3 |
| pathogenic_strong | 3 VUS -&gt; 3 VUS | 7349 |
| pathogenic_strong | 3 VUS -&gt; 4 Likely Pathogenic | 794 |
| pathogenic_strong | 3 VUS -&gt; 5 Pathogenic | 8 |
| pathogenic_strong | 3 VUS -&gt; 2 Likely Benign | 3 |
| benign_supporting | 3 VUS -&gt; 3 VUS | 6669 |
| benign_supporting | 3 VUS -&gt; 2 Likely Benign | 1477 |
| benign_supporting | 3 VUS -&gt; 4 Likely Pathogenic | 8 |
| benign_strong | 3 VUS -&gt; 2 Likely Benign | 7198 |
| benign_strong | 3 VUS -&gt; 3 VUS | 781 |
| benign_strong | 3 VUS -&gt; 1 Benign | 175 |

## VUS Combos Most Affected By Hypothetical Evidence

| scenario | original_combo | transition | count |
| --- | --- | --- | --- |
| pathogenic_supporting | PM2_Supporting | 3 VUS -&gt; 3 VUS | 5046 |
| benign_supporting | PM2_Supporting | 3 VUS -&gt; 3 VUS | 5046 |
| pathogenic_strong | PM2_Supporting | 3 VUS -&gt; 3 VUS | 5046 |
| benign_strong | PM2_Supporting | 3 VUS -&gt; 2 Likely Benign | 5046 |
| pathogenic_supporting | BP4+BP7+PM2_Supporting | 3 VUS -&gt; 3 VUS | 1276 |
| benign_supporting | BP4+BP7+PM2_Supporting | 3 VUS -&gt; 2 Likely Benign | 1276 |
| pathogenic_strong | BP4+BP7+PM2_Supporting | 3 VUS -&gt; 3 VUS | 1276 |
| benign_strong | BP4+BP7+PM2_Supporting | 3 VUS -&gt; 2 Likely Benign | 1276 |
| pathogenic_supporting | PM2_Supporting+PP3 | 3 VUS -&gt; 3 VUS | 791 |
| benign_supporting | PM2_Supporting+PP3 | 3 VUS -&gt; 3 VUS | 791 |
| pathogenic_strong | PM2_Supporting+PP3 | 3 VUS -&gt; 4 Likely Pathogenic | 791 |
| benign_strong | PM2_Supporting+PP3 | 3 VUS -&gt; 2 Likely Benign | 791 |
| pathogenic_supporting | PM2_Supporting+PS3 | 3 VUS -&gt; 4 Likely Pathogenic | 407 |
| benign_supporting | PM2_Supporting+PS3 | 3 VUS -&gt; 3 VUS | 407 |
| pathogenic_strong | PM2_Supporting+PS3 | 3 VUS -&gt; 3 VUS | 407 |
| benign_strong | PM2_Supporting+PS3 | 3 VUS -&gt; 3 VUS | 407 |
| pathogenic_supporting | none | 3 VUS -&gt; 3 VUS | 331 |
| benign_supporting | none | 3 VUS -&gt; 3 VUS | 331 |
| pathogenic_strong | none | 3 VUS -&gt; 3 VUS | 331 |
| benign_strong | none | 3 VUS -&gt; 3 VUS | 331 |
| pathogenic_supporting | BS3 | 3 VUS -&gt; 2 Likely Benign | 170 |
| benign_supporting | BS3 | 3 VUS -&gt; 2 Likely Benign | 170 |
| pathogenic_strong | BS3 | 3 VUS -&gt; 3 VUS | 170 |
| benign_strong | BS3 | 3 VUS -&gt; 1 Benign | 170 |
| pathogenic_supporting | PP3 | 3 VUS -&gt; 3 VUS | 51 |

## Example Variants Whose Class Changes After Criterion Removal

| gene | c_notation | p_notation | removed_criterion | original_class | new_class | original_points | new_points | original_combo |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | c.3418A&gt;G | p.(Ser1140Gly) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA2 | c.7469T&gt;C | p.(Ile2490Thr) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA2 | c.8460A&gt;C | p.(Val2820=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.2612C&gt;T | p.(Pro871Leu) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.4039A&gt;G | p.(Arg1347Gly) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.4535G&gt;T | p.(Ser1512Ile) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.4956G&gt;A | p.(Met1652Ile) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA2 | c.8830A&gt;T | p.(Ile2944Phe) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA2 | c.9738C&gt;T | p.(Ala3246=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.114G&gt;A | p.(Lys38=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.981A&gt;G | p.(Thr327=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.1067A&gt;G | p.(Gln356Arg) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.2077G&gt;A | p.(Asp693Asn) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.2082C&gt;T | p.(Ser694=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.2311T&gt;C | p.(Leu771=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.2458A&gt;G | p.(Lys820Glu) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.2733A&gt;G | p.(Gly911=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.3024G&gt;A | p.(Met1008Ile) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.3083G&gt;A | p.(Arg1028His) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.3113A&gt;G | p.(Glu1038Gly) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.3119G&gt;A | p.(Ser1040Asn) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.3548A&gt;G | p.(Lys1183Arg) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.4308T&gt;C | p.(Ser1436=) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.4837A&gt;G | p.(Ser1613Gly) | BA1 | 1 | 3 | -99 | 0 | BA1 |
| BRCA1 | c.5579A&gt;C | p.(His1860Pro) | BA1 | 1 | 3 | -99 | 0 | BA1 |

## Example VUS Variants Moved By Hypothetical Evidence

| scenario | gene | c_notation | p_notation | original_points | new_points | new_class | original_combo |
| --- | --- | --- | --- | --- | --- | --- | --- |
| pathogenic_supporting | BRCA1 | c.181T&gt;A | p.(Cys61Ser) | 9 | 10 | 5 | PM2_Supporting+PP4+PS3 |
| pathogenic_supporting | BRCA1 | c.5145C&gt;A | p.(Ser1715Arg) | 9 | 10 | 5 | PM2_Supporting+PS1+PS3 |
| pathogenic_supporting | BRCA2 | c.7940T&gt;C | p.(Leu2647Pro) | 9 | 10 | 5 | PM2_Supporting+PP4+PS3 |
| pathogenic_supporting | BRCA1 | c.131G&gt;C | p.(Cys44Ser) | 9 | 10 | 5 | PM2_Supporting+PS1+PS3 |
| pathogenic_supporting | BRCA1 | c.5207T&gt;G | p.(Val1736Gly) | 9 | 10 | 5 | PM2_Supporting+PP4+PS3 |
| pathogenic_supporting | BRCA1 | c.5212G&gt;C | p.(Gly1738Arg) | 9 | 10 | 5 | PM2_Supporting+PS1+PS3 |
| pathogenic_supporting | BRCA1 | c.5509T&gt;A | p.(Trp1837Arg) | 9 | 10 | 5 | PM2_Supporting+PS1+PS3 |
| pathogenic_supporting | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | 5 | 6 | 4 | PP3+PS3 |
| pathogenic_supporting | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | 5 | 6 | 4 | PP3+PS3 |
| pathogenic_supporting | BRCA1 | c.3G&gt;A | p.(Met1Ile) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | 5 | 6 | 4 | PP3+PS3 |
| pathogenic_supporting | BRCA2 | c.9218A&gt;G | p.(Asp3073Gly) | 5 | 6 | 4 | PP3+PS3 |
| pathogenic_supporting | BRCA1 | c.53T&gt;G | p.(Met18Arg) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.190T&gt;C | p.(Cys64Arg) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.5056C&gt;A | p.(His1686Asn) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.5096G&gt;T | p.(Arg1699Leu) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.5057A&gt;C | p.(His1686Pro) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.5060T&gt;A | p.(Val1687Asp) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA2 | c.7522G&gt;C | p.(Gly2508Arg) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA2 | c.7857G&gt;C | p.(Trp2619Cys) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.173C&gt;A | p.(Pro58His) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.218T&gt;G | p.(Leu73Arg) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.245T&gt;G | p.(Leu82Arg) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.5009G&gt;A | p.(Arg1670Lys) | 5 | 6 | 4 | PM2_Supporting+PS3 |
| pathogenic_supporting | BRCA1 | c.5056C&gt;G | p.(His1686Asp) | 5 | 6 | 4 | PM2_Supporting+PS3 |

## Interpretation

- Criteria with many class changes are not necessarily more important
  biologically; they are more important for the current automated decision
  boundary.
- `BP1`, `PVS1`, `PM5_PTC`, `BS3`, `PS3`, and frequency criteria are expected to
  be highly sensitive because they carry strong point weight or satisfy Table 3
  combinations.
- The VUS hypothetical section is useful for curation planning: it identifies
  which VUS groups are one accepted evidence item away from a different class.
- This analysis should not be used to add evidence that is not present in the
  real data.

## Outputs

- `tables/criterion_sensitivity/criterion_removal_summary.csv`
- `tables/criterion_sensitivity/criterion_removal_transitions.csv`
- `tables/criterion_sensitivity/criterion_removal_changed_variants.csv`
- `tables/criterion_sensitivity/vus_hypothetical_evidence_summary.csv`
- `tables/criterion_sensitivity/vus_hypothetical_combo_impact.csv`
- `tables/criterion_sensitivity/vus_hypothetical_changed_variants.csv`
- `plots/14_criterion_sensitivity/criterion_removal_changed_counts.svg`
- `plots/14_criterion_sensitivity/vus_hypothetical_evidence_impact.svg`
