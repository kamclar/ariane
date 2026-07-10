# Benign Counterexamples In Pathogenic Domains

Generated: 2026-06-20 17:24

Input: `Exploratory_analysis\precomputed_classification\tables\benign_structure_function\benign_structure_function_variants.csv`

## Purpose

This analysis finds benign/likely benign generated variants that fall inside
structurally or functionally important protein regions. These are
counterexamples to the tempting but unsafe shortcut: "variant is in a domain,
therefore it is pathogenic".

They are not errors by default. They are useful reminders that region/domain
context is not an ACMG/ENIGMA criterion by itself and must be interpreted with
the actual evidence, such as BS3, BP1, BP4, BP7, PS3, PP3, PS1, PVS1, or PM5.

## Counterexample Tiers

- `A_exact_or_known_residue_context`: benign variant overlaps a source-extracted
  functional annotation or a known pathogenic-residue context
- `B_pathogenic_enriched_domain`: benign variant lies in a domain that is much
  more frequent in pathogenic than benign variants after group normalization
- `C_region_level_interface_context`: benign variant lies in a UniProt
  region-level interaction annotation
- `D_broad_domain_context`: benign variant lies in a broad important domain

## Summary

Total benign counterexamples: 17955

| counterexample_tier | count |
| --- | --- |
| C_region_level_interface_context | 15199 |
| D_broad_domain_context | 2267 |
| B_pathogenic_enriched_domain | 283 |
| A_exact_or_known_residue_context | 206 |

## Feature Summary

| feature | count |
| --- | --- |
| BRC repeats / RAD51-binding region | 9242 |
| C-terminal RAD51-binding/nuclear localization context | 1834 |
| BRCT phosphopeptide-binding region | 1178 |
| RING zinc-binding/E3 ligase region | 576 |
| DNA-binding domain / helical-OB-DSS1 region | 256 |
| coiled-coil PALB2 interaction region | 22 |
| PALB2 binding N-terminal region | 8 |

## Gene Feature Summary

| gene | feature | count |
| --- | --- | --- |
| BRCA2 | BRC repeats / RAD51-binding region | 9242 |
| BRCA2 | C-terminal RAD51-binding/nuclear localization context | 1834 |
| BRCA1 | BRCT phosphopeptide-binding region | 1178 |
| BRCA1 | RING zinc-binding/E3 ligase region | 576 |
| BRCA2 | DNA-binding domain / helical-OB-DSS1 region | 256 |
| BRCA1 | coiled-coil PALB2 interaction region | 22 |
| BRCA2 | PALB2 binding N-terminal region | 8 |

## Active/Interface Summary

| curated_active_site_status | count |
| --- | --- |
| region_level:UniProt:Region:Interaction with RAD51 (1003-2082) [source_extracted_region] | 5450 |
| region_level:UniProt:Region:Interaction with RAD51 (1003-2082) [source_extracted_region];region_level:UniProt:Region:Interaction with POLH (1338-1781) [source_extracted_region] | 3783 |
| region_level:UniProt:Region:Interaction with NPM1 (639-1000) [source_extracted_region] | 3076 |
| not_curated_in_current_dataset | 2402 |
| region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 1037 |
| region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | 595 |
| region_level:UniProt:Region:Interaction with HSF2BP (2270-2337) [source_extracted_region] | 527 |
| region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 498 |
| region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 214 |
| region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 98 |
| region_level:UniProt:Region:Interaction with PALB2 (1-40) [source_extracted_region] | 76 |
| region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 38 |
| UniProt:Mutagenesis:Abolished interaction with DMC1. [source_extracted];region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 27 |
| region_level:UniProt:Region:Interaction with PALB2 (1397-1424) [source_extracted_region] | 19 |
| UniProt:Mutagenesis:Does not affect interaction with DMC1. [source_extracted];region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 16 |
| UniProt:Mutagenesis:Decreased interaction with DMC1. [source_extracted];region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 15 |
| UniProt:Mutagenesis:Does not abolish ABRAXAS1 binding, but abolishes formation of a heterotetramer with ABRAXAS1. [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 11 |
| UniProt:Mutagenesis:Impaired interaction with RAD51. [source_extracted] | 9 |
| UniProt:Mutagenesis:No effect on affinity for a BRIP1 phosphopeptide. [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 8 |
| UniProt:Mutagenesis:Does not abolish ABRAXAS1 binding, but impairs formation of a heterotetramer with ABRAXAS1. [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 6 |

## Criteria Summary

| criteria_combo | count |
| --- | --- |
| BP1+PM2_Supporting | 13933 |
| BS3+PM2_Supporting | 1703 |
| BP1 | 1652 |
| BP4+BP7 | 119 |
| BS3+PM2_Supporting+PP3 | 117 |
| BP1+BS3+PM2_Supporting | 81 |
| BP1+BP5 | 70 |
| BP1+BP5+PM2_Supporting | 44 |
| BP1+BS1_Supporting | 36 |
| BP5+BS3 | 29 |
| BP1+BS1_Strong | 24 |
| BA1 | 23 |
| BP1+BS3 | 22 |
| BP5+BS3+PM2_Supporting | 13 |
| BS1_Supporting+BS3 | 10 |
| BP1+BP5+BS3 | 10 |
| BP5 | 10 |
| BP1+BP5+BS1_Supporting | 7 |
| BP4+BP7+BS1_Supporting | 6 |
| BP1+BP5+BS1_Strong | 6 |

## Top Counterexamples

| counterexample_score | counterexample_tier | gene | c_notation | p_notation | variant_type | predicted_class | structure_features | curated_active_site_status | criteria_combo |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 264 | A_exact_or_known_residue_context | BRCA2 | c.9370A&gt;C | p.(Asn3124His) | missense | 2 | DNA-binding domain / helical-OB-DSS1 region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| 254 | A_exact_or_known_residue_context | BRCA2 | c.8187G&gt;T | p.(Lys2729Asn) | missense | 1 | DNA-binding domain / helical-OB-DSS1 region | UniProt:Natural variant:in BC; benign; no effect on homologous recombination-mediated DNA repair; no effect on interaction with SEM1; dbSNP:rs80359065 [source_extracted];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | BP5+BS3 |
| 244 | A_exact_or_known_residue_context | BRCA2 | c.7469T&gt;C | p.(Ile2490Thr) | missense | 1 | DNA-binding domain / helical-OB-DSS1 region | UniProt:Natural variant:no effect on homologous recombination-mediated DNA repair; no effect on interaction with SEM1; dbSNP:rs11571707 [source_extracted];region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | BA1 |
| 96 | A_exact_or_known_residue_context | BRCA1 | c.117T&gt;C | p.(Cys39=) | synonymous | 2 | RING zinc-binding/E3 ligase region | zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys39 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting |
| 96 | A_exact_or_known_residue_context | BRCA1 | c.123C&gt;T | p.(His41=) | synonymous | 2 | RING zinc-binding/E3 ligase region | zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue His41 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting |
| 96 | A_exact_or_known_residue_context | BRCA1 | c.192T&gt;C | p.(Cys64=) | synonymous | 2 | RING zinc-binding/E3 ligase region | UniProt:Natural variant:in BC; no interaction with BAP1; dbSNP:rs80357064 [source_extracted];zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys64 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5206G&gt;A | p.(Val1736Ile) | missense | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC and FANCS; pathogenic; decreased localization to DNA damage sites and reduced interaction with UIMC1/RAP80; dbSNP:rs45553935 [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5208C&gt;A | p.(Val1736=) | synonymous | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC and FANCS; pathogenic; decreased localization to DNA damage sites and reduced interaction with UIMC1/RAP80; dbSNP:rs45553935 [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5208C&gt;G | p.(Val1736=) | synonymous | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC and FANCS; pathogenic; decreased localization to DNA damage sites and reduced interaction with UIMC1/RAP80; dbSNP:rs45553935 [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5208C&gt;T | p.(Val1736=) | synonymous | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC and FANCS; pathogenic; decreased localization to DNA damage sites and reduced interaction with UIMC1/RAP80; dbSNP:rs45553935 [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5214A&gt;C | p.(Gly1738=) | synonymous | 2 | BRCT phosphopeptide-binding region | UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5214A&gt;G | p.(Gly1738=) | synonymous | 2 | BRCT phosphopeptide-binding region | UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5214A&gt;T | p.(Gly1738=) | synonymous | 2 | BRCT phosphopeptide-binding region | UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5323A&gt;C | p.(Met1775Leu) | missense | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5323A&gt;T | p.(Met1775Leu) | missense | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5324T&gt;C | p.(Met1775Thr) | missense | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 92 | A_exact_or_known_residue_context | BRCA1 | c.5325G&gt;A | p.(Met1775Ile) | missense | 2 | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 89 | A_exact_or_known_residue_context | BRCA1 | c.5117G&gt;C | p.(Gly1706Ala) | missense | 1 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BP5+BS3 |
| 88 | A_exact_or_known_residue_context | BRCA1 | c.109A&gt;C | p.(Thr37Pro) | missense | 2 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting |
| 88 | A_exact_or_known_residue_context | BRCA1 | c.109A&gt;G | p.(Thr37Ala) | missense | 2 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting |
| 88 | A_exact_or_known_residue_context | BRCA1 | c.109A&gt;T | p.(Thr37Ser) | missense | 2 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting |
| 88 | A_exact_or_known_residue_context | BRCA1 | c.111A&gt;C | p.(Thr37=) | synonymous | 2 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting |
| 88 | A_exact_or_known_residue_context | BRCA1 | c.111A&gt;G | p.(Thr37=) | synonymous | 2 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5053A&gt;T | p.(Thr1685Ser) | missense | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5055T&gt;A | p.(Thr1685=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5055T&gt;C | p.(Thr1685=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5055T&gt;G | p.(Thr1685=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5089T&gt;A | p.(Cys1697Ser) | missense | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5089T&gt;G | p.(Cys1697Gly) | missense | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5091T&gt;C | p.(Cys1697=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5118A&gt;C | p.(Gly1706=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5118A&gt;G | p.(Gly1706=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5118A&gt;T | p.(Gly1706=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5144G&gt;C | p.(Ser1715Thr) | missense | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5290C&gt;A | p.(Leu1764Ile) | missense | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5290C&gt;G | p.(Leu1764Val) | missense | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5292A&gt;C | p.(Leu1764=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5292A&gt;G | p.(Leu1764=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5292A&gt;T | p.(Leu1764=) | synonymous | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |
| 84 | A_exact_or_known_residue_context | BRCA1 | c.5296A&gt;C | p.(Ile1766Leu) | missense | 2 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting |

## Interpretation

The most important use of this table is negative control. If a future narrative
or model says that all variants in BRCT, RING, BRCA2 DBD, PALB2, RAD51/BRC, or
SEM1/DSS1 contexts should be treated as suspicious, these benign counterexamples
show why that is too simple.

The practical conclusion is that domain or interface context should be used as
a prioritization and interpretation layer, not as a classification rule. A VUS
in one of these regions becomes more interesting only when the local
neighborhood, ACMG/ENIGMA criteria, variant mechanism, and residue-level
evidence point in the same direction.

The counterexamples should be read together with the pathogenic-region and VUS
prioritization analyses. A VUS in a pathogenic region becomes more interesting
when its local neighborhood, criteria, and residue-level evidence point in the
same direction. Domain membership alone is not enough.

## Outputs

- `tables/benign_counterexamples/benign_counterexamples.csv`
- `tables/benign_counterexamples/tier_summary.csv`
- `tables/benign_counterexamples/feature_summary.csv`
- `tables/benign_counterexamples/gene_feature_summary.csv`
- `tables/benign_counterexamples/active_site_summary.csv`
- `tables/benign_counterexamples/criteria_summary.csv`
- `plots/20_benign_counterexamples/tier_summary.svg`
- `plots/20_benign_counterexamples/feature_summary.svg`
- `plots/20_benign_counterexamples/active_site_summary.svg`
