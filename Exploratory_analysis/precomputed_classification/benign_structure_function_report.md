# Benign Structure Function Mapping

Generated: 2026-06-20 17:16

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Purpose

This is the benign counterpart to the pathogenic structure/function mapping.
It asks where generated benign and likely benign variants, classes 1 and 2, lie
relative to the same structural domains and UniProt active/interface annotations.

This does not mean that a structure annotation is benign evidence by itself. It
is a landscape check: if benign variants also occur in a region, then location
alone is not sufficient for interpretation.

## Dataset

- Benign/likely benign variants analyzed: 36952
- Variants inside mapped structure/function features: 13116
- Variants with exact or region-level UniProt active/interface annotation: 15553

## Feature Summary

| feature | count |
| --- | --- |
| outside_mapped_structure_function_features | 23836 |
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
| BRCA1 | outside_mapped_structure_function_features | 12498 |
| BRCA2 | outside_mapped_structure_function_features | 11338 |
| BRCA2 | BRC repeats / RAD51-binding region | 9242 |
| BRCA2 | C-terminal RAD51-binding/nuclear localization context | 1834 |
| BRCA1 | BRCT phosphopeptide-binding region | 1178 |
| BRCA1 | RING zinc-binding/E3 ligase region | 576 |
| BRCA2 | DNA-binding domain / helical-OB-DSS1 region | 256 |
| BRCA1 | coiled-coil PALB2 interaction region | 22 |
| BRCA2 | PALB2 binding N-terminal region | 8 |

## Curated Active Site Or Interface Status

| curated_active_site_status | count |
| --- | --- |
| not_curated_in_current_dataset | 21399 |
| region_level:UniProt:Region:Interaction with RAD51 (1003-2082) [source_extracted_region] | 9233 |
| region_level:UniProt:Region:Interaction with POLH (1338-1781) [source_extracted_region] | 3783 |
| region_level:UniProt:Region:Interaction with NPM1 (639-1000) [source_extracted_region] | 3076 |
| region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 1134 |
| region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | 609 |
| region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 537 |
| region_level:UniProt:Region:Interaction with HSF2BP (2270-2337) [source_extracted_region] | 527 |
| region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 222 |
| region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 142 |
| region_level:UniProt:Region:Interaction with PALB2 (1-40) [source_extracted_region] | 76 |
| UniProt:Mutagenesis:Abolished interaction with DMC1. [source_extracted] | 27 |
| region_level:UniProt:Region:Interaction with PALB2 (1397-1424) [source_extracted_region] | 19 |
| UniProt:Mutagenesis:Does not affect interaction with DMC1. [source_extracted] | 16 |
| UniProt:Mutagenesis:Decreased interaction with DMC1. [source_extracted] | 15 |
| UniProt:Natural variant:in BC | 13 |
|  pathogenic | 12 |
| UniProt:Mutagenesis:Does not abolish ABRAXAS1 binding, but abolishes formation of a heterotetramer with ABRAXAS1. [source_extracted] | 11 |
| UniProt:Mutagenesis:Impaired interaction with RAD51. [source_extracted] | 9 |
| UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | 8 |

## Variant Types

| variant_type | count |
| --- | --- |
| missense | 28762 |
| synonymous | 8189 |
| nonsense | 1 |

## Criteria Combinations

| criteria_combo | count |
| --- | --- |
| BP1+PM2_Supporting | 31088 |
| BP1 | 2866 |
| BS3+PM2_Supporting | 1708 |
| BP1+BS3+PM2_Supporting | 261 |
| BP1+BP5 | 182 |
| BS3+PM2_Supporting+PP3 | 119 |
| BP4+BP7 | 119 |
| BP1+BP5+PM2_Supporting | 103 |
| BP1+BS3 | 93 |
| BP1+BS1_Supporting | 74 |
| BP1+BP5+BS3 | 61 |
| BP1+BS1_Strong | 50 |
| BA1 | 48 |
| BP5+BS3 | 30 |
| BP1+BP5+BS1_Strong | 23 |
| BP1+BP5+BS1_Supporting | 17 |
| BP5+BS3+PM2_Supporting | 13 |
| BS1_Supporting+BS3 | 11 |
| BP1+BP5+BS3+PM2_Supporting | 11 |
| BP1+BP5+BS1_Strong+BS3 | 11 |

## Benign Versus Pathogenic Feature Comparison

The percentages are within each group. This helps avoid reading raw counts as
enrichment when benign variants greatly outnumber pathogenic variants in the
generated snapshot.

| feature | benign_count | pathogenic_count | benign_percent | pathogenic_percent | pathogenic_to_benign_rate_ratio |
| --- | --- | --- | --- | --- | --- |
| outside_mapped_structure_function_features | 23836 | 1284 | 64.505 | 52.601 | 0.815 |
| BRC repeats / RAD51-binding region | 9242 | 487 | 25.011 | 19.951 | 0.798 |
| DNA-binding domain / helical-OB-DSS1 region | 256 | 362 | 0.693 | 14.830 | 21.406 |
| BRCT phosphopeptide-binding region | 1178 | 140 | 3.188 | 5.735 | 1.799 |
| RING zinc-binding/E3 ligase region | 576 | 80 | 1.559 | 3.277 | 2.103 |
| C-terminal RAD51-binding/nuclear localization context | 1834 | 57 | 4.963 | 2.335 | 0.470 |
| coiled-coil PALB2 interaction region | 22 | 17 | 0.060 | 0.696 | 11.698 |
| PALB2 binding N-terminal region | 8 | 14 | 0.022 | 0.574 | 26.492 |

## Benign Versus Pathogenic Active/Interface Comparison

| curated_active_site_status | benign_count | pathogenic_count | benign_percent | pathogenic_percent | pathogenic_to_benign_rate_ratio |
| --- | --- | --- | --- | --- | --- |
| not_curated_in_current_dataset | 21399 | 1305 | 57.910 | 53.462 | 0.923 |
| region_level:UniProt:Region:Interaction with RAD51 (1003-2082) [source_extracted_region] | 5450 | 274 | 14.749 | 11.225 | 0.761 |
| region_level:UniProt:Region:Interaction with RAD51 (1003-2082) [source_extracted_region];region_level:UniProt:Region:Interaction with POLH (1338-1781) [source_extracted_region] | 3783 | 213 | 10.238 | 8.726 | 0.852 |
| region_level:UniProt:Region:Interaction with NPM1 (639-1000) [source_extracted_region] | 3076 | 180 | 8.324 | 7.374 | 0.886 |
| region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 98 | 127 | 0.265 | 5.203 | 19.618 |
| region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 498 | 71 | 1.348 | 2.909 | 2.158 |
| region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 1037 | 54 | 2.806 | 2.212 | 0.788 |
| region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | 595 | 50 | 1.610 | 2.048 | 1.272 |
| region_level:UniProt:Region:Interaction with HSF2BP (2270-2337) [source_extracted_region] | 527 | 38 | 1.426 | 1.557 | 1.092 |
| region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 38 | 24 | 0.103 | 0.983 | 9.561 |
| region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 214 | 16 | 0.579 | 0.655 | 1.132 |
| region_level:UniProt:Region:Interaction with PALB2 (1-40) [source_extracted_region] | 76 | 15 | 0.206 | 0.615 | 2.988 |
| region_level:UniProt:Region:Interaction with PALB2 (1397-1424) [source_extracted_region] | 19 | 12 | 0.051 | 0.492 | 9.561 |
| UniProt:Mutagenesis:No effect on interaction with BAP1. [source_extracted] | 2 | 7 | 0.005 | 0.287 | 52.983 |
| UniProt:Natural variant:in BC; no interaction with BAP1; dbSNP:rs80357064 [source_extracted];zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys64 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 1 | 5 | 0.003 | 0.205 | 75.690 |
| region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region];region_level:UniProt:Motif:Nuclear export signal; masked by interaction with SEM1 (2682-2698) [source_extracted_region] | 4 | 5 | 0.011 | 0.205 | 18.923 |
| UniProt:Natural variant:in BC and OC; pathogenic; no interaction with BAP1; dbSNP:rs28897672 [source_extracted];zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys61 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 0 | 4 | 0.000 | 0.164 | NA |
| zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys27 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 0 | 4 | 0.000 | 0.164 | NA |
| zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys44 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 0 | 4 | 0.000 | 0.164 | NA |
| UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | 3 | 3 | 0.008 | 0.123 | 15.138 |
| UniProt:Mutagenesis:Does not abolish ABRAXAS1 binding, but abolishes formation of a heterotetramer with ABRAXAS1. [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 11 | 3 | 0.030 | 0.123 | 4.129 |
| UniProt:Mutagenesis:Decreased interaction with DMC1. [source_extracted];region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 15 | 2 | 0.041 | 0.082 | 2.018 |
| UniProt:Mutagenesis:Disrupts interaction with SEM1. [source_extracted];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 0 | 2 | 0.000 | 0.082 | NA |
| UniProt:Mutagenesis:Does not abolish ABRAXAS1 binding, but impairs formation of a heterotetramer with ABRAXAS1. [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 6 | 2 | 0.016 | 0.082 | 5.046 |
| UniProt:Mutagenesis:Does not affect interaction with DMC1. [source_extracted];region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region] | 16 | 2 | 0.043 | 0.082 | 1.892 |
| UniProt:Natural variant:in BC; abolishes interaction with PALB2; dbSNP:rs80359214 [source_extracted];UniProt:Natural variant:in BC; abolishes interaction with PALB2; dbSNP:rs80359182 [source_extracted];region_level:UniProt:Region:Interaction with PALB2 (1-40) [source_extracted_region] | 0 | 2 | 0.000 | 0.082 | NA |
| UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | 4 | 2 | 0.011 | 0.082 | 7.569 |
| UniProt:Natural variant:in FANCD1; hypersensitive to DNA damage; reduced homology-directed repair activity; no effect on interaction with SEM1; dbSNP:rs80359013 [source_extracted];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 0 | 2 | 0.000 | 0.082 | NA |
| zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys47 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 0 | 2 | 0.000 | 0.082 | NA |
| UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 5 | 1 | 0.014 | 0.041 | 3.028 |

## Top Benign Variants In Mapped Features

| gene | c_notation | p_notation | variant_type | predicted_class | structure_features | curated_active_site_status | criteria_combo |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | c.4G&gt;T | p.(Asp2Tyr) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.5A&gt;C | p.(Asp2Ala) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.6T&gt;A | p.(Asp2Glu) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.6T&gt;G | p.(Asp2Glu) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.7T&gt;A | p.(Leu3Ile) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.7T&gt;C | p.(Leu3=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.7T&gt;G | p.(Leu3Val) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.8T&gt;C | p.(Leu3Ser) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.9A&gt;C | p.(Leu3Phe) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.9A&gt;T | p.(Leu3Phe) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.10T&gt;A | p.(Ser4Thr) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.10T&gt;C | p.(Ser4Pro) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.10T&gt;G | p.(Ser4Ala) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.11C&gt;A | p.(Ser4Tyr) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.11C&gt;G | p.(Ser4Cys) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.11C&gt;T | p.(Ser4Phe) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.12T&gt;A | p.(Ser4=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.12T&gt;C | p.(Ser4=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.12T&gt;G | p.(Ser4=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.13G&gt;A | p.(Ala5Thr) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.13G&gt;T | p.(Ala5Ser) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.14C&gt;A | p.(Ala5Asp) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.14C&gt;G | p.(Ala5Gly) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.15T&gt;A | p.(Ala5=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.15T&gt;C | p.(Ala5=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.15T&gt;G | p.(Ala5=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.16C&gt;A | p.(Leu6Ile) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.16C&gt;G | p.(Leu6Val) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.17T&gt;A | p.(Leu6His) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.17T&gt;C | p.(Leu6Pro) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.17T&gt;G | p.(Leu6Arg) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.18T&gt;A | p.(Leu6=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.18T&gt;C | p.(Leu6=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.18T&gt;G | p.(Leu6=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.19C&gt;G | p.(Arg7Gly) | missense | 1 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BP5+BS3+PM2_Supporting |
| BRCA1 | c.20G&gt;A | p.(Arg7His) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS1_Supporting+BS3 |
| BRCA1 | c.20G&gt;C | p.(Arg7Pro) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.20G&gt;T | p.(Arg7Leu) | missense | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.21C&gt;A | p.(Arg7=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |
| BRCA1 | c.21C&gt;G | p.(Arg7=) | synonymous | 2 | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | BS3+PM2_Supporting |

## Interpretation

- Benign/likely benign variants are expected to be numerous because BP1, BP4,
  BP7, BS3 and related benign-direction rules classify many generated SNVs.
- Structural location alone is not enough: benign variants can occur inside
  broad domains such as BRCT, RING, BRC/RAD51 context, or BRCA2 DBD.
- Exact residue-level functional annotations are more useful than broad domain
  membership, but they still require expert interpretation and should not be
  treated as an automatic criterion outside ACMG/ENIGMA rules.

## Outputs

- `tables/benign_structure_function/benign_structure_function_variants.csv`
- `tables/benign_structure_function/benign_feature_summary.csv`
- `tables/benign_structure_function/benign_gene_feature_summary.csv`
- `tables/benign_structure_function/benign_active_site_summary.csv`
- `tables/benign_structure_function/benign_vs_pathogenic_feature_comparison.csv`
- `tables/benign_structure_function/benign_vs_pathogenic_active_site_comparison.csv`
- `plots/19_benign_structure_function/benign_feature_summary.svg`
- `plots/19_benign_structure_function/benign_active_site_summary.svg`
- `plots/19_benign_structure_function/benign_vs_pathogenic_feature_comparison.svg`
- `plots/19_benign_structure_function/benign_vs_pathogenic_active_site_comparison.svg`
- `plots/19_benign_structure_function/brca1_benign_lollipop.svg`
- `plots/19_benign_structure_function/brca2_benign_lollipop.svg`
