# BS3 In Functional Domains

Generated: 2026-06-20 17:52

Input: `Exploratory_analysis\precomputed_classification\tables\benign_structure_function\benign_structure_function_variants.csv`

## Purpose

This analysis focuses on benign/likely benign variants in functional domains or
active/interface contexts where `BS3` is present. It asks when functional
benign evidence is acting as the main counterweight to domain or interface
context, and where conflicts such as `BS3+PP3` or high SpliceAI appear.

This is a review aid only. It does not change classifications and it does not
create new ACMG/ENIGMA criteria.

## Summary

- Benign domain/interface variants analyzed: 17955
- With `BS3`: 2017
- With `BS3` plus `PP3` or SpliceAI >= 0.20: 126

## Review Buckets

| review_bucket | count |
| --- | --- |
| non_BS3_benign_domain_context | 15938 |
| BS3_in_region_level_interface_context | 1379 |
| BS3_in_broad_domain_context | 443 |
| BS3_with_PP3_and_high_SpliceAI | 126 |
| BS3_in_exact_or_source_active_context | 69 |

## BS3 By Domain

| feature | benign_domain_count | bs3_count | bs3_percent | bs3_pp3_count | bs3_high_spliceai_count | bs3_missense_count | bs3_synonymous_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BRCT phosphopeptide-binding region | 1178 | 1174 | 99.66 | 27 | 27 | 849 | 325 |
| RING zinc-binding/E3 ligase region | 576 | 571 | 99.13 | 92 | 92 | 423 | 148 |
| DNA-binding domain / helical-OB-DSS1 region | 256 | 127 | 49.61 | 3 | 3 | 127 | 0 |
| BRC repeats / RAD51-binding region | 9242 | 32 | 0.35 | 0 | 0 | 32 | 0 |
| C-terminal RAD51-binding/nuclear localization context | 1834 | 14 | 0.76 | 0 | 0 | 14 | 0 |
| coiled-coil PALB2 interaction region | 22 | 14 | 63.64 | 0 | 0 | 14 | 0 |
| PALB2 binding N-terminal region | 8 | 4 | 50.00 | 0 | 0 | 4 | 0 |

## BS3 Variant Mechanisms

| variant_mechanism_group | count |
| --- | --- |
| missense | 1530 |
| synonymous | 487 |

## BS3 Criteria Combinations

| criteria_combo | count |
| --- | --- |
| BS3+PM2_Supporting | 1703 |
| BS3+PM2_Supporting+PP3 | 117 |
| BP1+BS3+PM2_Supporting | 81 |
| BP5+BS3 | 29 |
| BP1+BS3 | 22 |
| BP5+BS3+PM2_Supporting | 13 |
| BS1_Supporting+BS3 | 10 |
| BP1+BP5+BS3 | 10 |
| BS3+PP3 | 5 |
| BP5+BS1_Strong+BS3 | 5 |
| BS1_Strong+BS3 | 5 |
| BP5+BS1_Supporting+BS3 | 4 |
| BP1+BP5+BS1_Strong+BS3 | 4 |
| BP5+BS3+PP3 | 2 |
| BP1+BP5+BS3+PM2_Supporting | 2 |
| BP1+BS1_Supporting+BS3 | 2 |
| BS1_Supporting+BS3+PP3 | 1 |
| BP5+BS1_Supporting+BS3+PP3 | 1 |
| BP1+BP5+BS1_Supporting+BS3 | 1 |

## Top BS3 Domain Conflicts

| bs3_domain_conflict_score | review_bucket | gene | c_notation | p_notation | variant_type | spliceai_score | structure_features | curated_active_site_status | criteria_combo |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 100 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.192T&gt;C | p.(Cys64=) | synonymous | 0.87 | RING zinc-binding/E3 ligase region | UniProt:Natural variant:in BC; no interaction with BAP1; dbSNP:rs80357064 [source_extracted];zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys64 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 100 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.213G&gt;C | p.(Arg71Ser) | missense | 0.39 | RING zinc-binding/E3 ligase region | UniProt:Mutagenesis:No effect on interaction with BAP1. [source_extracted] | BS3+PM2_Supporting+PP3 |
| 100 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.4985T&gt;G | p.(Phe1662Cys) | missense | 0.27 | BRCT phosphopeptide-binding region | UniProt:Mutagenesis:Does not abolish ABRAXAS1 binding, but abolishes formation of a heterotetramer with ABRAXAS1. [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5196T&gt;G | p.(His1732Gln) | missense | 0.99 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5470A&gt;G | p.(Ile1824Val) | missense | 0.93 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA2 | c.6842G&gt;T | p.(Gly2281Val) | missense | 0.92 | outside_mapped_structure_function_features | region_level:UniProt:Region:Interaction with HSF2BP (2270-2337) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.193A&gt;C | p.(Lys65Gln) | missense | 0.87 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.195G&gt;C | p.(Lys65Asn) | missense | 0.87 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.195G&gt;T | p.(Lys65Asn) | missense | 0.86 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.194A&gt;G | p.(Lys65Arg) | missense | 0.84 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.194A&gt;T | p.(Lys65Met) | missense | 0.8 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5196T&gt;A | p.(His1732Gln) | missense | 0.71 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.162G&gt;T | p.(Gln54His) | missense | 0.7 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.187T&gt;G | p.(Leu63Val) | missense | 0.7 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5270A&gt;T | p.(Asp1757Val) | missense | 0.69 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.194A&gt;C | p.(Lys65Thr) | missense | 0.67 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.83T&gt;A | p.(Leu28Gln) | missense | 0.57 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.165G&gt;T | p.(Lys55Asn) | missense | 0.55 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.4931A&gt;T | p.(Glu1644Val) | missense | 0.53 | outside_mapped_structure_function_features | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA2 | c.7447A&gt;G | p.(Ser2483Gly) | missense | 0.51 | DNA-binding domain / helical-OB-DSS1 region | region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | BS1_Supporting+BS3+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.193A&gt;G | p.(Lys65Glu) | missense | 0.49 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.4947A&gt;T | p.(Arg1649Ser) | missense | 0.41 | outside_mapped_structure_function_features | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5042C&gt;G | p.(Thr1681Ser) | missense | 0.38 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.135A&gt;T | p.(Lys45Asn) | missense | 0.37 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA2 | c.7448G&gt;A | p.(Ser2483Asn) | missense | 0.35 | DNA-binding domain / helical-OB-DSS1 region | region_level:UniProt:Region:Interaction with FANCD2 (2350-2545) [source_extracted_region];region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | BP5+BS3+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5140G&gt;A | p.(Val1714Ile) | missense | 0.34 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.167A&gt;T | p.(Lys56Ile) | missense | 0.33 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5450A&gt;G | p.(Glu1817Gly) | missense | 0.33 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.189A&gt;T | p.(Leu63Phe) | missense | 0.32 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5000A&gt;G | p.(Lys1667Arg) | missense | 0.29 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.162G&gt;C | p.(Gln54His) | missense | 0.28 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.165G&gt;C | p.(Lys55Asn) | missense | 0.28 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5277G&gt;T | p.(Lys1759Asn) | missense | 0.28 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5453A&gt;T | p.(Asp1818Val) | missense | 0.28 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.163A&gt;C | p.(Lys55Gln) | missense | 0.27 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.164A&gt;T | p.(Lys55Met) | missense | 0.27 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.163A&gt;G | p.(Lys55Glu) | missense | 0.26 | RING zinc-binding/E3 ligase region | region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5035C&gt;G | p.(Leu1679Val) | missense | 0.26 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5045A&gt;C | p.(Glu1682Ala) | missense | 0.26 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BS3+PM2_Supporting+PP3 |
| 85 | BS3_with_PP3_and_high_SpliceAI | BRCA1 | c.5471T&gt;G | p.(Ile1824Ser) | missense | 0.25 | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BS3+PM2_Supporting+PP3 |

## Interpretation

The main question is not whether domain context exists, but whether the BS3
assay is measuring the relevant mechanism for that region. A benign functional
result is strongest when the assay captures the disease-relevant mechanism for
that domain. It is more ambiguous when the same variant also has a PP3/splice
signal or when the assay may not detect splicing, NMD, or another mechanism.

Missense and synonymous counterexamples should be interpreted separately.
Missense variants are more directly relevant to protein-function hypotheses.
Synonymous variants in protein domains are often more about splice/RNA context
than protein active-site biology.

## Outputs

- `tables/bs3_domain_conflicts/bs3_domain_conflict_variants.csv`
- `tables/bs3_domain_conflicts/review_bucket_summary.csv`
- `tables/bs3_domain_conflicts/domain_conflict_summary.csv`
- `tables/bs3_domain_conflicts/variant_mechanism_summary.csv`
- `tables/bs3_domain_conflicts/criteria_summary.csv`
- `plots/21_bs3_domain_conflicts/review_bucket_summary.svg`
- `plots/21_bs3_domain_conflicts/domain_conflict_summary.svg`
- `plots/21_bs3_domain_conflicts/domain_mechanism_stacked.svg`
