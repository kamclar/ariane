# Positional Context Follow-up Analyses

Generated: 2026-06-21 11:30

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Purpose

These follow-up analyses ask why positional context can look informative in the
precomputed Module 1 map. The goal is to separate regional signals driven by
variant type, splice prediction, functional evidence, and mixed positions or
codons.

This remains a context analysis only. It does not create ACMG/ENIGMA evidence.

## 1. Variant Type Effect

| variant_type | total_count | benign_count | vus_count | pathogenic_count | pathogenic_fraction | truncation_driver_count | splice_driver_count | functional_pathogenic_driver_count | functional_benign_driver_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| missense | 35241 | 28762 | 6376 | 103 | 0.0029 | 0 | 845 | 546 | 1970 |
| synonymous | 9891 | 8189 | 1693 | 9 | 0.0009 | 0 | 215 | 11 | 557 |
| nonsense | 2397 | 1 | 70 | 2326 | 0.9704 | 2397 | 146 | 147 | 1 |
| initiation_codon | 18 | 0 | 15 | 3 | 0.1667 | 0 | 1 | 10 | 0 |

## 2. Splice-driven Regional Signal

Top bins by splice-driver fraction:

| gene | bin | cds_start | cds_end | pathogenic_fraction | splice_driver_fraction | pp3_count | high_spliceai_count | near_boundary_count | total_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | 3 | 201 | 300 | 0.0933 | 0.3100 | 84 | 93 | 96 | 300 |
| BRCA2 | 92 | 9101 | 9200 | 0.0633 | 0.2167 | 59 | 65 | 66 | 300 |
| BRCA1 | 44 | 4301 | 4400 | 0.0600 | 0.1967 | 50 | 59 | 66 | 300 |
| BRCA1 | 2 | 101 | 200 | 0.1000 | 0.1700 | 44 | 50 | 66 | 300 |
| BRCA1 | 51 | 5001 | 5100 | 0.1233 | 0.1667 | 43 | 50 | 66 | 300 |
| BRCA2 | 90 | 8901 | 9000 | 0.0867 | 0.1667 | 43 | 50 | 66 | 300 |
| BRCA2 | 91 | 9001 | 9100 | 0.0767 | 0.1633 | 43 | 49 | 0 | 300 |
| BRCA1 | 7 | 601 | 700 | 0.0100 | 0.1567 | 38 | 47 | 78 | 300 |
| BRCA2 | 78 | 7701 | 7800 | 0.0733 | 0.1100 | 31 | 33 | 18 | 300 |
| BRCA2 | 6 | 501 | 600 | 0.0433 | 0.1033 | 25 | 31 | 66 | 300 |
| BRCA1 | 41 | 4001 | 4100 | 0.0700 | 0.1000 | 23 | 30 | 45 | 300 |
| BRCA2 | 98 | 9701 | 9800 | 0.0600 | 0.0833 | 24 | 25 | 0 | 300 |
| BRCA2 | 76 | 7501 | 7600 | 0.0300 | 0.0833 | 24 | 25 | 0 | 300 |
| BRCA2 | 75 | 7401 | 7500 | 0.0567 | 0.0800 | 22 | 24 | 66 | 300 |
| BRCA2 | 97 | 9601 | 9700 | 0.0567 | 0.0800 | 21 | 24 | 66 | 300 |
| BRCA2 | 5 | 401 | 500 | 0.0367 | 0.0800 | 22 | 24 | 132 | 300 |
| BRCA1 | 6 | 501 | 600 | 0.0333 | 0.0733 | 20 | 22 | 120 | 300 |
| BRCA2 | 7 | 601 | 700 | 0.0233 | 0.0700 | 20 | 21 | 132 | 300 |
| BRCA2 | 79 | 7801 | 7900 | 0.0600 | 0.0667 | 19 | 20 | 48 | 300 |
| BRCA1 | 8 | 701 | 800 | 0.0433 | 0.0667 | 20 | 20 | 0 | 300 |

## 3. Functional-evidence Regional Signal

Top bins by functional evidence count:

| gene | bin | cds_start | cds_end | pathogenic_fraction | ps3_ps1_pp4_count | bs3_count | bs3_pp3_count | ps3_pp3_count | total_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | 52 | 5101 | 5200 | 0.0867 | 82 | 198 | 3 | 6 | 300 |
| BRCA1 | 54 | 5301 | 5400 | 0.0667 | 46 | 225 | 1 | 0 | 300 |
| BRCA1 | 50 | 4901 | 5000 | 0.0600 | 33 | 238 | 4 | 4 | 300 |
| BRCA1 | 1 | 1 | 100 | 0.0733 | 59 | 211 | 4 | 3 | 300 |
| BRCA1 | 2 | 101 | 200 | 0.1000 | 75 | 191 | 34 | 9 | 300 |
| BRCA1 | 55 | 5401 | 5500 | 0.0367 | 33 | 231 | 8 | 5 | 300 |
| BRCA1 | 53 | 5201 | 5300 | 0.0700 | 72 | 185 | 3 | 1 | 300 |
| BRCA1 | 51 | 5001 | 5100 | 0.1233 | 75 | 179 | 10 | 25 | 300 |
| BRCA1 | 3 | 201 | 300 | 0.0933 | 52 | 201 | 53 | 15 | 300 |
| BRCA1 | 56 | 5501 | 5592 | 0.0616 | 46 | 132 | 0 | 0 | 276 |
| BRCA1 | 49 | 4801 | 4900 | 0.0367 | 1 | 36 | 2 | 1 | 300 |
| BRCA2 | 79 | 7801 | 7900 | 0.0600 | 24 | 10 | 0 | 0 | 300 |
| BRCA2 | 82 | 8101 | 8200 | 0.0500 | 17 | 7 | 0 | 1 | 300 |
| BRCA2 | 84 | 8301 | 8400 | 0.0400 | 9 | 15 | 0 | 0 | 300 |
| BRCA2 | 80 | 7901 | 8000 | 0.0700 | 12 | 9 | 0 | 3 | 300 |
| BRCA2 | 81 | 8001 | 8100 | 0.0467 | 8 | 12 | 0 | 2 | 300 |
| BRCA1 | 48 | 4701 | 4800 | 0.0300 | 0 | 20 | 0 | 0 | 300 |
| BRCA2 | 93 | 9201 | 9300 | 0.0633 | 10 | 9 | 0 | 1 | 300 |
| BRCA2 | 85 | 8401 | 8500 | 0.0633 | 8 | 9 | 0 | 2 | 300 |
| BRCA2 | 94 | 9301 | 9400 | 0.0533 | 8 | 8 | 0 | 0 | 300 |

## 4. Mixed Position/Codon Context After Driver Exclusions

| filter | variant_count | positions_total | positions_mixed | positions_with_benign_and_pathogenic | codons_total | codons_mixed | codons_with_benign_and_pathogenic |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all_variants | 47547 | 15849 | 3104 | 1654 | 5283 | 2565 | 1737 |
| no_truncation | 45150 | 15849 | 1186 | 36 | 5283 | 906 | 56 |
| no_splice | 46340 | 15716 | 2762 | 1572 | 5273 | 2420 | 1636 |
| no_truncation_no_splice | 44089 | 15699 | 933 | 11 | 5270 | 755 | 20 |
| no_truncation_no_splice_no_functional | 41209 | 14966 | 462 | 0 | 5117 | 414 | 0 |

## Interpretation

- Nonsense/PTC variants strongly explain many pathogenic-enriched regions.
- Splice driver signal identifies bins where local context is tied to `PP3`,
  high SpliceAI, or proximity to CDS exon boundaries.
- Functional evidence explains BRCA1 RING and BRCT mixed regions, especially
  where `PS3` and `BS3` are both common.
- Mixed positions and codons remain useful for review interpretation because a
  coordinate can host multiple specific variant consequences.

The practical conclusion is unchanged: positional context is useful for triage
and explanation, but classification must remain variant-level and criterion
based.

## Outputs

- `tables/positional_context_followup/variant_type_driver_summary.csv`
- `tables/positional_context_followup/regional_driver_ablation_by_type.csv`
- `tables/positional_context_followup/splice_driver_regions.csv`
- `tables/positional_context_followup/functional_driver_regions.csv`
- `tables/positional_context_followup/mixed_context_driver_ablation.csv`
- `plots/26_positional_context_followup/variant_type_grouped_class_stacked.svg`
- `plots/26_positional_context_followup/splice_driver_vs_pathogenic_fraction.svg`
- `plots/26_positional_context_followup/functional_driver_vs_pathogenic_fraction.svg`
- `plots/26_positional_context_followup/mixed_position_codon_driver_ablation.svg`
