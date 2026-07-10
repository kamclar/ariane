# Regional Driver Decomposition

Generated: 2026-06-21 11:21

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Purpose

This analysis asks what drives local pathogenic enrichment in the precomputed
ARIANE Module 1 coding SNV map. It is not a classifier and it does not add new
ACMG/ENIGMA evidence. It decomposes regional signal into rule/mechanism
components so that local context can be interpreted more safely.

Core interpretation:

> Local pathogenic enrichment is a context signal, not an ACMG/ENIGMA evidence
> criterion.

## Method

The CDS of BRCA1 and BRCA2 was split into 100 bp bins. For each bin, the
analysis counted benign/VUS/pathogenic generated classes and several driver
signals:

- truncation driver: nonsense, `PVS1`, or `PM5_PTC`
- splice driver: `PP3` or SpliceAI >= 0.20
- pathogenic functional driver: `PS3`, `PS1`, or `PP4`
- benign functional driver: `BS3`
- benign rule driver: `BP1`, `BP4`, `BP5`, `BP7`, `BA1`, or `BS1`
- frequency/absence driver: `PM2`, `BS1`, or `BA1`
- low-information variants: `PM2_Supporting` only or no criteria

The analysis also recalculated regional pathogenic fraction after excluding
major drivers:

- all variants
- excluding truncation drivers
- excluding splice drivers
- excluding functional evidence drivers
- excluding truncation, splice, and functional evidence together

## Regional Categories

| region_category | count |
| --- | --- |
| truncation_driven | 146 |
| functional_pathogenic_driven | 10 |
| benign_evidence_dominant | 3 |

## Driver Ablation Summary

| gene | metric | mean | max | bins_ge_0_10 | bins_ge_0_25 |
| --- | --- | --- | --- | --- | --- |
| BRCA1 | pathogenic_fraction | 0.0523 | 0.1233 | 2 | 0 |
| BRCA1 | pathogenic_fraction_no_truncation | 0.0059 | 0.0931 | 0 | 0 |
| BRCA1 | pathogenic_fraction_no_splice | 0.0456 | 0.0805 | 0 | 0 |
| BRCA1 | pathogenic_fraction_no_functional | 0.0403 | 0.0808 | 0 | 0 |
| BRCA1 | residual_pathogenic_fraction | 0.0000 | 0.0000 | 0 | 0 |
| BRCA2 | pathogenic_fraction | 0.0506 | 0.0867 | 0 | 0 |
| BRCA2 | pathogenic_fraction_no_truncation | 0.0007 | 0.0176 | 0 | 0 |
| BRCA2 | pathogenic_fraction_no_splice | 0.0490 | 0.0836 | 0 | 0 |
| BRCA2 | pathogenic_fraction_no_functional | 0.0503 | 0.0887 | 0 | 0 |
| BRCA2 | residual_pathogenic_fraction | 0.0000 | 0.0000 | 0 | 0 |
| both | pathogenic_fraction | 0.0512 | 0.1233 | 2 | 0 |
| both | pathogenic_fraction_no_truncation | 0.0025 | 0.0931 | 0 | 0 |
| both | pathogenic_fraction_no_splice | 0.0478 | 0.0836 | 0 | 0 |
| both | pathogenic_fraction_no_functional | 0.0468 | 0.0887 | 0 | 0 |
| both | residual_pathogenic_fraction | 0.0000 | 0.0000 | 0 | 0 |

## Top Pathogenic-Enriched Bins

| gene | bin | cds_start | cds_end | aa_start | aa_end | pathogenic_fraction | pathogenic_count | total_count | region_category | top_structure_feature | truncation_pathogenic_fraction | splice_pathogenic_fraction | functional_pathogenic_fraction | residual_pathogenic_fraction | mixed_codon_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | 51 | 5001 | 5100 | 1667 | 1700 | 0.1233 | 37 | 300 | functional_pathogenic_driven | BRCT phosphopeptide-binding region | 0.2703 | 0.8378 | 1.0000 | 0.0000 | 33 |
| BRCA1 | 2 | 101 | 200 | 34 | 67 | 0.1000 | 30 | 300 | functional_pathogenic_driven | RING zinc-binding/E3 ligase region | 0.5667 | 0.5000 | 1.0000 | 0.0000 | 31 |
| BRCA1 | 3 | 201 | 300 | 67 | 100 | 0.0933 | 28 | 300 | functional_pathogenic_driven | RING zinc-binding/E3 ligase region | 0.4643 | 0.8571 | 0.9286 | 0.0000 | 30 |
| BRCA1 | 52 | 5101 | 5200 | 1701 | 1734 | 0.0867 | 26 | 300 | functional_pathogenic_driven | BRCT phosphopeptide-binding region | 0.7308 | 0.1538 | 1.0000 | 0.0000 | 30 |
| BRCA2 | 90 | 8901 | 9000 | 2967 | 3000 | 0.0867 | 26 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 1.0000 | 0.2692 | 0.0000 | 0.0000 | 23 |
| BRCA2 | 86 | 8501 | 8600 | 2834 | 2867 | 0.0833 | 25 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 25 |
| BRCA1 | 37 | 3601 | 3700 | 1201 | 1234 | 0.0800 | 24 | 300 | truncation_driven | none | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 22 |
| BRCA2 | 91 | 9001 | 9100 | 3001 | 3034 | 0.0767 | 23 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 0.9565 | 0.3043 | 0.0435 | 0.0000 | 23 |
| BRCA1 | 1 | 1 | 100 | 1 | 34 | 0.0733 | 22 | 300 | functional_pathogenic_driven | RING zinc-binding/E3 ligase region | 0.7273 | 0.3182 | 1.0000 | 0.0000 | 31 |
| BRCA2 | 53 | 5201 | 5300 | 1734 | 1767 | 0.0733 | 22 | 300 | truncation_driven | BRC repeats / RAD51-binding region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 15 |
| BRCA2 | 63 | 6201 | 6300 | 2067 | 2100 | 0.0733 | 22 | 300 | truncation_driven | BRC repeats / RAD51-binding region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 16 |
| BRCA2 | 69 | 6801 | 6900 | 2267 | 2300 | 0.0733 | 22 | 300 | truncation_driven | none | 1.0000 | 0.1364 | 0.0000 | 0.0000 | 23 |
| BRCA2 | 78 | 7701 | 7800 | 2567 | 2600 | 0.0733 | 22 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 1.0000 | 0.0909 | 0.0000 | 0.0000 | 21 |
| BRCA1 | 41 | 4001 | 4100 | 1334 | 1367 | 0.0700 | 21 | 300 | truncation_driven | none | 1.0000 | 0.3333 | 0.0000 | 0.0000 | 27 |
| BRCA1 | 53 | 5201 | 5300 | 1734 | 1767 | 0.0700 | 21 | 300 | functional_pathogenic_driven | BRCT phosphopeptide-binding region | 0.6667 | 0.0476 | 1.0000 | 0.0000 | 32 |
| BRCA2 | 80 | 7901 | 8000 | 2634 | 2667 | 0.0700 | 21 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 0.7619 | 0.2857 | 0.2381 | 0.0000 | 23 |
| BRCA1 | 54 | 5301 | 5400 | 1767 | 1800 | 0.0667 | 20 | 300 | functional_pathogenic_driven | BRCT phosphopeptide-binding region | 0.8000 | 0.0000 | 1.0000 | 0.0000 | 31 |
| BRCA2 | 26 | 2501 | 2600 | 834 | 867 | 0.0667 | 20 | 300 | truncation_driven | none | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 17 |
| BRCA2 | 42 | 4101 | 4200 | 1367 | 1400 | 0.0667 | 20 | 300 | truncation_driven | BRC repeats / RAD51-binding region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 18 |
| BRCA2 | 51 | 5001 | 5100 | 1667 | 1700 | 0.0667 | 20 | 300 | truncation_driven | BRC repeats / RAD51-binding region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 16 |
| BRCA2 | 64 | 6301 | 6400 | 2101 | 2134 | 0.0667 | 20 | 300 | truncation_driven | none | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 16 |
| BRCA2 | 87 | 8601 | 8700 | 2867 | 2900 | 0.0667 | 20 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 1.0000 | 0.0500 | 0.0000 | 0.0000 | 22 |
| BRCA1 | 34 | 3301 | 3400 | 1101 | 1134 | 0.0633 | 19 | 300 | truncation_driven | none | 1.0000 | 0.0526 | 0.0000 | 0.0000 | 20 |
| BRCA1 | 47 | 4601 | 4700 | 1534 | 1567 | 0.0633 | 19 | 300 | truncation_driven | none | 0.9474 | 0.1579 | 0.0526 | 0.0000 | 20 |
| BRCA2 | 13 | 1201 | 1300 | 401 | 434 | 0.0633 | 19 | 300 | truncation_driven | none | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 16 |
| BRCA2 | 32 | 3101 | 3200 | 1034 | 1067 | 0.0633 | 19 | 300 | truncation_driven | BRC repeats / RAD51-binding region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 15 |
| BRCA2 | 50 | 4901 | 5000 | 1634 | 1667 | 0.0633 | 19 | 300 | truncation_driven | BRC repeats / RAD51-binding region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 15 |
| BRCA2 | 85 | 8401 | 8500 | 2801 | 2834 | 0.0633 | 19 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 0.8947 | 0.2105 | 0.1579 | 0.0000 | 21 |
| BRCA2 | 89 | 8801 | 8900 | 2934 | 2967 | 0.0633 | 19 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 22 |
| BRCA2 | 92 | 9101 | 9200 | 3034 | 3067 | 0.0633 | 19 | 300 | truncation_driven | DNA-binding domain / helical-OB-DSS1 region | 0.9474 | 0.3684 | 0.0526 | 0.0000 | 22 |

## Interpretation

Regional pathogenic enrichment is not one biological phenomenon. In this
precomputed Module 1 map it can arise from several different sources:

- truncating possibilities and `PVS1`/`PM5_PTC`
- splice prediction and `PP3`
- functional evidence such as `PS3` or benign evidence such as `BS3`
- benign rule structure such as `BP1`, `BP4`, or `BP7`
- mixed codons and positions where different substitutions activate different
  criteria

The practical use is triage: a region with high pathogenic enrichment should
tell the curator what to inspect first, not how to classify the VUS.

## Outputs

- `tables/regional_driver_decomposition/regional_driver_bins.csv`
- `tables/regional_driver_decomposition/regional_signal_categories.csv`
- `tables/regional_driver_decomposition/regional_driver_ablation_summary.csv`
- `tables/regional_driver_decomposition/top_pathogenic_enriched_bins.csv`
- `plots/25_regional_driver_decomposition/brca1_regional_driver_tracks.svg`
- `plots/25_regional_driver_decomposition/brca2_regional_driver_tracks.svg`
- `plots/25_regional_driver_decomposition/regional_signal_categories.svg`
- `plots/25_regional_driver_decomposition/pathogenic_fraction_driver_ablation.svg`
