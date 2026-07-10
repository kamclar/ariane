# Domain And Region Context Comparison

Generated: 2026-06-23 17:49

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

This analysis replaces a direct BRCA1-versus-BRCA2 whole-gene comparison with a
domain and broad-region view. Direct whole-gene comparison is limited because
BRCA1 and BRCA2 are different proteins with different lengths, domain
architecture, exon structure, and structured evidence coverage.

The domain view is still descriptive. Domain membership is not an ACMG/ENIGMA
criterion by itself. It helps explain where generated Module 1 signals are
concentrated and where benign counterexamples exist.

## Outputs

- `tables/gene_comparison/gene_summary.csv`
- `tables/gene_comparison/domain_region_summary.csv`
- `tables/gene_comparison/domain_region_vus_examples.csv`
- `tables/gene_comparison/normalized_position_bins.csv`
- `tables/gene_comparison/brca1_vs_brca2_bin_differences.csv`
- `plots/07_gene_comparison/domain_region_signal_heatmap.svg`
- `plots/07_gene_comparison/domain_region_grouped_class_mix.svg`
- `plots/07_gene_comparison/domain_vs_domain_within_domain_percent_heatmap.svg`
- `plots/07_gene_comparison/domain_vs_domain_normalized_density_heatmap.svg`

## Whole-Gene Summary

| Metric | BRCA1 | BRCA2 |
| --- | ---: | ---: |
| Variants | 16776 | 30771 |
| Benign percent | 85.09% | 73.70% |
| VUS percent | 9.69% | 21.22% |
| Pathogenic percent | 5.23% | 5.08% |
| High SpliceAI percent | 3.31% | 2.11% |
| Missense percent | 74.34% | 74.00% |
| Nonsense percent | 4.79% | 5.18% |
| Synonymous percent | 20.82% | 20.79% |

## Domain And Region Summary

| gene | domain_context | domain_range | aa_length | n_variants | variants_per_100aa | benign_percent | vus_percent | pathogenic_percent | high_spliceai_percent | ps3_percent | bs3_percent |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | RING zinc-binding/E3 ligase | 2-101 | 100 | 900 | 900.00 | 64.00 | 27.11 | 8.89 | 17.44 | 19.67 | 67.22 |
| BRCA1 | coiled-coil PALB2 interaction | 1391-1424 | 34 | 306 | 900.00 | 7.19 | 87.25 | 5.56 | 5.23 | 0.00 | 4.58 |
| BRCA1 | BRCT phosphopeptide-binding | 1650-1857 | 208 | 1872 | 900.00 | 62.93 | 29.59 | 7.48 | 4.91 | 19.87 | 67.36 |
| BRCA1 | outside mapped domains | outside | 1522 | 13698 | 900.00 | 91.24 | 4.09 | 4.67 | 2.12 | 0.14 | 2.73 |
| BRCA2 | PALB2 binding N-terminal | 10-40 | 31 | 279 | 900.00 | 2.87 | 92.11 | 5.02 | 6.45 | 1.08 | 1.79 |
| BRCA2 | BRC repeats / RAD51-binding | 1000-2080 | 1081 | 9729 | 900.00 | 94.99 | 0.00 | 5.01 | 0.00 | 0.04 | 0.33 |
| BRCA2 | DNA-binding domain / helical-OB-DSS1 | 2481-3186 | 706 | 6354 | 900.00 | 4.03 | 90.27 | 5.70 | 5.93 | 1.86 | 2.71 |
| BRCA2 | C-terminal RAD51/NLS context | 3187-3418 | 232 | 2088 | 900.00 | 87.84 | 9.43 | 2.73 | 3.54 | 0.05 | 0.72 |
| BRCA2 | outside mapped domains | outside | 1369 | 12321 | 900.00 | 92.02 | 2.75 | 5.23 | 1.46 | 0.04 | 0.41 |

## BRCA1 Domains Versus BRCA2 Domains

The main comparison is domain-to-domain and region-to-region, not whole-gene
BRCA1 versus whole-gene BRCA2. Two normalizations are useful:

1. Within-domain percent asks: among generated SNVs in this domain, what
   fraction is benign, VUS, pathogenic, high-SpliceAI, PS3, BS3, and so on?
2. Density per 100 amino acids asks: how many generated signals occur per
   protein length unit in that domain?

The first view is best for class mix. The second view is useful when comparing
small regions such as BRCA1 RING or BRCA2 PALB2-binding context with larger
regions such as BRCA2 BRC repeats or the DNA-binding domain.

## VUS Examples Inside Mapped Domains

These are examples with a nonzero heuristic signal inside mapped domains. They
are not clinical recommendations and do not imply pathogenicity.

| gene | domain_context | score | reasons | c_notation | p_notation | spliceai_score |
| --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | BRCT phosphopeptide-binding | 70 | SpliceAI&gt;=0.20;PS3;PP3 | c.5044G&gt;A | p.(Glu1682Lys) | 0.57 |
| BRCA1 | BRCT phosphopeptide-binding | 70 | SpliceAI&gt;=0.20;PS3;PP3 | c.5123C&gt;A | p.(Ala1708Glu) | 0.22 |
| BRCA1 | BRCT phosphopeptide-binding | 70 | SpliceAI&gt;=0.20;PS3;PP3 | c.5123C&gt;G | p.(Ala1708Gly) | 0.4 |
| BRCA1 | BRCT phosphopeptide-binding | 70 | SpliceAI&gt;=0.20;PS3;PP3 | c.5453A&gt;G | p.(Asp1818Gly) | 0.35 |
| BRCA1 | RING zinc-binding/E3 ligase | 70 | SpliceAI&gt;=0.20;PS3;PP3 | c.190T&gt;G | p.(Cys64Gly) | 0.65 |
| BRCA2 | DNA-binding domain / helical-OB-DSS1 | 70 | SpliceAI&gt;=0.20;PS3;PP3 | c.9218A&gt;G | p.(Asp3073Gly) | 0.22 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.4987A&gt;C | p.(Met1663Leu) | 0.23 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.4992C&gt;A | p.(Leu1664=) | 0.43 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.4998C&gt;T | p.(Tyr1666=) | 0.23 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5007C&gt;T | p.(Ala1669=) | 0.21 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5027T&gt;C | p.(Leu1676Ser) | 0.21 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5044G&gt;C | p.(Glu1682Gln) | 0.27 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5046A&gt;T | p.(Glu1682Asp) | 0.43 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5048A&gt;G | p.(Glu1683Gly) | 0.38 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5051C&gt;G | p.(Thr1684Ser) | 0.32 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5074G&gt;T | p.(Asp1692Tyr) | 0.74 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5078C&gt;G | p.(Ala1693Gly) | 0.31 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5332G&gt;C | p.(Asp1778His) | 0.32 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5332G&gt;T | p.(Asp1778Tyr) | 0.58 |
| BRCA1 | BRCT phosphopeptide-binding | 45 | SpliceAI&gt;=0.20;PP3 | c.5339T&gt;A | p.(Leu1780Gln) | 0.25 |

## Legacy Relative-CDS Bin Comparison

The older normalized-bin comparison is retained as a descriptive quality check,
but it should not be interpreted as homologous domain comparison. Positive
difference means BRCA1 higher; negative difference means BRCA2 higher at the
same relative CDS percentage.

| bin | relative_start_percent | relative_end_percent | brca1_vus_percent | brca2_vus_percent | diff_vus_percent |
| --- | --- | --- | --- | --- | --- |
| 103 | 85.00 | 85.83 | 0.72 | 91.76 | -91.04 |
| 89 | 73.33 | 74.17 | 0.71 | 90.70 | -89.99 |
| 102 | 84.17 | 85.00 | 0.71 | 90.31 | -89.60 |
| 93 | 76.67 | 77.50 | 0.72 | 90.31 | -89.59 |
| 98 | 80.83 | 81.67 | 0.00 | 89.41 | -89.41 |
| 99 | 81.67 | 82.50 | 1.42 | 90.70 | -89.28 |
| 104 | 85.83 | 86.67 | 1.42 | 88.76 | -87.34 |
| 105 | 86.67 | 87.50 | 2.90 | 90.20 | -87.30 |

## Reading Notes

The most useful interpretation is within-region rather than gene-versus-gene.
For example, BRCA1 BRCT and BRCA2 DNA-binding/BRC regions are not equivalent
domains, but both can be inspected as functionally important contexts with
different Module 1 signal patterns. Domain context should always be combined
with exact variant consequence, accepted ACMG/ENIGMA evidence, functional/RNA
metadata, and benign counterexamples.
