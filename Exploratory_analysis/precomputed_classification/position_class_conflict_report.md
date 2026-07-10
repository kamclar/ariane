# Position And Codon Class Conflict Analysis

Generated: 2026-06-21 10:03

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Purpose

This analysis separates a coding position from a concrete variant. A single
coding position can have up to three possible SNV substitutions, and those
substitutions can receive different generated ARIANE classes. The same is even
more true at codon level, where up to nine SNVs can affect one amino-acid
position.

This matters for VUS interpretation because a local pathogenic-looking region
does not mean every possible substitution at that position is pathogenic. The
opposite is also true: a benign variant at a position does not prove that a
different substitution at the same position is benign.

This is an exploratory sanity-check of the generated Module 1 snapshot only. It
does not reclassify variants.

## Angles Tested

1. Same coding nucleotide position: do different SNVs at the same `c.` position
   fall into different grouped classes?
2. Same codon: do different SNVs in the same codon fall into different grouped
   classes?
3. VUS inside mixed contexts: are there VUS at positions or codons where both
   benign and pathogenic generated variants also exist?
4. Mechanistic hints: do mixed positions differ by criteria, variant type,
   points, or SpliceAI score?

## Position-Level Summary

| level | pattern | count |
| --- | --- | --- |
| cds_position | benign_only | 10793 |
| cds_position | vus_only | 1945 |
| cds_position | benign_and_pathogenic | 1553 |
| cds_position | benign_and_vus | 1064 |
| cds_position | pathogenic_and_vus | 386 |
| cds_position | all_three_groups | 101 |
| cds_position | pathogenic_only | 7 |

## Codon-Level Summary

| level | pattern | count |
| --- | --- | --- |
| codon | benign_only | 2401 |
| codon | benign_and_pathogenic | 1394 |
| codon | benign_and_vus | 581 |
| codon | all_three_groups | 343 |
| codon | vus_only | 317 |
| codon | pathogenic_and_vus | 247 |

## Same Coding Position With Both Benign And Pathogenic Generated Variants

| gene | cds_pos | codon | group_pattern | variant_count | benign_count | vus_count | pathogenic_count | class_set | variant_types | criteria_codes | max_spliceai | example_variants |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | 8 | 3 | benign_and_pathogenic | 3 | 1 | 0 | 2 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.090 | c.8T&gt;A p.(Leu3Ter) [5]; c.8T&gt;C p.(Leu3Ser) [2]; c.8T&gt;G p.(Leu3Ter) [5] |
| BRCA1 | 25 | 9 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.220 | c.25G&gt;A p.(Glu9Lys) [2]; c.25G&gt;C p.(Glu9Gln) [2]; c.25G&gt;T p.(Glu9Ter) [5] |
| BRCA1 | 28 | 10 | all_three_groups | 3 | 1 | 1 | 1 | 2+3+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.460 | c.28G&gt;A p.(Glu10Lys) [3]; c.28G&gt;C p.(Glu10Gln) [2]; c.28G&gt;T p.(Glu10Ter) [5] |
| BRCA1 | 34 | 12 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.130 | c.34C&gt;A p.(Gln12Lys) [2]; c.34C&gt;G p.(Gln12Glu) [2]; c.34C&gt;T p.(Gln12Ter) [5] |
| BRCA1 | 55 | 19 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.310 | c.55C&gt;A p.(Gln19Lys) [2]; c.55C&gt;G p.(Gln19Glu) [2]; c.55C&gt;T p.(Gln19Ter) [5] |
| BRCA1 | 58 | 20 | all_three_groups | 3 | 1 | 1 | 1 | 2+3+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.200 | c.58A&gt;C p.(Lys20Gln) [3]; c.58A&gt;G p.(Lys20Glu) [2]; c.58A&gt;T p.(Lys20Ter) [5] |
| BRCA1 | 67 | 23 | all_three_groups | 3 | 1 | 1 | 1 | 2+3+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.090 | c.67G&gt;A p.(Glu23Lys) [2]; c.67G&gt;C p.(Glu23Gln) [3]; c.67G&gt;T p.(Glu23Ter) [5] |
| BRCA1 | 72 | 24 | all_three_groups | 3 | 1 | 1 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.020 | c.72T&gt;A p.(Cys24Ter) [5]; c.72T&gt;C p.(Cys24=) [2]; c.72T&gt;G p.(Cys24Trp) [3] |
| BRCA1 | 85 | 29 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.150 | c.85G&gt;A p.(Glu29Lys) [2]; c.85G&gt;C p.(Glu29Gln) [2]; c.85G&gt;T p.(Glu29Ter) [5] |
| BRCA1 | 89 | 30 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.060 | c.89T&gt;A p.(Leu30Ter) [5]; c.89T&gt;C p.(Leu30Ser) [2]; c.89T&gt;G p.(Leu30Trp) [2] |
| BRCA1 | 94 | 32 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.150 | c.94A&gt;C p.(Lys32Gln) [2]; c.94A&gt;G p.(Lys32Glu) [2]; c.94A&gt;T p.(Lys32Ter) [5] |
| BRCA1 | 97 | 33 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.060 | c.97G&gt;A p.(Glu33Lys) [2]; c.97G&gt;C p.(Glu33Gln) [2]; c.97G&gt;T p.(Glu33Ter) [5] |
| BRCA1 | 112 | 38 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.070 | c.112A&gt;C p.(Lys38Gln) [2]; c.112A&gt;G p.(Lys38Glu) [2]; c.112A&gt;T p.(Lys38Ter) [5] |
| BRCA1 | 117 | 39 | all_three_groups | 3 | 1 | 1 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.020 | c.117T&gt;A p.(Cys39Ter) [5]; c.117T&gt;C p.(Cys39=) [2]; c.117T&gt;G p.(Cys39Trp) [3] |
| BRCA1 | 133 | 45 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 1+2+5 | missense+nonsense | BP5+BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.190 | c.133A&gt;C p.(Lys45Gln) [1]; c.133A&gt;G p.(Lys45Glu) [2]; c.133A&gt;T p.(Lys45Ter) [5] |
| BRCA1 | 148 | 50 | all_three_groups | 3 | 1 | 1 | 1 | 2+3+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.140 | c.148A&gt;C p.(Lys50Gln) [2]; c.148A&gt;G p.(Lys50Glu) [3]; c.148A&gt;T p.(Lys50Ter) [5] |
| BRCA1 | 160 | 54 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.130 | c.160C&gt;A p.(Gln54Lys) [2]; c.160C&gt;G p.(Gln54Glu) [2]; c.160C&gt;T p.(Gln54Ter) [5] |
| BRCA1 | 163 | 55 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.730 | c.163A&gt;C p.(Lys55Gln) [2]; c.163A&gt;G p.(Lys55Glu) [2]; c.163A&gt;T p.(Lys55Ter) [5] |
| BRCA1 | 166 | 56 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.520 | c.166A&gt;C p.(Lys56Gln) [2]; c.166A&gt;G p.(Lys56Glu) [2]; c.166A&gt;T p.(Lys56Ter) [5] |
| BRCA1 | 176 | 59 | benign_and_pathogenic | 3 | 1 | 0 | 2 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.410 | c.176C&gt;A p.(Ser59Ter) [5]; c.176C&gt;G p.(Ser59Ter) [5]; c.176C&gt;T p.(Ser59Leu) [2] |
| BRCA1 | 178 | 60 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.290 | c.178C&gt;A p.(Gln60Lys) [2]; c.178C&gt;G p.(Gln60Glu) [2]; c.178C&gt;T p.(Gln60Ter) [5] |
| BRCA1 | 188 | 63 | benign_and_pathogenic | 3 | 1 | 0 | 2 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.290 | c.188T&gt;A p.(Leu63Ter) [5]; c.188T&gt;C p.(Leu63Ser) [2]; c.188T&gt;G p.(Leu63Ter) [5] |
| BRCA1 | 192 | 64 | benign_and_pathogenic | 3 | 1 | 0 | 2 | 2+4+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.870 | c.192T&gt;A p.(Cys64Ter) [5]; c.192T&gt;C p.(Cys64=) [2]; c.192T&gt;G p.(Cys64Trp) [4] |
| BRCA1 | 193 | 65 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+5 | missense+nonsense | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.870 | c.193A&gt;C p.(Lys65Gln) [2]; c.193A&gt;G p.(Lys65Glu) [2]; c.193A&gt;T p.(Lys65Ter) [5] |
| BRCA1 | 196 | 66 | benign_and_pathogenic | 3 | 2 | 0 | 1 | 2+4 | missense | BS3+PM2_Supporting+PP3+PS3 | 0.650 | c.196A&gt;C p.(Asn66His) [2]; c.196A&gt;G p.(Asn66Asp) [2]; c.196A&gt;T p.(Asn66Tyr) [4] |

## Same Codon With Both Benign And Pathogenic Generated Variants

| gene | codon | cds_positions | group_pattern | variant_count | benign_count | vus_count | pathogenic_count | class_set | variant_types | criteria_codes | max_spliceai | example_variants |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | 3 | 7+8+9 | all_three_groups | 9 | 6 | 1 | 2 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.090 | c.7T&gt;A p.(Leu3Ile) [2]; c.7T&gt;C p.(Leu3=) [2]; c.7T&gt;G p.(Leu3Val) [2]; c.8T&gt;A p.(Leu3Ter) [5]; c.8T&gt;C p.(Leu3Ser) [2]; c.8T&gt;G p.(Leu3Ter) [5]; c.9A&gt;C p.(Leu3Phe) [2]; c.9A&gt;G p.(Leu3=) [3] |
| BRCA1 | 9 | 25+26+27 | benign_and_pathogenic | 9 | 8 | 0 | 1 | 2+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.220 | c.25G&gt;A p.(Glu9Lys) [2]; c.25G&gt;C p.(Glu9Gln) [2]; c.25G&gt;T p.(Glu9Ter) [5]; c.26A&gt;C p.(Glu9Ala) [2]; c.26A&gt;G p.(Glu9Gly) [2]; c.26A&gt;T p.(Glu9Val) [2]; c.27A&gt;C p.(Glu9Asp) [2]; c.27A&gt;G p.(Glu9=) [2] |
| BRCA1 | 10 | 28+29+30 | all_three_groups | 9 | 6 | 2 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.460 | c.28G&gt;A p.(Glu10Lys) [3]; c.28G&gt;C p.(Glu10Gln) [2]; c.28G&gt;T p.(Glu10Ter) [5]; c.29A&gt;C p.(Glu10Ala) [2]; c.29A&gt;G p.(Glu10Gly) [2]; c.29A&gt;T p.(Glu10Val) [3]; c.30A&gt;C p.(Glu10Asp) [2]; c.30A&gt;G p.(Glu10=) [2] |
| BRCA1 | 12 | 34+35+36 | all_three_groups | 9 | 6 | 2 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.130 | c.34C&gt;A p.(Gln12Lys) [2]; c.34C&gt;G p.(Gln12Glu) [2]; c.34C&gt;T p.(Gln12Ter) [5]; c.35A&gt;C p.(Gln12Pro) [3]; c.35A&gt;G p.(Gln12Arg) [2]; c.35A&gt;T p.(Gln12Leu) [2]; c.36A&gt;C p.(Gln12His) [2]; c.36A&gt;G p.(Gln12=) [3] |
| BRCA1 | 18 | 52+53+54 | all_three_groups | 9 | 5 | 3 | 1 | 2+3+5 | missense | BS3+PM2_Supporting+PP4+PS3 | 0.190 | c.52A&gt;C p.(Met18Leu) [2]; c.52A&gt;G p.(Met18Val) [3]; c.52A&gt;T p.(Met18Leu) [2]; c.53T&gt;A p.(Met18Lys) [3]; c.53T&gt;C p.(Met18Thr) [5]; c.53T&gt;G p.(Met18Arg) [3]; c.54G&gt;A p.(Met18Ile) [2]; c.54G&gt;C p.(Met18Ile) [2] |
| BRCA1 | 19 | 55+56+57 | all_three_groups | 9 | 7 | 1 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.310 | c.55C&gt;A p.(Gln19Lys) [2]; c.55C&gt;G p.(Gln19Glu) [2]; c.55C&gt;T p.(Gln19Ter) [5]; c.56A&gt;C p.(Gln19Pro) [3]; c.56A&gt;G p.(Gln19Arg) [2]; c.56A&gt;T p.(Gln19Leu) [2]; c.57G&gt;A p.(Gln19=) [2]; c.57G&gt;C p.(Gln19His) [2] |
| BRCA1 | 20 | 58+59+60 | all_three_groups | 9 | 6 | 2 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.200 | c.58A&gt;C p.(Lys20Gln) [3]; c.58A&gt;G p.(Lys20Glu) [2]; c.58A&gt;T p.(Lys20Ter) [5]; c.59A&gt;C p.(Lys20Thr) [2]; c.59A&gt;G p.(Lys20Arg) [2]; c.59A&gt;T p.(Lys20Ile) [2]; c.60A&gt;C p.(Lys20Asn) [3]; c.60A&gt;G p.(Lys20=) [2] |
| BRCA1 | 22 | 64+65+66 | all_three_groups | 9 | 6 | 1 | 2 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.070 | c.64T&gt;A p.(Leu22Ile) [2]; c.64T&gt;C p.(Leu22=) [2]; c.64T&gt;G p.(Leu22Val) [2]; c.65T&gt;A p.(Leu22Ter) [5]; c.65T&gt;C p.(Leu22Ser) [3]; c.65T&gt;G p.(Leu22Ter) [5]; c.66A&gt;C p.(Leu22Phe) [2]; c.66A&gt;G p.(Leu22=) [2] |
| BRCA1 | 23 | 67+68+69 | all_three_groups | 9 | 3 | 5 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.110 | c.67G&gt;A p.(Glu23Lys) [2]; c.67G&gt;C p.(Glu23Gln) [3]; c.67G&gt;T p.(Glu23Ter) [5]; c.68A&gt;C p.(Glu23Ala) [2]; c.68A&gt;G p.(Glu23Gly) [3]; c.68A&gt;T p.(Glu23Val) [2]; c.69G&gt;A p.(Glu23=) [3]; c.69G&gt;C p.(Glu23Asp) [3] |
| BRCA1 | 24 | 70+71+72 | all_three_groups | 9 | 1 | 7 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.040 | c.70T&gt;A p.(Cys24Ser) [3]; c.70T&gt;C p.(Cys24Arg) [3]; c.70T&gt;G p.(Cys24Gly) [3]; c.71G&gt;A p.(Cys24Tyr) [3]; c.71G&gt;C p.(Cys24Ser) [3]; c.71G&gt;T p.(Cys24Phe) [3]; c.72T&gt;A p.(Cys24Ter) [5]; c.72T&gt;C p.(Cys24=) [2] |
| BRCA1 | 29 | 85+86+87 | all_three_groups | 9 | 7 | 1 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.160 | c.85G&gt;A p.(Glu29Lys) [2]; c.85G&gt;C p.(Glu29Gln) [2]; c.85G&gt;T p.(Glu29Ter) [5]; c.86A&gt;C p.(Glu29Ala) [2]; c.86A&gt;G p.(Glu29Gly) [3]; c.86A&gt;T p.(Glu29Val) [2]; c.87G&gt;A p.(Glu29=) [2]; c.87G&gt;C p.(Glu29Asp) [2] |
| BRCA1 | 30 | 88+89+90 | all_three_groups | 9 | 7 | 1 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.060 | c.88T&gt;A p.(Leu30Met) [2]; c.88T&gt;C p.(Leu30=) [2]; c.88T&gt;G p.(Leu30Val) [2]; c.89T&gt;A p.(Leu30Ter) [5]; c.89T&gt;C p.(Leu30Ser) [2]; c.89T&gt;G p.(Leu30Trp) [2]; c.90G&gt;A p.(Leu30=) [2]; c.90G&gt;C p.(Leu30Phe) [2] |
| BRCA1 | 32 | 94+95+96 | benign_and_pathogenic | 9 | 8 | 0 | 1 | 2+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.160 | c.94A&gt;C p.(Lys32Gln) [2]; c.94A&gt;G p.(Lys32Glu) [2]; c.94A&gt;T p.(Lys32Ter) [5]; c.95A&gt;C p.(Lys32Thr) [2]; c.95A&gt;G p.(Lys32Arg) [2]; c.95A&gt;T p.(Lys32Met) [2]; c.96G&gt;A p.(Lys32=) [2]; c.96G&gt;C p.(Lys32Asn) [2] |
| BRCA1 | 33 | 97+98+99 | all_three_groups | 9 | 7 | 1 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.070 | c.97G&gt;A p.(Glu33Lys) [2]; c.97G&gt;C p.(Glu33Gln) [2]; c.97G&gt;T p.(Glu33Ter) [5]; c.98A&gt;C p.(Glu33Ala) [2]; c.98A&gt;G p.(Glu33Gly) [3]; c.98A&gt;T p.(Glu33Val) [2]; c.99A&gt;C p.(Glu33Asp) [2]; c.99A&gt;G p.(Glu33=) [2] |
| BRCA1 | 38 | 112+113+114 | benign_and_pathogenic | 9 | 8 | 0 | 1 | 1+2+5 | missense+nonsense+synonymous | BA1+BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.100 | c.112A&gt;C p.(Lys38Gln) [2]; c.112A&gt;G p.(Lys38Glu) [2]; c.112A&gt;T p.(Lys38Ter) [5]; c.113A&gt;C p.(Lys38Thr) [2]; c.113A&gt;G p.(Lys38Arg) [2]; c.113A&gt;T p.(Lys38Met) [2]; c.114G&gt;A p.(Lys38=) [1]; c.114G&gt;C p.(Lys38Asn) [2] |
| BRCA1 | 39 | 115+116+117 | all_three_groups | 9 | 1 | 7 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.020 | c.115T&gt;A p.(Cys39Ser) [3]; c.115T&gt;C p.(Cys39Arg) [3]; c.115T&gt;G p.(Cys39Gly) [3]; c.116G&gt;A p.(Cys39Tyr) [3]; c.116G&gt;C p.(Cys39Ser) [3]; c.116G&gt;T p.(Cys39Phe) [3]; c.117T&gt;A p.(Cys39Ter) [5]; c.117T&gt;C p.(Cys39=) [2] |
| BRCA1 | 41 | 121+122+123 | all_three_groups | 9 | 1 | 7 | 1 | 2+3+5 | missense+synonymous | BS3+PM2_Supporting+PP4+PS3 | 0.020 | c.121C&gt;A p.(His41Asn) [3]; c.121C&gt;G p.(His41Asp) [3]; c.121C&gt;T p.(His41Tyr) [3]; c.122A&gt;C p.(His41Pro) [3]; c.122A&gt;G p.(His41Arg) [5]; c.122A&gt;T p.(His41Leu) [3]; c.123C&gt;A p.(His41Gln) [3]; c.123C&gt;G p.(His41Gln) [3] |
| BRCA1 | 45 | 133+134+135 | all_three_groups | 9 | 7 | 1 | 1 | 1+2+3+5 | missense+nonsense+synonymous | BP5+BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.370 | c.133A&gt;C p.(Lys45Gln) [1]; c.133A&gt;G p.(Lys45Glu) [2]; c.133A&gt;T p.(Lys45Ter) [5]; c.134A&gt;C p.(Lys45Thr) [2]; c.134A&gt;G p.(Lys45Arg) [2]; c.134A&gt;T p.(Lys45Ile) [3]; c.135A&gt;C p.(Lys45Asn) [2]; c.135A&gt;G p.(Lys45=) [2] |
| BRCA1 | 50 | 148+149+150 | all_three_groups | 9 | 7 | 1 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.140 | c.148A&gt;C p.(Lys50Gln) [2]; c.148A&gt;G p.(Lys50Glu) [3]; c.148A&gt;T p.(Lys50Ter) [5]; c.149A&gt;C p.(Lys50Thr) [2]; c.149A&gt;G p.(Lys50Arg) [2]; c.149A&gt;T p.(Lys50Ile) [2]; c.150A&gt;C p.(Lys50Asn) [2]; c.150A&gt;G p.(Lys50=) [2] |
| BRCA1 | 54 | 160+161+162 | benign_and_pathogenic | 9 | 8 | 0 | 1 | 2+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.700 | c.160C&gt;A p.(Gln54Lys) [2]; c.160C&gt;G p.(Gln54Glu) [2]; c.160C&gt;T p.(Gln54Ter) [5]; c.161A&gt;C p.(Gln54Pro) [2]; c.161A&gt;G p.(Gln54Arg) [2]; c.161A&gt;T p.(Gln54Leu) [2]; c.162G&gt;A p.(Gln54=) [2]; c.162G&gt;C p.(Gln54His) [2] |
| BRCA1 | 55 | 163+164+165 | benign_and_pathogenic | 9 | 8 | 0 | 1 | 2+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.730 | c.163A&gt;C p.(Lys55Gln) [2]; c.163A&gt;G p.(Lys55Glu) [2]; c.163A&gt;T p.(Lys55Ter) [5]; c.164A&gt;C p.(Lys55Thr) [2]; c.164A&gt;G p.(Lys55Arg) [2]; c.164A&gt;T p.(Lys55Met) [2]; c.165G&gt;A p.(Lys55=) [2]; c.165G&gt;C p.(Lys55Asn) [2] |
| BRCA1 | 56 | 166+167+168 | benign_and_pathogenic | 9 | 8 | 0 | 1 | 2+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.520 | c.166A&gt;C p.(Lys56Gln) [2]; c.166A&gt;G p.(Lys56Glu) [2]; c.166A&gt;T p.(Lys56Ter) [5]; c.167A&gt;C p.(Lys56Thr) [2]; c.167A&gt;G p.(Lys56Arg) [2]; c.167A&gt;T p.(Lys56Ile) [2]; c.168A&gt;C p.(Lys56Asn) [2]; c.168A&gt;G p.(Lys56=) [2] |
| BRCA1 | 59 | 175+176+177 | benign_and_pathogenic | 9 | 7 | 0 | 2 | 2+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.410 | c.175T&gt;A p.(Ser59Thr) [2]; c.175T&gt;C p.(Ser59Pro) [2]; c.175T&gt;G p.(Ser59Ala) [2]; c.176C&gt;A p.(Ser59Ter) [5]; c.176C&gt;G p.(Ser59Ter) [5]; c.176C&gt;T p.(Ser59Leu) [2]; c.177A&gt;C p.(Ser59=) [2]; c.177A&gt;G p.(Ser59=) [2] |
| BRCA1 | 60 | 178+179+180 | all_three_groups | 9 | 6 | 2 | 1 | 2+3+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 0.290 | c.178C&gt;A p.(Gln60Lys) [2]; c.178C&gt;G p.(Gln60Glu) [2]; c.178C&gt;T p.(Gln60Ter) [5]; c.179A&gt;C p.(Gln60Pro) [3]; c.179A&gt;G p.(Gln60Arg) [3]; c.179A&gt;T p.(Gln60Leu) [2]; c.180G&gt;A p.(Gln60=) [2]; c.180G&gt;C p.(Gln60His) [2] |
| BRCA1 | 63 | 187+188+189 | benign_and_pathogenic | 9 | 7 | 0 | 2 | 2+5 | missense+nonsense+synonymous | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 0.700 | c.187T&gt;A p.(Leu63Ile) [2]; c.187T&gt;C p.(Leu63=) [2]; c.187T&gt;G p.(Leu63Val) [2]; c.188T&gt;A p.(Leu63Ter) [5]; c.188T&gt;C p.(Leu63Ser) [2]; c.188T&gt;G p.(Leu63Ter) [5]; c.189A&gt;C p.(Leu63Phe) [2]; c.189A&gt;G p.(Leu63=) [2] |

## VUS In Fully Mixed Contexts

Same-position contexts with VUS plus benign plus pathogenic generated variants:
101

Same-codon contexts with VUS plus benign plus pathogenic generated variants:
343

These are useful manual-review candidates because the local context is not
uniform. For such variants, the exact nucleotide change, protein consequence,
splice signal, and criterion combination matter more than the position label.

## What Drives The Mixed Contexts?

Same-position mixed contexts by variant-type combination:

| field | value | count |
| --- | --- | --- |
| variant_types | missense+nonsense | 1352 |
| variant_types | missense+nonsense+synonymous | 187 |
| variant_types | nonsense+synonymous | 82 |
| variant_types | missense | 23 |
| variant_types | missense+synonymous | 6 |
| variant_types | synonymous | 4 |

Same-codon mixed contexts by variant-type combination:

| field | value | count |
| --- | --- | --- |
| variant_types | missense+nonsense+synonymous | 1684 |
| variant_types | missense+synonymous | 34 |
| variant_types | missense+nonsense | 15 |
| variant_types | missense | 4 |

Same-position mixed contexts by criterion-code combination:

| field | value | count |
| --- | --- | --- |
| criteria_codes | BP1+PM2_Supporting+PM5_PTC+PVS1 | 1337 |
| criteria_codes | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 66 |
| criteria_codes | BP1+PM2_Supporting+PVS1 | 47 |
| criteria_codes | BP1+BP5+PM2_Supporting+PM5_PTC+PVS1 | 31 |
| criteria_codes | BP1+BS3+PM2_Supporting+PM5_PTC+PVS1 | 20 |
| criteria_codes | BP1+PM5_PTC+PVS1 | 20 |
| criteria_codes | BS3+PM2_Supporting+PP3+PS3 | 20 |
| criteria_codes | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 13 |
| criteria_codes | BS3+PM2_Supporting+PM5_PTC+PVS1 | 11 |
| criteria_codes | BS3+PM2_Supporting+PS3+PVS1 | 10 |
| criteria_codes | BP1+BP5+BS3+PM2_Supporting+PM5_PTC+PVS1 | 8 |
| criteria_codes | BP1+BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 8 |
| criteria_codes | BS3+PM2_Supporting+PP4+PS3 | 8 |
| criteria_codes | BP1+PM2_Supporting+PM5_PTC+PP3+PVS1 | 6 |
| criteria_codes | BP1+PM2_Supporting+PM5_PTC+PS3+PVS1 | 6 |
| criteria_codes | BP4+BP7+PM2_Supporting+PM5_PTC+PVS1 | 5 |
| criteria_codes | BP1+BS1_Supporting+PM2_Supporting+PM5_PTC+PVS1 | 4 |
| criteria_codes | BP5+BS3+PM2_Supporting+PM5_PTC+PVS1 | 3 |
| criteria_codes | BP5+PM2_Supporting+PM5_PTC+PVS1 | 3 |
| criteria_codes | BA1+BP1+PM2_Supporting+PM5_PTC+PVS1 | 2 |

Same-codon mixed contexts by criterion-code combination:

| field | value | count |
| --- | --- | --- |
| criteria_codes | BP1+PM2_Supporting+PM5_PTC+PVS1 | 1165 |
| criteria_codes | BP1+BP5+PM2_Supporting+PM5_PTC+PVS1 | 62 |
| criteria_codes | BP1+PM2_Supporting+PM5_PTC+PP3+PVS1 | 58 |
| criteria_codes | BP1+BS3+PM2_Supporting+PM5_PTC+PVS1 | 55 |
| criteria_codes | BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 52 |
| criteria_codes | BP1+PM2_Supporting+PVS1 | 37 |
| criteria_codes | BS3+PM2_Supporting+PM5_PTC+PP3+PS3+PVS1 | 27 |
| criteria_codes | BP4+BP7+BS3+PM2_Supporting+PM5_PTC+PVS1 | 21 |
| criteria_codes | BP1+PM5_PTC+PVS1 | 20 |
| criteria_codes | BP1+BS1_Supporting+PM2_Supporting+PM5_PTC+PVS1 | 17 |
| criteria_codes | BP1+BP5+BS3+PM2_Supporting+PM5_PTC+PVS1 | 16 |
| criteria_codes | BP4+BP7+PM2_Supporting+PM5_PTC+PVS1 | 16 |
| criteria_codes | BP1+PM2_Supporting+PP3+PVS1 | 14 |
| criteria_codes | BS3+PM2_Supporting+PP3+PS3 | 13 |
| criteria_codes | BP1+BS1_Strong+PM2_Supporting+PM5_PTC+PVS1 | 12 |
| criteria_codes | BA1+BP1+PM2_Supporting+PM5_PTC+PVS1 | 9 |
| criteria_codes | BS3+PM2_Supporting+PP4+PS3 | 8 |
| criteria_codes | BP1+BP5+BS1_Supporting+PM2_Supporting+PM5_PTC+PVS1 | 7 |
| criteria_codes | BS3+PM2_Supporting+PS3+PVS1 | 7 |
| criteria_codes | BP1+BS3+PM2_Supporting+PM5_PTC+PS3+PVS1 | 6 |

The dominant pattern should be read biologically and technically: many mixed
positions are expected because one possible nucleotide substitution creates a
premature stop codon while another creates a missense or synonymous consequence.
The more interesting subset is where the mixed context is driven by splice
signal, functional evidence, PS1/PS3, or other non-truncating mechanisms.

## Outputs

- `tables/position_class_conflicts/position_class_summary.csv`
- `tables/position_class_conflicts/codon_class_summary.csv`
- `tables/position_class_conflicts/same_position_benign_pathogenic.csv`
- `tables/position_class_conflicts/same_codon_benign_pathogenic.csv`
- `tables/position_class_conflicts/vus_with_mixed_position_context.csv`
- `tables/position_class_conflicts/vus_with_mixed_codon_context.csv`
- `tables/position_class_conflicts/position_mixed_variant_type_summary.csv`
- `tables/position_class_conflicts/codon_mixed_variant_type_summary.csv`
- `tables/position_class_conflicts/position_mixed_criteria_summary.csv`
- `tables/position_class_conflicts/codon_mixed_criteria_summary.csv`
- `plots/22_position_class_conflicts/position_group_patterns.svg`
- `plots/22_position_class_conflicts/codon_group_patterns.svg`
- `plots/22_position_class_conflicts/same_position_mixed_examples.svg`
- `plots/22_position_class_conflicts/same_codon_mixed_examples.svg`
