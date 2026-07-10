# VUS in Pathogenic Regions

Generated: 2026-06-21 10:30

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Purpose

This analysis asks which VUS lie in regions dominated by generated
pathogenic/likely pathogenic variants and which VUS are locally closer to the
pathogenic group than to the benign/likely benign group.

This is a prioritization and interpretation layer only. It does not change any
classification.

## Definitions

- pathogenic group: predicted classes 4 and 5
- benign group: predicted classes 1 and 2
- VUS group: predicted class 3
- in the local-neighborhood scatter plots, each dot is one VUS
- the VUS is the center of the window; `+/-20 bp` means 20 coding bases
  upstream and 20 coding bases downstream of that VUS position
- x-axis: number of benign-group variants within that window around the VUS
- y-axis: number of pathogenic-group variants within that window around the VUS
- broad pathogenic context: VUS is in a pathogenic hotspot bin, or has at least 5
  pathogenic-group variants within +/-20 coding bases, or at least 10 within
  +/-50 coding bases, or at least 20 within +/-100 coding bases
- strong pathogenic region: VUS is in a pathogenic hotspot bin, or the local
  pathogenic group clearly outnumbers the benign group by stricter thresholds
- closer to pathogenic: nearest class 4/5 variant is closer than nearest class
  1/2 variant in coding position

## Summary

| summary_type | label | count |
| --- | --- | --- |
| closer_to_group | closer_to_benign | 3853 |
| closer_to_group | closer_to_pathogenic | 3480 |
| closer_to_group | equal_distance | 821 |
| broad_pathogenic_context | broad_pathogenic_context | 7938 |
| broad_pathogenic_context | no_broad_pathogenic_context | 216 |
| strong_pathogenic_region | not_strong_pathogenic_region | 5853 |
| strong_pathogenic_region | strong_pathogenic_region | 2301 |
| intersection | strong_pathogenic_region_and_closer_to_pathogenic | 1451 |

## Top VUS In Pathogenic Regions And Closer To Pathogenic

| pathogenic_proximity_score | gene | c_notation | p_notation | variant_type | spliceai_score | strong_pathogenic_region | broad_pathogenic_context | closer_to_group | nearest_pathogenic_distance | nearest_benign_distance | pathogenic_count_20bp | benign_count_20bp | criteria |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 203 | BRCA2 | c.9007G&gt;C | p.(Gly3003Arg) | missense | 0.28 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 0 | 5 | 16 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 203 | BRCA2 | c.9007G&gt;A | p.(Gly3003Arg) | missense | 0.2 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 0 | 5 | 16 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 200 | BRCA2 | c.9001A&gt;G | p.(Thr3001Ala) | missense | 0.5 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 4 | 16 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 200 | BRCA2 | c.9001A&gt;T | p.(Thr3001Ser) | missense | 0.47 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 4 | 16 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 200 | BRCA2 | c.9001A&gt;C | p.(Thr3001Pro) | missense | 0.46 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 4 | 16 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 200 | BRCA2 | c.9000A&gt;T | p.(Leu3000Phe) | missense | 0.34 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 3 | 16 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 200 | BRCA2 | c.9000A&gt;C | p.(Leu3000Phe) | missense | 0.24 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 3 | 16 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 197 | BRCA2 | c.9009A&gt;T | p.(Gly3003=) | synonymous | 0.76 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 3 | 14 | 2 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 195 | BRCA2 | c.8993C&gt;A | p.(Ser2998Tyr) | missense | 0.32 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 4 | 15 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 195 | BRCA2 | c.8993C&gt;G | p.(Ser2998Cys) | missense | 0.27 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 4 | 15 | 3 | PP3:Supporting:1 |
| 195 | BRCA2 | c.8993C&gt;T | p.(Ser2998Phe) | missense | 0.2 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 4 | 15 | 3 | PP3:Supporting:1 |
| 193 | BRCA2 | c.9006A&gt;T | p.(Glu3002Asp) | missense | 0.35 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 6 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 193 | BRCA2 | c.9008G&gt;T | p.(Gly3003Val) | missense | 0.29 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 4 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 193 | BRCA2 | c.9006A&gt;C | p.(Glu3002Asp) | missense | 0.24 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 6 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 191 | BRCA2 | c.9005A&gt;G | p.(Glu3002Gly) | missense | 0.37 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 7 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 191 | BRCA2 | c.9005A&gt;T | p.(Glu3002Val) | missense | 0.31 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 7 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 191 | BRCA2 | c.9005A&gt;C | p.(Glu3002Ala) | missense | 0.25 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 7 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 190 | BRCA2 | c.8990A&gt;T | p.(Tyr2997Phe) | missense | 0.28 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 2 | 15 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 190 | BRCA2 | c.8991T&gt;C | p.(Tyr2997=) | synonymous | 0.24 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 0 | 3 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 189 | BRCA2 | c.8983G&gt;T | p.(Asp2995Tyr) | missense | 0.29 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 5 | 12 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 189 | BRCA2 | c.8983G&gt;A | p.(Asp2995Asn) | missense | 0.21 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 5 | 12 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 189 | BRCA2 | c.8983G&gt;C | p.(Asp2995His) | missense | 0.21 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 5 | 12 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 188 | BRCA2 | c.8981C&gt;T | p.(Ser2994Leu) | missense | 0.26 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 0 | 7 | 12 | 4 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 186 | BRCA2 | c.9002C&gt;A | p.(Thr3001Lys) | missense | 0.56 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 5 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 186 | BRCA2 | c.9002C&gt;G | p.(Thr3001Arg) | missense | 0.56 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 5 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 186 | BRCA2 | c.9003A&gt;C | p.(Thr3001=) | synonymous | 0.56 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 6 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 186 | BRCA2 | c.9003A&gt;G | p.(Thr3001=) | synonymous | 0.56 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 6 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 186 | BRCA2 | c.9003A&gt;T | p.(Thr3001=) | synonymous | 0.56 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 1 | 6 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 186 | BRCA2 | c.9004G&gt;C | p.(Glu3002Gln) | missense | 0.56 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 0 | 7 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| 186 | BRCA2 | c.9002C&gt;T | p.(Thr3001Ile) | missense | 0.48 | strong_pathogenic_dense_20bp;strong_pathogenic_dense_50bp | pathogenic_dense_20bp;pathogenic_dense_50bp;pathogenic_dense_100bp | closer_to_pathogenic | 2 | 5 | 14 | 3 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |

## Interpretation

These variants are useful manual-review candidates because their local generated
classification landscape is more pathogenic than benign. The signal can come
from dense pathogenic neighborhoods, pathogenic hotspot bins, or shorter
distance to class 4/5 variants than to class 1/2 variants.

This does not prove pathogenicity. It can reflect real biological constraint,
rule behavior in the automated snapshot, local splice effects, or clustering of
criteria such as PS3, PP3, PVS1, or PM5_PTC.

## Outputs

- `tables/vus_pathogenic_regions/vus_pathogenic_region_annotated.csv`
- `tables/vus_pathogenic_regions/vus_pathogenic_region_top_candidates.csv`
- `tables/vus_pathogenic_regions/vus_pathogenic_region_summary.csv`
- `plots/18_vus_pathogenic_regions/vus_pathogenic_region_summary.svg`
- `plots/18_vus_pathogenic_regions/vus_pathogenic_vs_benign_neighborhood_20bp.svg`
- `plots/18_vus_pathogenic_regions/vus_pathogenic_vs_benign_neighborhood_50bp.svg`
