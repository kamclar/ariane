# BRCA1 BRCT1 Mixed Local Neighborhood

Generated: 2026-06-21 10:31

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Why This Region Was Flagged

In the VUS local-neighborhood scatter plot, a visible green cluster appears in
the benign-enriched sector. Those dots are VUS with many neighboring generated
variants of both directions, but with more class 1/2 than class 4/5 neighbors
within +/-20 coding bp.

The cluster maps mainly to:

- gene: `BRCA1`
- coding range: `c.5044-c.5143`
- protein range: approximately `p.1682-p.1715`
- local feature context: BRCA1 BRCT phosphopeptide-binding region, inside the
  UniProt BRCT 1 subdomain

This is a prioritization/interpretation finding only. It does not reclassify
any VUS.

## What The Region Represents

This range is inside the BRCA1 C-terminal BRCT region. The generated landscape
is mixed because different substitutions in the same short region can trigger
different automated evidence patterns:

- many missense substitutions receive benign-direction evidence such as `BS3`
  or no strong pathogenic evidence
- some substitutions have `PP3` because of local SpliceAI signal
- some substitutions have `PS3`
- nearby nonsense substitutions receive truncation-related pathogenic evidence
  such as `PVS1` and `PM5_PTC`

So the region is not uniformly benign or pathogenic. It is a mixed BRCT1
neighborhood where the exact nucleotide/protein consequence matters.

## All Generated SNVs In The Region

Total generated SNVs in region: 300

Grouped classes:

| field | value | count |
| --- | --- | --- |
| class_group | benign | 139 |
| class_group | vus | 114 |
| class_group | pathogenic | 47 |

Exact classes:

| field | value | count |
| --- | --- | --- |
| predicted_class | 2 | 138 |
| predicted_class | 3 | 114 |
| predicted_class | 4 | 25 |
| predicted_class | 5 | 22 |
| predicted_class | 1 | 1 |

Variant types:

| field | value | count |
| --- | --- | --- |
| variant_type | missense | 216 |
| variant_type | synonymous | 69 |
| variant_type | nonsense | 15 |

SpliceAI bins:

| spliceai_bin | count |
| --- | --- |
| &gt;=0.50 | 14 |
| 0.20-0.49 | 32 |
| 0.10-0.19 | 77 |
| &lt;0.10 | 177 |

Most frequent criteria:

| criterion | count |
| --- | --- |
| PM2_Supporting | 273 |
| BS3 | 151 |
| PS3 | 119 |
| PP3 | 41 |
| PVS1 | 15 |
| PM5_PTC | 15 |
| PP4 | 7 |
| BP5 | 3 |
| BP4 | 3 |
| BP7 | 3 |

## VUS In The Region

VUS in region: 114

VUS SpliceAI bins:

| spliceai_bin | count |
| --- | --- |
| &gt;=0.50 | 2 |
| 0.20-0.49 | 7 |
| 0.10-0.19 | 37 |
| &lt;0.10 | 68 |

VUS criteria:

| criterion | count |
| --- | --- |
| PM2_Supporting | 93 |
| PS3 | 72 |
| BS3 | 12 |
| PP3 | 9 |
| BP4 | 3 |
| BP7 | 3 |
| BP5 | 1 |

Top VUS by local benign/pathogenic neighbor counts:

| gene | c_notation | p_notation | variant_type | spliceai_score | benign_count_20bp | pathogenic_count_20bp | criteria |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | c.5048A&gt;G | p.(Glu1683Gly) | missense | 0.38 | 75 | 15 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA1 | c.5046A&gt;T | p.(Glu1682Asp) | missense | 0.43 | 74 | 16 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | missense | 0.57 | 73 | 15 | BP5:Strong:-4;PM2_Supporting:Supporting:1;PP3:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5044G&gt;C | p.(Glu1682Gln) | missense | 0.27 | 73 | 15 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA1 | c.5050A&gt;C | p.(Thr1684Pro) | missense | 0.03 | 73 | 15 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5051C&gt;G | p.(Thr1684Ser) | missense | 0.32 | 72 | 15 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA1 | c.5051C&gt;A | p.(Thr1684Asn) | missense | 0.14 | 72 | 15 | PM2_Supporting:Supporting:1 |
| BRCA1 | c.5052T&gt;G | p.(Thr1684=) | synonymous | 0.18 | 70 | 17 | BS3:Strong:-4 |
| BRCA1 | c.5052T&gt;C | p.(Thr1684=) | synonymous | 0.01 | 70 | 17 | BS3:Strong:-4 |
| BRCA1 | c.5138T&gt;A | p.(Val1713Glu) | missense | 0.09 | 69 | 15 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5138T&gt;C | p.(Val1713Ala) | missense | 0.09 | 69 | 15 | PS3:Strong:4 |
| BRCA1 | c.5138T&gt;G | p.(Val1713Gly) | missense | 0.08 | 69 | 15 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5131A&gt;G | p.(Lys1711Glu) | missense | 0.05 | 69 | 15 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5053A&gt;C | p.(Thr1685Pro) | missense | 0.09 | 68 | 19 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5143A&gt;G | p.(Ser1715Gly) | missense | 0.08 | 66 | 16 | PM2_Supporting:Supporting:1 |
| BRCA1 | c.5143A&gt;T | p.(Ser1715Cys) | missense | 0.07 | 66 | 16 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5141T&gt;C | p.(Val1714Ala) | missense | 0.04 | 66 | 16 | PM2_Supporting:Supporting:1 |
| BRCA1 | c.5141T&gt;A | p.(Val1714Asp) | missense | 0.02 | 66 | 16 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5140G&gt;T | p.(Val1714Phe) | missense | 0.05 | 66 | 15 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5134T&gt;G | p.(Trp1712Gly) | missense | 0.03 | 65 | 17 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5142T&gt;G | p.(Val1714=) | synonymous | 0.08 | 65 | 16 | BS3:Strong:-4 |
| BRCA1 | c.5120T&gt;G | p.(Ile1707Ser) | missense | 0.03 | 64 | 13 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5128G&gt;A | p.(Gly1710Arg) | missense | 0.06 | 63 | 17 |  |
| BRCA1 | c.5119A&gt;T | p.(Ile1707Phe) | missense | 0.01 | 62 | 13 | PM2_Supporting:Supporting:1;PS3:Strong:4 |
| BRCA1 | c.5074G&gt;T | p.(Asp1692Tyr) | missense | 0.74 | 61 | 23 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |

## Boundary Context

Boundary annotation rows in this region: 300

Nearest boundary bin summary:

| field | value | count |
| --- | --- | --- |
| splice_boundary_distance_bin | 21-50 | 138 |
| splice_boundary_distance_bin | 11-20 | 90 |
| splice_boundary_distance_bin | 6-10 | 36 |
| splice_boundary_distance_bin | 3-5 | 18 |
| splice_boundary_distance_bin | 1-2 | 12 |
| splice_boundary_distance_bin | 0 | 6 |

Donor/acceptor-like context:

| field | value | count |
| --- | --- | --- |
| nearest_site_type | donor_like | 183 |
| nearest_site_type | acceptor_like | 117 |

## Interpretation

This cluster is interesting because it sits in an important functional BRCT
region but is locally benign-enriched in the generated landscape. The green
position in the scatter plot should therefore be interpreted as:

- not "these VUS are benign"
- but "this BRCT1 subregion has many generated benign/likely benign neighbors,
  despite also containing pathogenic-direction evidence nearby"

The practical manual-review question is whether each VUS is driven by:

- protein functional evidence (`PS3`/`BS3`)
- splice prediction (`PP3`/SpliceAI)
- nearby truncating consequences
- or a mixed rule pattern that needs curated functional assay metadata

## Outputs

- `tables/brca1_brct1_mixed_cluster/all_region_variants.csv`
- `tables/brca1_brct1_mixed_cluster/vus_region_variants.csv`
- `tables/brca1_brct1_mixed_cluster/boundary_region_variants.csv`
- `plots/23_brca1_brct1_mixed_cluster/brca1_brct1_mixed_cluster_class_profile.svg`
