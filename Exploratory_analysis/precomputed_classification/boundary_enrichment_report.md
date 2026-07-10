# Boundary Enrichment Analysis

Generated: 2026-06-19 21:14

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

This analysis refines the earlier boundary-distance view. Internal CDS exon starts are treated as acceptor-like boundaries and internal CDS exon ends are treated as donor-like boundaries in transcript coordinate space. The first CDS start and last CDS end are not treated as splice junctions.

Because the snapshot contains coding SNVs only, this is an exonic-side approximation. It does not include intronic variants.

## Outputs

- `tables/boundary_enrichment/variant_boundary_annotations.csv`
- `tables/boundary_enrichment/boundary_distance_group_summary.csv`
- `tables/boundary_enrichment/donor_acceptor_site_summary.csv`
- `tables/boundary_enrichment/relative_position_summary.csv`
- `tables/boundary_enrichment/threshold_sensitivity_summary.csv`
- `tables/boundary_enrichment/exon_boundary_summary.csv`
- `tables/boundary_enrichment/variant_type_boundary_summary.csv`
- `tables/boundary_enrichment/boundary_conflict_sanity_checks.csv`
- `tables/boundary_enrichment/strongest_boundary_vus_pathogenic_variants.csv`
- `plots/09_boundary_enrichment/donor_acceptor_near_boundary_rate.svg`
- `plots/09_boundary_enrichment/boundary_distance_rate_heatmap.svg`
- `plots/09_boundary_enrichment/relative_position_spliceai_rate.svg`
- `plots/09_boundary_enrichment/threshold_sensitivity_near_boundary.svg`
- `plots/09_boundary_enrichment/brca1_exon_boundary_heatmap.svg`
- `plots/09_boundary_enrichment/brca2_exon_boundary_heatmap.svg`
- `plots/09_boundary_enrichment/variant_type_boundary_rate_heatmap.svg`

## Donor-like versus Acceptor-like Summary

| gene | site_type | total_variants | spliceai_ge_0_20 | percent | benign_percent | vus_percent | pathogenic_percent |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | acceptor_like | 8469 | 250 | 2.95 | 0.65 | 16.88 | 11.30 |
| BRCA1 | donor_like | 8307 | 306 | 3.68 | 1.08 | 21.01 | 18.00 |
| BRCA2 | acceptor_like | 16215 | 424 | 2.61 | 0.03 | 10.71 | 6.54 |
| BRCA2 | donor_like | 14556 | 225 | 1.55 | 0.02 | 6.29 | 3.77 |

## Highest Near-boundary Rates

Near-boundary here means distance bin 0 or 1-2 bp.

| gene | site_type | group | distance_bin | total_variants | spliceai_ge_0_20 | percent |
| --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | donor_like | pathogenic | 0 | 16 | 16 | 100.00 |
| BRCA2 | donor_like | vus | 0 | 60 | 51 | 85.00 |
| BRCA2 | donor_like | pathogenic | 0 | 9 | 7 | 77.78 |
| BRCA1 | donor_like | vus | 0 | 30 | 21 | 70.00 |
| BRCA1 | acceptor_like | pathogenic | 0 | 8 | 5 | 62.50 |
| BRCA1 | donor_like | pathogenic | 1-2 | 18 | 11 | 61.11 |
| BRCA1 | acceptor_like | vus | 1-2 | 41 | 22 | 53.66 |
| BRCA2 | acceptor_like | pathogenic | 1-2 | 4 | 2 | 50.00 |
| BRCA1 | acceptor_like | vus | 0 | 29 | 14 | 48.28 |
| BRCA1 | acceptor_like | pathogenic | 1-2 | 12 | 5 | 41.67 |
| BRCA1 | donor_like | vus | 1-2 | 39 | 16 | 41.03 |
| BRCA2 | donor_like | vus | 1-2 | 77 | 31 | 40.26 |
| BRCA2 | acceptor_like | vus | 0 | 47 | 14 | 29.79 |
| BRCA2 | donor_like | pathogenic | 1-2 | 24 | 6 | 25.00 |
| BRCA1 | donor_like | benign | 0 | 17 | 3 | 17.65 |
| BRCA2 | acceptor_like | vus | 1-2 | 88 | 15 | 17.05 |
| BRCA1 | acceptor_like | benign | 0 | 26 | 2 | 7.69 |
| BRCA1 | donor_like | benign | 1-2 | 69 | 4 | 5.80 |

## Highest Exact Relative-position Rates

Only rows with at least 3 variants are shown here.

| gene | site_type | group | relative_boundary_position | total_variants | spliceai_ge_0_20 | percent |
| --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | donor_like | pathogenic | 0 | 16 | 16 | 100.00 |
| BRCA2 | donor_like | pathogenic | -1 | 4 | 4 | 100.00 |
| BRCA2 | donor_like | vus | 0 | 60 | 51 | 85.00 |
| BRCA2 | donor_like | pathogenic | 0 | 9 | 7 | 77.78 |
| BRCA1 | acceptor_like | pathogenic | 12 | 4 | 3 | 75.00 |
| BRCA1 | acceptor_like | pathogenic | 1 | 7 | 5 | 71.43 |
| BRCA1 | donor_like | pathogenic | -1 | 7 | 5 | 71.43 |
| BRCA1 | donor_like | vus | 0 | 30 | 21 | 70.00 |
| BRCA1 | acceptor_like | pathogenic | 5 | 3 | 2 | 66.67 |
| BRCA1 | acceptor_like | pathogenic | 6 | 3 | 2 | 66.67 |
| BRCA1 | donor_like | pathogenic | -18 | 3 | 2 | 66.67 |
| BRCA1 | donor_like | pathogenic | -11 | 3 | 2 | 66.67 |
| BRCA2 | acceptor_like | pathogenic | 18 | 3 | 2 | 66.67 |
| BRCA1 | acceptor_like | pathogenic | 0 | 8 | 5 | 62.50 |
| BRCA1 | donor_like | pathogenic | -20 | 8 | 5 | 62.50 |
| BRCA1 | acceptor_like | vus | 2 | 19 | 11 | 57.89 |
| BRCA1 | donor_like | vus | -1 | 21 | 12 | 57.14 |
| BRCA1 | donor_like | pathogenic | -2 | 11 | 6 | 54.55 |

## Boundary Conflict Sanity Checks

Conflict rows: 222

| boundary_conflict_reason | gene | c_notation | p_notation | group | spliceai_score | nearest_site_type | relative_boundary_position |
| --- | --- | --- | --- | --- | --- | --- | --- |
| benign_high_spliceai_near_boundary | BRCA1 | c.5196T&gt;G | p.(His1732Gln) | benign | 0.99 | acceptor_like | 2 |
| benign_high_spliceai_near_boundary | BRCA1 | c.5470A&gt;G | p.(Ile1824Val) | benign | 0.93 | acceptor_like | 2 |
| benign_high_spliceai_near_boundary | BRCA2 | c.6842G&gt;T | p.(Gly2281Val) | benign | 0.92 | acceptor_like | 0 |
| benign_high_spliceai_near_boundary | BRCA1 | c.5073A&gt;G | p.(Thr1691=) | benign | 0.71 | donor_like | -1 |
| benign_high_spliceai_near_boundary | BRCA1 | c.5196T&gt;A | p.(His1732Gln) | benign | 0.71 | acceptor_like | 2 |
| benign_high_spliceai_near_boundary | BRCA1 | c.83T&gt;A | p.(Leu28Gln) | benign | 0.57 | acceptor_like | 2 |
| benign_high_spliceai_near_boundary | BRCA1 | c.213G&gt;C | p.(Arg71Ser) | benign | 0.39 | acceptor_like | 0 |
| benign_high_spliceai_near_boundary | BRCA1 | c.135A&gt;T | p.(Lys45Asn) | benign | 0.37 | acceptor_like | 0 |
| benign_high_spliceai_near_boundary | BRCA1 | c.301T&gt;G | p.(Tyr101Asp) | benign | 0.32 | donor_like | 0 |
| benign_high_spliceai_near_boundary | BRCA1 | c.5277G&gt;T | p.(Lys1759Asn) | benign | 0.28 | donor_like | 0 |
| benign_high_spliceai_near_boundary | BRCA1 | c.4985T&gt;G | p.(Phe1662Cys) | benign | 0.27 | donor_like | -1 |
| benign_high_spliceai_near_boundary | BRCA1 | c.210A&gt;G | p.(Lys70=) | benign | 0.24 | donor_like | -2 |
| benign_high_spliceai_near_boundary | BRCA1 | c.5277G&gt;C | p.(Lys1759Asn) | benign | 0.24 | donor_like | 0 |
| benign_high_spliceai_near_boundary | BRCA1 | c.300G&gt;A | p.(Glu100=) | benign | 0.22 | donor_like | -1 |
| benign_high_spliceai_near_boundary | BRCA2 | c.6935A&gt;T | p.(Asp2312Val) | benign | 0.22 | donor_like | -2 |
| pathogenic_low_spliceai_near_boundary | BRCA1 | c.4099G&gt;T | p.(Glu1367Ter) | pathogenic | 0.07 | acceptor_like | 2 |
| pathogenic_low_spliceai_near_boundary | BRCA2 | c.473C&gt;G | p.(Ser158Ter) | pathogenic | 0.07 | donor_like | -2 |
| pathogenic_low_spliceai_near_boundary | BRCA2 | c.7615C&gt;T | p.(Gln2539Ter) | pathogenic | 0.07 | donor_like | -2 |

## Strongest Boundary VUS/Pathogenic Variants

Variants are VUS or pathogenic-group, within 2 bp of an internal CDS splice boundary, and SpliceAI >=0.20.

| gene | c_notation | p_notation | group | spliceai_score | nearest_site_type | relative_boundary_position | criteria |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA2 | c.8754G&gt;C | p.(Glu2918Asp) | vus | 0.99 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.8754G&gt;A | p.(Glu2918=) | vus | 0.98 | donor_like | 0 | PP3:Supporting:1 |
| BRCA2 | c.425G&gt;T | p.(Ser142Ile) | vus | 0.97 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.1909G&gt;T | p.(Gly637Cys) | vus | 0.97 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA1 | c.4185G&gt;T | p.(Gln1395His) | vus | 0.96 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.65C&gt;T | p.(Ala22Val) | vus | 0.96 | donor_like | -2 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA1 | c.4185G&gt;C | p.(Gln1395His) | vus | 0.95 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.1909G&gt;C | p.(Gly637Arg) | vus | 0.95 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.6937G&gt;T | p.(Gly2313Cys) | vus | 0.95 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.7435G&gt;T | p.(Asp2479Tyr) | vus | 0.95 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.8487G&gt;C | p.(Gln2829His) | vus | 0.95 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.8487G&gt;T | p.(Gln2829His) | vus | 0.95 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.8754G&gt;T | p.(Glu2918Asp) | vus | 0.95 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.425G&gt;C | p.(Ser142Thr) | vus | 0.94 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.475G&gt;T | p.(Val159Leu) | vus | 0.94 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.477G&gt;A | p.(Val159=) | vus | 0.94 | acceptor_like | 1 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.516G&gt;T | p.(Lys172Asn) | vus | 0.94 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |
| BRCA2 | c.6937G&gt;A | p.(Gly2313Ser) | vus | 0.94 | donor_like | 0 | PM2_Supporting:Supporting:1;PP3:Supporting:1 |

## Reading Notes

If donor-like and acceptor-like rates are similar, the main signal is general proximity to a splice boundary. If one side is consistently higher, that suggests side-specific splice disruption patterns worth reviewing against RNA or curated rules.
