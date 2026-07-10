# Hotspot and Coldspot Report

Generated: 2026-06-19 18:36

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

This is an exploratory bin-level scan over BRCA1 and BRCA2 coding SNVs. Each gene is split into 120 CDS bins. Hotspots are bins with at least 5 signal variants and at least 2x enrichment over the full coding-SNV background. Coldspots are bins with expected signal count at least 5 but observed count 0.

Signals:

- pathogenic: predicted classes 4 and 5
- vus: predicted class 3
- high_spliceai: SpliceAI >= 0.20
- pvs1, pm5_ptc, bs3, bp1, pp3: criterion present

Interpretation caveat: this is an automated Module 1 landscape, not final expert classification. Bin-level enrichment is descriptive and should be followed by variant-level review.

## Output Files

- `tables/hotspot_coldspot_bins_all_signals.csv`
- `tables/hotspots_<signal>.csv`
- `tables/coldspots_<signal>.csv`
- `plots/05_hotspots/hotspot_heatmap_all_signals.svg`
- `plots/05_hotspots/hotspot_heatmap_<signal>.svg`
- `plots/05_hotspots/top_hotspots_<signal>.svg`

## pathogenic

Total signal count: 2441

Hotspot bins found: 3
Coldspot bins found: 5

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA1 | 109 | 5033-5079 | 141 | 28 | 19.86% | 7.24 | 3.87x | 3.68e-09 |
| BRCA1 | 5 | 187-233 | 138 | 21 | 15.22% | 7.08 | 2.96x | 1.724e-05 |
| BRCA1 | 111 | 5127-5172 | 141 | 17 | 12.06% | 7.24 | 2.35x | 0.001362 |

Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 118 | 10001-10086 | 258 | 0 | 0.00% | 13.25 | 0.00x | 1 |
| BRCA2 | 119 | 10087-10171 | 255 | 0 | 0.00% | 13.09 | 0.00x | 1 |
| BRCA2 | 120 | 10172-10257 | 255 | 0 | 0.00% | 13.09 | 0.00x | 1 |
| BRCA1 | 14 | 606-652 | 141 | 0 | 0.00% | 7.24 | 0.00x | 0.9993 |
| BRCA1 | 13 | 560-605 | 138 | 0 | 0.00% | 7.08 | 0.00x | 0.9992 |

## vus

Total signal count: 8154

Hotspot bins found: 35
Coldspot bins found: 108

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 95 | 8035-8120 | 258 | 241 | 93.41% | 44.25 | 5.45x | 1.11e-16 |
| BRCA2 | 111 | 9403-9487 | 255 | 238 | 93.33% | 43.73 | 5.44x | 8.882e-16 |
| BRCA2 | 107 | 9061-9145 | 255 | 235 | 92.16% | 43.73 | 5.37x | 8.882e-16 |
| BRCA2 | 91 | 7693-7778 | 258 | 237 | 91.86% | 44.25 | 5.36x | 1.11e-16 |
| BRCA2 | 103 | 8719-8803 | 255 | 234 | 91.76% | 43.73 | 5.35x | 8.882e-16 |
| BRCA2 | 92 | 7779-7863 | 255 | 233 | 91.37% | 43.73 | 5.33x | 8.882e-16 |
| BRCA2 | 97 | 8206-8291 | 258 | 235 | 91.09% | 44.25 | 5.31x | 1.11e-16 |
| BRCA2 | 96 | 8121-8205 | 255 | 232 | 90.98% | 43.73 | 5.31x | 8.882e-16 |
| BRCA2 | 89 | 7522-7607 | 258 | 234 | 90.70% | 44.25 | 5.29x | 1.11e-16 |
| BRCA2 | 99 | 8377-8462 | 258 | 234 | 90.70% | 44.25 | 5.29x | 1.11e-16 |

Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 3 | 171-256 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 13 | 1026-1111 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 24 | 1966-2051 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 26 | 2137-2222 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 28 | 2308-2393 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 30 | 2479-2564 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 32 | 2650-2735 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 34 | 2821-2906 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 36 | 2992-3077 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |
| BRCA2 | 38 | 3163-3248 | 258 | 0 | 0.00% | 44.25 | 0.00x | 1 |

## high_spliceai

Total signal count: 1205

Hotspot bins found: 37
Coldspot bins found: 63

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA1 | 94 | 4334-4380 | 141 | 56 | 39.72% | 3.57 | 15.67x | 0 |
| BRCA1 | 5 | 187-233 | 138 | 46 | 33.33% | 3.50 | 13.15x | 0 |
| BRCA2 | 106 | 8975-9060 | 258 | 81 | 31.40% | 6.54 | 12.39x | 3.331e-16 |
| BRCA2 | 107 | 9061-9145 | 255 | 74 | 29.02% | 6.46 | 11.45x | 2.22e-16 |
| BRCA1 | 109 | 5033-5079 | 141 | 40 | 28.37% | 3.57 | 11.19x | 0 |
| BRCA1 | 6 | 234-279 | 141 | 39 | 27.66% | 3.57 | 10.91x | 0 |
| BRCA1 | 7 | 280-326 | 141 | 33 | 23.40% | 3.57 | 9.23x | 0 |
| BRCA2 | 114 | 9659-9744 | 258 | 45 | 17.44% | 6.54 | 6.88x | 3.331e-16 |
| BRCA1 | 88 | 4055-4100 | 138 | 24 | 17.39% | 3.50 | 6.86x | 6.358e-13 |
| BRCA1 | 4 | 140-186 | 141 | 24 | 17.02% | 3.57 | 6.72x | 9.905e-13 |

Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 13 | 1026-1111 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 17 | 1368-1453 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 19 | 1539-1624 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 26 | 2137-2222 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 28 | 2308-2393 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 30 | 2479-2564 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 32 | 2650-2735 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 34 | 2821-2906 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 36 | 2992-3077 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |
| BRCA2 | 38 | 3163-3248 | 258 | 0 | 0.00% | 6.54 | 0.00x | 0.9986 |

## pvs1

Total signal count: 2327

Hotspot bins found: 0
Coldspot bins found: 5

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| none | | | | | | | | |


Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 118 | 10001-10086 | 258 | 0 | 0.00% | 12.63 | 0.00x | 1 |
| BRCA2 | 119 | 10087-10171 | 255 | 0 | 0.00% | 12.48 | 0.00x | 1 |
| BRCA2 | 120 | 10172-10257 | 255 | 0 | 0.00% | 12.48 | 0.00x | 1 |
| BRCA1 | 14 | 606-652 | 141 | 0 | 0.00% | 6.90 | 0.00x | 0.999 |
| BRCA1 | 13 | 560-605 | 138 | 0 | 0.00% | 6.75 | 0.00x | 0.9988 |

## pm5_ptc

Total signal count: 2241

Hotspot bins found: 1
Coldspot bins found: 12

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 106 | 8975-9060 | 258 | 25 | 9.69% | 12.16 | 2.06x | 0.0008223 |

Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 114 | 9659-9744 | 258 | 0 | 0.00% | 12.16 | 0.00x | 1 |
| BRCA2 | 116 | 9830-9915 | 258 | 0 | 0.00% | 12.16 | 0.00x | 1 |
| BRCA2 | 118 | 10001-10086 | 258 | 0 | 0.00% | 12.16 | 0.00x | 1 |
| BRCA2 | 115 | 9745-9829 | 255 | 0 | 0.00% | 12.02 | 0.00x | 1 |
| BRCA2 | 117 | 9916-10000 | 255 | 0 | 0.00% | 12.02 | 0.00x | 1 |
| BRCA2 | 119 | 10087-10171 | 255 | 0 | 0.00% | 12.02 | 0.00x | 1 |
| BRCA2 | 120 | 10172-10257 | 255 | 0 | 0.00% | 12.02 | 0.00x | 1 |
| BRCA1 | 14 | 606-652 | 141 | 0 | 0.00% | 6.65 | 0.00x | 0.9987 |
| BRCA1 | 119 | 5499-5545 | 141 | 0 | 0.00% | 6.65 | 0.00x | 0.9987 |
| BRCA1 | 13 | 560-605 | 138 | 0 | 0.00% | 6.50 | 0.00x | 0.9985 |

## bs3

Total signal count: 2528

Hotspot bins found: 22
Coldspot bins found: 56

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA1 | 106 | 4894-4939 | 141 | 122 | 86.52% | 7.50 | 16.27x | 0 |
| BRCA1 | 116 | 5360-5405 | 141 | 117 | 82.98% | 7.50 | 15.61x | 0 |
| BRCA1 | 118 | 5453-5498 | 138 | 111 | 80.43% | 7.34 | 15.13x | 0 |
| BRCA1 | 107 | 4940-4986 | 141 | 107 | 75.89% | 7.50 | 14.27x | 0 |
| BRCA1 | 115 | 5313-5359 | 138 | 104 | 75.36% | 7.34 | 14.17x | 0 |
| BRCA1 | 1 | 1-46 | 141 | 105 | 74.47% | 7.50 | 14.01x | 0 |
| BRCA1 | 4 | 140-186 | 141 | 104 | 73.76% | 7.50 | 13.87x | 0 |
| BRCA1 | 5 | 187-233 | 138 | 100 | 72.46% | 7.34 | 13.63x | 0 |
| BRCA1 | 6 | 234-279 | 141 | 99 | 70.21% | 7.50 | 13.21x | 0 |
| BRCA1 | 112 | 5173-5219 | 141 | 99 | 70.21% | 7.50 | 13.21x | 0 |

Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 5 | 342-427 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 11 | 855-940 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 13 | 1026-1111 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 15 | 1197-1282 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 19 | 1539-1624 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 24 | 1966-2051 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 26 | 2137-2222 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 28 | 2308-2393 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 32 | 2650-2735 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |
| BRCA2 | 34 | 2821-2906 | 258 | 0 | 0.00% | 13.72 | 0.00x | 1 |

## bp1

Total signal count: 34855

Hotspot bins found: 0
Coldspot bins found: 42

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| none | | | | | | | | |


Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 89 | 7522-7607 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 91 | 7693-7778 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 93 | 7864-7949 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 95 | 8035-8120 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 97 | 8206-8291 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 99 | 8377-8462 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 102 | 8633-8718 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 104 | 8804-8889 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 106 | 8975-9060 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |
| BRCA2 | 108 | 9146-9231 | 258 | 0 | 0.00% | 189.13 | 0.00x | 1 |

## pp3

Total signal count: 1060

Hotspot bins found: 32
Coldspot bins found: 66

Top hotspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA1 | 94 | 4334-4380 | 141 | 50 | 35.46% | 3.14 | 15.91x | 1.11e-16 |
| BRCA1 | 5 | 187-233 | 138 | 41 | 29.71% | 3.08 | 13.33x | 0 |
| BRCA2 | 106 | 8975-9060 | 258 | 70 | 27.13% | 5.75 | 12.17x | 0 |
| BRCA2 | 107 | 9061-9145 | 255 | 67 | 26.27% | 5.68 | 11.79x | 0 |
| BRCA1 | 109 | 5033-5079 | 141 | 36 | 25.53% | 3.14 | 11.45x | 1.11e-16 |
| BRCA1 | 6 | 234-279 | 141 | 34 | 24.11% | 3.14 | 10.82x | 1.11e-16 |
| BRCA1 | 7 | 280-326 | 141 | 30 | 21.28% | 3.14 | 9.54x | 1.11e-16 |
| BRCA2 | 114 | 9659-9744 | 258 | 41 | 15.89% | 5.75 | 7.13x | 0 |
| BRCA1 | 4 | 140-186 | 141 | 21 | 14.89% | 3.14 | 6.68x | 2.744e-11 |
| BRCA1 | 17 | 746-792 | 141 | 20 | 14.18% | 3.14 | 6.36x | 1.847e-10 |

Top coldspots:

| Gene | Bin | CDS range | n | Count | Rate | Expected | Enrichment | Poisson tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA2 | 3 | 171-256 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 13 | 1026-1111 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 17 | 1368-1453 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 19 | 1539-1624 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 24 | 1966-2051 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 26 | 2137-2222 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 28 | 2308-2393 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 30 | 2479-2564 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 32 | 2650-2735 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |
| BRCA2 | 34 | 2821-2906 | 258 | 0 | 0.00% | 5.75 | 0.00x | 0.9968 |

