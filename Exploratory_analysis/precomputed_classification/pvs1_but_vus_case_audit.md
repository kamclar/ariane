# PVS1 But VUS Case Audit

Generated: 2026-06-20

Variant: `BRCA2 c.9925G>T p.(Glu3309Ter)`

Input snapshot: `variant_space_scan/outputs/brca_module1_full_snv_classification.csv`

## Why This Case Was Reviewed

The criteria sanity audit found one variant with `PVS1` that still classified as VUS:

`BRCA2 c.9925G>T p.(Glu3309Ter)`

This is worth checking because `PVS1` is very strong pathogenic evidence. A VUS result can still be correct when there is contradictory benign evidence, but it should be explicit why the classification does not become likely pathogenic.

## Observed Automated Result

| field | value |
| --- | --- |
| gene | BRCA2 |
| c_notation | c.9925G>T |
| p_notation | p.(Glu3309Ter) |
| variant_type | nonsense |
| SpliceAI reference transcript score | 0.01 |
| GRCh37 | 13:32972575:G>T |
| GRCh38 | 13:32398438:G>T |
| applied criteria | BS3 Strong, PM2 Supporting, PVS1 Very Strong |
| total points | 5 |
| predicted class | 3 VUS |
| classification note | Contradictory benign and pathogenic evidence, point-based classification, expert review required |

## Evidence Sources In The Application

### BS3

`BS3 Strong` comes from `backend/data/enigma_table9.json`.

The local Table 9 entry for `BRCA2:c.9925G>T` states that this variant was reported by one calibrated study to exhibit protein function similar to benign control variants, with PMID `29988080`, and marks `BS3` as met.

This gives `-4` points.

### PVS1

`PVS1 Very Strong` comes from the Table 4 PTC rule in `backend/data/enigma_table4.json` and `backend/modules/table4.py`.

The relevant BRCA2 Table 4 critical boundary is:

| field | value |
| --- | --- |
| exon | E27 |
| boundary amino acid | 3309 |
| rule text in local JSON | PTC < p.Thr3310 -> PVS1, PTC > p.Glu3309 -> PVS1_N/A |

`p.(Glu3309Ter)` is at amino acid 3309. The implementation treats PTC at or before p.3309 as PVS1, because the truncation removes the critical C-terminal region. Therefore `PVS1 Very Strong` is applied.

This gives `+8` points.

### PM2 Supporting

In the full precomputed snapshot, the variant is absent from the local gnomAD lookup with adequate local coverage context, so `PM2_Supporting` is applied.

This gives `+1` point.

### SpliceAI

The precomputed reference-transcript SpliceAI cache gives a score of `0.01` for `BRCA2 c.9925G>T`, selected transcript `ENST00000380152.8` / `NM_000059.4`. No splice criterion is driving this result.

## Why The Final Class Is VUS

The evidence is contradictory:

| criterion | direction | strength | points |
| --- | --- | --- | --- |
| PVS1 | pathogenic | Very Strong | +8 |
| PM2_Supporting | pathogenic/supporting | Supporting | +1 |
| BS3 | benign | Strong | -4 |

Total: `8 + 1 - 4 = 5`.

Because there is both benign and pathogenic evidence, the application does not use a simple pathogenic ENIGMA combination such as `PVS1 + PM2_Supporting`. It switches to the point-based contradictory-evidence path and adds an expert-review note.

In that point system:

| total points | class |
| --- | --- |
| 6 to 9 | class 4 Likely Pathogenic |
| -1 to 5 | class 3 VUS |

The total is 5, so the automated result remains VUS.

## Interpretation

This looks like an expected edge case, not a parsing or implementation bug.

The variant sits exactly at the BRCA2 C-terminal PVS1 boundary, so the current PVS1 result is consistent with the local Table 4 rule. At the same time, Table 9 provides calibrated functional benign evidence (`BS3 Strong`). That conflict is strong enough to prevent automatic class 4.

The result should remain flagged for expert review. The curator question is not whether the code can add the points; it is whether `BS3 Strong` from the functional assay should be allowed to counterbalance `PVS1 Very Strong` for this specific terminal nonsense variant under the intended ENIGMA interpretation.

## Practical Outcome

Recommended handling in the analysis dataset:

- keep `PVS1_but_VUS` as a high-priority manual review flag
- do not treat this as a failed automated classification
- keep the classification as class 3 VUS unless curator review decides that the functional evidence or PVS1 boundary interpretation should be handled differently
- optionally add this variant as a regression/example case for contradictory evidence handling

