# BS3 Functional Source Audit

Generated: 2026-06-22 10:42

## Sources

- Local ENIGMA Table 9 JSON: `backend\data\enigma_table9.json`
- Downloaded source: Findlay et al. 2018 Nature Supplementary Table 1,
  `external_sources/brca1_bs3_functional_evidence_audit_sources/findlay_2018_brca1_sge/Findlay_2018_Nature_Supplementary_Table_1_SGE_scores.xlsx`
- Downloaded source guide: Findlay et al. 2018 Nature Supplementary Information,
  `external_sources/brca1_bs3_functional_evidence_audit_sources/findlay_2018_brca1_sge/Findlay_2018_Nature_Supplementary_Information.pdf`
- Downloaded update source: Dace & Findlay 2023 interim report cited by ENIGMA Table 9,
  `external_sources/brca1_bs3_functional_evidence_audit_sources/dace_findlay_2023_interim/Dace_Findlay_2023_BRCA1_interim_report.pdf`

## Purpose

This audit checks whether BRCA1 `BS3` entries in the local ENIGMA Table 9
resource are supported variant-by-variant by Findlay et al. 2018 saturation
genome editing source rows, rather than being inferred only from broad domain
membership.

## Match Summary

| region | brca1_bs3_table9_count | matched_findlay_sge_count | matched_percent |
| --- | --- | --- | --- |
| BRCA1_BRCT_region | 1865 | 1861 | 99.79 |
| BRCA1_N_terminal_RING_region | 871 | 868 | 99.66 |
| other_BRCA1_region | 338 | 85 | 25.15 |

## Function Class Summary

| region | findlay_func_class | count |
| --- | --- | --- |
| BRCA1_BRCT_region | FUNC | 1861 |
| BRCA1_N_terminal_RING_region | FUNC | 844 |
| other_BRCA1_region | not_in_findlay | 253 |
| other_BRCA1_region | FUNC | 85 |
| BRCA1_N_terminal_RING_region | LOF | 24 |
| BRCA1_BRCT_region | not_in_findlay | 4 |
| BRCA1_N_terminal_RING_region | not_in_findlay | 3 |

## Preliminary Interpretation

Most BRCA1 `BS3` records in Table 9 appear to cite PMID:30209399, which is the
Findlay et al. BRCA1 saturation genome editing study. The audit therefore
supports treating these as variant-level functional evidence from a calibrated
assay source, not merely as a domain-level rule. However, this does not by
itself prove that the local application strength is always correct. The next
step is to inspect discordant or unmatched records, and to verify how ENIGMA
Table 9 converts Findlay functional classes and RNA scores into `BS3`.

## Outputs

- `tables/bs3_functional_source_audit/brca1_bs3_table9_vs_findlay_sge.csv`
- `tables/bs3_functional_source_audit/brca1_bs3_findlay_match_summary.csv`
- `tables/bs3_functional_source_audit/brca1_bs3_findlay_function_class_summary.csv`
- `external_sources/brca1_bs3_functional_evidence_audit_sources/findlay_2018_brca1_sge/Findlay_2018_Nature_Supplementary_Table_1_SGE_scores.normalized.csv`
