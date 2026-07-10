# External Sources For Functional Evidence Audit

This directory stores external source files used to audit structured functional
evidence in the ARIANE Module 1 exploratory analyses.

The first audit target is `BS3` in BRCA1 RING and BRCT regions, because the
generated map contains very high `BS3` coverage in these domains. The goal is
to determine whether `BS3` is supported variant-by-variant by source functional
assay data, or whether any local mapping is too broad.

Files should be kept separate from the ARIANE application backend data until
they have been reviewed and explicitly promoted.

## Directory Map

### `brca1_bs3_functional_evidence_audit_sources`

External source package for auditing BRCA1 `BS3` evidence in RING/BRCT regions.

Contains:

- `findlay_2018_brca1_sge`
  - Findlay et al. 2018 Nature supplementary information PDF.
  - Findlay et al. 2018 Supplementary Table 1 XLSX with BRCA1 saturation genome editing scores.
  - A normalized CSV copy derived from the XLSX for local audit scripts.
- `dace_findlay_2023_interim`
  - Dace & Findlay 2023 interim report cited in ENIGMA Table 9 for updated interpretation of selected variants.

These files are audit sources only. They are not application runtime data.
