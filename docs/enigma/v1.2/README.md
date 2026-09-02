# ClinGen ENIGMA BRCA1/2 VCEP v1.2 source bundle

This directory contains the official files attached to ClinGen CSpec record
`GN092`, version `1.2.0`, released 2025-01-09.

The files under `source/` are stored without modification. Their source URLs,
sizes and SHA-256 checksums are pinned in `manifest.json`. Before using these
files to change classification logic, verify the bundle with:

```powershell
python scripts/verify_enigma_source_bundle.py
```

Primary registry record:

https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2.0

The current ClinGen registry record was checked on 2026-09-01. Its latest
approved specification is still version 1.2, accepted on 2025-01-09. The
corresponding Zenodo record was created on 2026-07-18 as a distribution record
for the same approved version. It is not a new VCEP rules release:

https://zenodo.org/records/21434315

All five current CSpec attachments were downloaded again on 2026-09-01 and
matched the sizes and SHA-256 checksums in `manifest.json`. The UCSC ENIGMA
track update dated 2026-08-18 is a separate variant-level data release. It does
not change the version of the VCEP specification stored in this directory.

The local bundle is evidence for implementation review. It does not permit the
application to silently choose between contradictory statements in different
official files. Such conflicts must be documented and resolved explicitly.
