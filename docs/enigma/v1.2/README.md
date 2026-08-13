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

The local bundle is evidence for implementation review. It does not permit the
application to silently choose between contradictory statements in different
official files. Such conflicts must be documented and resolved explicitly.
