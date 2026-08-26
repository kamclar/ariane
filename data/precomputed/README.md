# ARIANE Precomputed Snapshots

This directory stores versioned precomputed data artifacts for fast lookup and
analysis. These snapshots are not the source of truth for expert variant
classification unless their rule and data versions are explicitly accepted.

## BRCA Module 1 SNV Classification Snapshot

Files:

- `brca_module1_snv_classification_snapshot.index.json`
- `brca_module1_snv_classification_snapshot.metadata.json`

Source:

- `variant_space_scan/outputs/brca_module1_full_snv_classification.csv`
- `variant_space_scan/outputs/brca_module1_full_snv_classification.summary.json`

Scope:

- BRCA1/BRCA2 coding SNVs
- ARIANE Module 1 automated criteria
- Precomputed coordinates and historical reference-transcript SpliceAI scores

Intended use:

- Fast coordinate and notation lookup by `GENE:c.notation`
- Landscape analysis
- Baseline comparison between future ARIANE/data versions
- Optional display as a precomputed snapshot alongside live/current
  classification

Runtime does not use the embedded SpliceAI score or precomputed classification.
SpliceAI is computed on demand through the configured profile-pinned service.

Rebuild command:

```powershell
python scripts\build_classification_snapshot.py
```
