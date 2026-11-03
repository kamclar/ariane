# Experiment artifacts

This directory belongs only to the SpliceAI search-distance experiment. It is
not read by the ARIANE classifier at runtime.

- `distance_comparison_20260915.json` is the versioned pilot result.
- `spliceai_distance_experiment_report.pdf` is the generated report.
- `build/` contains local resumable checkpoints and raw arrays. It is ignored
  by Git because the contents are large and reproducible from the pinned
  experiment inputs.
