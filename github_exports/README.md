# GitHub Exports

This directory contains clean, shareable exports generated from the ARIANE
working analysis.

## Main Export

- `vus_exploration_module1/`: GitHub-ready directory for the BRCA1/2 Module 1
  VUS exploration.
- `vus_exploration_module1.zip`: optional archive for sending or storing a
  snapshot. The directory is the preferred source for GitHub.

## Refresh Workflow

From the ARIANE workspace root:

```powershell
python scripts\prepare_vus_github_export.py
```

This rebuilds `github_exports/vus_exploration_module1` from the current
analysis outputs. Commit or upload that directory to GitHub.

Only create the zip when a separate archive is needed:

```powershell
Compress-Archive -Path github_exports\vus_exploration_module1\* -DestinationPath github_exports\vus_exploration_module1.zip -Force
```
