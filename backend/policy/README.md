# Policy

This package is the shared policy layer for the backend.

- `gene.py` loads and validates the versioned gene and VCEP policy manifest.
- `classification.py` combines evidence into an ENIGMA classification.
- `spliceai.py` defines rule-facing SpliceAI applicability, completeness, and
  provenance decisions.

The package must not import classification DAG nodes, evidence lookups, review
builders, or presentation code. Those layers depend on policy, not the other
way around.
