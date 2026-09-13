# Variant processing

This package owns HGVS parsing and normalization, local sequence providers,
reference validation, variant typing, and shared variant-coordinate helpers.
It consumes checksum-validated local snapshots through `backend/reference_data`.
It does not acquire classification evidence, call live lookups, or classify
evidence.
