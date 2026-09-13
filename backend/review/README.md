# Review

This package builds and validates structured expert-review workflows. Review
modules may propose or validate evidence but do not acquire remote data.

- `definitions.py` owns form text, presentation metadata and resource links.
- `strength.py` derives criterion strength from structured evidence and policy.
- `validation.py` reports form completeness without assigning evidence.
- `service.py` builds the amended working classification.
- `manual_evidence.py` is a compatibility import facade with no rule logic.
