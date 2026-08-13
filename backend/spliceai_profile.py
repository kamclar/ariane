"""Immutable ENIGMA v1.2 SpliceAI scoring profile.

Classification code and cache builders import the same profile so that a
cache created with different SpliceAI parameters cannot be accepted silently.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROFILE_PATH = PROJECT_ROOT / "data" / "spliceai" / "enigma_v1_2_spliceai_profile.json"

_EXPECTED = {
    "profile_id": "enigma-brca-v1.2-appendix-j-spliceai-raw-10kb-v1",
    "genome_assembly": "GRCh38",
    "max_distance": 10000,
    "mask": 0,
    "annotation_subset": "basic",
    "transcript_policy": "reference_transcript",
    "aggregation": "maximum_raw_delta",
    "delta_score_fields": ["DS_AG", "DS_AL", "DS_DG", "DS_DL"],
    "reference_score_fields": [
        "DS_AG_REF", "DS_AL_REF", "DS_DG_REF", "DS_DL_REF",
    ],
    "alternate_score_fields": [
        "DS_AG_ALT", "DS_AL_ALT", "DS_DG_ALT", "DS_DL_ALT",
    ],
}


def _load_profile() -> dict[str, Any]:
    try:
        value = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    except Exception as exc:
        raise RuntimeError(
            f"ENIGMA SpliceAI scoring profile could not be loaded: {type(exc).__name__}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise RuntimeError("ENIGMA SpliceAI scoring profile must be a JSON object")
    mismatches = [
        key for key, expected in _EXPECTED.items() if value.get(key) != expected
    ]
    thresholds = value.get("thresholds") or {}
    if thresholds != {
        "bp4_max_inclusive": 0.1,
        "uninformative_above": 0.1,
        "uninformative_below": 0.2,
        "pp3_min_inclusive": 0.2,
    }:
        mismatches.append("thresholds")
    if mismatches:
        raise RuntimeError(
            "ENIGMA SpliceAI scoring profile does not match Appendix J: "
            + ", ".join(sorted(set(mismatches)))
        )
    return value


SPLICEAI_PROFILE = _load_profile()
SPLICEAI_PROFILE_SHA256 = hashlib.sha256(PROFILE_PATH.read_bytes()).hexdigest()
SPLICEAI_PROFILE_ID = str(SPLICEAI_PROFILE["profile_id"])
SPLICEAI_MAX_DISTANCE = int(SPLICEAI_PROFILE["max_distance"])
SPLICEAI_MASK = int(SPLICEAI_PROFILE["mask"])
SPLICEAI_ANNOTATION_SUBSET = str(SPLICEAI_PROFILE["annotation_subset"])
SPLICEAI_GENOME_ASSEMBLY = str(SPLICEAI_PROFILE["genome_assembly"])
SPLICEAI_TRANSCRIPT_POLICY_REQUIRED = str(SPLICEAI_PROFILE["transcript_policy"])
SPLICEAI_AGGREGATION = str(SPLICEAI_PROFILE["aggregation"])
SPLICEAI_DELTA_FIELDS = tuple(SPLICEAI_PROFILE["delta_score_fields"])
SPLICEAI_REFERENCE_FIELDS = tuple(SPLICEAI_PROFILE["reference_score_fields"])
SPLICEAI_ALTERNATE_FIELDS = tuple(SPLICEAI_PROFILE["alternate_score_fields"])
SPLICEAI_LOW_THRESHOLD = float(
    SPLICEAI_PROFILE["thresholds"]["bp4_max_inclusive"]
)
SPLICEAI_HIGH_THRESHOLD = float(
    SPLICEAI_PROFILE["thresholds"]["pp3_min_inclusive"]
)


def scoring_profile_metadata() -> dict[str, Any]:
    """Return the fields every immutable cache must repeat verbatim."""
    return {
        "scoring_profile_id": SPLICEAI_PROFILE_ID,
        "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
        "genome_assembly": SPLICEAI_GENOME_ASSEMBLY,
        "distance": SPLICEAI_MAX_DISTANCE,
        "mask": SPLICEAI_MASK,
        "annotation_subset": SPLICEAI_ANNOTATION_SUBSET,
        "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY_REQUIRED,
        "aggregation": SPLICEAI_AGGREGATION,
        "delta_score_fields": list(SPLICEAI_DELTA_FIELDS),
        "reference_score_fields": list(SPLICEAI_REFERENCE_FIELDS),
        "alternate_score_fields": list(SPLICEAI_ALTERNATE_FIELDS),
    }


def validate_scoring_metadata(metadata: dict[str, Any]) -> list[str]:
    """Return human-readable mismatches; an empty list means exact match."""
    required = scoring_profile_metadata()
    return [
        f"{key}={metadata.get(key)!r}, expected {expected!r}"
        for key, expected in required.items()
        if metadata.get(key) != expected
    ]
