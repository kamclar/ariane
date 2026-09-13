"""Immutable ENIGMA v1.2 SpliceAI scoring profile.

Classification code and cache builders import the same profile so that a
cache created with different SpliceAI parameters cannot be accepted silently.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from backend.domain.lazy_mapping import LazyJsonObject


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROFILE_PATH = PROJECT_ROOT / "data" / "spliceai" / "enigma_v1_2_spliceai_profile.json"

_EXPECTED = {
    "profile_id": "enigma-brca-v1.2-appendix-j-spliceai-raw-10kb-v3",
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
    "reference_transcripts": {
        "BRCA1": {"refseq": "NM_007294.4", "ensembl": "ENST00000357654.9"},
        "BRCA2": {"refseq": "NM_000059.4", "ensembl": "ENST00000380152.8"},
    },
}
_EXPECTED_SHA256 = "dd0113cbee03f8847f263d138c98389c67b1dfa782c984ca6220e6d7f8513906"


@lru_cache(maxsize=1)
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
    approved_engine = value.get("approved_engine") or {}
    if approved_engine != {
        "software": "SpliceAI",
        "python_package_version": "1.3.4",
        "model_commit": "7f36ca847e1b1885167dab79681dbb75c09c6743",
        "lookup_server": "Broad Institute SpliceAI Lookup, self-hosted",
        "lookup_annotation": "GENCODE v49 basic, restricted to the pinned reference transcripts",
        "reference_transcript_annotation": (
            "data/spliceai/gencode_v49_basic_reference_transcripts.tsv"
        ),
        "reference_transcript_annotation_sha256": (
            "932290cb981fad001ba03b08aecdb7d961e88101f2d6bd686ab5afff3bc14faf"
        ),
        "docker_image": (
            "docker.io/weisburd/spliceai-38@sha256:"
            "c2bb3d5c65eee01087b5dd130d72233fc812d158dce9779ae63117e826ec8399"
        ),
    }:
        mismatches.append("approved_engine")
    if mismatches:
        raise RuntimeError(
            "ENIGMA SpliceAI scoring profile does not match Appendix J: "
            + ", ".join(sorted(set(mismatches)))
        )
    digest = hashlib.sha256(PROFILE_PATH.read_bytes()).hexdigest()
    if digest != _EXPECTED_SHA256:
        raise RuntimeError(
            "ENIGMA SpliceAI scoring profile checksum does not match the pinned release"
        )
    return value


def validate_spliceai_profile() -> dict[str, Any]:
    """Load and validate the pinned profile during application startup."""
    return _load_profile()


SPLICEAI_PROFILE = LazyJsonObject(_load_profile)
SPLICEAI_PROFILE_SHA256 = _EXPECTED_SHA256
SPLICEAI_PROFILE_ID = str(_EXPECTED["profile_id"])
SPLICEAI_MAX_DISTANCE = int(_EXPECTED["max_distance"])
SPLICEAI_MASK = int(_EXPECTED["mask"])
SPLICEAI_ANNOTATION_SUBSET = str(_EXPECTED["annotation_subset"])
SPLICEAI_GENOME_ASSEMBLY = str(_EXPECTED["genome_assembly"])
SPLICEAI_TRANSCRIPT_POLICY_REQUIRED = str(_EXPECTED["transcript_policy"])
SPLICEAI_AGGREGATION = str(_EXPECTED["aggregation"])
SPLICEAI_DELTA_FIELDS = tuple(_EXPECTED["delta_score_fields"])
SPLICEAI_REFERENCE_FIELDS = tuple(_EXPECTED["reference_score_fields"])
SPLICEAI_ALTERNATE_FIELDS = tuple(_EXPECTED["alternate_score_fields"])
SPLICEAI_REFERENCE_TRANSCRIPTS = dict(_EXPECTED["reference_transcripts"])
SPLICEAI_APPROVED_DOCKER_IMAGE = (
    "docker.io/weisburd/spliceai-38@sha256:"
    "c2bb3d5c65eee01087b5dd130d72233fc812d158dce9779ae63117e826ec8399"
)
SPLICEAI_LOW_THRESHOLD = 0.1
SPLICEAI_HIGH_THRESHOLD = 0.2


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
