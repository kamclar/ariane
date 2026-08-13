import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from backend.data_health import clear_issue, get_data_issues
from backend.lookups import spliceai
from backend.spliceai_profile import (
    SPLICEAI_PROFILE_ID,
    SPLICEAI_PROFILE_SHA256,
    scoring_profile_metadata,
)


def complete_entry():
    return {
        "status": "ok",
        "score": 0.2,
        "max_delta_field": "DS_AG",
        "delta_scores": {
            "DS_AG": 0.2, "DS_AL": 0.1, "DS_DG": 0.0, "DS_DL": 0.0,
        },
        "reference_scores": {
            "DS_AG_REF": 0.3, "DS_AL_REF": 0.2,
            "DS_DG_REF": 0.1, "DS_DL_REF": 0.0,
        },
        "alternate_scores": {
            "DS_AG_ALT": 0.1, "DS_AL_ALT": 0.1,
            "DS_DG_ALT": 0.1, "DS_DL_ALT": 0.0,
        },
    }


def write_cache(path: Path, entries: dict, **metadata_overrides):
    path.write_text(json.dumps(entries), encoding="utf-8")
    metadata = {
        **scoring_profile_metadata(),
        "expected_variants": len(entries),
        "cache_entries": len(entries),
        "status_ok": len(entries),
        "status_error": 0,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        **metadata_overrides,
    }
    path.with_name(path.stem + ".metadata.json").write_text(
        json.dumps(metadata), encoding="utf-8"
    )


class SpliceAIAppendixJProfileTests(unittest.TestCase):
    def tearDown(self):
        clear_issue("SpliceAI coding SNV cache")
        clear_issue("SpliceAI intronic cache")

    def test_profile_is_fixed_to_appendix_j(self):
        self.assertEqual(SPLICEAI_PROFILE_ID, "enigma-brca-v1.2-appendix-j-spliceai-raw-10kb-v1")
        self.assertEqual(spliceai.SPLICEAI_MAX_DISTANCE, 10000)
        self.assertEqual(spliceai.SPLICEAI_MASK, 0)
        self.assertEqual(spliceai.SPLICEAI_LOW_THRESHOLD, 0.1)
        self.assertEqual(spliceai.SPLICEAI_HIGH_THRESHOLD, 0.2)

    def test_complete_cache_with_exact_profile_is_loaded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coding = root / "coding.json"
            intronic = root / "intronic.json"
            write_cache(coding, {"BRCA1:c.1A>G": complete_entry()})
            write_cache(intronic, {})
            old = spliceai.SPLICEAI_PRECOMPUTED_CACHE
            old_enabled = spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE
            try:
                spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = True
                spliceai.SPLICEAI_PRECOMPUTED_CACHE = None
                with patch.object(spliceai, "SPLICEAI_PRECOMPUTED_CACHE_PATH", coding), patch.object(
                    spliceai, "SPLICEAI_INTRONIC_CACHE_PATH", intronic
                ):
                    loaded = spliceai._load_precomputed_cache()
                self.assertIn("BRCA1:c.1A>G", loaded)
            finally:
                spliceai.SPLICEAI_PRECOMPUTED_CACHE = old
                spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = old_enabled

    def test_distance_50_cache_is_rejected_without_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coding = root / "coding.json"
            intronic = root / "intronic.json"
            write_cache(coding, {"BRCA1:c.1A>G": complete_entry()}, distance=50)
            write_cache(intronic, {})
            old = spliceai.SPLICEAI_PRECOMPUTED_CACHE
            old_enabled = spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE
            try:
                spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = True
                spliceai.SPLICEAI_PRECOMPUTED_CACHE = None
                with patch.object(spliceai, "SPLICEAI_PRECOMPUTED_CACHE_PATH", coding), patch.object(
                    spliceai, "SPLICEAI_INTRONIC_CACHE_PATH", intronic
                ):
                    loaded = spliceai._load_precomputed_cache()
                self.assertNotIn("BRCA1:c.1A>G", loaded)
                issue = next(
                    item for item in get_data_issues()
                    if item["component"] == "SpliceAI coding SNV cache"
                )
                self.assertIn("distance=50", issue["reason"])
                self.assertIn("expected 10000", issue["reason"])
            finally:
                spliceai.SPLICEAI_PRECOMPUTED_CACHE = old
                spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = old_enabled

    def test_precomputed_caches_are_disabled_in_current_api_primary_mode(self):
        self.assertFalse(spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE)

    def test_runtime_record_requires_profile_and_ref_alt_audit(self):
        entry = complete_entry()
        entry.update({
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
            "genome_assembly": "GRCh38",
            "distance": 10000,
            "mask": 0,
            "annotation_subset": "basic",
            "aggregation": "maximum_raw_delta",
            "transcript_policy": "reference_transcript",
        })
        self.assertTrue(spliceai._runtime_entry_matches_profile(entry))
        entry.pop("alternate_scores")
        self.assertFalse(spliceai._runtime_entry_matches_profile(entry))


if __name__ == "__main__":
    unittest.main()
