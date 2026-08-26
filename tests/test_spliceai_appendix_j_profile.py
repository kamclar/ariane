import unittest

from backend.lookups import spliceai
from backend.spliceai_profile import (
    SPLICEAI_PROFILE_ID,
    SPLICEAI_PROFILE_SHA256,
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


class SpliceAIAppendixJProfileTests(unittest.TestCase):
    def test_profile_is_fixed_to_appendix_j(self):
        self.assertEqual(SPLICEAI_PROFILE_ID, "enigma-brca-v1.2-appendix-j-spliceai-raw-10kb-v1")
        self.assertEqual(spliceai.SPLICEAI_MAX_DISTANCE, 10000)
        self.assertEqual(spliceai.SPLICEAI_MASK, 0)
        self.assertEqual(spliceai.SPLICEAI_LOW_THRESHOLD, 0.1)
        self.assertEqual(spliceai.SPLICEAI_HIGH_THRESHOLD, 0.2)

    def test_precomputed_variant_space_is_not_supported_by_runtime(self):
        self.assertFalse(hasattr(spliceai, "SPLICEAI_USE_PRECOMPUTED_CACHE"))
        self.assertFalse(hasattr(spliceai, "_load_precomputed_cache"))

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
