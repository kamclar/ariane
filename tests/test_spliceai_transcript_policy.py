import unittest
import json
from unittest.mock import MagicMock, patch

from backend.lookups import spliceai


def score_row(transcript, refseq, *, ds_al, ds_dl=0.0):
    return {
        "t_id": transcript,
        "t_refseq_ids": [refseq],
        "DS_AG": "0.00",
        "DS_AL": str(ds_al),
        "DS_DG": "0.00",
        "DS_DL": str(ds_dl),
        "DS_AG_REF": "0.10",
        "DS_AL_REF": "0.20",
        "DS_DG_REF": "0.30",
        "DS_DL_REF": "0.40",
        "DS_AG_ALT": "0.10",
        "DS_AL_ALT": "0.20",
        "DS_DG_ALT": "0.30",
        "DS_DL_ALT": "0.40",
    }


def cache_entry(score=0.01, max_field="DS_AG"):
    delta = {"DS_AG": 0.01, "DS_AL": 0.0, "DS_DG": 0.0, "DS_DL": 0.0}
    return {
        "status": "ok",
        "score": score,
        "max_delta_field": max_field,
        "delta_scores": delta,
        "reference_scores": {
            "DS_AG_REF": 0.1, "DS_AL_REF": 0.2,
            "DS_DG_REF": 0.3, "DS_DL_REF": 0.4,
        },
        "alternate_scores": {
            "DS_AG_ALT": 0.11, "DS_AL_ALT": 0.2,
            "DS_DG_ALT": 0.3, "DS_DL_ALT": 0.4,
        },
    }


class SpliceAITranscriptPolicyTests(unittest.TestCase):
    def setUp(self):
        self._old_policy = spliceai.SPLICEAI_TRANSCRIPT_POLICY
        self._old_use_precomputed = spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE
        self._old_precomputed_cache = spliceai.SPLICEAI_PRECOMPUTED_CACHE
        self._old_score_cache = dict(spliceai.SPLICEAI_CACHE)
        self._old_status_cache = dict(spliceai.SPLICEAI_STATUS_CACHE)

    def tearDown(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = self._old_policy
        spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = self._old_use_precomputed
        spliceai.SPLICEAI_PRECOMPUTED_CACHE = self._old_precomputed_cache
        spliceai.SPLICEAI_CACHE.clear()
        spliceai.SPLICEAI_CACHE.update(self._old_score_cache)
        spliceai.SPLICEAI_STATUS_CACHE.clear()
        spliceai.SPLICEAI_STATUS_CACHE.update(self._old_status_cache)

    def test_reference_transcript_policy_uses_brca1_reference_record(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = "reference_transcript"
        records = [
            score_row("ENST00000634433.2", "NM_001407648.1", ds_al=0.31),
            score_row("ENST00000357654.9", "NM_007294.4", ds_al=0.23),
        ]

        selected = spliceai._select_spliceai_score("BRCA1", records)

        self.assertEqual(selected["score"], 0.23)
        self.assertEqual(selected["selected_transcript"], "ENST00000357654.9")
        self.assertEqual(selected["max_any_transcript_score"], 0.31)
        self.assertEqual(selected["max_any_transcript"], "ENST00000634433.2")

    def test_max_any_transcript_policy_uses_highest_record(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = "max_any_transcript"
        records = [
            score_row("ENST00000634433.2", "NM_001407648.1", ds_al=0.31),
            score_row("ENST00000357654.9", "NM_007294.4", ds_al=0.23),
        ]

        selected = spliceai._select_spliceai_score("BRCA1", records)

        self.assertEqual(selected["score"], 0.31)
        self.assertEqual(selected["selected_transcript"], "ENST00000634433.2")
        self.assertEqual(selected["reference_transcript_score"], 0.23)
        self.assertEqual(selected["reference_transcript"], "ENST00000357654.9")

    def test_reference_policy_uses_precomputed_cache_before_api(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = "reference_transcript"
        spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = True
        spliceai.SPLICEAI_PRECOMPUTED_CACHE = {
            "BRCA1:c.6T>A": {
                **cache_entry(),
                "source": "test precomputed cache",
                "grch38": "17:43124091:A>T",
            }
        }
        spliceai.SPLICEAI_CACHE.clear()
        spliceai.SPLICEAI_STATUS_CACHE.clear()

        score = spliceai.get_spliceai_score("BRCA1", "c.6T>A")

        self.assertEqual(score, 0.01)
        status = spliceai.SPLICEAI_STATUS_CACHE["BRCA1:c.6T>A"]
        self.assertEqual(status["status"], "ok")
        self.assertEqual(status["source"], "test precomputed cache")
        self.assertIn("precomputed", status["reason"])
        self.assertEqual(status["transcript_policy"], "reference_transcript")
        self.assertEqual(status["selected_transcript"], "ENST00000357654.9")
        self.assertEqual(status["max_delta_field"], "DS_AG")
        self.assertEqual(status["grch38"], "17:43124091:A>T")
        self.assertEqual(status["cache_key"], "BRCA1:c.6T>A")

    def test_max_any_policy_skips_reference_precomputed_cache(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = "max_any_transcript"
        spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = True
        spliceai.SPLICEAI_PRECOMPUTED_CACHE = {
            "BRCA1:c.6T>A": cache_entry()
        }
        spliceai.SPLICEAI_CACHE.clear()
        spliceai.SPLICEAI_STATUS_CACHE.clear()

        self.assertIsNone(spliceai._lookup_precomputed_score("BRCA1", "c.6T>A"))

    def test_runtime_cache_key_contains_immutable_scoring_profile(self):
        key = spliceai._cache_key("BRCA1", "c.6T>A")
        self.assertIn(spliceai.SPLICEAI_PROFILE_ID, key)
        self.assertNotEqual(key, "reference_transcript:BRCA1:c.6T>A")

    def test_api_request_uses_appendix_j_parameters_and_keeps_ref_alt_scores(self):
        payload = {
            "hg": "38",
            "genomeVersion": "38",
            "distance": 10000,
            "mask": 0,
            "source": "test SpliceAI",
            "scores": [score_row("ENST00000357654.9", "NM_007294.4", ds_al=0.23)],
        }
        response = MagicMock()
        response.__enter__.return_value.read.return_value = json.dumps(payload).encode()
        with patch.object(spliceai.urllib.request, "urlopen", return_value=response) as mocked:
            result = spliceai._query_spliceai_api("BRCA1", "17", 1, "A", "G")
        requested_url = mocked.call_args.args[0].full_url
        self.assertIn("distance=10000", requested_url)
        self.assertIn("mask=0", requested_url)
        self.assertIn("bc=basic", requested_url)
        self.assertEqual(result["score"], 0.23)
        self.assertEqual(result["reference_scores"]["DS_AL_REF"], 0.2)
        self.assertEqual(result["alternate_scores"]["DS_AL_ALT"], 0.2)

    def test_api_response_with_wrong_distance_is_rejected(self):
        payload = {
            "hg": "38", "distance": 50, "mask": 0,
            "scores": [score_row("ENST00000357654.9", "NM_007294.4", ds_al=0.23)],
        }
        response = MagicMock()
        response.__enter__.return_value.read.return_value = json.dumps(payload).encode()
        with patch.object(spliceai.urllib.request, "urlopen", return_value=response):
            result = spliceai._query_spliceai_api("BRCA1", "17", 1, "A", "G")
        self.assertIsNone(result["score"])
        self.assertIn("profile mismatch", result["error"])


if __name__ == "__main__":
    unittest.main()
