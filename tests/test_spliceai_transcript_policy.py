import unittest
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import tempfile
import urllib.error
from unittest.mock import MagicMock, patch

from backend.lookups import spliceai
from backend.policy.spliceai import spliceai_failure_is_retryable


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


class SpliceAITranscriptPolicyTests(unittest.TestCase):
    def test_concurrent_cache_entry_updates_are_merged(self):
        with tempfile.TemporaryDirectory() as directory:
            runtime_dir = Path(directory)
            cache_path = runtime_dir / "spliceai_api_cache.json"
            keys = [
                spliceai._cache_key("BRCA1", "c.1A>G"),
                spliceai._cache_key("BRCA1", "c.2A>G"),
            ]
            with patch.object(spliceai, "RUNTIME_CACHE_DIR", runtime_dir), patch.object(
                spliceai, "SPLICEAI_API_CACHE_PATH", cache_path
            ):
                with ThreadPoolExecutor(max_workers=2) as executor:
                    saved = list(executor.map(
                        lambda item: spliceai._persist_api_cache_entry(
                            item[0], {"score": item[1]}
                        ),
                        zip(keys, (0.1, 0.2)),
                    ))

                self.assertEqual(saved, [True, True])
                persisted = json.loads(cache_path.read_text(encoding="utf-8"))
                self.assertEqual(set(persisted), set(keys))

    def test_runtime_cache_discards_retired_profile_entries(self):
        current_key = (
            f"{spliceai.SPLICEAI_PROFILE_ID}:reference_transcript:BRCA1:c.1A>G"
        )
        cache = spliceai._current_profile_cache_entries({
            current_key: {"score": 0.2},
            "retired-profile:reference_transcript:BRCA1:c.2A>G": {"score": 0.3},
        })
        self.assertEqual(cache, {current_key: {"score": 0.2}})

    def setUp(self):
        self._old_policy = spliceai.SPLICEAI_TRANSCRIPT_POLICY
        self._old_score_cache = dict(spliceai.SPLICEAI_CACHE)
        self._old_status_cache = dict(spliceai.SPLICEAI_STATUS_CACHE)

    def tearDown(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = self._old_policy
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

    def test_reference_transcript_policy_rejects_a_different_transcript_version(self):
        records = [
            score_row("ENST00000357654.10", "NM_007294.4", ds_al=0.23),
        ]

        selected = spliceai._select_spliceai_score("BRCA1", records)

        self.assertIsNone(selected["score"])
        self.assertIn("exact reference transcript", selected["error"])
        self.assertIn("ENST00000357654.9", selected["error"])

    def test_runtime_cache_requires_the_exact_profile_transcript(self):
        row = score_row("ENST00000357654.9", "NM_007294.4", ds_al=0.23)
        audit = spliceai._row_score_audit(row)
        entry = {
            **audit,
            "scoring_profile_id": spliceai.SPLICEAI_PROFILE_ID,
            "scoring_profile_sha256": spliceai.SPLICEAI_PROFILE_SHA256,
            "genome_assembly": spliceai.SPLICEAI_GENOME_ASSEMBLY,
            "distance": spliceai.SPLICEAI_MAX_DISTANCE,
            "mask": spliceai.SPLICEAI_MASK,
            "annotation_subset": spliceai.SPLICEAI_ANNOTATION_SUBSET,
            "aggregation": spliceai.SPLICEAI_AGGREGATION,
            "transcript_policy": "reference_transcript",
            "selected_transcript": "ENST00000357654.9",
        }

        self.assertTrue(spliceai._runtime_entry_matches_profile(entry, "BRCA1"))
        entry["selected_transcript"] = "ENST00000357654.10"
        self.assertFalse(spliceai._runtime_entry_matches_profile(entry, "BRCA1"))

    def test_brca2_7805_plus_9_current_api_components_select_0206(self):
        records = [
            score_row(
                "ENST00000380152.8",
                "NM_000059.4",
                ds_al=0.016,
                ds_dl=0.049,
            ),
            score_row(
                "ENST00000700202.2",
                "NM_001406720.1",
                ds_al=0.90,
            ),
        ]
        records[0]["DS_AG"] = "0.030"
        records[0]["DS_DG"] = "0.206"

        selected = spliceai._select_spliceai_score("BRCA2", records)

        self.assertEqual(selected["selected_transcript"], "ENST00000380152.8")
        self.assertEqual(selected["reference_transcript"], "ENST00000380152.8")
        self.assertEqual(selected["max_delta_field"], "DS_DG")
        self.assertEqual(selected["score"], 0.206)
        self.assertEqual(selected["max_any_transcript_score"], 0.9)

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

    def test_api_request_preserves_multi_nucleotide_ref_and_alt_alleles(self):
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
            result = spliceai._query_spliceai_api("BRCA1", "17", 100, "TG", "GC")

        requested_url = mocked.call_args.args[0].full_url
        self.assertIn("variant=chr17-100-TG-GC", requested_url)
        self.assertEqual(result["score"], 0.23)

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

    def test_transient_timeout_is_retried_once(self):
        payload = {
            "hg": "38",
            "genomeVersion": "38",
            "distance": 10000,
            "mask": 0,
            "source": "test SpliceAI",
            "scores": [score_row(
                "ENST00000357654.9",
                "NM_007294.4",
                ds_al=0.23,
            )],
        }
        response = MagicMock()
        response.__enter__.return_value.read.return_value = json.dumps(payload).encode()
        with patch.object(
            spliceai.urllib.request,
            "urlopen",
            side_effect=[TimeoutError("temporary timeout"), response],
        ) as mocked, patch.object(spliceai.time, "sleep"), patch.object(
            spliceai, "SPLICEAI_API_ATTEMPTS", 2
        ):
            result = spliceai._query_spliceai_api("BRCA1", "17", 1, "A", "G")

        self.assertEqual(result["score"], 0.23)
        self.assertEqual(mocked.call_count, 2)

    def test_http_400_is_not_retried(self):
        error = urllib.error.HTTPError(
            "https://example.test",
            400,
            "Bad Request",
            None,
            None,
        )
        with patch.object(
            spliceai.urllib.request,
            "urlopen",
            side_effect=error,
        ) as mocked, patch.object(spliceai.time, "sleep"):
            result = spliceai._query_spliceai_api("BRCA1", "17", 1, "A", "G")

        self.assertIsNone(result["score"])
        self.assertFalse(result["retryable"])
        self.assertEqual(mocked.call_count, 1)

    def test_only_transient_api_failures_are_retryable(self):
        self.assertTrue(spliceai_failure_is_retryable(
            "api_error", "TimeoutError: The read operation timed out"
        ))
        self.assertTrue(spliceai_failure_is_retryable(
            "api_error", "HTTP Error 503: Service Unavailable"
        ))
        self.assertFalse(spliceai_failure_is_retryable(
            "api_error", "HTTP Error 404: Not Found"
        ))
        self.assertFalse(spliceai_failure_is_retryable(
            "api_error", "RuntimeError: unexpected response"
        ))


if __name__ == "__main__":
    unittest.main()
