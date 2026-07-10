import unittest

from backend.lookups import spliceai


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
            {
                "t_id": "ENST00000634433.2",
                "t_refseq_ids": ["NM_001407648.1"],
                "DS_AG": "0.00",
                "DS_AL": "0.31",
                "DS_DG": "0.00",
                "DS_DL": "0.00",
            },
            {
                "t_id": "ENST00000357654.9",
                "t_refseq_ids": ["NM_007294.4"],
                "DS_AG": "0.00",
                "DS_AL": "0.23",
                "DS_DG": "0.00",
                "DS_DL": "0.00",
            },
        ]

        selected = spliceai._select_spliceai_score("BRCA1", records)

        self.assertEqual(selected["score"], 0.23)
        self.assertEqual(selected["selected_transcript"], "ENST00000357654.9")
        self.assertEqual(selected["max_any_transcript_score"], 0.31)
        self.assertEqual(selected["max_any_transcript"], "ENST00000634433.2")

    def test_max_any_transcript_policy_uses_highest_record(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = "max_any_transcript"
        records = [
            {
                "t_id": "ENST00000634433.2",
                "t_refseq_ids": ["NM_001407648.1"],
                "DS_AG": "0.00",
                "DS_AL": "0.31",
                "DS_DG": "0.00",
                "DS_DL": "0.00",
            },
            {
                "t_id": "ENST00000357654.9",
                "t_refseq_ids": ["NM_007294.4"],
                "DS_AG": "0.00",
                "DS_AL": "0.23",
                "DS_DG": "0.00",
                "DS_DL": "0.00",
            },
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
                "status": "ok",
                "score": 0.01,
                "max_delta_field": "DS_AG",
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

    def test_max_any_policy_skips_reference_precomputed_cache(self):
        spliceai.SPLICEAI_TRANSCRIPT_POLICY = "max_any_transcript"
        spliceai.SPLICEAI_USE_PRECOMPUTED_CACHE = True
        spliceai.SPLICEAI_PRECOMPUTED_CACHE = {
            "BRCA1:c.6T>A": {"status": "ok", "score": 0.01}
        }
        spliceai.SPLICEAI_CACHE.clear()
        spliceai.SPLICEAI_STATUS_CACHE.clear()

        self.assertIsNone(spliceai._lookup_precomputed_score("BRCA1", "c.6T>A"))


if __name__ == "__main__":
    unittest.main()
