import copy
import unittest

from backend.modules.exon_cnv_evidence import (
    load_exon_cnv_evidence_snapshot,
    lookup_exon_cnv_evidence,
    validate_exon_cnv_evidence_snapshot,
)
from scripts.build_exon_cnv_evidence_snapshot import spans_all_exons


class ExonCnvEvidenceTests(unittest.TestCase):
    def test_population_match_requires_deletion_to_span_the_complete_exon(self):
        intervals = [{
            "chrom": "13",
            "start_1_based_inclusive": 32906409,
            "end_1_based_inclusive": 32907524,
        }]
        complete = {
            "#chrom": "13", "start": "32906408", "end": "32907524",
            "svtype": "DEL",
        }
        misses_first_base = {**complete, "start": "32906409"}
        wrong_type = {**complete, "svtype": "DUP"}

        self.assertTrue(spans_all_exons(complete, intervals))
        self.assertFalse(spans_all_exons(misses_first_base, intervals))
        self.assertFalse(spans_all_exons(wrong_type, intervals))

    def test_snapshot_is_exon_wide_and_contains_no_variant_allowlist(self):
        payload = load_exon_cnv_evidence_snapshot()

        self.assertEqual(len(payload["exons"]), 50)
        self.assertIn("BRCA1:E10(11)", payload["exons"])
        self.assertIn("BRCA2:E10", payload["exons"])
        self.assertNotIn("records", payload)
        self.assertNotIn("c_notation", payload["exons"]["BRCA2:E10"])
        self.assertNotIn("criteria", payload["exons"]["BRCA2:E10"])

    def test_brca2_exon10_deletion_runs_the_appendix_g_graph(self):
        result = lookup_exon_cnv_evidence(
            "BRCA2", "c.(793+1_794-1)_(1909+1_1910-1)del"
        )

        self.assertTrue(result["found"])
        self.assertEqual(result["exon"], "E10")
        criteria = {item["code"]: item for item in result["criteria"]}
        self.assertEqual(criteria["PM2_Supporting"]["points"], 1)
        self.assertEqual(
            [step["step"] for step in result["decision_trace"]],
            [
                "variant_type",
                "table4_exon_mapping",
                "grch37_exon_interval",
                "appendix_g_size",
                "gnomad_sv_exon_match",
            ],
        )
        self.assertTrue(all(step["status"] == "pass" for step in result["decision_trace"]))
        self.assertEqual(
            result["source_identity"]["sha256"],
            "c843ff53b4bf36c7f733cb08565860065b3b0189375d135e33db0886381598d8",
        )

    def test_non_exon_exact_deletion_fails_closed(self):
        result = lookup_exon_cnv_evidence("BRCA2", "c.794_1908del")

        self.assertFalse(result["found"])
        self.assertEqual(result["criteria"], [])
        self.assertIn("unambiguously", result["reason"])

    def test_tampered_snapshot_fails_closed(self):
        payload = copy.deepcopy(load_exon_cnv_evidence_snapshot())
        payload["exons"]["BRCA2:E10"]["coding_length_bp"] = 1

        with self.assertRaisesRegex(RuntimeError, "exon checksum mismatch"):
            validate_exon_cnv_evidence_snapshot(payload)


if __name__ == "__main__":
    unittest.main()
