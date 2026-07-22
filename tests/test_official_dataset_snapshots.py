import json
import unittest
from pathlib import Path


DATA = Path(__file__).resolve().parents[1] / "backend" / "data"


class OfficialDatasetSnapshotTests(unittest.TestCase):
    def test_table9_is_lossless_and_preserves_splicing_evidence(self):
        data = json.loads((DATA / "enigma_table9.json").read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], 2)
        self.assertEqual(data["row_count"], 4731)
        self.assertEqual(len(data["variants"]), 4731)
        entry = data["variants"]["BRCA1:c.3891_3893del"]
        self.assertEqual(len(entry), 14)
        self.assertEqual(entry["splice_result_published"], "no aberration (PMID: 18273839)")
        self.assertEqual(entry["spliceai_prediction"], 0)
        self.assertEqual(entry["predicted_or_observed_splicing"], "N, no aberration")
        self.assertEqual(
            sum(entry["code"] == "None" for entry in data["variants"].values()),
            437,
        )

    def test_table4_is_lossless_and_preserves_warnings(self):
        data = json.loads((DATA / "enigma_table4.json").read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], 2)
        self.assertEqual(data["source_columns"], 20)
        self.assertEqual(len(data["source_rows"]), 493)
        self.assertIn(
            "GC splice site",
            data["deletion_rules"]["BRCA2"]["E15"]["notes"],
        )
        self.assertEqual(len(data["splice_rules"]["BRCA1"]), 264)
        self.assertEqual(len(data["splice_rules"]["BRCA2"]), 311)

    def test_st7_is_lossless(self):
        data = json.loads((DATA / "st7_reference_set.json").read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], 2)
        self.assertEqual(data["source_columns"], 28)
        self.assertEqual(data["total_variants"], 773)
        self.assertTrue(all(len(record) == 28 for record in data["variants"]))

    def test_informational_reference_sets_have_official_provenance(self):
        residues = json.loads(
            (DATA / "clinically_important_residues.json").read_text(encoding="utf-8")
        )
        brca1 = sum(
            len(domain.get("pathogenic_residues", []))
            for domain in residues["domains"]["BRCA1"].values()
        )
        brca2 = sum(
            len(domain.get("pathogenic_residues", []))
            for domain in residues["domains"]["BRCA2"].values()
        )
        self.assertEqual((brca1, brca2), (36, 8))

        splice_ps1 = json.loads(
            (DATA / "splice_ps1_reference_set.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            splice_ps1["curation_status"],
            "pilot_unreviewed_not_for_automatic_scoring",
        )
        self.assertEqual(len(splice_ps1["variants"]), 75)
        self.assertEqual(
            len({record["source_row"] for record in splice_ps1["variants"]}), 75
        )
        self.assertNotIn("GN097", json.dumps(splice_ps1))


if __name__ == "__main__":
    unittest.main()
