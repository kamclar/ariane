import unittest

from backend.lookups.precomputed import (
    load_classification_snapshot_metadata,
    lookup_classification_snapshot,
)


class PrecomputedSnapshotTests(unittest.TestCase):
    def test_metadata_is_available(self):
        metadata = load_classification_snapshot_metadata()
        self.assertEqual(metadata["n_records"], 47547)
        self.assertEqual(metadata["status"], "snapshot_not_authoritative")
        self.assertIn("index_sha256", metadata)

    def test_lookup_existing_variant(self):
        result = lookup_classification_snapshot("brca1", "c.1A>G")
        self.assertIsNotNone(result)
        self.assertEqual(result["snapshot_status"], "snapshot_not_authoritative")
        self.assertEqual(result["record"]["predicted_class"], 5)
        self.assertEqual(result["record"]["predicted_label"], "Pathogenic")

    def test_lookup_missing_variant(self):
        self.assertIsNone(lookup_classification_snapshot("BRCA1", "c.999999A>G"))


if __name__ == "__main__":
    unittest.main()
