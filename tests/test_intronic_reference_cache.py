import hashlib
import json
import re
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
COORDINATES = ROOT / "data/coordinates/brca_intronic_snv_coordinates.json"
COORDINATE_METADATA = ROOT / "data/coordinates/brca_intronic_snv_coordinates.metadata.json"
SPLICEAI = ROOT / "data/spliceai/spliceai_brca_intronic_snv_reference_cache.json"
SPLICEAI_METADATA = ROOT / "data/spliceai/spliceai_brca_intronic_snv_reference_cache.metadata.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class IntronicCoordinateDatasetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.coordinates = json.loads(COORDINATES.read_text(encoding="utf-8"))
        cls.metadata = json.loads(COORDINATE_METADATA.read_text(encoding="utf-8"))

    def test_dataset_is_complete_and_checksum_matches(self):
        self.assertEqual(self.metadata["window_bp_per_exon_boundary"], 50)
        self.assertEqual(self.metadata["variants"], 13800)
        self.assertEqual(len(self.coordinates), 13800)
        self.assertEqual(self.metadata["sha256"], sha256(COORDINATES))

    def test_every_site_has_three_alternatives_and_both_builds(self):
        sites = {}
        for key, entry in self.coordinates.items():
            match = re.fullmatch(r"(BRCA[12]):(c\.\d+[+-]\d+)([ACGT])>([ACGT])", key)
            self.assertIsNotNone(match, key)
            site = match.group(1), match.group(2), match.group(3)
            sites.setdefault(site, set()).add(match.group(4))
            self.assertEqual(entry["status"], "ok")
            self.assertTrue(entry["grch37"])
            self.assertTrue(entry["grch38"])
        self.assertEqual(len(sites), 4600)
        self.assertTrue(all(len(alternatives) == 3 for alternatives in sites.values()))

    def test_runtime_uses_generated_map_without_network(self):
        from backend.lookups import coordinates

        key = next(key for key in self.coordinates if key.startswith("BRCA2:c.67+50"))
        c_notation = key.split(":", 1)[1]
        coordinates._RESOLVER_CACHE.pop(key, None)
        coordinates._load_coords_cache()
        with patch.object(coordinates, "_resolve_variantvalidator", side_effect=AssertionError("network called")):
            result = coordinates.resolve_variant("BRCA2", c_notation)
        self.assertEqual(result.status, "ok")
        self.assertEqual(result.source, "versioned_intronic_coordinate_map")


class IntronicSpliceAIDatasetTests(unittest.TestCase):
    @unittest.skipUnless(SPLICEAI.exists() and SPLICEAI_METADATA.exists(), "SpliceAI build is not complete")
    def test_spliceai_dataset_is_complete_and_checksum_matches(self):
        cache = json.loads(SPLICEAI.read_text(encoding="utf-8"))
        metadata = json.loads(SPLICEAI_METADATA.read_text(encoding="utf-8"))
        self.assertEqual(metadata["coordinate_variants"], 13800)
        self.assertEqual(metadata["status_ok"], 13800)
        self.assertEqual(metadata["status_error"], 0)
        self.assertEqual(len(cache), 13800)
        self.assertTrue(all(entry.get("status") == "ok" for entry in cache.values()))
        self.assertEqual(metadata["sha256"], sha256(SPLICEAI))


if __name__ == "__main__":
    unittest.main()
