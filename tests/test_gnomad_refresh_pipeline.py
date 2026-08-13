import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

from scripts import refresh_gnomad_panel_snapshot as refresh


class GnomadRefreshPipelineTests(unittest.TestCase):
    def test_manifest_rejects_interval_only_gene_activation(self):
        manifest = json.loads(
            refresh.DEFAULT_MANIFEST.read_text(encoding="utf-8")
        )
        target = deepcopy(manifest["targets"][0])
        target["gene"] = "GENE3"
        target.pop("classification_policy_id")
        manifest["targets"].append(target)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(
                RuntimeError, "GENE3 has no active gene-specific policy"
            ):
                refresh._load_manifest(path)

    def test_manifest_rejects_new_gene_pointing_to_brca_policy(self):
        manifest = json.loads(
            refresh.DEFAULT_MANIFEST.read_text(encoding="utf-8")
        )
        target = deepcopy(manifest["targets"][0])
        target["gene"] = "GENE3"
        manifest["targets"].append(target)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(
                RuntimeError, "not approved for target GENE3"
            ):
                refresh._load_manifest(path)

    def test_manifest_pins_enigma_brca_frequency_policy(self):
        manifest = refresh._load_manifest(refresh.DEFAULT_MANIFEST)
        policy = manifest["classification_policies"]["enigma_brca_v1_2"]
        frequency = policy["frequency_criteria"]
        self.assertEqual(policy["policy_scope"], "gnomad_frequency_criteria_only")
        self.assertEqual(set(policy["applicable_genes"]), {"BRCA1", "BRCA2"})
        self.assertFalse(
            policy["automation_scope"]["inherit_other_vcep_criteria"]
        )
        self.assertEqual(frequency["ba1"]["threshold"], 0.001)
        self.assertEqual(frequency["bs1"]["strong"]["threshold"], 0.0001)
        self.assertEqual(
            frequency["bs1"]["supporting"]["lower_threshold"], 0.00002
        )
        self.assertEqual(frequency["pm2"]["minimum_mean_depth"], 25.0)
        self.assertTrue(all(
            target["classification_policy_id"] == "enigma_brca_v1_2"
            for target in manifest["targets"]
        ))

    def test_overlapping_future_gene_targets_do_not_make_coverage_ambiguous(self):
        manifest = {
            "region_padding_bp": 0,
            "targets": [
                {
                    "gene": "GENE1",
                    "intervals": {
                        "GRCh38": {"chrom": "1", "start": 100, "end": 200}
                    },
                },
                {
                    "gene": "GENE2",
                    "intervals": {
                        "GRCh38": {"chrom": "1", "start": 150, "end": 250}
                    },
                },
            ],
        }
        self.assertEqual(
            refresh._gene_for_position(manifest, "GRCh38", "chr1", 175),
            "GENE1,GENE2",
        )

    def test_candidate_url_checks_the_equivalent_sites_product(self):
        dataset = {
            "dataset_key": "gnomad_v3_1_2_genomes_grch38",
            "release": "3.1.2",
            "sites_metadata_url": (
                "https://storage.googleapis.com/bucket/release/3.1.2/ht/genomes/"
                "gnomad.genomes.v3.1.2.sites.ht/metadata.json.gz"
            ),
        }
        self.assertEqual(
            refresh._candidate_sites_metadata_url(dataset, "3.1.3"),
            "https://storage.googleapis.com/bucket/release/3.1.3/ht/genomes/"
            "gnomad.genomes.v3.1.3.sites.ht/metadata.json.gz",
        )

    def test_update_checker_does_not_treat_other_release_product_as_update(self):
        manifest = {
            "datasets": [
                {
                    "dataset_key": "gnomad_v3_1_2_genomes_grch38",
                    "release": "3.1.2",
                    "sites_metadata_url": (
                        "https://storage.googleapis.com/bucket/release/3.1.2/ht/"
                        "genomes/gnomad.genomes.v3.1.2.sites.ht/metadata.json.gz"
                    ),
                }
            ]
        }
        with patch.object(refresh, "_load_manifest", return_value=manifest), patch.object(
            refresh, "_verify_sources", return_value={"verified": True}
        ), patch.object(
            refresh, "_discover_releases", return_value=["3.1.2", "3.1.3"]
        ), patch.object(
            refresh, "_source_exists", return_value=(False, 404)
        ):
            report = refresh.check_updates(Path("unused.json"))

        self.assertEqual(
            report["newer_same_major_release_directories"],
            {"3.1.2": ["3.1.3"]},
        )
        self.assertEqual(
            report["newer_equivalent_product_releases"], {"3.1.2": []}
        )
        self.assertFalse(report["automatic_activation"])
        self.assertEqual(report["activation_status"], "up_to_date")
        self.assertEqual(
            report["equivalent_product_checks"][0]["http_status"], 404
        )


if __name__ == "__main__":
    unittest.main()
