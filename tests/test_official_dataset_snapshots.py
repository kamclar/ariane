import json
import hashlib
import unittest
from collections import Counter
from pathlib import Path


DATA = Path(__file__).resolve().parents[1] / "backend" / "data"


class OfficialDatasetSnapshotTests(unittest.TestCase):
    def test_pathogenic_founder_snapshot_has_integrity_and_provenance(self):
        data = json.loads(
            (DATA / "brca_pathogenic_founder_variants.json").read_text(encoding="utf-8")
        )
        canonical = json.dumps(
            data["variants"],
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        self.assertEqual(
            hashlib.sha256(canonical).hexdigest(),
            data["metadata"]["records_sha256"],
        )
        self.assertEqual(len(data["variants"]), 8)
        self.assertTrue(all(
            source.get("content_sha256_at_access")
            for source in data["metadata"]["sources"]
        ))
        self.assertIn(
            ("BRCA1", "c.181T>G"),
            {
                (record["gene"], record["canonical_c_notation"])
                for record in data["variants"]
            },
        )

    def test_hail_filtering_allele_frequency_reference_values(self):
        from scripts.build_gnomad_v3_brca_snapshot import filtering_allele_frequency

        self.assertEqual(filtering_allele_frequency(1, 133_976), 0.0)
        self.assertEqual(filtering_allele_frequency(2, 133_976), 0.00000248)
        self.assertEqual(filtering_allele_frequency(3, 100_000), 0.00000797)

    def test_gnomad_cache_has_complete_non_cancer_faf95(self):
        path = DATA / "gnomad" / "gnomad_brca_frequency_snapshot.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        manifest = json.loads(
            (DATA / "gnomad" / "gnomad_panel_manifest.json").read_text(encoding="utf-8")
        )
        records = [
            record
            for entries in data["variants"].values()
            for record in entries
        ]
        self.assertEqual(
            Counter(record["dataset"] for record in records),
            Counter({
                "gnomad_v2_1_1_exomes_grch37": 5602,
                "gnomad_v3_1_2_genomes_grch38": 37148,
            }),
        )
        self.assertTrue(all(record["faf95_max"] is not None for record in records))
        self.assertTrue(all(
            record["faf95_scope"] == "non_cancer_non_founder_ancestries"
            for record in records
        ))
        self.assertTrue(all(record["faf95_method"] for record in records))
        self.assertEqual(
            data["metadata"]["manifest_sha256"],
            hashlib.sha256(
                json.dumps(
                    manifest,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
            ).hexdigest(),
        )
        self.assertEqual(
            data["metadata"]["records_sha256"],
            hashlib.sha256(
                json.dumps(
                    data["variants"],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
            ).hexdigest(),
        )
        self.assertFalse(data["metadata"]["automatic_release_activation"])
        self.assertEqual(
            data["metadata"]["classification_policies"],
            manifest["classification_policies"],
        )
        self.assertFalse(data["metadata"]["v3_faf95"]["api_faf95_used"])
        self.assertEqual(
            {
                record["dataset"]: record["faf95_method"]
                for record in records
            },
            {
                "gnomad_v2_1_1_exomes_grch37": (
                    "official_gnomad_hail_table_non_cancer_faf95"
                ),
                "gnomad_v3_1_2_genomes_grch38": (
                    "hail.experimental.filtering_allele_frequency_from_official_non_cancer_ac_an"
                ),
            },
        )
        self.assertTrue(all(
            item["source_identity"]["etag"] and item["source_identity"]["x_goog_hash"]
            for item in data["metadata"]["extraction_log"]
        ))
        scored = set(
            manifest["classification_policies"]["enigma_brca_v1_2"]
            ["frequency_criteria"]["scored_non_founder_ancestries"]
        )
        context_by_dataset = {
            item["dataset_key"]: {
                population["code"]
                for population in item["excluded_population_context"]
            }
            for item in manifest["datasets"]
        }
        self.assertTrue(all(
            set(record["non_founder_ac_by_ancestry"]) == scored
            and record["non_founder_observed"]
            == any(
                (value or 0) > 0
                for value in record["non_founder_ac_by_ancestry"].values()
            )
            and set(record["excluded_population_context"])
            == context_by_dataset[record["dataset"]]
            for record in records
        ))
        self.assertEqual(
            sum(
                not record["non_founder_observed"]
                and any(
                    (population.get("ac") or 0) > 0
                    for population in record["excluded_population_context"].values()
                )
                for record in records
            ),
            1644,
        )
        self.assertTrue(all(
            population["used_for_ba1_bs1"] is False
            and population["used_for_pm2_presence"] is False
            for record in records
            for population in record["excluded_population_context"].values()
        ))

        for variant_id in ("13-32972540-C-CATTT", "13-32398403-C-CATTT"):
            record = data["variants"][variant_id][0]
            self.assertEqual(record["faf95_max"], 0.0)
            self.assertGreater(record["popmax_af"], 0.00002)

        # The gnomAD API exposes the main-dataset FAF95 even when a subset is
        # selected.  The scored value must be the Hail-native non-cancer value.
        self.assertEqual(
            data["variants"]["17-43068553-TAA-T"][0]["faf95_max"],
            0.00010019,
        )

    def test_gnomad_coverage_snapshot_has_official_sources_and_integrity(self):
        path = DATA / "gnomad" / "gnomad_brca_coverage_snapshot.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        mapping = data["coverage_by_position"]
        self.assertEqual(data["metadata"]["schema_version"], 2)
        self.assertEqual(data["metadata"]["records"], len(mapping))
        self.assertEqual(
            data["metadata"]["records_sha256"],
            hashlib.sha256(
                json.dumps(
                    mapping,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
            ).hexdigest(),
        )
        self.assertEqual(
            {item["dataset"] for item in data["metadata"]["datasets"]},
            {
                "gnomad_v2_1_1_exomes_grch37",
                "gnomad_v3_1_2_genomes_grch38",
            },
        )
        self.assertTrue(all(
            item["coverage_source_identity"]["etag"]
            and item["coverage_source_identity"]["x_goog_hash"]
            for item in data["metadata"]["datasets"]
        ))
        self.assertTrue(all(
            "threshold" not in item and "passes" not in item
            for item in mapping.values()
        ))

    def test_gnomad_panel_manifest_is_gene_extensible_and_pins_activation(self):
        manifest = json.loads(
            (DATA / "gnomad" / "gnomad_panel_manifest.json").read_text(encoding="utf-8")
        )
        self.assertFalse(
            manifest["classification_policies"]["enigma_brca_v1_2"]
            ["automatic_release_activation"]
        )
        self.assertEqual(
            {target["gene"] for target in manifest["targets"]},
            {"BRCA1", "BRCA2"},
        )
        self.assertTrue(all(
            set(target["intervals"]) == {"GRCh37", "GRCh38"}
            for target in manifest["targets"]
        ))
        self.assertTrue(all(
            target["activation_status"] == "active"
            and target["classification_policy_id"] == "enigma_brca_v1_2"
            and target["reference_transcript"]
            and target["vcep_specification"]["version"] == "1.2.0"
            for target in manifest["targets"]
        ))
        self.assertEqual(
            set(
                manifest["classification_policies"]["enigma_brca_v1_2"]
                ["frequency_criteria"]["scored_non_founder_ancestries"]
            ),
            {"afr", "amr", "eas", "nfe", "sas"},
        )
        self.assertFalse(
            manifest["future_gene_activation"]["inherit_existing_policy_by_default"]
        )
        self.assertFalse(
            manifest["future_gene_activation"]["activate_from_interval_only"]
        )
        self.assertEqual(
            {
                item["code"]
                for dataset in manifest["datasets"]
                for item in dataset["excluded_population_context"]
                if item["category"] == "founder_population"
            },
            {"ami", "asj", "fin"},
        )

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
        self.assertEqual(
            data["critical_boundaries"]["BRCA1"]["at_or_before"],
            {"pvs1_code": "PVS1", "pm5_code": "PM5_Strong (PTC)"},
        )
        self.assertEqual(
            data["critical_boundaries"]["BRCA1"]["after"],
            {"pvs1_code": "PVS1_N/A", "pm5_code": "PM5_N/A"},
        )
        self.assertEqual(
            data["critical_boundaries"]["BRCA2"]["at_or_before"],
            {"pvs1_code": "PVS1", "pm5_code": "PM5_Strong (PTC)"},
        )

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

        protein_ps1 = json.loads(
            (DATA / "ps1_protein_reference_registry.json").read_text(encoding="utf-8")
        )
        self.assertEqual(protein_ps1["schema_version"], 3)
        self.assertEqual(
            protein_ps1["candidate_source"]["usage"],
            "enigma_reference_set_with_per_record_ps1_eligibility",
        )
        self.assertEqual(protein_ps1["reference_count"], 60)
        self.assertEqual(protein_ps1["status_counts"], {"eligible": 40, "excluded": 20})
        self.assertTrue(all(
            record["reference_splice_evidence"]["prediction_policy"]
            == "runtime_required"
            for record in protein_ps1["references"]
        ))
        self.assertTrue(all(
            "spliceai_score" not in record["reference_splice_evidence"]
            for record in protein_ps1["references"]
        ))
        self.assertEqual(len(protein_ps1["references"]), 60)
        self.assertEqual(
            {
                item["id"]
                for item in protein_ps1["reference_source_policy"]["accepted_classification_bases"]
            },
            {
                "enigma_st7_v1_2_reference_set",
                "external_vcep_assertion",
                "locally_recurated_under_enigma_vcep",
            },
        )
        self.assertTrue(
            all(len(value) == 64 for value in protein_ps1["source_checksums"].values())
        )
        extensions = json.loads(
            (DATA / "ps1_protein_reference_extensions.json").read_text(encoding="utf-8")
        )
        self.assertEqual(extensions["schema_version"], 1)
        self.assertEqual(extensions["records"], [])

        st2 = json.loads(
            (DATA / "enigma_st2_splice_evidence.json").read_text(encoding="utf-8")
        )
        self.assertEqual(st2["schema_version"], 1)
        self.assertEqual(st2["source_columns"], 11)
        self.assertEqual(st2["total_variants"], 220)
        self.assertEqual(len(st2["source_file_sha256"]), 64)
        self.assertTrue(all(len(record) == 12 for record in st2["variants"]))


if __name__ == "__main__":
    unittest.main()
