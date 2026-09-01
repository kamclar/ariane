import unittest
import asyncio
import json
import hashlib
import tempfile
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from fastapi import HTTPException

from backend.lookups.precomputed import (
    load_classification_snapshot_metadata,
    lookup_classification_snapshot,
)
from backend.modules.reference_validation import validate_reference_allele


class PrecomputedSnapshotTests(unittest.TestCase):
    def test_pp4_snapshot_json_writer_uses_lf_on_every_platform(self):
        from scripts.build_pp4_clinical_lr_snapshot import _write_json_lf

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            _write_json_lf(path, {"one": 1, "two": 2})
            payload = path.read_bytes()

        self.assertNotIn(b"\r\n", payload)
        self.assertTrue(payload.endswith(b"\n"))

    def test_pp4_snapshot_metadata_matches_index(self):
        root = Path(__file__).resolve().parents[1]
        index_path = root / "data/precomputed/brca_pp4_clinical_lr_snapshot.index.json"
        metadata_path = root / "data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json"
        source_manifest_path = root / "data/sources/enigma/clinical_lr_sources.manifest.json"
        indel_index_path = root / "data/precomputed/brca_normalized_indel_snapshot.index.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        records = json.loads(index_path.read_text(encoding="utf-8"))
        self.assertEqual(metadata["status"], "validated_derived_snapshot")
        self.assertEqual(metadata["source_manifest"]["schema_version"], 3)
        self.assertEqual(
            metadata["automatic_application_statuses"],
            {"eligible": 12656, "review_required": 294},
        )
        self.assertEqual(metadata["rows_seen"], 13481)
        self.assertEqual(metadata["records"], 12950)
        self.assertEqual(
            metadata["criteria"],
            {"BP5": 5497, "PP4": 1720, "not_informative": 5733},
        )
        self.assertEqual(metadata["records"], len(records))
        self.assertEqual(metadata["index_sha256"], hashlib.sha256(index_path.read_bytes()).hexdigest())
        self.assertEqual(
            metadata["source_manifest_sha256"],
            hashlib.sha256(source_manifest_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(metadata["excluded"]["no_reference_transcript_hgvsc"], 210)
        self.assertEqual(
            metadata["excluded"]["reference_transcript_hgvs_not_validated"], 26
        )
        self.assertEqual(
            metadata["normalization"]["normalized_indel_dependency"]["index_sha256"],
            hashlib.sha256(indel_index_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            metadata["normalization"]["provenance"]["normalization_engine"],
            "biocommons.hgvs",
        )

    def test_conflicting_normalized_c5266_alias_requires_source_review(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        canonical = evaluate_pp4_bp5("BRCA1", "c.5266dup")
        alias = evaluate_pp4_bp5("BRCA1", "c.5266dupC")
        self.assertFalse(canonical["applies"])
        self.assertEqual(canonical["application_status"], "review_required")
        self.assertTrue(canonical["double_counting_risk"])
        self.assertFalse(canonical["automatic_combination_allowed"])
        self.assertIsNone(canonical["likelihood_ratio"])
        self.assertIsNone(canonical["candidate_likelihood_ratio"])
        self.assertEqual(
            alias["candidate_likelihood_ratio"],
            canonical["candidate_likelihood_ratio"],
        )
        self.assertEqual(canonical["overlap_status"], "conflicting_normalized_source_rows")
        self.assertEqual(
            {pmid for item in canonical["source_components"] for pmid in item["pmids"]},
            {"31853058", "40413188"},
        )

    def test_c509_uses_publisher_combined_lr_as_bp5_strong(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.509G>A")
        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "BP5")
        self.assertEqual(result["strength"], "Strong")
        self.assertEqual(result["application_status"], "applied")
        self.assertEqual(result["overlap_status"], "source_curated_combination")
        self.assertFalse(result["double_counting_risk"])
        self.assertAlmostEqual(result["likelihood_ratio"], 0.03947)
        self.assertTrue(result["source_reported_overlap_caveat"])
        self.assertEqual(
            set(result["source_components"][0]["pmids"]),
            {"31131967", "31853058", "40413188"},
        )
        self.assertNotIn("local", result["reason"].lower())
        self.assertNotIn("snapshot", result["reason"].lower())
        self.assertTrue(result["single_strong_likely_benign_eligible"])
        self.assertGreaterEqual(result["likelihood_ratio_contribution_count"], 2)
        self.assertGreaterEqual(len(result["clinical_evidence_types"]), 2)
        comparison = result["threshold_comparison"]
        self.assertEqual(comparison["status"], "match")
        self.assertEqual(comparison["source_label"], "BP5 - Benign - Strong")
        self.assertEqual(comparison["vcep_label"], "BP5 Strong")
        self.assertEqual(comparison["threshold_operator"], "<=")
        self.assertEqual(comparison["threshold_value"], 0.05)

    def test_source_label_difference_is_compared_with_vcep_threshold_generically(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.4675+3A>T")

        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "PP4")
        self.assertEqual(result["strength"], "Moderate")
        comparison = result["threshold_comparison"]
        self.assertEqual(comparison["status"], "different")
        self.assertEqual(
            comparison["source_label"],
            "PP4 - Pathogenic - Supporting",
        )
        self.assertEqual(comparison["vcep_label"], "PP4 Moderate")
        self.assertEqual(comparison["threshold_operator"], ">=")
        self.assertEqual(comparison["threshold_value"], 4.3)
        self.assertIn("combined LR=4.31493", comparison["reason"])
        self.assertIn("VCEP threshold result", result["reason"])

    def test_unknown_source_acmg_label_fails_closed(self):
        from backend.modules.pp4_bp5 import _parse_source_acmg_label

        with self.assertRaisesRegex(RuntimeError, "Unsupported PP4/BP5 source ACMG label"):
            _parse_source_acmg_label("Probably moderate")

    def test_every_eligible_snapshot_record_uses_the_same_label_comparison_rule(self):
        from backend.modules.pp4_bp5 import (
            _threshold_comparison,
            load_pp4_bp5_snapshot,
            lr_to_bp5_strength,
            lr_to_pp4_strength,
        )

        snapshot, _aliases = load_pp4_bp5_snapshot()
        comparisons = []
        for entry in snapshot.values():
            lr = entry.get("combined_lr")
            if lr is None:
                continue
            gene = entry["gene"]
            pp4_strength = lr_to_pp4_strength(gene, lr)
            bp5_strength = lr_to_bp5_strength(gene, lr)
            code = "PP4" if pp4_strength else "BP5" if bp5_strength else None
            strength = pp4_strength or bp5_strength
            comparisons.append(_threshold_comparison(
                gene,
                lr,
                code,
                strength,
                entry["source_acmg_label"],
            ))

        self.assertTrue(comparisons)
        self.assertEqual(
            {item["status"] for item in comparisons},
            {"match", "different"},
        )
        self.assertEqual(
            sum(item["status"] == "different" for item in comparisons),
            40,
        )

    def test_single_lr_bp5_strong_is_not_single_strong_class2_eligible(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.1005C>A")
        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "BP5")
        self.assertEqual(result["strength"], "Strong")
        self.assertEqual(result["likelihood_ratio_contribution_count"], 1)
        self.assertEqual(len(result["clinical_evidence_types"]), 1)
        self.assertFalse(result["single_strong_likely_benign_eligible"])

    def test_pp4_unavailable_reason_does_not_expose_storage_terminology(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.999999A>G")
        self.assertFalse(result["applies"])
        self.assertIn("ENIGMA v1.2", result["reason"])
        self.assertNotIn("local", result["reason"].lower())
        self.assertNotIn("snapshot", result["reason"].lower())

    def test_bp5_snapshot_resolves_multibase_duplication_alias(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        canonical = evaluate_pp4_bp5("BRCA2", "c.9891_9894dup")
        source_spelling = evaluate_pp4_bp5("BRCA2", "c.9891_9894dupATTT")
        self.assertTrue(canonical["applies"])
        self.assertEqual(canonical["code"], "BP5")
        self.assertEqual(canonical["strength"], "Supporting")
        self.assertEqual(canonical["points"], -1)
        self.assertAlmostEqual(canonical["likelihood_ratio"], 0.41018)
        self.assertEqual(source_spelling["likelihood_ratio"], canonical["likelihood_ratio"])
        self.assertEqual(set(canonical["source_components"][0]["pmids"]), {"31853058"})

    def test_conflicting_normalized_c475_del_rows_require_source_review(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA2", "c.475+4del")
        self.assertFalse(result["applies"])
        self.assertEqual(result["application_status"], "review_required")
        self.assertEqual(result["overlap_status"], "conflicting_normalized_source_rows")
        self.assertIsNone(result["likelihood_ratio"])
        self.assertEqual(
            {pmid for item in result["source_components"] for pmid in item["pmids"]},
            {"31131967", "31853058", "40413188"},
        )

    def test_c3891_del_uses_publisher_combined_lr_as_bp5_strong(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.3891_3893del")
        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "BP5")
        self.assertEqual(result["strength"], "Strong")
        self.assertEqual(result["application_status"], "applied")
        self.assertAlmostEqual(result["likelihood_ratio"], 0.02896)
        self.assertEqual(
            set(result["source_components"][0]["pmids"]),
            {"31131967", "34597585", "40413188"},
        )

    def test_pp4_uses_publisher_combined_clinical_lr(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.4185G>A")
        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "PP4")
        self.assertEqual(result["strength"], "Strong")
        self.assertEqual(result["points"], 4)
        self.assertAlmostEqual(result["likelihood_ratio"], 328.18363)
        self.assertEqual(
            set(result["source_components"][0]["pmids"]),
            {"31131967", "31853058"},
        )

    def test_c4185_full_path_routes_unquantified_rna_to_review(self):
        from backend.classification_dag import ClassifierEngineMode
        from backend.main import CLASSIFICATION_ORCHESTRATION, _classify_one
        from backend.services import EvidenceOrchestrationService, ExternalEvidenceDependencies

        pm2 = {
            "PM2_Supporting": {
                "applies": True,
                "strength": "Supporting",
                "points": 1,
                "reason": "Absent with sufficient coverage in both required gnomAD datasets.",
            }
        }
        coordinates = {
            "chrom": "17", "pos": 43090944, "ref": "C", "alt": "T",
            "assembly": "GRCh38",
        }
        dependencies = replace(
            CLASSIFICATION_ORCHESTRATION.provider_dependencies,
            resolve_variant=lambda *_args, **_kwargs: None,
            get_grch37=lambda *_args, **_kwargs: coordinates,
            get_grch38=lambda *_args, **_kwargs: coordinates,
            gnomad_lookup=lambda **_kwargs: {"available": True},
            spliceai_lookup=lambda *_args, **_kwargs: 0.95,
            spliceai_status=lambda *_args, **_kwargs: {
                "status": "ok", "score": 0.95, "source": "test"
            },
            bayesdel_lookup=lambda *_args, **_kwargs: (None, None),
            bayesdel_status=lambda *_args, **_kwargs: {"status": "not_found"},
        )
        orchestration = EvidenceOrchestrationService(
            engine_mode=ClassifierEngineMode.DAG,
            provider_dependencies=dependencies,
            external_dependencies=ExternalEvidenceDependencies(
                clinvar_lookup=lambda *_args, **_kwargs: {"status": "not_found"},
                clingen_lookup=lambda *_args, **_kwargs: {"status": "not_found"},
            ),
        )
        with patch(
            "backend.classification_dag.nodes.population.evaluate_frequency_criteria",
            return_value=pm2,
        ), patch("backend.main.CLASSIFICATION_ORCHESTRATION", orchestration):
            result = asyncio.run(
                _classify_one("BRCA1", "c.4185G>A", "p.(Gln1395=)")
            )

        criteria = {criterion.name: criterion for criterion in result.criteria}
        self.assertEqual(set(criteria), {"PP3", "PM2_Supporting", "PP4"})
        self.assertEqual(criteria["PP3"].strength, "Supporting")
        self.assertEqual(criteria["PP4"].strength, "Strong")
        self.assertEqual(result.total_points, 6)
        self.assertEqual(result.predicted_class, 4)
        self.assertTrue(result.rna_review.recommended)
        self.assertEqual(result.rna_review.priority, "high")
        self.assertEqual(
            result.rna_review.manual_review_prefill["transcript_accession"],
            "NM_007294.4",
        )

    def test_pp4_snapshot_missing_metadata_fails_closed(self):
        from backend.modules import pp4_bp5

        original_snapshot, original_aliases = pp4_bp5._SNAPSHOT, pp4_bp5._ALIASES
        try:
            pp4_bp5._SNAPSHOT = None
            pp4_bp5._ALIASES = None
            with patch.object(pp4_bp5, "METADATA_PATH", Path("missing-pp4-metadata.json")):
                with self.assertRaisesRegex(RuntimeError, "metadata is missing"):
                    pp4_bp5.load_pp4_bp5_snapshot()
        finally:
            pp4_bp5._SNAPSHOT, pp4_bp5._ALIASES = original_snapshot, original_aliases

    def test_pp4_snapshot_missing_source_manifest_fails_closed(self):
        from backend.modules import pp4_bp5

        original_snapshot, original_aliases = pp4_bp5._SNAPSHOT, pp4_bp5._ALIASES
        try:
            pp4_bp5._SNAPSHOT = None
            pp4_bp5._ALIASES = None
            with patch.object(
                pp4_bp5, "SOURCE_MANIFEST_PATH", Path("missing-clinical-lr-manifest.json")
            ):
                with self.assertRaisesRegex(RuntimeError, "source manifest is missing"):
                    pp4_bp5.load_pp4_bp5_snapshot()
        finally:
            pp4_bp5._SNAPSHOT, pp4_bp5._ALIASES = original_snapshot, original_aliases

    def test_normalized_indel_snapshot_resolves_alias_and_protein(self):
        from backend.lookups.indels import lookup_indel_snapshot

        record = lookup_indel_snapshot("BRCA1", "c.5266dupC")
        self.assertIsNotNone(record)
        self.assertEqual(record["canonical_c_notation"], "c.5266dup")
        self.assertEqual(record["p_notation"], "p.(Gln1756ProfsTer74)")
        self.assertEqual(record["reference_transcript"], "NM_007294.4")
        self.assertEqual(record["grch37"]["ref"], "T")
        self.assertEqual(record["grch38"]["alt"], "TG")

    def test_indel_snapshot_records_and_excludes_ambiguous_aliases(self):
        from backend.lookups.indels import load_indel_snapshot, lookup_indel_snapshot

        index, aliases = load_indel_snapshot()
        self.assertEqual(len(index), 16511)
        self.assertEqual(aliases["BRCA2:c.3975_3978dup"], "BRCA2:c.3975_3978dup")
        self.assertNotIn("c.3975_3978dup", index["BRCA2:c.3975dup"]["input_c_notations"])
        self.assertIsNotNone(lookup_indel_snapshot("BRCA2", "c.3975_3978dup"))

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

    def test_reference_allele_accepts_real_brca1_change(self):
        validate_reference_allele("BRCA1", "c.181T>G")

    def test_reference_allele_rejects_wrong_brca1_change(self):
        with self.assertRaisesRegex(ValueError, r"is T .* not A"):
            validate_reference_allele("BRCA1", "c.181A>C")


class ClassificationInputIntegrationTests(unittest.TestCase):
    def test_c4676_acceptor_variant_with_two_very_strong_criteria_is_pathogenic(self):
        from backend.main import _classify_one

        with patch(
            "backend.lookups.spliceai.get_spliceai_score", return_value=0.996
        ), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense",
            return_value=(None, None),
        ), patch(
            "backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}
        ), patch(
            "backend.lookups.clingen.clingen_erepo_lookup",
            return_value={"status": "not_found"},
        ):
            result = asyncio.run(_classify_one("BRCA1", "c.4676-1G>A", "p.?"))

        criteria = {criterion.name: criterion for criterion in result.criteria}
        self.assertEqual(criteria["PVS1"].strength, "Very Strong")
        self.assertEqual(criteria["PP4"].strength, "Very Strong")
        self.assertNotIn("PM2_Supporting", criteria)
        self.assertTrue(
            any("coverage-region method" in warning for warning in result.warnings)
        )
        self.assertEqual(result.total_points, 16)
        self.assertEqual(result.predicted_class, 5)
        self.assertEqual(result.predicted_label, "Pathogenic")

    def test_abbreviated_frameshift_is_classified_with_canonical_output(self):
        from backend.main import _classify_one

        with patch("backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}), patch(
            "backend.lookups.clingen.clingen_erepo_lookup", return_value={"status": "not_found"}
        ):
            result = asyncio.run(
                _classify_one("BRCA1", "c.3668_3671dup", "p.(Cys1225fs)")
            )

        self.assertEqual(result.p_notation, "p.(Cys1225SerfsTer10)")
        self.assertEqual(result.predicted_class, 5)

    def test_terminal_frameshift_receives_paired_table4_codes(self):
        from backend.main import _classify_one

        with patch("backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}), patch(
            "backend.lookups.clingen.clingen_erepo_lookup", return_value={"status": "not_found"}
        ):
            result = asyncio.run(
                _classify_one("BRCA1", "c.5556_5560del", "p.(Tyr1853AspfsTer25)")
            )

        criteria = {criterion.name: criterion for criterion in result.criteria}
        self.assertEqual(criteria["PVS1"].strength, "Very Strong")
        self.assertEqual(criteria["PM5_PTC"].strength, "Strong")
        self.assertEqual(result.total_points, 12)
        self.assertEqual(result.predicted_class, 5)

    def test_exon_cnv_skips_small_variant_lookups_and_hides_provider_errors(self):
        from backend.main import _classify_one

        notation = "c.(793+1_794-1)_(1909+1_1910-1)del"
        with patch("backend.lookups.coordinates.resolve_variant") as coordinates, patch(
            "backend.lookups.spliceai.get_spliceai_score"
        ) as spliceai, patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense"
        ) as bayesdel, patch(
            "backend.lookups.clinvar.clinvar_lookup",
            return_value={"status": "api_timeout", "error": "raw provider error"},
        ), patch(
            "backend.lookups.clingen.clingen_erepo_lookup",
            return_value={"status": "api_timeout", "error": "raw provider error"},
        ):
            result = asyncio.run(_classify_one("BRCA2", notation, "p.(?)"))

        coordinates.assert_not_called()
        spliceai.assert_not_called()
        bayesdel.assert_not_called()
        warnings = "\n".join(result.warnings)
        self.assertIn(
            "PVS1 was not applied: ENIGMA Table 4 marks the BRCA2 E10 deletion as PVS1 N/A.",
            warnings,
        )
        self.assertIn("uncertain genomic breakpoints", warnings)
        self.assertIn("ClinVar comparison is temporarily unavailable", warnings)
        self.assertIn("ClinGen ERepo comparison is temporarily unavailable", warnings)
        self.assertNotIn("SpliceAI not available", warnings)
        self.assertNotIn("raw provider error", warnings)
        self.assertNotIn("HTTP 422", warnings)
        criteria = {criterion.name: criterion for criterion in result.criteria}
        self.assertEqual(criteria["PM2_Supporting"].strength, "Supporting")
        self.assertEqual(criteria["PM2_Supporting"].points, 1)
        self.assertNotIn("BS3", criteria)
        not_applicable = {
            criterion.name: criterion
            for criterion in result.not_applicable_criteria
        }
        self.assertEqual(not_applicable["PVS1"].strength, "N/A")
        self.assertNotIn(
            "PVS1",
            {criterion.name for criterion in result.excluded_criteria},
        )
        self.assertEqual(result.total_points, 1)
        self.assertEqual(result.predicted_class, 3)
        self.assertFalse(result.mixed_evidence)
        self.assertIn(
            "no applicable variant-specific recommendation in ENIGMA Table 9",
            "\n".join(result.warnings),
        )

    def test_c5266_conflicting_clinical_lr_rows_require_review(self):
        from backend.main import _classify_one

        with patch("backend.lookups.spliceai.get_spliceai_score", return_value=None), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense", return_value=(None, None)
        ), patch("backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}), patch(
            "backend.lookups.clingen.clingen_erepo_lookup", return_value={"status": "not_found"}
        ):
            result = asyncio.run(
                _classify_one("BRCA1", "c.5266dup", "p.(Gln1756ProfsTer74)")
            )

        self.assertNotIn("PP4", {criterion.name for criterion in result.criteria})
        self.assertEqual(result.clinical_lr_audit.application_status, "review_required")
        self.assertTrue(result.clinical_lr_audit.double_counting_risk)
        self.assertTrue(
            any("expert source review is required" in item.lower() for item in result.warnings)
        )

    def test_c509_full_classification_scores_publisher_combined_bp5(self):
        from backend.main import _classify_one

        with patch("backend.lookups.coordinates.resolve_variant", return_value=None), patch(
            "backend.lookups.spliceai.get_spliceai_score", return_value=0.02
        ), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense",
            return_value=(None, None),
        ), patch(
            "backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}
        ), patch(
            "backend.lookups.clingen.clingen_erepo_lookup",
            return_value={"status": "not_found"},
        ):
            result = asyncio.run(
                _classify_one("BRCA1", "c.509G>A", "p.(Arg170Gln)")
            )

        criteria = {criterion.name: criterion for criterion in result.criteria}
        self.assertNotIn("PP4", criteria)
        self.assertEqual(criteria["BP5"].strength, "Strong")
        self.assertEqual(criteria["BP5"].points, -4)
        self.assertEqual(result.clinical_lr_audit.application_status, "applied")
        self.assertAlmostEqual(
            result.clinical_lr_audit.likelihood_ratio,
            0.03947,
        )

    def test_c3247a_to_c_receives_bp1_and_variant_specific_pp4(self):
        from backend.main import _classify_one

        with patch("backend.lookups.coordinates.resolve_variant", return_value=None), patch(
            "backend.lookups.spliceai.get_spliceai_score", return_value=0.0
        ), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense",
            return_value=(None, None),
        ), patch(
            "backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}
        ), patch(
            "backend.lookups.clingen.clingen_erepo_lookup",
            return_value={"status": "not_found"},
        ):
            result = asyncio.run(
                _classify_one("BRCA1", "c.3247A>C", "p.(Met1083Leu)")
            )

        criteria = {criterion.name: criterion for criterion in result.criteria}
        self.assertEqual(criteria["BP1"].strength, "Strong")
        self.assertEqual(criteria["BP1"].points, -4)
        self.assertEqual(criteria["PP4"].strength, "Moderate")
        self.assertEqual(criteria["PP4"].points, 2)
        self.assertIn("combined LR=8.9313", criteria["PP4"].reason)
        self.assertIn("PMID 31853058, 40413188", criteria["PP4"].reason)
        self.assertEqual(result.total_points, -2)
        self.assertEqual(result.predicted_class, 2)

    def test_brca2_multibase_duplication_receives_bp5_supporting(self):
        from backend.main import _classify_one

        with patch("backend.lookups.spliceai.get_spliceai_score", return_value=None), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense", return_value=(None, None)
        ), patch("backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}), patch(
            "backend.lookups.clingen.clingen_erepo_lookup", return_value={"status": "not_found"}
        ):
            result = asyncio.run(
                _classify_one("BRCA2", "c.9891_9894dup", "p.(Gln3299fs)")
            )

        criteria = {criterion.name: criterion for criterion in result.criteria}
        self.assertEqual(criteria["PVS1"].strength, "Very Strong")
        self.assertEqual(criteria["PM5_PTC"].strength, "Strong")
        self.assertEqual(criteria["BP5"].strength, "Supporting")
        self.assertEqual(criteria["BP5"].points, -1)
        self.assertEqual(result.total_points, 11)
        self.assertEqual(result.predicted_class, 5)
        self.assertTrue(result.mixed_evidence)

    def test_general_indel_snapshot_rejects_random_protein_notation(self):
        from backend.main import _classify_one

        with self.assertRaises(HTTPException) as raised:
            asyncio.run(_classify_one("BRCA1", "c.3668_3671dup", "p.(Arg100Gly)"))
        self.assertEqual(raised.exception.status_code, 422)
        self.assertIn("p.(Cys1225SerfsTer10)", raised.exception.detail)

    def test_table9_indel_rejects_random_protein_notation(self):
        from backend.main import _classify_one

        with self.assertRaises(HTTPException) as raised:
            asyncio.run(
                _classify_one("BRCA1", "c.5266dup", "p.(Arg100Gly)")
            )
        self.assertEqual(raised.exception.status_code, 422)
        self.assertIn("p.(Gln1756ProfsTer74)", raised.exception.detail)

    def test_table9_asterisk_and_ter_protein_notations_are_equivalent(self):
        from backend.modules.hgvs import normalize_protein_notation
        from backend.modules.table9 import table9_protein_notation

        reviewed = table9_protein_notation("BRCA1", "c.5266dup")
        self.assertEqual(
            normalize_protein_notation(reviewed),
            normalize_protein_notation("p.(Gln1756ProfsTer74)"),
        )

    def test_c_only_nonsense_derives_protein_consequence(self):
        from backend.main import _classify_one

        with patch("backend.lookups.coordinates.resolve_variant", return_value=None), patch(
            "backend.lookups.spliceai.get_spliceai_score", return_value=None
        ), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense", return_value=(None, None)
        ), patch(
            "backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}
        ), patch(
            "backend.lookups.clingen.clingen_erepo_lookup", return_value={"status": "not_found"}
        ):
            result = asyncio.run(_classify_one("BRCA1", "c.303T>G", ""))
        self.assertEqual(result.p_notation, "p.(Tyr101Ter)")
        self.assertEqual(result.consequence_status, "sequence_derived")

    def test_nonsense_with_protein_notation_is_classified(self):
        from backend.main import _classify_one

        with patch("backend.lookups.spliceai.get_spliceai_score", return_value=None), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense", return_value=(None, None)
        ), patch("backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}), patch(
            "backend.lookups.clingen.clingen_erepo_lookup", return_value={"status": "not_found"}
        ):
            result = asyncio.run(
                _classify_one("BRCA1", "c.303T>G", "p.(Tyr101Ter)")
            )

        self.assertEqual(result.p_notation, "p.(Tyr101Ter)")
        pvs1 = next(criterion for criterion in result.criteria if criterion.name == "PVS1")
        self.assertTrue(pvs1.applies)
        self.assertEqual(result.predicted_class, 5)

    def test_wrong_reference_stops_before_classification(self):
        from backend.main import _classify_one

        with self.assertRaises(HTTPException) as raised:
            asyncio.run(_classify_one("BRCA1", "c.181A>C", ""))
        self.assertEqual(raised.exception.status_code, 422)
        self.assertIn("is T", raised.exception.detail)


if __name__ == "__main__":
    unittest.main()
