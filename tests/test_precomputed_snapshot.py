import unittest
import asyncio
import json
import hashlib
from pathlib import Path
from unittest.mock import patch

from fastapi import HTTPException

from backend.lookups.precomputed import (
    load_classification_snapshot_metadata,
    lookup_classification_snapshot,
)
from backend.modules.reference_validation import validate_reference_allele


class PrecomputedSnapshotTests(unittest.TestCase):
    def test_pp4_snapshot_metadata_matches_index(self):
        root = Path(__file__).resolve().parents[1]
        index_path = root / "data/precomputed/brca_pp4_clinical_lr_snapshot.index.json"
        metadata_path = root / "data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json"
        indel_index_path = root / "data/precomputed/brca_normalized_indel_snapshot.index.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        records = json.loads(index_path.read_text(encoding="utf-8"))
        self.assertEqual(metadata["status"], "validated_derived_snapshot")
        self.assertEqual(metadata["records"], len(records))
        self.assertEqual(metadata["index_sha256"], hashlib.sha256(index_path.read_bytes()).hexdigest())
        self.assertEqual(
            metadata["normalization"]["normalized_indel_dependency"]["index_sha256"],
            hashlib.sha256(indel_index_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            metadata["normalization"]["provenance"]["normalization_engine"],
            "biocommons.hgvs",
        )

    def test_pp4_snapshot_resolves_c5266_alias_as_very_strong(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        canonical = evaluate_pp4_bp5("BRCA1", "c.5266dup")
        alias = evaluate_pp4_bp5("BRCA1", "c.5266dupC")
        self.assertTrue(canonical["applies"])
        self.assertEqual(canonical["code"], "PP4")
        self.assertEqual(canonical["strength"], "Very Strong")
        self.assertEqual(canonical["points"], 8)
        self.assertAlmostEqual(canonical["likelihood_ratio"], 6.89647e45)
        self.assertEqual(alias["likelihood_ratio"], canonical["likelihood_ratio"])
        self.assertEqual({item["pmid"] for item in canonical["source_components"]}, {"31853058"})

    def test_pp4_public_reason_uses_enigma_method_not_storage_terminology(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.509G>A")
        self.assertTrue(result["applies"])
        self.assertEqual(result["strength"], "Moderate")
        self.assertIn("ENIGMA v1.2 Appendix B", result["reason"])
        self.assertIn("combined LR=6.1764", result["reason"])
        self.assertNotIn("local", result["reason"].lower())
        self.assertNotIn("snapshot", result["reason"].lower())

    def test_pp4_unavailable_reason_does_not_expose_storage_terminology(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.999999A>G")
        self.assertFalse(result["applies"])
        self.assertIn("ENIGMA v1.2 Appendix B", result["reason"])
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
        self.assertEqual(
            {item["pmid"] for item in canonical["source_components"]}, {"31853058"}
        )

    def test_pp4_snapshot_combines_independent_rows_after_hgvs_normalization(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA2", "c.475+4del")
        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "PP4")
        self.assertEqual(result["strength"], "Strong")
        self.assertAlmostEqual(result["likelihood_ratio"], 73.38448758438403)
        self.assertEqual(
            {item["pmid"] for item in result["source_components"]},
            {"31131967", "31853058"},
        )

    def test_pp4_uses_combined_clinical_lr_not_historical_posterior_probability(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.3891_3893del")
        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "PP4")
        self.assertEqual(result["strength"], "Strong")
        self.assertEqual(result["points"], 4)
        self.assertAlmostEqual(result["likelihood_ratio"], 28.554690389559802)
        self.assertEqual(
            {item["pmid"] for item in result["source_components"]}, {"31131967"}
        )

    def test_pp4_combines_available_appendix_b_clinical_lr_components(self):
        from backend.modules.pp4_bp5 import evaluate_pp4_bp5

        result = evaluate_pp4_bp5("BRCA1", "c.4185G>A")
        self.assertTrue(result["applies"])
        self.assertEqual(result["code"], "PP4")
        self.assertEqual(result["strength"], "Strong")
        self.assertEqual(result["points"], 4)
        self.assertAlmostEqual(result["likelihood_ratio"], 328.183625413548)
        self.assertEqual(
            {item["pmid"] for item in result["source_components"]},
            {"31131967", "31853058"},
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

    def test_c5266_automatically_receives_pp4_very_strong(self):
        from backend.main import _classify_one

        with patch("backend.lookups.spliceai.get_spliceai_score", return_value=None), patch(
            "backend.lookups.bayesdel.get_bayesdel_and_alphamissense", return_value=(None, None)
        ), patch("backend.lookups.clinvar.clinvar_lookup", return_value={"status": "not_found"}), patch(
            "backend.lookups.clingen.clingen_erepo_lookup", return_value={"status": "not_found"}
        ):
            result = asyncio.run(
                _classify_one("BRCA1", "c.5266dup", "p.(Gln1756ProfsTer74)")
            )

        pp4 = next(criterion for criterion in result.criteria if criterion.name == "PP4")
        self.assertTrue(pp4.applies)
        self.assertEqual(pp4.strength, "Very Strong")
        self.assertEqual(pp4.points, 8)

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
