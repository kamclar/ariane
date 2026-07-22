import asyncio
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from backend.data_validation import validate_required_datasets
from backend.data_health import clear_issue, get_data_issues, get_user_warnings, register_issue
from backend.lookup_execution import lookup_or_unavailable


class RequiredDatasetValidationTests(unittest.TestCase):
    def test_missing_required_dataset_stops_startup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "Table 4.*missing"):
                validate_required_datasets({
                    "table4": root / "missing-table4.json",
                    "table9": root / "missing-table9.json",
                    "st7": root / "missing-st7.json",
                })

    def test_invalid_required_dataset_stops_startup_with_reason(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table4 = root / "table4.json"
            table4.write_text("not-json", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "cannot be loaded"):
                validate_required_datasets({
                    "table4": table4,
                    "table9": root / "table9.json",
                    "st7": root / "st7.json",
                })


class GnomadFailClosedTests(unittest.TestCase):
    def test_fixture_is_never_selected_automatically(self):
        from backend.modules import frequency

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "fixture.json"
            fixture.write_text("{}", encoding="utf-8")
            with patch.object(frequency, "GNOMAD_CACHE_WITH_REAL_COVERAGE", root / "missing.json"), patch.object(
                frequency, "GNOMAD_CACHE_FIXTURE", fixture
            ):
                self.assertIsNone(frequency.choose_gnomad_cache_file())

    def test_untrusted_cache_cannot_produce_frequency_evidence(self):
        from backend.modules import frequency

        old_cache = frequency.GNOMAD_CACHE
        old_mode = frequency.GNOMAD_CACHE_MODE
        try:
            frequency.GNOMAD_CACHE = {"17-1-A-G": [{"dataset": "fixture"}]}
            frequency.GNOMAD_CACHE_MODE = "fixture_or_no_real_coverage"
            result = frequency.query_gnomad_dataset_local(
                "17-1-A-G",
                {"chrom": "17", "pos": 1, "ref": "A", "alt": "G"},
                frequency.GNOMAD_LOCAL_DATASET_CONFIG["v2_1_non_cancer"],
            )
            self.assertEqual(result["status"], "cache_untrusted")
            self.assertIn("not approved", result["errors"][0])
        finally:
            frequency.GNOMAD_CACHE = old_cache
            frequency.GNOMAD_CACHE_MODE = old_mode


class LookupDiagnosticsTests(unittest.IsolatedAsyncioTestCase):
    async def test_exception_is_logged_explained_and_returns_unavailable_default(self):
        diagnostics = []

        def failing_lookup():
            raise ConnectionError("service refused connection")

        with self.assertLogs("backend.lookup_execution", level="ERROR") as logs:
            result = await lookup_or_unavailable(
                failing_lookup, None, "Example service", diagnostics
            )
        self.assertIsNone(result)
        self.assertIn("ConnectionError", diagnostics[0])
        self.assertIn("service refused connection", diagnostics[0])
        self.assertIn("Example service", " ".join(logs.output))

    def test_bayesdel_preserves_api_failure_reason(self):
        from backend.lookups import bayesdel

        key = "BRCA1:c.999999A>G"
        bayesdel.BAYESDEL_CACHE.pop(key, None)
        bayesdel.BAYESDEL_STATUS_CACHE.pop(key, None)
        failed_coordinates = MagicMock()
        failed_coordinates.has_grch37.return_value = False
        with patch("backend.lookups.coordinates.resolve_variant", return_value=failed_coordinates), patch.object(
            bayesdel, "_save_cache"
        ):
            score, _ = bayesdel.get_bayesdel_and_alphamissense("BRCA1", "c.999999A>G")
        self.assertIsNone(score)
        self.assertEqual(bayesdel.BAYESDEL_STATUS_CACHE[key]["status"], "no_grch37_coords")
        self.assertIn("No GRCh37", bayesdel.BAYESDEL_STATUS_CACHE[key]["reason"])
        self.assertNotIn(key, bayesdel.BAYESDEL_CACHE)
        bayesdel.BAYESDEL_CACHE.pop(key, None)


class DataHealthTests(unittest.TestCase):
    def test_spliceai_runtime_cache_prefers_explicit_directory(self):
        from backend.lookups import spliceai

        with patch.dict(
            "os.environ",
            {
                "ARIANE_RUNTIME_CACHE_DIR": "/persistent/ariane",
                "RAILWAY_VOLUME_MOUNT_PATH": "/railway-volume",
            },
        ):
            self.assertEqual(
                spliceai.choose_runtime_cache_dir(), Path("/persistent/ariane")
            )

    def test_spliceai_runtime_cache_uses_railway_volume(self):
        from backend.lookups import spliceai

        with patch.dict(
            "os.environ",
            {"RAILWAY_VOLUME_MOUNT_PATH": "/railway-volume"},
            clear=True,
        ):
            self.assertEqual(
                spliceai.choose_runtime_cache_dir(),
                Path("/railway-volume/ariane-runtime-cache"),
            )

    def test_spliceai_cache_write_failure_explains_current_score_is_usable(self):
        from backend.lookups import spliceai

        clear_issue("SpliceAI API cache")
        with patch.object(
            spliceai.tempfile,
            "NamedTemporaryFile",
            side_effect=OSError(30, "Read-only file system"),
        ):
            saved = spliceai._save_api_cache({"example": {"score": 0.1}})
        self.assertFalse(saved)
        warning = next(
            item for item in get_user_warnings() if "SpliceAI API cache" in item
        )
        self.assertIn("score was obtained and used", warning)
        self.assertIn("this request is unaffected", warning)
        self.assertTrue(warning.startswith("Runtime cache persistence warning:"))
        self.assertNotIn("Data source degraded", warning)
        clear_issue("SpliceAI API cache")

    def test_registered_degradation_is_visible_to_user(self):
        clear_issue("test cache")
        register_issue("test cache", "checksum mismatch")
        self.assertIn(
            {"component": "test cache", "reason": "checksum mismatch"},
            get_data_issues(),
        )
        self.assertTrue(any("test cache" in warning and "checksum mismatch" in warning for warning in get_user_warnings()))
        clear_issue("test cache")

    def test_degradation_messages_hide_linux_deployment_path(self):
        clear_issue("test cache")
        register_issue(
            "test cache",
            "metadata is missing: /home/ubuntu/ariane/data/spliceai/cache.metadata.json",
        )
        issue = next(item for item in get_data_issues() if item["component"] == "test cache")
        self.assertEqual(
            issue["reason"],
            "metadata is missing: …ariane/data/spliceai/cache.metadata.json",
        )
        self.assertNotIn("/home/ubuntu", get_user_warnings()[0])
        clear_issue("test cache")

    def test_degradation_messages_hide_windows_deployment_path(self):
        clear_issue("test cache")
        register_issue(
            "test cache",
            r"cache is missing: F:\UOCHB\Enigma\ARIANE_app\ariane\data\cache.json",
        )
        issue = next(item for item in get_data_issues() if item["component"] == "test cache")
        self.assertEqual(issue["reason"], "cache is missing: …ariane/data/cache.json")
        clear_issue("test cache")

    def test_incomplete_intronic_spliceai_cache_is_not_loaded(self):
        from backend.lookups import spliceai

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coding = root / "coding.json"
            coding.write_text("{}", encoding="utf-8")
            intronic = root / "intronic.json"
            intronic.write_text(json.dumps({"BRCA1:c.1+1G>A": {"status": "ok", "score": 1.0}}), encoding="utf-8")
            old_cache = spliceai.SPLICEAI_PRECOMPUTED_CACHE
            try:
                spliceai.SPLICEAI_PRECOMPUTED_CACHE = None
                with patch.object(spliceai, "SPLICEAI_PRECOMPUTED_CACHE_PATH", coding), patch.object(
                    spliceai, "SPLICEAI_INTRONIC_CACHE_PATH", intronic
                ):
                    loaded = spliceai._load_precomputed_cache()
                self.assertNotIn("BRCA1:c.1+1G>A", loaded)
                self.assertTrue(any(
                    issue["component"] == "SpliceAI intronic cache" and "metadata is missing" in issue["reason"]
                    for issue in get_data_issues()
                ))
            finally:
                spliceai.SPLICEAI_PRECOMPUTED_CACHE = old_cache
                clear_issue("SpliceAI intronic cache")


class ClinVarAmbiguityTests(unittest.TestCase):
    def test_multiple_nonmatching_hgvs_results_are_ambiguous(self):
        from backend.lookups import clinvar

        failed_coordinates = MagicMock(status="failed")
        search_response = MagicMock()
        search_response.__enter__.return_value.read.return_value = json.dumps({
            "esearchresult": {"idlist": ["111", "222"]}
        }).encode()
        summary_response = MagicMock()
        summary_response.__enter__.return_value.read.return_value = json.dumps({
            "result": {
                "111": {"title": "NM_007294.4:c.100A>G"},
                "222": {"title": "NM_007294.4:c.200A>G"},
            }
        }).encode()

        with patch.object(clinvar, "resolve_variant", return_value=failed_coordinates), patch.object(
            clinvar.urllib.request, "urlopen", side_effect=[search_response, summary_response]
        ):
            result = clinvar.clinvar_search_variation_id("BRCA1", "c.300A>G")

        self.assertEqual(result["status"], "ambiguous")
        self.assertEqual(result["candidate_ids"], ["111", "222"])

    def test_multiple_exact_hgvs_results_are_ambiguous(self):
        from backend.lookups import clinvar

        failed_coordinates = MagicMock(status="failed")
        search_response = MagicMock()
        search_response.__enter__.return_value.read.return_value = json.dumps({
            "esearchresult": {"idlist": ["111", "222"]}
        }).encode()
        summary_response = MagicMock()
        summary_response.__enter__.return_value.read.return_value = json.dumps({
            "result": {
                "111": {"title": "NM_007294.4:c.300A>G"},
                "222": {"title": "NM_007294.4:c.300A>G alternate condition"},
            }
        }).encode()
        with patch.object(clinvar, "resolve_variant", return_value=failed_coordinates), patch.object(
            clinvar.urllib.request, "urlopen", side_effect=[search_response, summary_response]
        ):
            result = clinvar.clinvar_search_variation_id("BRCA1", "c.300A>G")
        self.assertEqual(result["status"], "ambiguous")
        self.assertEqual(result["candidate_ids"], ["111", "222"])


class RemainingFallbackTests(unittest.TestCase):
    def test_unknown_pvs1_code_has_no_implicit_weight(self):
        from backend.modules.table4 import parse_pvs1_code_strength

        self.assertEqual(parse_pvs1_code_strength("PVS1_UNRECOGNISED"), (None, 0, False))

    def test_exon_cnv_requires_exact_boundaries(self):
        from backend.modules.table4 import parse_exon_from_duplication_notation

        exact = "c.(80+1_81-1)_(134+1_135-1)dup"
        shifted = "c.(80+1_82-1)_(133+1_135-1)dup"
        self.assertEqual(parse_exon_from_duplication_notation(exact, "BRCA1"), "E3")
        self.assertIsNone(parse_exon_from_duplication_notation(shifted, "BRCA1"))

    def test_total_coordinate_failure_remains_retryable(self):
        from backend.lookups import coordinates

        key = "BRCA1:c.99998A>G"
        coordinates._RESOLVER_CACHE.pop(key, None)
        with patch.object(coordinates, "_resolve_precomputed_snapshot", return_value=None), patch.object(
            coordinates, "_resolve_variantvalidator", return_value=None
        ), patch.object(coordinates, "_resolve_mutalyzer", return_value=None), patch.object(
            coordinates.time, "sleep"
        ):
            result = coordinates.resolve_variant("BRCA1", "c.99998A>G")
        self.assertEqual(result.status, "failed")
        self.assertNotIn(key, coordinates._RESOLVER_CACHE)

    def test_clingen_api_failure_remains_retryable(self):
        from backend.lookups import clingen

        key = "BRCA1:c.301A>G"
        clingen.EREPO_CACHE.pop(key, None)
        with patch.object(clingen.urllib.request, "urlopen", side_effect=OSError("offline")):
            result = clingen.clingen_erepo_lookup("BRCA1", "c.301A>G")
        self.assertEqual(result["status"], "api_error")
        self.assertNotIn(key, clingen.EREPO_CACHE)

    def test_duplication_rule_never_substitutes_another_arrangement(self):
        from backend.modules import table4

        old_rules = table4.TABLE4_DATA["duplication_rules"]
        try:
            table4.TABLE4_DATA["duplication_rules"] = {
                "BRCA1": {"E2": {"Tandem": {"pvs1_code": "PVS1_Strong"}}},
                "BRCA2": {},
            }
            result = table4.table4_lookup_duplication("BRCA1", "E2", "Unknown")
        finally:
            table4.TABLE4_DATA["duplication_rules"] = old_rules
        self.assertFalse(result["found"])
        self.assertEqual(result["pvs1_points"], 0)
        self.assertIn("No exact", result["reason"])

    def test_clingen_multiple_interpretations_are_ambiguous(self):
        from backend.lookups import clingen

        response = MagicMock()
        response.__enter__.return_value.read.return_value = json.dumps({
            "variantInterpretations": [{"caid": "CA1"}, {"caid": "CA2"}]
        }).encode()
        clingen.EREPO_CACHE.clear()
        with patch.object(clingen.urllib.request, "urlopen", return_value=response):
            result = clingen.clingen_erepo_lookup("BRCA1", "c.300A>G")
        self.assertEqual(result["status"], "ambiguous")
        self.assertEqual(result["candidate_caids"], ["CA1", "CA2"])

    def test_failed_pvs1_evaluation_is_visible_for_frameshift(self):
        from backend.modules.classifier import evaluate_variant

        result = evaluate_variant(
            gene="BRCA1",
            variant_type="frameshift",
            p_notation="p.(Gly9999ValfsTer2)",
            c_notation="c.99999dup",
        )
        self.assertTrue(any("not found in Table 4 exon ranges" in warning for warning in result["warnings"]))


if __name__ == "__main__":
    unittest.main()
