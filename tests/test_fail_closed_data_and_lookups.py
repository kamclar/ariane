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

    def test_incomplete_coding_snv_snapshot_stops_startup(self):
        from backend.lookups import precomputed

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            index = root / "index.json"
            metadata = root / "metadata.json"
            index.write_text("{}", encoding="utf-8")
            metadata.write_text(
                json.dumps({"index_sha256": "not-the-file-checksum"}), encoding="utf-8"
            )
            with patch.object(precomputed, "CLASSIFICATION_SNAPSHOT_INDEX", index), patch.object(
                precomputed, "CLASSIFICATION_SNAPSHOT_METADATA", metadata
            ):
                with self.assertRaisesRegex(RuntimeError, "47,547 records"):
                    precomputed.validate_classification_snapshot()

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
    def test_population_service_attaches_founder_evidence_before_rules_run(self):
        from backend.population_frequency.models import GnomadRepository
        from backend.population_frequency.service import PopulationFrequencyService

        repository = GnomadRepository(
            variants={}, metadata={}, coverage_by_position={},
            frequency_path=None, frequency_status="missing",
            coverage_path=None, coverage_status="missing",
        )
        service = PopulationFrequencyService(
            repository,
            founder_lookup=lambda gene, c_notation: {
                "status": "reviewed_not_found",
                "is_pathogenic_founder": False,
                "reason": f"checked {gene}:{c_notation}",
            },
        )
        result = service.get_frequencies(
            gene="BRCA1",
            c_notation="c.509G>A",
            grch37={"chrom": "17", "pos": 1, "ref": "A", "alt": "G"},
        )
        self.assertFalse(result["founder_exception"]["is_pathogenic_founder"])
        self.assertEqual(
            result["founder_exception"]["reason"],
            "checked BRCA1:c.509G>A",
        )
        self.assertEqual(
            result["population_frequency_audit"]["founder_exception"],
            result["founder_exception"],
        )

    def test_population_service_reload_replaces_repository_for_bound_provider(self):
        from backend.population_frequency.models import GnomadRepository
        from backend.population_frequency.service import PopulationFrequencyService

        def repository(status):
            return GnomadRepository(
                variants={}, metadata={}, coverage_by_position={},
                frequency_path=None, frequency_status=status,
                coverage_path=None, coverage_status=status,
            )

        service = PopulationFrequencyService(
            repository("missing"),
            founder_lookup=lambda *_args: {
                "status": "reviewed_not_found",
                "is_pathogenic_founder": False,
                "reason": "test",
            },
        )
        provider = service.get_frequencies
        with patch(
            "backend.population_frequency.service.load_gnomad_repository",
            return_value=repository("approved_snapshot"),
        ):
            service.reload()
        result = provider(gene="BRCA1", c_notation="c.509G>A")
        self.assertEqual(service.repository.frequency_status, "approved_snapshot")
        self.assertEqual(result["cache_mode"], "approved_snapshot")

    def test_legacy_source_label_cannot_replace_extraction_provenance(self):
        from backend.population_frequency.lookup import dataset_extraction_ok
        from backend.population_frequency.models import GnomadRepository

        legacy_metadata = {
            "source": "gnomad_v2_1_1_exomes_grch37",
            "regions": {
                "GRCh37": {
                    "BRCA1": {"chrom": "17", "start": 1, "end": 100}
                }
            },
        }
        repository = GnomadRepository(
            variants={}, metadata=legacy_metadata, coverage_by_position={},
            frequency_path=None, frequency_status="approved_snapshot",
            coverage_path=None, coverage_status="approved_snapshot",
        )
        self.assertFalse(
            dataset_extraction_ok(
                repository,
                ["gnomad_v2_1_1_exomes_grch37"],
                {"chrom": "17", "pos": 50},
                "GRCh37",
            )
        )

    def test_founder_only_database_record_is_absent_for_outbred_presence(self):
        from backend.population_frequency.lookup import query_gnomad_dataset
        from backend.population_frequency.models import GnomadRepository
        from backend.population_frequency.policy import GNOMAD_LOCAL_DATASET_CONFIG

        variant_id = "13-100-A-G"
        record = {
            "variant_id": variant_id,
            "dataset": "gnomad_v2_1_1_exomes_grch37",
            "build": "GRCh37",
            "filter": "PASS",
            "faf95_max": 0.0,
            "faf95_pop": "afr",
            "faf95_scope": "non_cancer_non_founder_ancestries",
            "faf95_method": "official_gnomad_hail_table_non_cancer_faf95",
            "non_founder_observed": False,
            "non_founder_ac_by_ancestry": {
                "afr": 0, "amr": 0, "eas": 0, "nfe": 0, "sas": 0,
            },
            "excluded_population_context": {
                "asj": {
                    "label": "Ashkenazi Jewish",
                    "category": "founder_population",
                    "ac": 2,
                    "an": 1000,
                    "af": 0.002,
                    "faf95": 0.0004,
                    "used_for_ba1_bs1": False,
                    "used_for_pm2_presence": False,
                }
            },
        }
        dataset = "gnomad_v2_1_1_exomes_grch37"
        repository = GnomadRepository(
            variants={variant_id: (record,)},
            metadata={
                "regions": {
                    "GRCh37": {
                        "test": {"chrom": "13", "start": 1, "end": 200}
                    }
                },
                "region_padding_bp": 0,
                "extraction_log": [
                    {"dataset": dataset, "status": "ok", "chrom": "13"}
                ],
            },
            coverage_by_position={
                f"{dataset}|GRCh37|13|100": {
                    "mean_depth": 30.0,
                    "source": "test",
                }
            },
            frequency_path=None,
            frequency_status="approved_snapshot",
            coverage_path=None,
            coverage_status="approved_snapshot",
        )
        result = query_gnomad_dataset(
            repository,
            variant_id,
            {"chrom": "13", "pos": 100, "ref": "A", "alt": "G"},
            GNOMAD_LOCAL_DATASET_CONFIG["v2_1_non_cancer"],
            25.0,
            ["afr", "amr", "eas", "nfe", "sas"],
        )

        self.assertEqual(result["status"], "absent_in_non_founder_populations")
        self.assertFalse(result["found"])
        self.assertTrue(result["database_record_found"])
        self.assertEqual(
            result["exomes"]["excluded_population_context"]["asj"]["ac"], 2
        )

    def test_frequency_cache_requires_non_cancer_faf95_provenance(self):
        from backend.population_frequency.policy import (
            approved_manifest,
            approved_manifest_sha256,
            canonical_sha256,
        )
        from backend.population_frequency.snapshot_repository import validate_frequency_snapshot

        payload = {
            "metadata": {
                "schema_version": 2,
                "manifest_sha256": approved_manifest_sha256(),
                "automatic_release_activation": False,
                "classification_policies": approved_manifest()[
                    "classification_policies"
                ],
                "v2_faf95": {"raw_af_fallback_allowed": False},
                "v3_faf95": {"raw_af_fallback_allowed": False},
                "extraction_log": [
                    {
                        "dataset": dataset,
                        "status": "ok",
                        "source_identity": {"etag": "test", "x_goog_hash": "test"},
                    }
                    for dataset in (
                        "gnomad_v2_1_1_exomes_grch37",
                        "gnomad_v3_1_2_genomes_grch38",
                    )
                ],
            },
            "variants": {
                "17-1-A-G": [{
                    "dataset": "gnomad_v2_1_1_exomes_grch37",
                    "faf95_max": None,
                    "popmax_af": 0.001,
                }],
                "17-2-A-G": [{
                    "dataset": "gnomad_v3_1_2_genomes_grch38",
                    "faf95_max": None,
                    "popmax_af": 0.001,
                }],
            },
        }
        payload["metadata"]["records_sha256"] = canonical_sha256(payload["variants"])
        reason = validate_frequency_snapshot(payload)
        self.assertIn("lack ENIGMA-compatible non-cancer FAF95", reason)

    def test_coverage_cache_rejects_manifest_mismatch(self):
        from backend.population_frequency.snapshot_repository import validate_coverage_snapshot

        payload = {
            "metadata": {
                "schema_version": 2,
                "manifest_sha256": "wrong",
                "records": 1,
                "records_sha256": "wrong",
            },
            "coverage_by_position": {
                "gnomad_v2_1_1_exomes_grch37|GRCh37|17|1": {
                    "dataset_key": "gnomad_v2_1_1_exomes_grch37"
                }
            },
        }
        self.assertIn(
            "different panel/source manifest",
            validate_coverage_snapshot(payload),
        )

    def test_missing_approved_snapshot_has_no_alternate_selection(self):
        from backend.population_frequency.snapshot_repository import choose_frequency_snapshot

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.assertIsNone(choose_frequency_snapshot(root / "missing.json"))

    def test_untrusted_cache_cannot_produce_frequency_evidence(self):
        from backend.population_frequency.lookup import query_gnomad_dataset
        from backend.population_frequency.models import GnomadRepository
        from backend.population_frequency.policy import GNOMAD_LOCAL_DATASET_CONFIG

        repository = GnomadRepository(
            variants={"17-1-A-G": ({"dataset": "fixture"},)},
            metadata={},
            coverage_by_position={},
            frequency_path=None,
            frequency_status="unapproved",
            coverage_path=None,
            coverage_status="unapproved",
        )
        result = query_gnomad_dataset(
            repository,
            "17-1-A-G",
            {"chrom": "17", "pos": 1, "ref": "A", "alt": "G"},
            GNOMAD_LOCAL_DATASET_CONFIG["v2_1_non_cancer"],
            25.0,
            ["afr", "amr", "eas", "nfe", "sas"],
        )
        self.assertEqual(result["status"], "cache_untrusted")
        self.assertIn("not approved", result["errors"][0])


class LookupDiagnosticsTests(unittest.IsolatedAsyncioTestCase):
    def test_spliceai_deadline_finishes_before_nginx_timeout(self):
        from backend.lookup_execution import (
            EXTERNAL_LOOKUP_TIMEOUT,
            SERVICE_LOOKUP_TIMEOUTS,
        )

        self.assertEqual(EXTERNAL_LOOKUP_TIMEOUT, 12)
        self.assertEqual(SERVICE_LOOKUP_TIMEOUTS, {"SpliceAI": 30})
        self.assertLess(SERVICE_LOOKUP_TIMEOUTS["SpliceAI"], 60)

    async def test_spliceai_timeout_returns_unavailable_instead_of_hanging(self):
        from backend import lookup_execution

        diagnostics = []

        def slow_lookup():
            import time
            time.sleep(0.1)
            return 0.5

        with patch.dict(
            lookup_execution.SERVICE_LOOKUP_TIMEOUTS,
            {"SpliceAI": 0.01},
            clear=False,
        ):
            result = await lookup_or_unavailable(
                slow_lookup, None, "SpliceAI", diagnostics
            )
        self.assertIsNone(result)
        self.assertIn("timed out", diagnostics[0])

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

    def test_bayesdel_accepts_only_unambiguous_values(self):
        from backend.lookups import bayesdel

        self.assertEqual(
            bayesdel._select_unambiguous_bayesdel(0.31),
            (0.31, "single_value", None),
        )
        self.assertEqual(
            bayesdel._select_unambiguous_bayesdel([0.31, 0.31]),
            (0.31, "all_returned_values_identical", None),
        )
        score, basis, error = bayesdel._select_unambiguous_bayesdel(
            [0.31, 0.52]
        )
        self.assertIsNone(score)
        self.assertEqual(basis, "ambiguous_multiple_values")
        self.assertIn("without an explicit value-to-transcript mapping", error)

    def test_bayesdel_divergent_api_values_fail_closed(self):
        from backend.lookups import bayesdel

        response = MagicMock()
        response.__enter__.return_value = response
        response.read.return_value = json.dumps({
            "dbnsfp": {
                "bayesdel": {"no_af": {"score": [0.31, 0.52]}},
                "ensembl": {
                    "transcriptid": ["ENST00000357654", "ENST00000471181"]
                },
            }
        }).encode("utf-8")
        with patch("backend.lookups.bayesdel.urllib.request.urlopen", return_value=response):
            result = bayesdel.fetch_variant_data_myvariant(
                "BRCA1",
                "c.1A>G",
                {"chrom": "17", "pos": 1, "ref": "A", "alt": "G"},
            )

        self.assertIsNone(result["bayesdel"])
        self.assertEqual(result["status"], "ambiguous_transcript")
        self.assertEqual(result["reference_transcript"], "NM_007294.4")
        self.assertEqual(result["selection_basis"], "ambiguous_multiple_values")

    def test_bayesdel_does_not_cache_transient_no_score_response(self):
        from backend.lookups import bayesdel

        key = "BRCA1:c.5145C>A"
        bayesdel.BAYESDEL_CACHE.pop(key, None)
        bayesdel.BAYESDEL_STATUS_CACHE.pop(key, None)
        resolved = MagicMock()
        resolved.has_grch37.return_value = True
        resolved.grch37.chrom = "17"
        resolved.grch37.pos = 41215898
        resolved.grch37.ref = "G"
        resolved.grch37.alt = "T"
        no_score = {
            "bayesdel": None,
            "am_score": None,
            "am_class": None,
            "status": "no_score",
            "error": None,
        }
        scored = {
            "bayesdel": 0.333422,
            "am_score": None,
            "am_class": None,
            "status": "ok",
            "error": None,
            "selection_policy": bayesdel.BAYESDEL_SELECTION_POLICY,
            "selection_basis": "single_value",
            "reference_transcript": "NM_007294.4",
            "returned_transcripts": [],
        }

        with patch(
            "backend.lookups.coordinates.resolve_variant", return_value=resolved
        ), patch.object(
            bayesdel, "fetch_variant_data_myvariant", side_effect=[no_score, scored]
        ) as fetch, patch.object(bayesdel, "_save_cache") as save:
            first_score, _ = bayesdel.get_bayesdel_and_alphamissense(
                "BRCA1", "c.5145C>A"
            )
            self.assertIsNone(first_score)
            self.assertNotIn(key, bayesdel.BAYESDEL_CACHE)
            save.assert_not_called()

            second_score, _ = bayesdel.get_bayesdel_and_alphamissense(
                "BRCA1", "c.5145C>A"
            )

        self.assertEqual(second_score, 0.333422)
        self.assertEqual(fetch.call_count, 2)
        self.assertEqual(bayesdel.BAYESDEL_CACHE[key]["status"], "ok")
        save.assert_called_once()
        bayesdel.BAYESDEL_CACHE.pop(key, None)
        bayesdel.BAYESDEL_STATUS_CACHE.pop(key, None)

    def test_bayesdel_loader_ignores_persisted_empty_no_score_entries(self):
        from backend.lookups import bayesdel

        previous_cache = dict(bayesdel.BAYESDEL_CACHE)
        try:
            with tempfile.TemporaryDirectory() as directory:
                cache_path = Path(directory) / "bayesdel_cache.json"
                cache_path.write_text(
                    json.dumps({
                        "BRCA1:c.1A>G": {
                            "bayesdel": None,
                            "am_score": None,
                            "am_class": None,
                            "status": "no_score",
                        },
                        "BRCA1:c.2A>G": None,
                        "BRCA1:c.3A>G": {
                            "bayesdel": 0.42,
                            "am_score": None,
                            "am_class": None,
                            "status": "ok",
                            "selection_policy": bayesdel.BAYESDEL_SELECTION_POLICY,
                        },
                        "BRCA1:c.4A>G": {
                            "bayesdel": 0.55,
                            "am_score": None,
                            "am_class": None,
                            "status": "ok",
                        },
                    }),
                    encoding="utf-8",
                )
                bayesdel.BAYESDEL_CACHE.clear()
                with patch.object(bayesdel, "_CACHE_PATH", cache_path):
                    bayesdel._load_cache()

            self.assertNotIn("BRCA1:c.1A>G", bayesdel.BAYESDEL_CACHE)
            self.assertNotIn("BRCA1:c.2A>G", bayesdel.BAYESDEL_CACHE)
            self.assertEqual(
                bayesdel.BAYESDEL_CACHE["BRCA1:c.3A>G"]["bayesdel"], 0.42
            )
            self.assertNotIn("BRCA1:c.4A>G", bayesdel.BAYESDEL_CACHE)
        finally:
            bayesdel.BAYESDEL_CACHE.clear()
            bayesdel.BAYESDEL_CACHE.update(previous_cache)

    def test_bayesdel_cache_hit_does_not_call_api(self):
        from backend.lookups import bayesdel

        key = "BRCA1:c.999A>G"
        previous = bayesdel.BAYESDEL_CACHE.get(key)
        bayesdel.BAYESDEL_CACHE[key] = {
            "bayesdel": 0.41,
            "am_score": None,
            "am_class": None,
            "status": "ok",
            "reason": "Loaded from persistent runtime cache",
            "selection_policy": bayesdel.BAYESDEL_SELECTION_POLICY,
        }
        try:
            with patch.object(
                bayesdel,
                "fetch_variant_data_myvariant",
                side_effect=AssertionError("API must not be called on a cache hit"),
            ):
                score, _ = bayesdel.get_bayesdel_and_alphamissense(
                    "BRCA1", "c.999A>G"
                )
            self.assertEqual(score, 0.41)
        finally:
            if previous is None:
                bayesdel.BAYESDEL_CACHE.pop(key, None)
            else:
                bayesdel.BAYESDEL_CACHE[key] = previous


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

    def test_local_runtime_cache_is_outside_versioned_data_directories(self):
        from backend.runtime_cache import PROJECT_ROOT, choose_runtime_cache_dir

        with patch.dict("os.environ", {}, clear=True):
            self.assertEqual(
                choose_runtime_cache_dir(),
                PROJECT_ROOT / ".runtime-cache",
            )

    def test_mutable_cache_files_are_ignored_and_removed_from_data_directories(self):
        project_root = Path(__file__).resolve().parents[1]
        ignored = (project_root / ".gitignore").read_text(encoding="utf-8")
        for relative in (
            "backend/data/bayesdel_cache.json",
            "backend/data/coordinates_cache.json",
            "data/spliceai/spliceai_api_cache.json",
        ):
            self.assertFalse((project_root / relative).exists())
            self.assertIn(f"/{relative}", ignored)

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

    def test_precomputed_spliceai_variant_space_is_not_a_runtime_source(self):
        from backend.lookups import spliceai

        self.assertFalse(hasattr(spliceai, "SPLICEAI_PRECOMPUTED_CACHE"))
        self.assertFalse(hasattr(spliceai, "SPLICEAI_USE_PRECOMPUTED_CACHE"))
        self.assertFalse(hasattr(spliceai, "_load_precomputed_cache"))


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

    def test_missing_local_coordinates_fail_closed(self):
        from backend.lookups import coordinates

        key = "BRCA1:c.99998A>G"
        coordinates._RESOLVER_CACHE.pop(key, None)
        result = coordinates.resolve_variant("BRCA1", "c.99998A>G")
        self.assertEqual(result.status, "failed")
        self.assertEqual(result.source, "validated_local_sources")
        self.assertIn("Coordinate-dependent evidence was not evaluated", result.warnings[0])
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
        from tests.dag_test_support import classify_with_dag

        result = classify_with_dag(
            gene="BRCA1",
            variant_type="frameshift",
            p_notation="p.(Gly9999ValfsTer2)",
            c_notation="c.99999dup",
        )
        self.assertTrue(any("not found in Table 4 exon ranges" in warning for warning in result["warnings"]))


if __name__ == "__main__":
    unittest.main()
