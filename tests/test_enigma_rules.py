import unittest
from copy import deepcopy
from unittest.mock import patch
from pathlib import Path

from backend.population_frequency.coverage import (
    aggregate_coverage as _aggregate_coverage_from_dataset_results,
    lookup_coverage_by_position as _lookup_coverage_by_position,
)
from backend.population_frequency.criteria import evaluate_frequency_criteria
from backend.population_frequency.models import GnomadRepository
from backend.population_frequency.policy import (
    GNOMAD_LOCAL_DATASET_CONFIG,
    classification_policy_for_gene,
)
from backend.population_frequency.service import PopulationFrequencyService
from backend.classification_dag.policy import classify_by_enigma_combination
from tests.dag_test_support import classify_with_dag as evaluate_variant
from backend.modules.evidence_interactions import clinical_functional_risk_interactions
from backend.modules.exon_cnv_evidence import lookup_exon_cnv_evidence
from backend.modules.table9 import table9_lookup_ps3_bs3
from backend.modules.table4 import (
    TABLE4_DATA,
    parse_pvs1_code_strength,
    table4_lookup_splice,
)
from backend.modules.bp7 import evaluate_bp7
from backend.modules.pp3_bp4 import evaluate_pp3_bp4
from backend.modules.pvs1 import evaluate_pvs1
from backend.modules.pvs1_rna import evaluate_pvs1_rna
from backend.modules.ps1 import (
    compute_approval_basis_checksum,
    evaluate_ps1,
    select_vua_spliceai_for_ps1,
    validate_ps1_reference_registry,
)
import backend.modules.ps1 as ps1_module
from backend.lookups.founder_variants import lookup_pathogenic_founder_variant
from backend.modules.ps1_splice_evidence import evaluate_defined_splice_sources
from backend.modules.spliceai_policy import compare_table9_spliceai
from backend.modules.utils import is_in_functional_domain
from backend.modules.variant_type import infer_variant_type
from backend.modules.hgvs import split_combined_hgvs
from backend.modules.vus_explanation import explain_vus
from backend.modules.narrative import generate_narrative
try:
    from backend.models import VariantRequest
except ImportError:
    VariantRequest = None
try:
    import hgvs  # noqa: F401
    HGVS_RUNTIME_AVAILABLE = True
except ImportError:
    HGVS_RUNTIME_AVAILABLE = False


FREQUENCY_SERVICE = PopulationFrequencyService.load_default()


def get_gnomad_frequencies(gene, grch37=None, grch38=None, c_notation=""):
    return FREQUENCY_SERVICE.get_frequencies(
        gene=gene,
        c_notation=c_notation,
        grch37=grch37,
        grch38=grch38,
    )


def gnomad_data(
    *,
    max_af=None,
    found=False,
    v2_status="absent",
    v3_status="absent",
    v2_depth=30.0,
    v3_depth=30.0,
    pm2_absence_established=False,
    frequency_metric="faf95",
    gene="BRCA1",
    non_founder_allele_count=2,
):
    def dataset(status, depth, dataset_max_af=None):
        return {
            "status": status,
            "max_af": dataset_max_af,
            "coverage": {
                "mean_depth": depth,
                "classification_compatible": True,
            },
            "quality_filter_passed": True if status == "found" else None,
            "non_founder_allele_count": (
                non_founder_allele_count if status == "found" else 0
            ),
        }

    policy = classification_policy_for_gene(gene)
    return {
        "policy_id": policy["policy_id"],
        "classification_policy": policy,
        "frequency_policy": policy["frequency_criteria"],
        "status": "found" if found else "absent_with_coverage",
        "found": found,
        "max_af": max_af,
        "frequency_metric": frequency_metric,
        "pm2_absence_established": pm2_absence_established,
        "pm2_coverage_method": {
            "status": "approved_test_fixture",
            "automatic_assignment_allowed": True,
            "reason": "explicit approved method fixture",
        },
        "founder_exception": {
            "status": "reviewed_not_found",
            "is_pathogenic_founder": False,
            "reason": "test fixture records an authoritative negative review",
            "snapshot_version": "test",
        },
        "datasets": {
            "v2_1_non_cancer": dataset(v2_status, v2_depth, max_af if v2_status == "found" else None),
            "v3_1_non_cancer": dataset(
                v3_status,
                v3_depth,
                max_af if v3_status == "found" else None,
            ),
        },
    }


class VariantTypeTests(unittest.TestCase):
    @unittest.skipIf(
        VariantRequest is None or not HGVS_RUNTIME_AVAILABLE,
        "pydantic/hgvs runtime is not installed",
    )
    def test_protein_notation_is_derived_from_reference_snapshot(self):
        request = VariantRequest(gene="BRCA1", c_notation="c.4185G>A")
        self.assertEqual(request.p_notation, "p.(Gln1395=)")

        request = VariantRequest(
            gene="BRCA1", c_notation="c.4185G>A", p_notation=""
        )
        self.assertEqual(request.p_notation, "p.(Gln1395=)")

    @unittest.skipIf(
        VariantRequest is None or not HGVS_RUNTIME_AVAILABLE,
        "pydantic/hgvs runtime is not installed",
    )
    def test_tutorial_transcript_prefix_is_validated_and_removed(self):
        request = VariantRequest(
            gene="BRCA1",
            c_notation="NM_007294.4:c.509G>A",
            p_notation="p.Arg170Gln",
        )
        self.assertEqual(request.c_notation, "c.509G>A")
        self.assertEqual(request.p_notation, "p.(Arg170Gln)")
        with self.assertRaises(ValueError):
            VariantRequest(
                gene="BRCA1",
                c_notation="NM_000059.4:c.509G>A",
                p_notation="p.Arg170Gln",
            )

    @unittest.skipIf(
        VariantRequest is None or not HGVS_RUNTIME_AVAILABLE,
        "pydantic/hgvs runtime is not installed",
    )
    def test_common_c_hgvs_copy_paste_prefixes_are_normalized_safely(self):
        accepted = (
            ":c.5551_5552insT",
            ":  c.5551_5552insT",
            "BRCA1:c.5551_5552insT",
            "brca1   c.5551_5552insT",
            "BRCA1 NM_007294.4:c.5551_5552insT",
        )
        for notation in accepted:
            with self.subTest(notation=notation):
                request = VariantRequest(gene="BRCA1", c_notation=notation)
                self.assertEqual(request.c_notation, "c.5551_5552insT")
                self.assertEqual(request.p_notation, "p.(Asp1851ValfsTer29)")

        explicit_gene = VariantRequest(
            gene="BRCA1",
            c_notation="BRCA2:c.5551_5552insT",
        )
        self.assertEqual(explicit_gene.gene, "BRCA2")
        self.assertEqual(explicit_gene.reference_transcript, "NM_000059.4")
        with self.assertRaisesRegex(ValueError, "No reference transcript is configured"):
            VariantRequest(gene="BRCA1", c_notation="note:c.5551_5552insT")

    @unittest.skipIf(
        VariantRequest is None or not HGVS_RUNTIME_AVAILABLE,
        "pydantic/hgvs runtime is not installed",
    )
    def test_nonsense_protein_notation_without_parentheses_is_normalized(self):
        request = VariantRequest(
            gene="BRCA1",
            c_notation="c.5542C>T",
            p_notation="p.Gln1848Ter",
        )
        self.assertEqual(request.p_notation, "p.(Gln1848Ter)")

    def test_inframe_variants_use_normalized_types(self):
        self.assertEqual(infer_variant_type("c.123_125del", "p.(Val41del)"), "inframe_deletion")
        self.assertEqual(infer_variant_type("c.123_125insAAA", "p.(Val41_Gly42insLys)"), "inframe_insertion")
        self.assertEqual(infer_variant_type("c.123_125delinsAAA", "p.(Val41delinsLys)"), "inframe_delins")
        self.assertEqual(infer_variant_type("c.123dup", "p.(Val41dup)"), "inframe_insertion")

    def test_initiation_codon_is_not_treated_as_missense(self):
        self.assertEqual(infer_variant_type("c.1A>G", "p.(Met1Val)"), "initiation_codon")

    def test_methionine_at_position_1xxx_is_not_initiation_codon(self):
        # Met1083, Met1121 etc. must NOT be treated as initiation codon
        self.assertEqual(infer_variant_type("c.3247A>C", "p.(Met1083Leu)"), "missense")
        self.assertEqual(infer_variant_type("c.1121A>C", "p.(Met1121Leu)"), "missense")

    def test_met1_is_initiation_codon(self):
        self.assertEqual(infer_variant_type("c.1A>G", "p.(Met1Val)"), "initiation_codon")
        self.assertEqual(infer_variant_type("c.3A>G", "p.(Met1Ala)"), "initiation_codon")

    def test_5utr_not_classified_as_intronic(self):
        self.assertEqual(infer_variant_type("c.-10A>G", ""), "5utr")
        self.assertEqual(infer_variant_type("c.-1A>G", ""), "5utr")

    def test_3utr_not_classified_as_intronic(self):
        self.assertEqual(infer_variant_type("c.*10A>G", ""), "3utr")

    def test_combined_batch_hgvs_input_is_split_and_normalized(self):
        c_notation, p_notation = split_combined_hgvs(
            "c.6147_6149del (p.Val2050del)"
        )
        self.assertEqual(c_notation, "c.6147_6149del")
        self.assertEqual(p_notation, "p.(Val2050del)")

    def test_combined_nonsense_hgvs_without_protein_parentheses_is_normalized(self):
        c_notation, p_notation = split_combined_hgvs(
            "c.5542C>T p.Gln1848Ter"
        )
        self.assertEqual(c_notation, "c.5542C>T")
        self.assertEqual(p_notation, "p.(Gln1848Ter)")

    def test_exon_cnv_boundaries_are_not_misclassified_as_splice_sites(self):
        self.assertEqual(
            infer_variant_type("c.(80+1_81-1)_(134+1_135-1)dup", ""),
            "exon_duplication",
        )
        self.assertEqual(
            infer_variant_type("c.(80+1_81-1)_(134+1_135-1)del", ""),
            "exon_deletion",
        )


class CoordinateResolverTests(unittest.TestCase):
    def test_coding_snv_uses_validated_local_coordinates(self):
        from backend.lookups import coordinates

        key = "BRCA1:c.5366C>T"
        coordinates._RESOLVER_CACHE.pop(key, None)
        result = coordinates.resolve_variant("BRCA1", "c.5366C>T")

        self.assertEqual(result.status, "ok")
        self.assertEqual(result.source, "precomputed_snapshot")
        self.assertEqual(result.grch37.variant_id(), "17-41201178-G-A")
        self.assertEqual(result.grch38.variant_id(), "17-43049161-G-A")

    def test_runtime_coordinate_policy_disables_network_resolution(self):
        from backend.lookups import coordinates

        manifest = coordinates.validate_coordinate_source_manifest()
        self.assertFalse(manifest["runtime_policy"]["network_resolution_allowed"])
        source = Path(coordinates.__file__).read_text(encoding="utf-8")
        self.assertNotIn("urllib.request", source)
        self.assertNotIn("VariantValidator", source)
        self.assertNotIn("Mutalyzer", source)

    def test_reviewed_intronic_variants_use_versioned_coordinate_cache(self):
        from backend.lookups import coordinates

        expected = {
            "c.548-9A>G": ("17-41249315-T-C", "17-43097298-T-C"),
            "c.4987-6T>G": ("17-41219718-A-C", "17-43067701-A-C"),
        }
        for c_notation, (grch37, grch38) in expected.items():
            with self.subTest(c_notation=c_notation):
                result = coordinates.resolve_variant("BRCA1", c_notation)
                self.assertEqual(result.status, "ok")
                self.assertEqual(result.source, "versioned_intronic_coordinate_map")
                self.assertEqual(result.grch37.variant_id(), grch37)
                self.assertEqual(result.grch38.variant_id(), grch38)


class FrequencyTests(unittest.TestCase):
    def test_unknown_gene_has_no_inherited_brca_frequency_policy(self):
        data = get_gnomad_frequencies("TP53")
        self.assertEqual(data["status"], "policy_unavailable")
        self.assertIsNone(data["policy_id"])

        result = evaluate_frequency_criteria(
            {"status": "found", "max_af": 0.01, "policy_id": None},
            "missense",
            gene="TP53",
        )
        self.assertFalse(result["_gnomad_info"]["applies"])
        self.assertIn("gene-specific", result["_gnomad_info"]["reason"])

    def test_frequency_thresholds_are_read_from_gene_policy(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        policy = deepcopy(classification_policy_for_gene("BRCA1"))
        policy["frequency_criteria"]["ba1"]["threshold"] = 0.003
        data["classification_policy"] = policy
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BA1", result)
        self.assertIn("BS1_Strong", result)

    def test_tutorial_snvs_keep_coverage_auditable_but_do_not_establish_pm2(self):
        from backend.lookups.coordinates import resolve_variant

        for c_notation in ("c.4185G>A", "c.5217T>A"):
            with self.subTest(c_notation=c_notation):
                coords = resolve_variant("BRCA1", c_notation)
                self.assertEqual(coords.status, "ok")
                result = get_gnomad_frequencies(
                    "BRCA1", coords.grch37, coords.grch38
                )
                self.assertEqual(result["datasets"]["v2_1_non_cancer"]["status"], "absent")
                self.assertEqual(result["datasets"]["v3_1_non_cancer"]["status"], "absent")
                self.assertFalse(result["pm2_coverage_ok"])
                self.assertFalse(result["pm2_absence_established"])
                self.assertEqual(
                    result["status"], "pm2_coverage_method_unresolved"
                )
                self.assertFalse(
                    result["pm2_coverage_method"]["automatic_assignment_allowed"]
                )
                for coverage in result["coverage"]["datasets"].values():
                    self.assertEqual(coverage["coverage_scope"], "variant_reference_span")
                    self.assertEqual(coverage["positions_expected"], 1)
                    self.assertEqual(coverage["positions_available"], 1)
                self.assertTrue(
                    result["coverage"]["datasets"]["v2_1_non_cancer"]
                    ["classification_compatible"]
                )
                self.assertFalse(
                    result["coverage"]["datasets"]["v3_1_non_cancer"]
                    ["classification_compatible"]
                )
                self.assertEqual(
                    result["coverage"]["datasets"]["v3_1_non_cancer"]
                    ["compatibility_status"],
                    "unresolved_release_mismatch",
                )

    def test_pm2_configuration_requires_both_gnomad_versions(self):
        policy = classification_policy_for_gene("BRCA1")
        self.assertEqual(
            policy["frequency_criteria"]["pm2"]
            ["required_absence_dataset_runtime_keys"],
            ["v2_1_non_cancer", "v3_1_non_cancer"],
        )

        coverage = _aggregate_coverage_from_dataset_results(
            {
                "v2_1_non_cancer": {"coverage": {"mean_depth": 30.0}},
                "v3_1_non_cancer": {"coverage": {"mean_depth": 24.0}},
            },
            ["v2_1_non_cancer", "v3_1_non_cancer"],
            25.0,
        )
        self.assertFalse(coverage["passes_pm2"])

    def test_coverage_is_averaged_across_complete_reference_span(self):
        dataset = "gnomad_v2_1_1_exomes_grch37"
        coverage = {
            f"{dataset}|GRCh37|17|100": {"mean_depth": 24.0, "source": "test"},
            f"{dataset}|GRCh37|17|101": {"mean_depth": 28.0, "source": "test"},
        }
        repository = GnomadRepository(
            variants={}, metadata={}, coverage_by_position=coverage,
            frequency_path=None, frequency_status="test",
            coverage_path=None, coverage_status="test",
        )
        result = _lookup_coverage_by_position(
            repository,
            {"chrom": "17", "pos": 100, "ref": "AC", "alt": "GT"},
            dataset,
            "GRCh37",
            25.0,
            classification_compatible=True,
            compatibility_status="approved_test_fixture",
            compatibility_reason="test",
        )
        self.assertEqual(result["coverage_scope"], "variant_reference_span")
        self.assertEqual(result["positions_expected"], 2)
        self.assertEqual(result["positions_available"], 2)
        self.assertEqual(result["mean_depth"], 26.0)
        self.assertTrue(result["passes"])

    def test_incomplete_reference_span_coverage_fails_closed(self):
        dataset = "gnomad_v2_1_1_exomes_grch37"
        coverage = {
            f"{dataset}|GRCh37|17|100": {"mean_depth": 100.0, "source": "test"},
        }
        repository = GnomadRepository(
            variants={}, metadata={}, coverage_by_position=coverage,
            frequency_path=None, frequency_status="test",
            coverage_path=None, coverage_status="test",
        )
        result = _lookup_coverage_by_position(
            repository,
            {"chrom": "17", "pos": 100, "ref": "AC", "alt": "GT"},
            dataset,
            "GRCh37",
            25.0,
            classification_compatible=True,
            compatibility_status="approved_test_fixture",
            compatibility_reason="test",
        )
        self.assertIsNone(result["mean_depth"])
        self.assertFalse(result["passes"])
        self.assertEqual(result["missing_positions"], [101])

    def test_pm2_requires_absence_in_both_gnomad_versions(self):
        data = gnomad_data(pm2_absence_established=False)
        self.assertNotIn(
            "PM2_Supporting",
            evaluate_frequency_criteria(data, "missense", gene="BRCA1"),
        )

        data["pm2_absence_established"] = True
        self.assertIn(
            "PM2_Supporting",
            evaluate_frequency_criteria(data, "missense", gene="BRCA1"),
        )

    def test_pm2_is_not_applicable_for_small_indels(self):
        data = gnomad_data(pm2_absence_established=True)
        appendix_g = lookup_exon_cnv_evidence("BRCA1", "c.3891_3893del")
        self.assertFalse(
            evaluate_frequency_criteria(
                data,
                "inframe_deletion",
                gene="BRCA1",
                c_notation="c.3891_3893del",
                appendix_g_evidence=appendix_g,
            )["PM2"]["applies"]
        )

    def test_pm2_is_not_used_for_ptc_producing_indel(self):
        data = gnomad_data(pm2_absence_established=True)
        result = evaluate_frequency_criteria(
            data,
            "nonsense",
            gene="BRCA1",
            c_notation="c.5533_5534insG",
            appendix_g_evidence=lookup_exon_cnv_evidence(
                "BRCA1", "c.5533_5534insG"
            ),
        )
        self.assertFalse(result["PM2"]["applies"])
        self.assertNotIn("PM2_Supporting", result)
        self.assertIn("1 bp", result["PM2"]["reason"])
        self.assertIn(">50 bp", result["PM2"]["reason"])

    def test_large_indel_never_uses_small_variant_gnomad_pm2(self):
        data = gnomad_data(pm2_absence_established=True)
        appendix_g = lookup_exon_cnv_evidence("BRCA1", "c.100_150del")

        result = evaluate_frequency_criteria(
            data,
            "deletion",
            gene="BRCA1",
            c_notation="c.100_150del",
            appendix_g_evidence=appendix_g,
        )

        self.assertNotIn("PM2_Supporting", result)
        self.assertFalse(result["PM2"]["applies"])
        self.assertIn("unavailable", result["PM2"]["reason"].lower())

    def test_ba1_requires_depth_20(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found", v2_depth=19.0)
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BA1", result)
        self.assertIn("_gnomad_info", result)

        data["datasets"]["v2_1_non_cancer"]["coverage"]["mean_depth"] = 20.0
        self.assertIn(
            "BA1",
            evaluate_frequency_criteria(
                data, "missense", gene="BRCA1", c_notation="c.509G>A"
            ),
        )

    def test_ba1_bs1_are_suppressed_for_pathogenic_founder_variant(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        data["founder_exception"] = lookup_pathogenic_founder_variant(
            "BRCA1", "c.181T>G"
        )
        self.assertEqual(
            data["founder_exception"]["status"], "pathogenic_founder"
        )
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.181T>G"
        )
        self.assertNotIn("BA1", result)
        self.assertNotIn("BS1_Strong", result)
        self.assertIn("pathogenic founder", result["_gnomad_info"]["reason"])
        self.assertEqual(result["_excluded_criteria"]["BA1"]["points"], 0)
        self.assertFalse(result["_excluded_criteria"]["BA1"]["applies"])

    def test_bs1_supporting_threshold_is_reported_but_not_scored_for_c181t_g(self):
        data = gnomad_data(
            max_af=0.0000204,
            found=True,
            v3_status="found",
            non_founder_allele_count=4,
        )
        data["founder_exception"] = lookup_pathogenic_founder_variant(
            "BRCA1", "c.181T>G"
        )
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.181T>G"
        )
        self.assertNotIn("BS1_Supporting", result)
        excluded = result["_excluded_criteria"]["BS1_Supporting"]
        self.assertFalse(excluded["applies"])
        self.assertEqual(excluded["strength"], "Supporting")
        self.assertEqual(excluded["points"], 0)
        self.assertIn("0.00204%", excluded["reason"])
        self.assertIn("pathogenic founder", excluded["reason"])

    def test_ba1_reason_documents_excluded_founder_populations(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertIn("BA1", result)
        self.assertIn("founder and other non-scoring", result["BA1"]["reason"])
        self.assertIn("ENIGMA BRCA1/2 VCEP v1.2 Appendix G", result["BA1"]["reason"])

    def test_pm2_explains_founder_only_observation_policy(self):
        data = gnomad_data(max_af=0.0, pm2_absence_established=True)
        data["founder_context_only_observed"] = True
        result = evaluate_frequency_criteria(data, "missense", gene="BRCA1")
        self.assertIn("PM2_Supporting", result)
        self.assertIn(
            "observations confined to excluded founder/non-scoring populations",
            result["PM2_Supporting"]["reason"],
        )

    def test_population_frequency_founder_regression_matrix(self):
        single_outbred = gnomad_data(
            max_af=0.00001,
            found=True,
            v2_status="found",
            pm2_absence_established=False,
            non_founder_allele_count=1,
        )
        single_result = evaluate_frequency_criteria(
            single_outbred,
            "missense",
            gene="BRCA1",
            c_notation="c.509G>A",
        )
        self.assertNotIn("PM2_Supporting", single_result)
        self.assertFalse(any(code.startswith(("BA1", "BS1")) for code in single_result))

        founder_only = gnomad_data(
            max_af=0.0,
            found=False,
            pm2_absence_established=True,
        )
        founder_only["founder_context_only_observed"] = True
        founder_only_result = evaluate_frequency_criteria(
            founder_only,
            "missense",
            gene="BRCA1",
            c_notation="c.509G>A",
        )
        self.assertIn("PM2_Supporting", founder_only_result)
        self.assertIn(
            "excluded founder/non-scoring populations",
            founder_only_result["PM2_Supporting"]["reason"],
        )

        pathogenic_founder = gnomad_data(
            max_af=0.002,
            found=True,
            v2_status="found",
        )
        pathogenic_founder["founder_exception"] = lookup_pathogenic_founder_variant(
            "BRCA1", "c.181T>G"
        )
        pathogenic_founder_result = evaluate_frequency_criteria(
            pathogenic_founder,
            "missense",
            gene="BRCA1",
            c_notation="c.181T>G",
        )
        self.assertNotIn("BA1", pathogenic_founder_result)
        self.assertFalse(
            any(code.startswith("BS1") for code in pathogenic_founder_result)
        )
        self.assertIn(
            "pathogenic founder",
            pathogenic_founder_result["_gnomad_info"]["reason"],
        )

    def test_single_outbred_observation_cannot_create_bs1(self):
        data = gnomad_data(
            max_af=0.00006,
            found=True,
            v2_status="found",
            non_founder_allele_count=1,
        )
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BS1_Supporting", result)
        self.assertIn("single observation", result["_gnomad_info"]["reason"])

    def test_real_c5266dup_faf95_does_not_create_bs1(self):
        from backend.lookups.coordinates import resolve_variant

        coords = resolve_variant("BRCA1", "c.5266dup")
        data = get_gnomad_frequencies(
            "BRCA1", coords.grch37, coords.grch38, c_notation="c.5266dup"
        )
        self.assertGreater(data["max_af"], 0.00002)
        result = evaluate_frequency_criteria(
            data, "frameshift", gene="BRCA1", c_notation="c.5266dup"
        )
        self.assertNotIn("BS1_Supporting", result)
        self.assertIn("Ashkenazi Jewish founder", result["_gnomad_info"]["reason"])

    def test_ba1_bs1_founder_check_fails_closed_if_snapshot_is_unavailable(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        data["founder_exception"] = {
            "status": "unavailable",
            "is_pathogenic_founder": None,
            "reason": "test snapshot is invalid",
        }
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BA1", result)
        self.assertIn("authoritative review", result["_gnomad_info"]["reason"])

    def test_ba1_bs1_fail_closed_when_provider_omits_founder_evidence(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        data.pop("founder_exception")
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BA1", result)
        self.assertIn("does not contain a founder-exception result", result["_gnomad_info"]["reason"])

    def test_non_exhaustive_founder_snapshot_absence_is_unresolved(self):
        lookup = lookup_pathogenic_founder_variant("BRCA1", "c.509G>A")
        self.assertEqual(lookup["status"], "unresolved")
        self.assertIsNone(lookup["is_pathogenic_founder"])
        self.assertTrue(lookup["review_required"])

        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        data["founder_exception"] = lookup
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BA1", result)
        self.assertIn("authoritative review", result["_gnomad_info"]["reason"])

    def test_legacy_not_found_founder_status_does_not_authorize_ba1(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        data["founder_exception"] = {
            "status": "not_found",
            "is_pathogenic_founder": False,
            "reason": "legacy non-exhaustive lookup result",
        }
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BA1", result)
        self.assertIn("authoritative review", result["_gnomad_info"]["reason"])

    def test_incompatible_coverage_is_measured_but_not_classification_eligible(self):
        dataset = "gnomad_v3_1_2_genomes_grch38"
        repository = GnomadRepository(
            variants={},
            metadata={},
            coverage_by_position={
                f"{dataset}|GRCh38|17|100": {
                    "mean_depth": 35.0,
                    "source": "test r3.0.1 coverage",
                }
            },
            frequency_path=None,
            frequency_status="test",
            coverage_path=None,
            coverage_status="test",
        )
        result = _lookup_coverage_by_position(
            repository,
            {"chrom": "17", "pos": 100, "ref": "A", "alt": "G"},
            dataset,
            "GRCh38",
            25.0,
            classification_compatible=False,
            compatibility_status="unresolved_release_mismatch",
            compatibility_reason="test mismatch",
        )
        self.assertEqual(result["mean_depth"], 35.0)
        self.assertTrue(result["measurement_passes_threshold"])
        self.assertFalse(result["passes"])
        self.assertFalse(result["classification_compatible"])

    def test_ba1_bs1_require_passing_gnomad_record_qc(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found")
        data["datasets"]["v2_1_non_cancer"]["quality_filter_passed"] = False
        result = evaluate_frequency_criteria(
            data, "missense", gene="BRCA1", c_notation="c.509G>A"
        )
        self.assertNotIn("BA1", result)
        self.assertIn("did not pass dataset QC", result["_gnomad_info"]["reason"])

    def test_ba1_bs1_never_use_popmax_or_raw_af_fallback(self):
        for metric in ("popmax_af", "raw_af", "faf"):
            with self.subTest(metric=metric):
                data = gnomad_data(
                    max_af=0.002,
                    found=True,
                    v2_status="found",
                    frequency_metric=metric,
                )
                result = evaluate_frequency_criteria(
                    data, "missense", gene="BRCA1"
                )
                self.assertNotIn("BA1", result)
                self.assertNotIn("BS1_Strong", result)
                self.assertNotIn("BS1_Supporting", result)
                self.assertIn("FAF95", result["_gnomad_info"]["reason"])
                self.assertIn("cannot be used as a fallback", result["_gnomad_info"]["reason"])

    def test_found_variant_without_faf95_is_reported(self):
        data = gnomad_data(
            max_af=None,
            found=True,
            v2_status="found",
            frequency_metric=None,
        )
        result = evaluate_frequency_criteria(data, "frameshift", gene="BRCA1")
        self.assertFalse(result["PM2"]["applies"])
        self.assertIn("FAF95 is unavailable", result["_gnomad_info"]["reason"])

    def test_popmax_does_not_create_bs1_for_brca2_frameshift(self):
        from backend.lookups.coordinates import resolve_variant

        coords = resolve_variant("BRCA2", "c.9891_9894dup")
        data = get_gnomad_frequencies("BRCA2", coords.grch37, coords.grch38)
        self.assertEqual(data["max_af"], 0.0)
        self.assertEqual(data["frequency_metric"], "faf95")
        self.assertGreater(
            data["datasets"]["v3_1_non_cancer"]["genomes"]["popmax_af"],
            0.00002,
        )
        criteria = evaluate_frequency_criteria(
            data, "frameshift", gene="BRCA2"
        )
        self.assertNotIn("BS1_Supporting", criteria)

        result = evaluate_variant(
            gene="BRCA2",
            variant_type="frameshift",
            c_notation="c.9891_9894dup",
            p_notation="p.(Gln3299IlefsTer29)",
            gnomad_data=data,
        )
        self.assertEqual(result["criteria"]["PM5_PTC"]["strength"], "Strong")
        self.assertEqual(result["total_points"], 12)
        self.assertEqual(result["predicted_class"], 5)


class CriticalPtcBoundaryTests(unittest.TestCase):
    def test_brca1_tyr1845ter_insertion_gets_pvs1_pm5_but_not_pm2(self):
        c_notation = "c.5533_5534insG"
        p_notation = "p.(Tyr1845Ter)"
        result = evaluate_variant(
            gene="BRCA1",
            variant_type=infer_variant_type(c_notation, p_notation),
            p_notation=p_notation,
            c_notation=c_notation,
            gnomad_data=gnomad_data(pm2_absence_established=True),
        )
        self.assertEqual(result["criteria"]["PVS1"]["strength"], "Very Strong")
        self.assertEqual(result["criteria"]["PM5_PTC"]["strength"], "Strong")
        self.assertNotIn("PM2_Supporting", result["criteria"])
        self.assertEqual(result["total_points"], 12)
        self.assertEqual(result["predicted_class"], 5)

    def test_brca1_gln1848ter_gets_pvs1_and_pm5_strong(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="nonsense",
            p_notation="p.(Gln1848Ter)",
            c_notation="c.5542C>T",
        )
        self.assertEqual(result["criteria"]["PVS1"]["strength"], "Very Strong")
        self.assertEqual(result["criteria"]["PM5_PTC"]["strength"], "Strong")
        self.assertEqual(result["criteria"]["PS3"]["strength"], "Strong")
        self.assertEqual(result["total_points"], 16)
        self.assertEqual(result["predicted_class"], 5)

    def test_brca1_ptc_at_boundary_keeps_pm5_strong(self):
        result = evaluate_pvs1(
            "BRCA1", "nonsense", "p.(Leu1854Ter)", "c.5560C>T"
        )
        self.assertTrue(result["applies"])
        self.assertEqual(result["pm5_code"], "PM5_Strong (PTC)")
        self.assertEqual(result["pm5_points"], 4)

    def test_brca1_ptc_after_boundary_gets_neither_code(self):
        result = evaluate_pvs1(
            "BRCA1", "nonsense", "p.(Ile1855Ter)", "c.5563C>T"
        )
        self.assertFalse(result["applies"])
        self.assertIsNone(result["pm5_code"])
        self.assertEqual(result["pm5_points"], 0)


class SpliceTests(unittest.TestCase):
    def test_table9_spliceai_comparison_respects_enigma_band_boundaries(self):
        table9 = {"reviewed": True, "spliceai_prediction": 0.1}

        score, warnings = compare_table9_spliceai("BRCA1", 0.09, table9)
        self.assertEqual(score, 0.1)
        self.assertFalse(any("different ENIGMA prediction bands" in item for item in warnings))

        _, warnings = compare_table9_spliceai("BRCA1", 0.100001, table9)
        self.assertTrue(any("different ENIGMA prediction bands" in item for item in warnings))

        table9["spliceai_prediction"] = 0.2
        _, warnings = compare_table9_spliceai("BRCA1", 0.199999, table9)
        self.assertTrue(any("different ENIGMA prediction bands" in item for item in warnings))

    def test_table9_spliceai_comparison_never_supplies_missing_prediction(self):
        score, warnings = compare_table9_spliceai(
            "BRCA1", None, {"reviewed": True, "spliceai_prediction": 0.0}
        )
        self.assertEqual(score, 0.0)
        self.assertTrue(any("does not replace" in item for item in warnings))

    def test_unquantified_official_st2_patient_rna_requires_manual_review(self):
        tutorial_variant = evaluate_pvs1_rna("BRCA1", "c.4185G>A")
        another_st2_variant = evaluate_pvs1_rna("BRCA1", "c.80+5G>A")

        self.assertFalse(tutorial_variant["applies"])
        self.assertEqual(tutorial_variant["application_status"], "review_required")
        self.assertTrue(tutorial_variant["review_required"])
        self.assertEqual(tutorial_variant["points"], 0)
        self.assertEqual(tutorial_variant["table4_exon"], "E11(12)")
        self.assertEqual(
            tutorial_variant["manual_review_prefill"]["transcript_accession"],
            "NM_007294.4",
        )
        self.assertNotIn(
            "curated_strength", tutorial_variant["manual_review_prefill"]
        )
        self.assertFalse(another_st2_variant["applies"])
        self.assertEqual(another_st2_variant["application_status"], "review_required")

    def test_complex_st2_transcript_result_is_not_guessed(self):
        result = evaluate_pvs1_rna("BRCA1", "c.212+1G>T")
        self.assertIsNotNone(result.get("source_record"))
        self.assertFalse(result["applies"])
        self.assertIn("complex or partial", result["reason"])

    def test_c4185_keeps_prediction_and_routes_unquantified_rna_to_review(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="synonymous",
            p_notation="p.(Gln1395=)",
            c_notation="c.4185G>A",
            spliceai_score=0.95,
            gnomad_data=gnomad_data(pm2_absence_established=True),
            pp4_bp5_result={
                "applies": True,
                "code": "PP4",
                "strength": "Strong",
                "points": 4,
                "reason": "Combined clinical LR=328.184; PP4 Strong.",
            },
        )

        self.assertEqual(
            set(result["criteria"]),
            {"PP3", "PM2_Supporting", "PP4"},
        )
        self.assertEqual(result["criteria"]["PP3"]["strength"], "Supporting")
        self.assertEqual(result["total_points"], 6)
        self.assertEqual(result["predicted_class"], 4)
        self.assertTrue(result["rna_review"]["recommended"])
        self.assertEqual(result["rna_review"]["priority"], "high")
        self.assertEqual(
            result["rna_review"]["manual_review_prefill"]["transcript_accession"],
            "NM_007294.4",
        )

    def test_reviewed_intronic_variants_use_local_spliceai_and_apply_pp3(self):
        from backend.lookups.spliceai import get_spliceai_score

        expected_scores = {
            "c.548-9A>G": 0.86,
            # ENIGMA Table 9 and the Appendix J raw-delta profile both use
            # the maximum component. DS_DL=0.74 is larger than DS_AL=0.73.
            "c.4987-6T>G": 0.74,
        }
        for c_notation, expected_score in expected_scores.items():
            with self.subTest(c_notation=c_notation):
                score = get_spliceai_score("BRCA1", c_notation)
                self.assertEqual(score, expected_score)
                result = evaluate_variant(
                    gene="BRCA1",
                    variant_type=infer_variant_type(c_notation, ""),
                    p_notation="",
                    c_notation=c_notation,
                    spliceai_score=score,
                )
                self.assertEqual(result["criteria"]["PP3"]["strength"], "Supporting")
                self.assertEqual(result["criteria"]["PP3"]["points"], 1)

    def test_noncanonical_splice_prediction_does_not_create_pvs1(self):
        result = evaluate_pvs1(
            "BRCA1",
            "splice_site",
            "p.(?)",
            "c.100+3A>G",
            spliceai_score=0.9,
        )
        self.assertFalse(result["applies"])
        self.assertEqual(result["points"], 0)

    def test_intronic_splice_prediction_creates_pp3(self):
        result = evaluate_pp3_bp4("BRCA1", "intronic", "p.(?)", spliceai_score=0.9)
        self.assertEqual(result["PP3"]["points"], 1)

    def test_table4_rna_rule_is_outside_automated_scope(self):
        result = evaluate_pvs1("BRCA1", "splice_site", "p.(?)", "c.80+1G>A")
        self.assertTrue(result["requires_rna"])
        self.assertFalse(result["applies"])

    def test_intronic_bp7_applies_outside_conserved_motif(self):
        result = evaluate_bp7(
            "intronic",
            spliceai_score=0.05,
            bp4_met=True,
            c_notation="c.100+7A>G",
            gene="BRCA1",
        )
        self.assertTrue(result["applies"])
        self.assertEqual(result["points"], -1)

    def test_intronic_bp7_does_not_apply_inside_conserved_motif(self):
        result = evaluate_bp7(
            "intronic",
            spliceai_score=0.05,
            bp4_met=True,
            c_notation="c.100+6A>G",
            gene="BRCA1",
        )
        self.assertFalse(result["applies"])

    def test_ps1_requires_confirmed_low_spliceai(self):
        result = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.123A>G",
            p_notation="p.(His41Arg)",
            variant_type="missense",
            spliceai_score=None,
            vua_splice_evidence_status="none_identified",
            vua_splice_sources_checked=["ENIGMA Table 9", "ENIGMA ST2"],
        )
        self.assertFalse(result["applies"])
        self.assertTrue(result["review_required"])
        self.assertIn("SpliceAI is unavailable", result["reason"])

    def test_st7_reference_prefills_review_but_does_not_score_ps1(self):
        result = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.123A>G",
            p_notation="p.(His41Arg)",
            variant_type="missense",
            spliceai_score=0.01,
            vua_splice_evidence_status="none_identified",
            vua_splice_sources_checked=["ENIGMA Table 9", "ENIGMA ST2"],
            reference_spliceai_scores={"c.122A>G": 0.01},
        )
        self.assertFalse(result["applies"])
        self.assertEqual(result["points"], 0)
        self.assertTrue(result["review_required"])
        self.assertEqual(result["application_status"], "manual_review_required")
        self.assertEqual(result["candidates"][0]["c_notation"], "c.122A>G")
        self.assertEqual(result["candidates"][0]["reference_status"], "review_required")
        self.assertIn("no separate ENIGMA/ClinGen VCEP assertion", result["reason"])
        from backend.modules.protein_ps1_review import evaluate_protein_ps1_review

        display = evaluate_protein_ps1_review(result, gene="BRCA1")
        prefill = display["manual_review_prefill"]
        self.assertEqual(prefill["reference_variant"], "BRCA1 c.122A>G")
        self.assertEqual(prefill["reference_p_notation"], "p.(His41Arg)")
        self.assertEqual(prefill["reference_classification"], "Pathogenic")
        self.assertEqual(
            prefill["classification_verification"],
            "historical_classification_only",
        )
        self.assertTrue(prefill["same_missense_confirmed"])
        self.assertTrue(prefill["different_nucleotide_change_confirmed"])
        self.assertEqual(prefill["vua_spliceai_score"], 0.01)
        self.assertEqual(prefill["reference_spliceai_score"], 0.01)

    def test_ps1_fails_closed_when_reference_spliceai_is_unavailable(self):
        result = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.123A>G",
            p_notation="p.(His41Arg)",
            variant_type="missense",
            spliceai_score=0.01,
            vua_splice_evidence_status="none_identified",
            vua_splice_sources_checked=["ENIGMA Table 9", "ENIGMA ST2"],
            reference_spliceai_scores={"c.122A>G": None},
        )
        self.assertFalse(result["applies"])
        self.assertTrue(result["review_required"])
        self.assertIn("reference variant", result["reason"])

    def test_ps1_is_not_applied_when_reference_spliceai_exceeds_threshold(self):
        result = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.123A>G",
            p_notation="p.(His41Arg)",
            variant_type="missense",
            spliceai_score=0.01,
            vua_splice_evidence_status="none_identified",
            vua_splice_sources_checked=["ENIGMA Table 9", "ENIGMA ST2"],
            reference_spliceai_scores={"c.122A>G": 0.101},
        )
        self.assertFalse(result["applies"])
        self.assertFalse(result["review_required"])
        self.assertEqual(result["application_status"], "not_applicable")
        self.assertIn("0.101", result["reason"])

    def test_table9_score_does_not_replace_unavailable_ps1_prediction(self):
        table9 = table9_lookup_ps3_bs3("BRCA1", "c.131G>C")
        score, source = select_vua_spliceai_for_ps1(None, table9)
        self.assertIsNone(score)
        self.assertEqual(source, "configured SpliceAI source")
        splice_evidence = evaluate_defined_splice_sources(
            "BRCA1", "c.131G>C", table9
        )
        result = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.131G>C",
            p_notation="p.(Cys44Ser)",
            variant_type="missense",
            spliceai_score=score,
            vua_spliceai_source=source,
            vua_splice_evidence_status=splice_evidence["status"],
            vua_splice_sources_checked=splice_evidence["sources_checked"],
        )
        self.assertFalse(result["applies"])
        self.assertTrue(result["review_required"])
        self.assertIn("SpliceAI is unavailable", result["reason"])

    def test_st7_reference_with_splice_effect_is_excluded(self):
        result = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.139G>T",
            p_notation="p.(Cys47Tyr)",
            variant_type="missense",
            spliceai_score=0.01,
            vua_splice_evidence_status="none_identified",
            vua_splice_sources_checked=[
                "ENIGMA Specifications Table 9 v1.2",
                "ENIGMA Supplementary Table 2 v1.2",
            ],
        )
        self.assertFalse(result["applies"])
        self.assertFalse(result["review_required"])
        self.assertEqual(result["application_status"], "reference_ineligible")
        self.assertEqual(result["candidates"][0]["c_notation"], "c.140G>A")
        self.assertEqual(result["candidates"][0]["reference_status"], "excluded")
        from backend.modules.protein_ps1_review import evaluate_protein_ps1_review

        display = evaluate_protein_ps1_review(result)
        self.assertTrue(display["display"])
        self.assertFalse(display["recommended"])
        self.assertEqual(display["title"], "Protein PS1 not applicable")
        self.assertIn("different P/LP reference variant", display["summary"])
        self.assertIn("does not classify", display["reasons"][0])

    def test_approved_vcep_reference_can_automatically_score_ps1(self):
        approved = {
            "reference_id": "TEST_REF",
            "gene": "BRCA1",
            "transcript": "NM_007294.4",
            "c_notation": "c.122A>G",
            "p_notation": "p.(His41Arg)",
            "classification": "Pathogenic",
            "classification_verification": "external_vcep_assertion",
            "classification_source": "ClinGen ENIGMA expert panel test assertion",
            "candidate_source": "test",
        }
        with patch.object(ps1_module, "_LOADED", True), patch.object(
            ps1_module, "_ST7_LOOKUP", {}
        ), patch.object(
            ps1_module,
            "_APPROVED_LOOKUP",
            {"BRCA1": {"His41Arg": [approved]}},
        ):
            result = evaluate_ps1(
                gene="BRCA1",
                c_notation="c.123A>G",
                p_notation="p.(His41Arg)",
                variant_type="missense",
                spliceai_score=0.1,
                vua_splice_evidence_status="none_identified",
                vua_splice_sources_checked=["ENIGMA Table 9", "ENIGMA ST2"],
                reference_spliceai_scores={"c.122A>G": 0.01},
            )
        self.assertTrue(result["applies"])
        self.assertEqual(result["strength"], "Strong")
        self.assertEqual(result["points"], 4)
        self.assertEqual(result["application_status"], "auto_applied")

    def test_confirmed_or_predicted_splice_effect_blocks_protein_ps1(self):
        predicted = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.123A>G",
            p_notation="p.(His41Arg)",
            variant_type="missense",
            spliceai_score=0.101,
            vua_splice_evidence_status="none_identified",
            vua_splice_sources_checked=["ENIGMA Table 9", "ENIGMA ST2"],
        )
        confirmed = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.123A>G",
            p_notation="p.(His41Arg)",
            variant_type="missense",
            spliceai_score=0.01,
            vua_splice_evidence_status="abnormal",
            vua_splice_sources_checked=["ENIGMA Table 9", "ENIGMA ST2"],
        )
        self.assertFalse(predicted["applies"])
        self.assertFalse(predicted["review_required"])
        self.assertFalse(confirmed["applies"])
        self.assertFalse(confirmed["review_required"])

    def test_defined_splice_sources_detect_table9_splice_effect(self):
        table9 = table9_lookup_ps3_bs3("BRCA1", "c.181T>G")
        result = evaluate_defined_splice_sources("BRCA1", "c.181T>G", table9)
        self.assertEqual(result["status"], "abnormal")
        self.assertIn("ENIGMA Specifications Table 9 v1.2", result["sources_checked"])

    def test_defined_splice_sources_accept_published_no_aberration(self):
        table9 = table9_lookup_ps3_bs3("BRCA1", "c.110C>G")
        self.assertEqual(
            table9["predicted_or_observed_splicing"],
            "N, no aberration",
        )
        self.assertIn("no aberration", table9["splice_result_published"].lower())
        result = evaluate_defined_splice_sources("BRCA1", "c.110C>G", table9)
        self.assertEqual(result["status"], "normal")

    def test_ps1_registry_rejects_known_circular_dependency(self):
        def record(reference_id, c_notation, dependency):
            value = {
                "reference_id": reference_id,
                "gene": "BRCA1",
                "transcript": "NM_007294.4",
                "c_notation": c_notation,
                "p_notation": "p.(His41Arg)",
                "classification": "Pathogenic",
                "classification_verification": "external_vcep_assertion",
                "classification_source": "ClinGen ENIGMA test assertion",
                "classification_assertion": {
                    "organization": "ENIGMA BRCA1/2 VCEP",
                    "assertion_id": reference_id,
                    "ruleset_version": "1.2.0",
                    "accessed_at": "2026-08-12",
                },
                "status": "eligible",
                "status_reason": "Synthetic eligible test reference.",
                "protein_branch": "missense_runtime_spliceai_check_required",
                "protein_mechanism_evidence": {
                    "basis": "curated_protein_mechanism_assessment"
                },
                "reference_splice_evidence": {
                    "threshold": 0.1,
                    "prediction_policy": "runtime_required",
                    "confirmed_status": "none_identified",
                    "sources_checked": ["ENIGMA Table 9", "ENIGMA ST2"],
                    "checked_at": "2026-08-12",
                    "provenance": {"snapshot": "test"},
                },
                "classification_ps1_dependency": {
                    "used": True,
                    "reference_ids": [dependency],
                },
            }
            value["approval_basis_checksum"] = compute_approval_basis_checksum(value)
            return value

        registry = {
            "schema_version": 4,
            "registry_version": "test",
            "status": "active",
            "defined_splice_sources": ["ENIGMA Table 9", "ENIGMA ST2"],
            "reference_source_policy": {
                "accepted_classification_bases": [
                    {"id": "enigma_st7_v1_2_reference_set"},
                    {"id": "external_vcep_assertion"},
                    {"id": "locally_recurated_under_enigma_vcep"},
                ]
            },
            "source_checksums": {
                "st7_sha256": "0" * 64,
                "table9_sha256": "0" * 64,
                "st2_sha256": "0" * 64,
                "curated_extensions_sha256": "0" * 64,
            },
            "reference_count": 2,
            "status_counts": {"eligible": 2},
            "references": [
                record("A", "c.122A>G", "B"),
                record("B", "c.123A>G", "A"),
            ],
        }
        with self.assertRaisesRegex(RuntimeError, "circular dependency"):
            validate_ps1_reference_registry(registry)

    def test_splice_ps1_candidate_discovery_comes_directly_from_complete_st2(self):
        from backend.modules.ps1_splice_evidence import (
            list_splice_ps1_candidate_discovery,
        )

        data = list_splice_ps1_candidate_discovery()
        self.assertEqual(data["status"], "candidate_discovery_only")
        self.assertEqual(data["candidate_count"], 75)
        self.assertEqual(len(data["candidates"]), 75)
        self.assertTrue(all(
            item["eligibility_status"] == "candidate_discovery_only"
            for item in data["candidates"]
        ))
        self.assertTrue(all("suggested_strength" not in item for item in data["candidates"]))
        candidate = next(
            item
            for item in data["candidates"]
            if item["gene"] == "BRCA1" and item["reference_variant"] == "c.4185G>A"
        )
        self.assertEqual(candidate["classification"], "Pathogenic")
        self.assertTrue(candidate["reference_splice_event"])
        self.assertEqual(len(candidate["source_file_sha256"]), 64)


class ClassifierIntegrationTests(unittest.TestCase):
    def test_founder_frequency_exception_is_structured_and_not_scored(self):
        data = gnomad_data(
            max_af=0.0000204,
            found=True,
            v3_status="found",
            non_founder_allele_count=4,
        )
        data["founder_exception"] = lookup_pathogenic_founder_variant(
            "BRCA1", "c.181T>G"
        )
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="missense",
            p_notation="p.(Cys61Gly)",
            c_notation="c.181T>G",
            spliceai_score=0.0,
            gnomad_data=data,
        )
        self.assertNotIn("BS1_Supporting", result["criteria"])
        self.assertEqual(
            result["excluded_criteria"]["BS1_Supporting"]["points"], 0
        )

    def test_unapproved_protein_ps1_candidate_is_review_only(self):
        ps1_result = {
            "applies": False,
            "review_required": True,
            "application_status": "manual_review_required",
            "blocking_reasons": ["Reference VCEP verification is incomplete"],
            "vua_splice_sources_checked": ["ENIGMA Table 9", "ENIGMA ST2"],
            "vua_splice_evidence_status": "none_identified",
            "candidates": [{
                "key": "ST7|BRCA1|c.122A>G",
                "reference_id": "",
                "gene": "BRCA1",
                "transcript": "NM_007294.4",
                "c_notation": "c.122A>G",
                "p_notation": "p.(His41Arg)",
                "classification": "Pathogenic",
                "iarc_class": 5,
                "classification_basis": "enigma_multifactorial_likelihood_reference_set",
                "classification_source": "ST7 test source",
                "reference_status": "review_required",
                "status_reason": "VCEP verification absent",
                "source_dataset": "ENIGMA ST7",
            }],
        }
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="missense",
            p_notation="p.(His41Arg)",
            c_notation="c.123A>G",
            spliceai_score=0.01,
            ps1_result=ps1_result,
        )
        self.assertNotIn("PS1", result["criteria"])
        self.assertTrue(result["protein_ps1_review"]["recommended"])
        self.assertFalse(result["protein_ps1_review"]["is_evidence_criterion"])

    def test_frameshift_does_not_warn_about_inapplicable_bayesdel(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="frameshift",
            p_notation="p.(Gln1756ProfsTer74)",
            c_notation="c.5266dup",
            spliceai_score=0.0,
            bayesdel_score=None,
        )
        self.assertFalse(any("BayesDel" in warning for warning in result["warnings"]))

    def test_enigma_pathogenic_combinations_use_correct_strength_counts(self):
        very_strong_plus_one_moderate = {
            "PVS1": {"points": 8, "strength": "Very Strong"},
            "PM5": {"points": 2, "strength": "Moderate"},
        }
        self.assertEqual(
            classify_by_enigma_combination(very_strong_plus_one_moderate, 10)[0], 5
        )
        two_very_strong = {
            "PVS1": {"points": 8, "strength": "Very Strong"},
            "PP4": {"points": 8, "strength": "Very Strong"},
        }
        self.assertEqual(
            classify_by_enigma_combination(two_very_strong, 16)[:2],
            (5, "Pathogenic"),
        )
        two_strong = {
            "PS3": {"points": 4, "strength": "Strong"},
            "PS4": {"points": 4, "strength": "Strong"},
        }
        self.assertEqual(classify_by_enigma_combination(two_strong, 8)[0], 4)

    def test_table3_pathogenic_and_benign_edge_combinations(self):
        two_strong_plus_moderate = {
            "PS3": {"points": 4, "strength": "Strong"},
            "PS4": {"points": 4, "strength": "Strong"},
            "PM3": {"points": 2, "strength": "Moderate"},
        }
        self.assertEqual(
            classify_by_enigma_combination(two_strong_plus_moderate, 10)[:2],
            (5, "Pathogenic"),
        )

        benign_strong_plus_two_moderate = {
            "BS3": {"points": -4, "strength": "Strong"},
            "BS2": {"points": -2, "strength": "Moderate"},
            "BS4": {"points": -2, "strength": "Moderate"},
        }
        self.assertEqual(
            classify_by_enigma_combination(benign_strong_plus_two_moderate, -8)[:2],
            (1, "Benign"),
        )

        benign_strong_plus_supporting = {
            "BS3": {"points": -4, "strength": "Strong"},
            "BP7": {"points": -1, "strength": "Supporting"},
        }
        self.assertEqual(
            classify_by_enigma_combination(benign_strong_plus_supporting, -5)[:2],
            (2, "Likely Benign"),
        )

    def test_table9_functional_evidence_does_not_suppress_bp4(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="missense",
            p_notation="p.(Cys61Gly)",
            c_notation="c.181T>G",
            spliceai_score=0.05,
            bayesdel_score=0.10,
        )
        self.assertIn("PS3", result["criteria"])
        self.assertIn("BP4", result["criteria"])
        self.assertEqual(result["evidence_interactions"][0]["status"], "conflict")
        self.assertTrue(result["evidence_interactions"][0]["review_required"])

    def test_ps3_bp5_conflict_flags_possible_reduced_penetrance_without_asserting_it(self):
        interactions = clinical_functional_risk_interactions({
            "PS3": {"points": 4, "strength": "Strong"},
            "BP5": {"points": -4, "strength": "Strong"},
        })
        self.assertEqual(len(interactions), 1)
        warning = interactions[0]
        self.assertEqual(warning["status"], "conflict")
        self.assertTrue(warning["review_required"])
        self.assertIn("does not establish reduced penetrance", warning["reason"])
        self.assertEqual(warning["retained"], ["PS3", "BP5"])

    def test_reduced_penetrance_warning_requires_both_ps3_and_bp5(self):
        self.assertEqual(
            clinical_functional_risk_interactions({"PS3": {"points": 4}}),
            [],
        )
    def test_custom_donor_guard_is_not_part_of_active_scoring(self):
        project_root = Path(__file__).resolve().parents[1]
        self.assertFalse((project_root / "backend/modules/donor_guard.py").exists())
        self.assertNotIn("donor_guard", (project_root / "backend/main.py").read_text(encoding="utf-8"))

    def test_rna_dependent_pvs1_does_not_leak_into_score(self):
        unconfirmed = evaluate_variant(
            gene="BRCA1",
            variant_type="splice_site",
            p_notation="p.(?)",
            c_notation="c.212+1G>T",
        )
        self.assertNotIn("PVS1", unconfirmed["criteria"])
        self.assertEqual(unconfirmed["total_points"], 0)
        self.assertTrue(unconfirmed["rna_review"]["recommended"])
        self.assertEqual(unconfirmed["rna_review"]["priority"], "high")
        self.assertFalse(unconfirmed["rna_review"]["is_evidence_criterion"])
        self.assertIn("PVS1 (RNA)", unconfirmed["rna_review"]["potential_branches"])

    def test_splice_ps1_review_flags_predicted_splice_effect_without_scoring(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="missense",
            p_notation="p.(Cys64Gly)",
            c_notation="c.190T>G",
            spliceai_score=0.65,
        )
        self.assertTrue(result["splice_ps1_review"]["recommended"])
        self.assertFalse(result["splice_ps1_review"]["is_evidence_criterion"])
        self.assertIn("PS1 (splice)", result["splice_ps1_review"]["potential_branches"])
        self.assertNotIn("PS1", result["criteria"])

    def test_splice_ps1_review_does_not_flag_low_splice_intronic_variant(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="intronic",
            p_notation="p.(?)",
            c_notation="c.100+7A>G",
            spliceai_score=0.05,
        )
        self.assertFalse(result["splice_ps1_review"]["recommended"])

    def test_intronic_bp7_is_added_by_classifier(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="intronic",
            p_notation="p.(?)",
            c_notation="c.100+7A>G",
            spliceai_score=0.05,
        )
        self.assertEqual(result["criteria"]["BP4"]["points"], -1)
        self.assertEqual(result["criteria"]["BP7"]["points"], -1)

    def test_initiation_codon_is_flagged_without_automatic_points(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="initiation_codon",
            p_notation="p.(Met1Val)",
            c_notation="c.1A>G",
        )
        self.assertNotIn("PVS1", result["criteria"])
        self.assertTrue(any("Initiation codon" in warning for warning in result["warnings"]))
        self.assertTrue(result["initiation_review"]["recommended"])
        self.assertEqual(result["initiation_review"]["priority"], "high")
        self.assertIn("PVS1_INIT", result["initiation_review"]["potential_branches"])
        self.assertFalse(result["initiation_review"]["is_evidence_criterion"])

    def test_exon_duplication_uses_confirmed_tandem_input(self):
        notation = "c.(80+1_81-1)_(134+1_135-1)dup"
        variant_type = infer_variant_type(notation, "")
        unknown = evaluate_variant(
            gene="BRCA1",
            variant_type=variant_type,
            p_notation="p.(?)",
            c_notation=notation,
        )
        tandem = evaluate_variant(
            gene="BRCA1",
            variant_type=variant_type,
            p_notation="p.(?)",
            c_notation=notation,
            dup_type="Tandem",
        )
        self.assertEqual(variant_type, "exon_duplication")
        self.assertEqual(unknown["criteria"]["PVS1"]["strength"], "Moderate")
        self.assertEqual(tandem["criteria"]["PVS1"]["strength"], "Strong")

    def test_table9_splice_score_does_not_override_configured_score(self):
        table9_result = table9_lookup_ps3_bs3("BRCA1", "c.3891_3893del")
        self.assertEqual(table9_result["splice_result_published"], "no aberration (PMID: 18273839)")
        self.assertEqual(table9_result["spliceai_prediction"], 0)
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="inframe_deletion",
            p_notation="p.(Ser1298del)",
            c_notation="c.3891_3893del",
            spliceai_score=0.15,
            table9_result=table9_result,
        )
        self.assertEqual(result["criteria"]["BS3"]["points"], -4)
        self.assertNotIn("BP1", result["criteria"])
        self.assertEqual(result["total_points"], -4)
        self.assertEqual(result["predicted_class"], 3)
        self.assertEqual(result["predicted_label"], "VUS")
        self.assertTrue(any("Table 9=0.000" in warning for warning in result["warnings"]))
        self.assertTrue(any("different ENIGMA prediction bands" in warning for warning in result["warnings"]))

    def test_pvs1_very_strong_alone_remains_vus_with_explanation(self):
        result = evaluate_variant(
            gene="BRCA2",
            variant_type="nonsense",
            p_notation="p.(Cys161Ter)",
            c_notation="c.483T>A",
        )
        self.assertEqual(result["criteria"]["PVS1"]["strength"], "Very Strong")
        self.assertNotIn("PM5_PTC", result["criteria"])
        self.assertEqual(result["total_points"], 8)
        self.assertEqual(result["predicted_class"], 3)
        self.assertIn(
            "requires at least one additional Supporting",
            result["classification_note"],
        )

    def test_terminal_frameshift_uses_one_table4_row_for_pvs1_and_pm5(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="frameshift",
            p_notation="p.(Tyr1853AspfsTer25)",
            c_notation="c.5556_5560del",
        )

        self.assertEqual(result["criteria"]["PVS1"]["strength"], "Very Strong")
        self.assertEqual(result["criteria"]["PM5_PTC"]["strength"], "Strong")
        self.assertEqual(result["total_points"], 12)
        self.assertEqual(result["predicted_class"], 5)
        reason = result["criteria"]["PVS1"]["reason"]
        self.assertIn("First altered residue p.1853 <= p.1854", reason)
        self.assertIn("PVS1, PM5_Strong (PTC)", reason)
        self.assertIn("PM5 follows the selected Table 4 row", reason)

    def test_brca2_terminal_frameshift_uses_the_same_boundary_policy(self):
        result = evaluate_variant(
            gene="BRCA2",
            variant_type="frameshift",
            p_notation="p.(Glu3309AlafsTer20)",
            c_notation="c.9925del",
        )

        self.assertEqual(result["criteria"]["PVS1"]["strength"], "Very Strong")
        self.assertEqual(result["criteria"]["PM5_PTC"]["strength"], "Strong")
        self.assertEqual(result["total_points"], 12)
        self.assertEqual(result["predicted_class"], 5)
        reason = result["criteria"]["PVS1"]["reason"]
        self.assertIn("First altered residue p.3309 <= p.3309", reason)
        self.assertIn("PVS1, PM5_Strong (PTC)", reason)
        self.assertIn("PM5 follows the selected Table 4 row", reason)

    def test_cross_exon_frameshift_keeps_pm5_from_nucleotide_change_exon(self):
        result = evaluate_variant(
            gene="BRCA2",
            variant_type="frameshift",
            p_notation="p.(Leu167AlafsTer20)",
            c_notation="c.500del",
        )

        # c.500 is in E6, while the predicted stop p.186 is in E7. Appendix D
        # assigns the paired PM5 code by the exon containing the nucleotide
        # change. E6 has PM5_N/A, whereas E7 has PM5_Strong.
        self.assertEqual(result["criteria"]["PVS1"]["strength"], "Very Strong")
        self.assertNotIn("PM5_PTC", result["criteria"])
        self.assertEqual(result["total_points"], 8)
        self.assertEqual(result["predicted_class"], 3)
        reason = result["criteria"]["PVS1"]["reason"]
        self.assertIn("predicted termination p.186", reason)
        self.assertIn("nucleotide-change exon E6", reason)

    def test_bp1_strong_alone_is_likely_benign_exception(self):
        result = evaluate_variant(
            gene="BRCA2",
            variant_type="inframe_deletion",
            p_notation="p.(Val2050del)",
            c_notation="c.6147_6149del",
            spliceai_score=0.05,
        )
        self.assertEqual(result["criteria"]["BP1"]["strength"], "Strong")
        self.assertEqual(result["predicted_class"], 2)
        self.assertEqual(result["predicted_label"], "Likely Benign")

    def test_combined_batch_hgvs_input_enables_bp1_strong(self):
        c_notation, p_notation = split_combined_hgvs(
            "c.6147_6149del (p.Val2050del)"
        )
        result = evaluate_variant(
            gene="BRCA2",
            variant_type=infer_variant_type(c_notation, p_notation),
            p_notation=p_notation,
            c_notation=c_notation,
            spliceai_score=0.05,
        )
        self.assertEqual(result["criteria"]["BP1"]["strength"], "Strong")
        self.assertEqual(result["predicted_class"], 2)
        self.assertEqual(result["predicted_label"], "Likely Benign")


class Table4RareSpliceBranchTests(unittest.TestCase):
    def _assert_branch(
        self,
        *,
        gene,
        c_notation,
        source_code,
        source_count,
        strength,
        points,
        requires_rna,
    ):
        observed_count = sum(
            entry.get("pvs1_code") == source_code
            for gene_rules in TABLE4_DATA["splice_rules"].values()
            for entry in gene_rules.values()
        )
        self.assertEqual(observed_count, source_count)

        lookup = table4_lookup_splice(gene, c_notation)
        self.assertTrue(lookup["found"])
        self.assertEqual(lookup["pvs1_code"], source_code)
        self.assertEqual(
            parse_pvs1_code_strength(source_code),
            (strength, points, requires_rna),
        )

        module_result = evaluate_pvs1(
            gene,
            "splice_site",
            "p.(?)",
            c_notation=c_notation,
        )
        self.assertEqual(module_result["pvs1_code"], source_code)

        result = evaluate_variant(
            gene=gene,
            variant_type="splice_site",
            p_notation="p.(?)",
            c_notation=c_notation,
        )
        if strength is None:
            self.assertFalse(module_result["applies"])
            self.assertNotIn("PVS1", result["criteria"])
            self.assertEqual(
                result["not_applicable_criteria"]["PVS1"]["strength"],
                "N/A",
            )
        elif requires_rna:
            self.assertFalse(module_result["applies"])
            self.assertTrue(module_result["requires_rna"])
            self.assertNotIn("PVS1", result["criteria"])
            self.assertNotIn("PVS1_RNA", result["criteria"])
            self.assertTrue(result["rna_review"]["recommended"])
            self.assertIn("PVS1 (RNA)", result["rna_review"]["potential_branches"])
        else:
            self.assertTrue(module_result["applies"])
            self.assertEqual(module_result["strength"], strength)
            self.assertEqual(module_result["points"], points)
            self.assertEqual(result["criteria"]["PVS1"]["strength"], strength)
            self.assertEqual(result["criteria"]["PVS1"]["points"], points)

    def test_pvs1_moderate_splice_branch(self):
        self._assert_branch(
            gene="BRCA1",
            c_notation="c.671-2A>G",
            source_code="PVS1_Moderate",
            source_count=10,
            strength="Moderate",
            points=2,
            requires_rna=False,
        )

    def test_pvs1_supporting_rna_splice_branch(self):
        self._assert_branch(
            gene="BRCA1",
            c_notation="c.301+1G>C",
            source_code="PVS1_Supporting (RNA)",
            source_count=11,
            strength="Supporting",
            points=1,
            requires_rna=True,
        )

    def test_pvs1_strong_splice_branch(self):
        self._assert_branch(
            gene="BRCA2",
            c_notation="c.8331+1G>C",
            source_code="PVS1_Strong",
            source_count=3,
            strength="Strong",
            points=4,
            requires_rna=False,
        )

    def test_pvs1_strong_rna_splice_branch(self):
        self._assert_branch(
            gene="BRCA2",
            c_notation="c.8331+1G>A",
            source_code="PVS1_Strong (RNA)",
            source_count=3,
            strength="Strong",
            points=4,
            requires_rna=True,
        )

    def test_pvs1_not_applicable_rna_splice_branch(self):
        self._assert_branch(
            gene="BRCA1",
            c_notation="c.594-2A>C",
            source_code="PVS1_N/A (RNA)",
            source_count=2,
            strength=None,
            points=0,
            requires_rna=True,
        )


class VusExplanationTests(unittest.TestCase):
    def result_with(self, criteria, predicted_class=3, total_points=0):
        return {
            "predicted_class": predicted_class,
            "total_points": total_points,
            "criteria": criteria,
        }

    def criterion(self, strength, points):
        return {"applies": True, "strength": strength, "points": points, "reason": "test"}

    def test_pm2_only_explanation(self):
        explanation = explain_vus(
            self.result_with({"PM2_Supporting": self.criterion("Supporting", 1)}, total_points=1)
        )
        self.assertEqual(explanation["category"], "pm2_only")
        self.assertEqual(explanation["review_priority"], "low")

    def test_ps3_pm2_one_step_short_explanation(self):
        explanation = explain_vus(
            self.result_with(
                {
                    "PM2_Supporting": self.criterion("Supporting", 1),
                    "PS3": self.criterion("Strong", 4),
                },
                total_points=5,
            )
        )
        self.assertEqual(explanation["category"], "strong_pathogenic_one_step_short")
        self.assertIn("PP1", explanation["what_to_check"])

    def test_ps3_pp3_is_tier_a_explanation(self):
        explanation = explain_vus(
            self.result_with(
                {
                    "PP3": self.criterion("Supporting", 1),
                    "PS3": self.criterion("Strong", 4),
                },
                total_points=5,
            )
        )
        self.assertEqual(explanation["category"], "ps3_pp3_one_step_short")
        self.assertEqual(explanation["tier"], "A")

    def test_pvs1_explanation_says_vus_is_expected_not_erroneous(self):
        explanation = explain_vus(
            self.result_with({"PVS1": self.criterion("Very Strong", 8)}, total_points=8)
        )
        self.assertEqual(explanation["category"], "pvs1_needs_support")
        self.assertIn("does not mean that the classification", explanation["summary"])
        self.assertIn("VUS is therefore the expected result", explanation["summary"])
        self.assertIn("not that an error was detected", explanation["summary"])
        self.assertIn("Do not add evidence solely", explanation["what_to_check"])

    def test_bp4_bp7_pm2_explanation(self):
        explanation = explain_vus(
            self.result_with(
                {
                    "BP4": self.criterion("Supporting", -1),
                    "BP7": self.criterion("Supporting", -1),
                    "PM2_Supporting": self.criterion("Supporting", 1),
                },
                total_points=-1,
            )
        )
        self.assertEqual(explanation["category"], "bp4_bp7_pm2_benign_leaning")

    def test_pvs1_bs3_conflict_explanation(self):
        explanation = explain_vus(
            self.result_with(
                {
                    "PVS1": self.criterion("Very Strong", 8),
                    "BS3": self.criterion("Strong", -4),
                    "PM2_Supporting": self.criterion("Supporting", 1),
                },
                total_points=5,
            )
        )
        self.assertEqual(explanation["category"], "conflicting_pvs1_bs3")
        self.assertEqual(explanation["tier"], "C")

    def test_pp3_bs3_conflict_explanation(self):
        explanation = explain_vus(
            self.result_with(
                {
                    "PP3": self.criterion("Supporting", 1),
                    "BS3": self.criterion("Strong", -4),
                },
                total_points=-3,
            )
        )
        self.assertEqual(explanation["category"], "conflicting_pp3_bs3")

    def test_non_vus_has_no_explanation(self):
        explanation = explain_vus(
            self.result_with({"BP1": self.criterion("Strong", -4)}, predicted_class=2, total_points=-4)
        )
        self.assertIsNone(explanation)

    def test_narrative_includes_vus_explanation(self):
        result = self.result_with(
            {
                "PM2_Supporting": self.criterion("Supporting", 1),
                "PS3": self.criterion("Strong", 4),
            },
            total_points=5,
        )
        narrative = generate_narrative(
            gene="BRCA1",
            c_notation="c.3G>A",
            p_notation="p.(Met1Ile)",
            variant_type="initiation_codon",
            result=result,
            spliceai_score=0.23,
        )
        self.assertIn("VUS explanation", narrative)
        self.assertIn("Strong pathogenic evidence", narrative)

    def test_narrative_distinguishes_same_residue_reference_from_submitted_variant(self):
        result = self.result_with(
            {"PS3": self.criterion("Strong", 4)},
            total_points=4,
        )
        result["residue_info"] = {
            "domain": "BRCT",
            "known_pathogenic_at_position": [
                {"c": "c.5143A>C", "aa": "Ser1715Arg"},
            ],
        }
        narrative = generate_narrative(
            gene="BRCA1",
            c_notation="c.5145C>A",
            p_notation="p.(Ser1715Arg)",
            variant_type="missense",
            result=result,
        )
        self.assertIn("c.5143A>C p.(Ser1715Arg)", narrative)
        self.assertIn("not necessarily the submitted nucleotide variant", narrative)


class GoldenCaseRegressionTests(unittest.TestCase):
    def applied_codes(self, result):
        return {
            name
            for name, criterion in result["criteria"].items()
            if criterion.get("applies", True)
        }

    def assert_golden_case(
        self,
        result,
        *,
        predicted_class,
        total_points,
        criteria,
        vus_category=None,
    ):
        self.assertEqual(result["predicted_class"], predicted_class)
        self.assertEqual(result["total_points"], total_points)
        self.assertEqual(self.applied_codes(result), set(criteria))
        explanation = explain_vus(result)
        if predicted_class == 3:
            self.assertIsNotNone(explanation)
            self.assertEqual(explanation["category"], vus_category)
        else:
            self.assertIsNone(explanation)

    def test_pm2_only_vus_golden_case(self):
        result = evaluate_variant(
            gene="BRCA2",
            variant_type="initiation_codon",
            p_notation="p.(Met1Val)",
            c_notation="c.1A>G",
            spliceai_score=0.01,
            gnomad_data=gnomad_data(pm2_absence_established=True),
        )
        self.assert_golden_case(
            result,
            predicted_class=3,
            total_points=1,
            criteria={"PM2_Supporting"},
            vus_category="pm2_only",
        )

    def test_ps3_pm2_one_step_short_golden_case(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="initiation_codon",
            p_notation="p.(Met1Ile)",
            c_notation="c.3G>A",
            spliceai_score=0.23,
            gnomad_data=gnomad_data(pm2_absence_established=True),
            table9_result=table9_lookup_ps3_bs3("BRCA1", "c.3G>A"),
        )
        self.assert_golden_case(
            result,
            predicted_class=3,
            total_points=5,
            criteria={"PM2_Supporting", "PS3"},
            vus_category="strong_pathogenic_one_step_short",
        )

    def test_ps3_pp3_one_step_short_golden_case(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="missense",
            p_notation="p.(Cys64Gly)",
            c_notation="c.190T>G",
            spliceai_score=0.65,
            table9_result=table9_lookup_ps3_bs3("BRCA1", "c.190T>G"),
        )
        self.assert_golden_case(
            result,
            predicted_class=3,
            total_points=5,
            criteria={"PP3", "PS3"},
            vus_category="ps3_pp3_one_step_short",
        )

    def test_bp4_bp7_pm2_benign_leaning_vus_golden_case(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="synonymous",
            p_notation="p.(Gln74=)",
            c_notation="c.222A>G",
            spliceai_score=0.10,
            gnomad_data=gnomad_data(pm2_absence_established=True),
        )
        self.assert_golden_case(
            result,
            predicted_class=3,
            total_points=-1,
            criteria={"BP4", "BP7", "PM2_Supporting"},
            vus_category="bp4_bp7_pm2_benign_leaning",
        )

    def test_pvs1_bs3_conflict_golden_case(self):
        result = evaluate_variant(
            gene="BRCA2",
            variant_type="nonsense",
            p_notation="p.(Glu3309Ter)",
            c_notation="c.9925G>T",
            spliceai_score=0.01,
            gnomad_data=gnomad_data(pm2_absence_established=True),
            table9_result=table9_lookup_ps3_bs3("BRCA2", "c.9925G>T"),
        )
        self.assert_golden_case(
            result,
            predicted_class=4,
            total_points=9,
            criteria={"BS3", "PM2_Supporting", "PVS1", "PM5_PTC"},
        )
        self.assertTrue(result["mixed_evidence"])
        self.assertEqual(result["evidence_direction"], "mixed")
        self.assertEqual(result["pathogenic_points"], 13)
        self.assertEqual(result["benign_points"], -4)

    def test_pp3_bs3_conflict_golden_case(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="missense",
            p_notation="p.(His1732Gln)",
            c_notation="c.5196T>G",
            spliceai_score=0.99,
        )
        self.assert_golden_case(
            result,
            predicted_class=2,
            total_points=-3,
            criteria={"BS3", "PP3"},
        )
        self.assertTrue(result["mixed_evidence"])
        self.assertEqual(result["pathogenic_points"], 1)
        self.assertEqual(result["benign_points"], -4)

    def test_likely_benign_has_no_vus_explanation_golden_case(self):
        result = evaluate_variant(
            gene="BRCA2",
            variant_type=infer_variant_type("c.6147_6149del", "p.(Val2050del)"),
            p_notation="p.(Val2050del)",
            c_notation="c.6147_6149del",
            spliceai_score=0.05,
        )
        self.assert_golden_case(
            result,
            predicted_class=2,
            total_points=-4,
            criteria={"BP1"},
        )

    def test_pathogenic_has_no_vus_explanation_golden_case(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="initiation_codon",
            p_notation="p.(Met1Val)",
            c_notation="c.1A>G",
            spliceai_score=0.05,
            gnomad_data=gnomad_data(pm2_absence_established=True),
            table9_result=table9_lookup_ps3_bs3("BRCA1", "c.1A>G"),
            pp4_bp5_result={
                "applies": True,
                "code": "PP4",
                "strength": "Very Strong",
                "points": 8,
                "reason": "synthetic multifactorial likelihood evidence",
            },
        )
        self.assert_golden_case(
            result,
            predicted_class=5,
            total_points=13,
            criteria={"PM2_Supporting", "PP4", "PS3"},
        )


class DomainTests(unittest.TestCase):
    def test_brca1_ring_starts_at_residue_two(self):
        self.assertFalse(is_in_functional_domain("BRCA1", 1)[0])
        self.assertTrue(is_in_functional_domain("BRCA1", 2)[0])


if __name__ == "__main__":
    unittest.main()
