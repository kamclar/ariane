import unittest
from pathlib import Path

from backend.modules.frequency import (
    GNOMAD_LOCAL_DATASET_CONFIG,
    _aggregate_coverage_from_dataset_results,
    evaluate_frequency_criteria,
)
from backend.modules.classifier import evaluate_variant
from backend.modules.table9 import table9_lookup_ps3_bs3
from backend.modules.bp7 import evaluate_bp7
from backend.modules.pp3_bp4 import evaluate_pp3_bp4
from backend.modules.pvs1 import evaluate_pvs1
from backend.modules.ps1 import evaluate_ps1
from backend.modules.utils import is_in_functional_domain
from backend.modules.variant_type import infer_variant_type
from backend.modules.hgvs import split_combined_hgvs
from backend.modules.vus_explanation import explain_vus
from backend.modules.narrative import generate_narrative


def gnomad_data(
    *,
    max_af=None,
    found=False,
    v2_status="absent",
    v3_status="absent",
    v2_depth=30.0,
    v3_depth=30.0,
    pm2_absence_established=False,
):
    def dataset(status, depth, dataset_max_af=None):
        return {
            "status": status,
            "max_af": dataset_max_af,
            "coverage": {"mean_depth": depth},
        }

    return {
        "status": "found" if found else "absent_with_coverage",
        "found": found,
        "max_af": max_af,
        "frequency_metric": "faf95",
        "pm2_absence_established": pm2_absence_established,
        "datasets": {
            "v2_1_non_cancer": dataset(v2_status, v2_depth, max_af if v2_status == "found" else None),
            "v3_1_non_cancer": dataset(v3_status, v3_depth),
        },
    }


class VariantTypeTests(unittest.TestCase):
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


class FrequencyTests(unittest.TestCase):
    def test_pm2_configuration_requires_both_gnomad_versions(self):
        self.assertTrue(GNOMAD_LOCAL_DATASET_CONFIG["v2_1_non_cancer"]["required_for_pm2"])
        self.assertTrue(GNOMAD_LOCAL_DATASET_CONFIG["v3_1_non_cancer"]["required_for_pm2"])

        coverage = _aggregate_coverage_from_dataset_results(
            {
                "v2_1_non_cancer": {"coverage": {"mean_depth": 30.0}},
                "v3_1_non_cancer": {"coverage": {"mean_depth": 24.0}},
            }
        )
        self.assertFalse(coverage["passes_pm2"])

    def test_pm2_requires_absence_in_both_gnomad_versions(self):
        data = gnomad_data(pm2_absence_established=False)
        self.assertNotIn("PM2_Supporting", evaluate_frequency_criteria(data, "missense"))

        data["pm2_absence_established"] = True
        self.assertIn("PM2_Supporting", evaluate_frequency_criteria(data, "missense"))

    def test_pm2_is_not_used_for_indels(self):
        data = gnomad_data(pm2_absence_established=True)
        self.assertFalse(evaluate_frequency_criteria(data, "inframe_deletion")["PM2"]["applies"])

    def test_ba1_requires_depth_20(self):
        data = gnomad_data(max_af=0.002, found=True, v2_status="found", v2_depth=19.0)
        result = evaluate_frequency_criteria(data, "missense")
        self.assertNotIn("BA1", result)
        self.assertIn("_gnomad_info", result)

        data["datasets"]["v2_1_non_cancer"]["coverage"]["mean_depth"] = 20.0
        self.assertIn("BA1", evaluate_frequency_criteria(data, "missense"))

class SpliceTests(unittest.TestCase):
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
        )
        self.assertTrue(result["applies"])
        self.assertEqual(result["points"], -1)

    def test_intronic_bp7_does_not_apply_inside_conserved_motif(self):
        result = evaluate_bp7(
            "intronic",
            spliceai_score=0.05,
            bp4_met=True,
            c_notation="c.100+6A>G",
        )
        self.assertFalse(result["applies"])

    def test_ps1_requires_confirmed_low_spliceai(self):
        result = evaluate_ps1(
            gene="BRCA1",
            c_notation="c.509G>A",
            p_notation="p.(Arg170Gln)",
            variant_type="missense",
            spliceai_score=None,
        )
        self.assertFalse(result["applies"])
        self.assertIn("SpliceAI score not available", result["reason"])

    def test_splice_ps1_reference_pilot_is_unreviewed_seed_only(self):
        project_root = Path(__file__).resolve().parents[1]
        path = project_root / "backend/data/splice_ps1_reference_set.json"
        self.assertTrue(path.exists())
        import json
        data = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(
            data["curation_status"],
            "pilot_unreviewed_not_for_automatic_scoring",
        )
        self.assertGreaterEqual(len(data["variants"]), 1)
        self.assertTrue(
            any(
                variant["gene"] == "BRCA1"
                and variant["reference_variant"] == "c.4185G>A"
                for variant in data["variants"]
            )
        )

    def test_splice_ps1_reference_candidates_are_available_for_ui_prefill(self):
        from backend.modules.splice_ps1_reference import (
            load_splice_ps1_reference_candidates,
        )

        data = load_splice_ps1_reference_candidates("BRCA1")
        self.assertEqual(
            data["curation_status"],
            "pilot_unreviewed_not_for_automatic_scoring",
        )
        self.assertGreaterEqual(len(data["candidates"]), 1)
        candidate = next(
            item
            for item in data["candidates"]
            if item["reference_variant"] == "c.4185G>A"
        )
        self.assertEqual(candidate["gene"], "BRCA1")
        self.assertEqual(candidate["classification"], "Pathogenic")
        self.assertEqual(candidate["prefill_strength_suggestion"], "Strong")
        self.assertIn("exon 12", candidate["splice_event_label"])
        self.assertIn("row", candidate["source_label"])


class ClassifierIntegrationTests(unittest.TestCase):
    def test_custom_donor_guard_is_not_part_of_active_scoring(self):
        project_root = Path(__file__).resolve().parents[1]
        self.assertFalse((project_root / "backend/modules/donor_guard.py").exists())
        self.assertNotIn("donor_guard", (project_root / "backend/main.py").read_text(encoding="utf-8"))

    def test_rna_dependent_pvs1_does_not_leak_into_score(self):
        unconfirmed = evaluate_variant(
            gene="BRCA1",
            variant_type="splice_site",
            p_notation="p.(?)",
            c_notation="c.80+1G>A",
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
        unknown = evaluate_variant(
            gene="BRCA1",
            variant_type="exon_duplication",
            p_notation="p.(?)",
            c_notation=notation,
        )
        tandem = evaluate_variant(
            gene="BRCA1",
            variant_type="exon_duplication",
            p_notation="p.(?)",
            c_notation=notation,
            dup_type="Tandem",
        )
        self.assertEqual(unknown["criteria"]["PVS1"]["strength"], "Moderate")
        self.assertEqual(tandem["criteria"]["PVS1"]["strength"], "Strong")

    def test_bs3_strong_alone_remains_vus(self):
        table9_result = table9_lookup_ps3_bs3("BRCA1", "c.3891_3893del")
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="inframe_deletion",
            p_notation="p.(Ser1298del)",
            c_notation="c.3891_3893del",
            table9_result=table9_result,
        )
        self.assertEqual(result["criteria"]["BS3"]["points"], -4)
        self.assertEqual(result["total_points"], -4)
        self.assertEqual(result["predicted_class"], 3)
        self.assertEqual(result["predicted_label"], "VUS")

    def test_pvs1_very_strong_alone_remains_vus_with_explanation(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="frameshift",
            p_notation="p.(Asp1851ValfsTer29)",
            c_notation="c.5551_5552insT",
        )
        self.assertEqual(result["criteria"]["PVS1"]["strength"], "Very Strong")
        self.assertEqual(result["total_points"], 8)
        self.assertEqual(result["predicted_class"], 3)
        self.assertIn(
            "requires at least one additional Supporting",
            result["classification_note"],
        )

    def test_bp1_strong_alone_is_likely_benign_exception(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="inframe_deletion",
            p_notation="p.(Ser1298del)",
            c_notation="c.3891_3893del",
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
            predicted_class=3,
            total_points=5,
            criteria={"BS3", "PM2_Supporting", "PVS1"},
            vus_category="conflicting_pvs1_bs3",
        )

    def test_pp3_bs3_conflict_golden_case(self):
        result = evaluate_variant(
            gene="BRCA1",
            variant_type="missense",
            p_notation="p.(His1732Gln)",
            c_notation="c.5196T>G",
            spliceai_score=0.99,
            table9_result={
                "applies": True,
                "code": "BS3",
                "strength": "Strong",
                "points": -4,
                "reason": "synthetic calibrated functional benign evidence",
            },
        )
        self.assert_golden_case(
            result,
            predicted_class=2,
            total_points=-3,
            criteria={"BS3", "PP3"},
        )

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
