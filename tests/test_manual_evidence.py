import unittest

from backend.lookups.clinvar import clinvar_review_stars
from backend.modules.manual_evidence import (
    evaluate_manual_evidence,
    suggest_strength,
)


class ManualStrengthSuggestionTests(unittest.TestCase):
    def test_ps4_requires_all_enigma_thresholds(self):
        self.assertEqual(
            suggest_strength(
                "PS4",
                {
                    "p_value": 0.05,
                    "odds_ratio": 4,
                    "lower_ci": 2.01,
                    "case_control_country_matched": True,
                    "case_control_ethnicity_matched": True,
                },
            ),
            "Strong",
        )
        self.assertIsNone(
            suggest_strength(
                "PS4",
                {
                    "p_value": 0.05,
                    "odds_ratio": 4,
                    "lower_ci": 2.0,
                    "case_control_country_matched": True,
                    "case_control_ethnicity_matched": True,
                },
            )
        )

    def test_ps4_requires_country_and_ethnicity_matching(self):
        common = {
            "p_value": 0.01,
            "odds_ratio": 5,
            "lower_ci": 2.1,
            "case_control_country_matched": True,
            "case_control_ethnicity_matched": True,
        }
        for missing_field in (
            "case_control_country_matched",
            "case_control_ethnicity_matched",
        ):
            evidence = {**common, missing_field: False}
            self.assertIsNone(suggest_strength("PS4", evidence))

    def test_pm3_and_bs2_point_thresholds(self):
        for code in ("PM3", "BS2"):
            required = {
                "cooccurring_variant_classification_basis": "vcep_specifications",
            }
            if code == "PM3":
                required["vua_benign_population_review"] = "does_not_meet"
            self.assertEqual(
                suggest_strength(code, {**required, "evidence_points": 1}), "Supporting"
            )
            self.assertEqual(
                suggest_strength(code, {**required, "evidence_points": 2}), "Moderate"
            )
            self.assertEqual(
                suggest_strength(code, {**required, "evidence_points": 4}), "Strong"
            )

    def test_pm3_and_bs2_require_vcep_classified_cooccurring_variant(self):
        for code in ("PM3", "BS2"):
            evidence = {"evidence_points": 4}
            if code == "PM3":
                evidence["vua_benign_population_review"] = "does_not_meet"
            self.assertIsNone(suggest_strength(code, evidence))

    def test_pm3_requires_no_benign_population_criterion(self):
        common = {
            "evidence_points": 4,
            "cooccurring_variant_classification_basis": "vcep_specifications",
        }
        self.assertIsNone(suggest_strength("PM3", common))
        self.assertIsNone(suggest_strength("PM3", {
            **common,
            "vua_benign_population_review": "meets",
        }))

    def test_pp1_likelihood_ratio_thresholds(self):
        self.assertEqual(
            suggest_strength("PP1", {"likelihood_ratio": 2.08}), "Supporting"
        )
        self.assertEqual(
            suggest_strength("PP1", {"likelihood_ratio": 4.3}), "Moderate"
        )
        self.assertEqual(
            suggest_strength("PP1", {"likelihood_ratio": 18.7}), "Strong"
        )
        self.assertEqual(suggest_strength("PP1", {"likelihood_ratio": 350}), "Strong")
        self.assertEqual(suggest_strength("PP1", {
            "likelihood_ratio": 350,
            "very_strong_effect_basis": "predicted_protein",
        }), "Very Strong")

    def test_pp4_combined_clinical_lr_thresholds_require_provenance(self):
        def evidence(lr):
            return {
                "combined_clinical_lr": lr,
                "source_review_status": "appendix_b",
                "source_pmid": "31853058",
                "clinical_data_summary": "Family history and case-control data; overlap reviewed.",
            }

        self.assertEqual(suggest_strength("PP4", evidence(2.08)), "Supporting")
        self.assertEqual(suggest_strength("PP4", evidence(4.3)), "Moderate")
        self.assertEqual(suggest_strength("PP4", evidence(18.7)), "Strong")
        self.assertEqual(suggest_strength("PP4", evidence(350)), "Very Strong")
        incomplete = evidence(350)
        incomplete["source_pmid"] = "99999999"
        self.assertIsNone(suggest_strength("PP4", incomplete))

    def test_pp4_accepts_zanti_enigma_case_control_source(self):
        self.assertEqual(
            suggest_strength("PP4", {
                "combined_clinical_lr": 18.7,
                "source_review_status": "appendix_b",
                "source_pmid": "40413188",
                "clinical_data_summary": (
                    "Zanti 2025 BRIDGES, CARRIERS and UK Biobank case-control LR."
                ),
            }),
            "Strong",
        )

    def test_pp4_accepts_canvar_log10_lr_and_acmg_point_scales(self):
        common = {
            "source_review_status": "appendix_b",
            "source_pmid": "31853058",
            "clinical_data_summary": "Personal and family history model; one published study.",
        }
        self.assertEqual(
            suggest_strength("PP4", {
                **common,
                "clinical_lr_value": 45.84,
                "clinical_lr_scale": "log10_lr",
            }),
            "Very Strong",
        )
        self.assertEqual(suggest_strength("PP4", {
            **common,
            "clinical_lr_value": 8,
            "clinical_lr_scale": "acmg_points",
        }), "Very Strong")
        self.assertEqual(suggest_strength("PP4", {
            **common,
            "clinical_lr_value": 350,
            "clinical_lr_scale": "lr",
        }), "Very Strong")

    def test_pp4_other_source_requires_review_and_unreviewed_does_not_score(self):
        common = {
            "clinical_lr_value": 350,
            "clinical_lr_scale": "lr",
            "source_citation": "PMID:99999999",
            "clinical_data_summary": "Variant-specific clinical LR; overlap assessed.",
        }
        self.assertIsNone(suggest_strength("PP4", {
            **common,
            "source_review_status": "unreviewed",
        }))
        self.assertEqual(suggest_strength("PP4", {
            **common,
            "source_review_status": "other_reviewed",
            "source_reviewed_by": "Expert reviewer",
            "source_review_rationale": "Compatible multifactorial model and independent cohort.",
        }), "Very Strong")

    def test_bp5_uses_benign_clinical_lr_thresholds_and_reviewed_provenance(self):
        common = {
            "clinical_lr_scale": "lr",
            "source_review_status": "appendix_b",
            "source_pmid": "31853058",
            "clinical_data_summary": "Segregation and pathology contributions; overlap reviewed.",
        }
        self.assertEqual(
            suggest_strength("BP5", {**common, "clinical_lr_value": 0.48}),
            "Supporting",
        )
        self.assertEqual(
            suggest_strength("BP5", {**common, "clinical_lr_value": 0.05}),
            "Strong",
        )
        self.assertIsNone(suggest_strength("BP5", {
            **common,
            "clinical_lr_value": 0.05,
            "source_review_status": "unreviewed",
            "source_citation": "PMID:31853058",
        }))

    def test_manual_ps3_bs3_require_complete_enigma_calibration_review(self):
        common = {
            "assay_scope": "protein_only",
            "calibration_status": "reviewed_under_enigma_vcep",
            "pathogenic_and_benign_controls_confirmed": True,
            "curated_strength": "Strong",
            "assay_name": "Saturation functional assay",
            "source_citation": "PMID:99999999",
            "calibration_summary": "OddsPath calibration reviewed against pathogenic and benign controls.",
            "variant_result_summary": "Three replicates support the reported functional category.",
            "functional_reviewed_by": "Expert reviewer",
        }
        self.assertEqual(
            suggest_strength("PS3", {**common, "functional_conclusion": "abnormal"}),
            "Strong",
        )
        self.assertEqual(
            suggest_strength("BS3", {**common, "functional_conclusion": "normal"}),
            "Strong",
        )
        incomplete = {**common, "functional_conclusion": "abnormal"}
        incomplete["pathogenic_and_benign_controls_confirmed"] = False
        self.assertIsNone(suggest_strength("PS3", incomplete))
        unsupported_strength = {
            **common,
            "functional_conclusion": "abnormal",
            "curated_strength": "Moderate",
        }
        self.assertIsNone(suggest_strength("PS3", unsupported_strength))

    def test_bs4_likelihood_ratio_thresholds(self):
        self.assertEqual(
            suggest_strength("BS4", {"likelihood_ratio": 0.48}), "Supporting"
        )
        self.assertEqual(
            suggest_strength("BS4", {"likelihood_ratio": 0.23}), "Moderate"
        )
        self.assertEqual(
            suggest_strength("BS4", {"likelihood_ratio": 0.05}), "Strong"
        )
        self.assertEqual(
            suggest_strength("BS4", {"likelihood_ratio": 0.00285}),
            "Very Strong",
        )

    def test_pvs1_rna_requires_curated_mrna_only_record(self):
        evidence = {
            "assay_scope": "mrna_only",
            "rna_conclusion": "damaging",
            "functional_transcript_remaining": "absent_or_minimal",
            "curated_strength": "Strong",
            "transcript_accession": "NM_007294.4",
            "tissue_or_cell_type": "lymphocytes",
            "nmd_assessed": "yes",
        }
        self.assertEqual(suggest_strength("PVS1_RNA", evidence), "Strong")
        evidence["functional_transcript_remaining"] = "uncertain"
        self.assertIsNone(suggest_strength("PVS1_RNA", evidence))
        evidence["functional_transcript_remaining"] = "absent_or_minimal"
        evidence["assay_scope"] = "combined_mrna_protein"
        self.assertIsNone(suggest_strength("PVS1_RNA", evidence))

    def test_bp7_rna_is_strong_only_with_variant_context_eligibility(self):
        evidence = {
            "assay_scope": "mrna_only",
            "rna_conclusion": "no_damaging_effect",
            "transcript_accession": "NM_000059.4",
            "tissue_or_cell_type": "blood",
            "nmd_assessed": "not_applicable",
        }
        intronic_context = {
            "gene": "BRCA2",
            "c_notation": "c.8953+3T>C",
            "p_notation": "p.?",
        }
        self.assertEqual(
            suggest_strength(
                "BP7_RNA", evidence, variant_context=intronic_context
            ),
            "Strong",
        )
        self.assertIsNone(suggest_strength("BP7_RNA", evidence))

    def test_bp7_rna_domain_missense_requires_table9_bs3(self):
        evidence = {
            "assay_scope": "mrna_only",
            "rna_conclusion": "no_damaging_effect",
            "transcript_accession": "NM_007294.4",
            "tissue_or_cell_type": "blood",
            "nmd_assessed": "not_applicable",
        }
        context = {
            "gene": "BRCA1",
            "c_notation": "c.5123C>T",
            "p_notation": "p.(Ala1708Val)",
        }
        self.assertIsNone(
            suggest_strength("BP7_RNA", evidence, variant_context=context)
        )
        bs3 = [{
            "name": "BS3",
            "applies": True,
            "strength": "Strong",
            "points": -4,
            "decision_path": {
                "sources": [{"source_id": "enigma-v1.2-table9"}],
            },
        }]
        self.assertEqual(
            suggest_strength(
                "BP7_RNA",
                evidence,
                variant_context=context,
                base_criteria=bs3,
            ),
            "Strong",
        )

    def test_complete_manual_bs3_can_satisfy_bp7_rna_domain_stipulation(self):
        manual_bs3 = {
            "code": "BS3",
            "enabled": True,
            "evidence": {
                "assay_scope": "protein_only",
                "functional_conclusion": "normal",
                "calibration_status": "reviewed_under_enigma_vcep",
                "pathogenic_and_benign_controls_confirmed": True,
                "curated_strength": "Strong",
                "assay_name": "Calibrated functional assay",
                "source_citation": "PMID:99999999",
                "calibration_summary": "Reviewed against pathogenic and benign controls.",
                "variant_result_summary": "Variant result is similar to benign controls.",
                "functional_reviewed_by": "Expert reviewer",
            },
        }
        manual_bp7 = {
            "code": "BP7_RNA",
            "enabled": True,
            "evidence": {
                "assay_scope": "mrna_only",
                "rna_conclusion": "no_damaging_effect",
                "transcript_accession": "NM_007294.4",
                "tissue_or_cell_type": "blood",
                "nmd_assessed": "not_applicable",
            },
        }
        result = evaluate_manual_evidence(
            [],
            [manual_bp7, manual_bs3],
            {
                "gene": "BRCA1",
                "c_notation": "c.5123C>T",
                "p_notation": "p.(Ala1708Val)",
            },
        )
        applied = {
            item["code"] for item in result["manual_criteria"] if item["applies"]
        }
        self.assertEqual(applied, {"BS3", "BP7_RNA"})

    def test_bp7_rna_does_not_accept_unproven_bs3_for_domain_missense(self):
        evidence = {
            "assay_scope": "mrna_only",
            "rna_conclusion": "no_damaging_effect",
            "transcript_accession": "NM_007294.4",
            "tissue_or_cell_type": "blood",
            "nmd_assessed": "not_applicable",
        }
        context = {
            "gene": "BRCA1",
            "c_notation": "c.5123C>T",
            "p_notation": "p.(Ala1708Val)",
        }
        self.assertIsNone(suggest_strength(
            "BP7_RNA",
            evidence,
            variant_context=context,
            base_criteria=[{
                "name": "BS3",
                "applies": True,
                "strength": "Strong",
                "points": -4,
            }],
        ))

    def test_pvs1_init_requires_curated_flowchart_record(self):
        evidence = {
            "met1_loss_confirmed": True,
            "alternative_start_assessed": "yes",
            "nearest_alternative_start": "p.Met16",
            "upstream_pathogenic_evidence": "yes",
            "functional_domain_impact": "yes",
            "curated_strength": "Moderate",
            "initiation_flowchart_rationale": "start-loss flowchart supports PVS1_Moderate",
        }
        self.assertEqual(suggest_strength("PVS1_INIT", evidence), "Moderate")
        evidence["met1_loss_confirmed"] = False
        self.assertIsNone(suggest_strength("PVS1_INIT", evidence))

    def test_ps1_splice_requires_curated_same_event_record(self):
        evidence = {
            "reference_variant": "BRCA1 c.4185G>A",
            "reference_classification": "Pathogenic",
            "reference_classification_source": "ENIGMA ST2/ST7",
            "same_splice_event_confirmed": True,
            "vua_splice_event": "exon 12 deletion",
            "reference_splice_event": "exon 12 deletion",
            "prediction_strength_comparison": "similar",
            "curated_strength": "Moderate",
            "ps1_splice_rationale": "same predicted exon skipping event; Table 17 supports PS1_Moderate",
        }
        self.assertEqual(suggest_strength("PS1_SPLICE", evidence), "Moderate")
        evidence["prediction_strength_comparison"] = "weaker"
        self.assertIsNone(suggest_strength("PS1_SPLICE", evidence))

    def test_ps1_protein_strength_is_derived_from_verified_reference_class(self):
        evidence = {
            "reference_variant": "BRCA1 c.122A>G",
            "reference_p_notation": "p.(His41Arg)",
            "reference_classification": "Pathogenic",
            "classification_verification": "external_vcep_assertion",
            "classification_source": "ClinGen ENIGMA expert panel assertion",
            "same_missense_confirmed": True,
            "different_nucleotide_change_confirmed": True,
            "vua_spliceai_score": 0.01,
            "reference_spliceai_score": 0.02,
            "splice_source_check_completed": True,
            "splice_sources_checked": [
                "ENIGMA Specifications Table 9 v1.2",
                "ENIGMA Supplementary Table 2 v1.2",
            ],
            "vua_confirmed_splice_status": "none_identified",
            "reference_confirmed_splice_status": "normal",
            "reference_classification_used_ps1": "unknown",
            "ps1_protein_rationale": "VCEP classification and both splice checks reviewed.",
        }
        self.assertEqual(suggest_strength("PS1_PROTEIN", evidence), "Strong")
        evidence["reference_classification"] = "Likely Pathogenic"
        self.assertEqual(suggest_strength("PS1_PROTEIN", evidence), "Moderate")
        evidence["classification_verification"] = "historical_classification_only"
        self.assertIsNone(suggest_strength("PS1_PROTEIN", evidence))

    def test_ps1_protein_known_dependency_requires_reciprocal_check(self):
        evidence = {
            "reference_variant": "BRCA1 c.122A>G",
            "reference_p_notation": "p.(His41Arg)",
            "reference_classification": "Pathogenic",
            "classification_verification": "external_vcep_assertion",
            "classification_source": "ClinGen ENIGMA expert panel assertion",
            "same_missense_confirmed": True,
            "different_nucleotide_change_confirmed": True,
            "vua_spliceai_score": 0.01,
            "reference_spliceai_score": 0.02,
            "splice_source_check_completed": True,
            "splice_sources_checked": [
                "ENIGMA Specifications Table 9 v1.2",
                "ENIGMA Supplementary Table 2 v1.2",
            ],
            "vua_confirmed_splice_status": "none_identified",
            "reference_confirmed_splice_status": "none_identified",
            "reference_classification_used_ps1": "yes",
            "ps1_protein_rationale": "Reference classification used PS1.",
        }
        self.assertIsNone(suggest_strength("PS1_PROTEIN", evidence))
        evidence["reference_ps1_dependency_reference"] = "BRCA1 c.130T>A"
        evidence["direct_reciprocal_dependency_excluded"] = True
        self.assertEqual(suggest_strength("PS1_PROTEIN", evidence), "Strong")


class ManualEvidenceClassificationTests(unittest.TestCase):
    def test_pp4_unreviewed_source_is_audited_without_points(self):
        result = evaluate_manual_evidence(
            [],
            [{
                "code": "PP4",
                "enabled": True,
                "evidence": {
                    "clinical_lr_value": 350,
                    "clinical_lr_scale": "lr",
                    "source_review_status": "unreviewed",
                    "source_citation": "PMID:99999999",
                    "clinical_data_summary": "Candidate clinical LR pending review.",
                },
            }],
        )
        criterion = result["manual_criteria"][0]
        self.assertFalse(criterion["applies"])
        self.assertEqual(criterion["points"], 0)

    def test_pp4_adds_strength_from_combined_clinical_lr(self):
        result = evaluate_manual_evidence(
            [],
            [{
                "code": "PP4",
                "enabled": True,
                "evidence": {
                    "combined_clinical_lr": 350,
                    "source_review_status": "appendix_b",
                    "source_pmid": "31853058",
                    "clinical_data_summary": "Variant-specific combined clinical evidence; overlap reviewed.",
                },
                "references": ["PMID:31853058"],
            }],
        )
        criterion = result["manual_criteria"][0]
        self.assertTrue(criterion["applies"])
        self.assertEqual(criterion["selected_strength"], "Very Strong")
        self.assertEqual(criterion["points"], 8)

    def test_pp4_cannot_be_enabled_without_provenance(self):
        with self.assertRaisesRegex(ValueError, "recorded source"):
            evaluate_manual_evidence(
                [],
                [{
                    "code": "PP4",
                    "enabled": True,
                    "evidence": {"combined_clinical_lr": 350},
                }],
            )

    def test_bp5_strong_single_criterion_route_requires_multiple_reviewed_types(self):
        evidence = {
            "clinical_lr_value": 0.05,
            "clinical_lr_scale": "lr",
            "source_review_status": "appendix_b",
            "source_pmid": "31853058",
            "clinical_data_summary": "Segregation and pathology LR; overlap reviewed.",
            "clinical_evidence_types": ["segregation", "pathology"],
            "independence_review_confirmed": True,
        }
        result = evaluate_manual_evidence([], [{
            "code": "BP5",
            "enabled": True,
            "evidence": evidence,
            "references": ["PMID:31853058"],
        }])
        criterion = result["manual_criteria"][0]
        self.assertEqual(criterion["selected_strength"], "Strong")
        self.assertTrue(criterion["single_strong_likely_benign_eligible"])
        self.assertEqual(result["predicted_class"], 2)

        without_independence = evaluate_manual_evidence([], [{
            "code": "BP5",
            "enabled": True,
            "evidence": {**evidence, "independence_review_confirmed": False},
            "references": ["PMID:31853058"],
        }])
        self.assertFalse(
            without_independence["manual_criteria"][0][
                "single_strong_likely_benign_eligible"
            ]
        )
        self.assertEqual(without_independence["predicted_class"], 3)

    def test_pp1_with_automatic_pp4_requires_independence_review(self):
        base = [{
            "name": "PP4", "applies": True, "strength": "Strong", "points": 4,
            "reason": "Local clinical LR snapshot",
        }]
        manual = [{
            "code": "PP1", "enabled": True,
            "evidence": {"likelihood_ratio": 4.3},
        }]
        with self.assertRaisesRegex(ValueError, "independent observations"):
            evaluate_manual_evidence(base, manual)

    def test_pp1_with_automatic_pp4_is_allowed_after_independence_review(self):
        base = [{
            "name": "PP4", "applies": True, "strength": "Strong", "points": 4,
            "reason": "Local clinical LR snapshot",
        }]
        manual = [{
            "code": "PP1", "enabled": True,
            "evidence": {
                "likelihood_ratio": 4.3,
                "independent_from_pp4_bp5": True,
                "independence_rationale": "Pedigree observations are absent from the PP4 source cohorts.",
            },
        }]
        result = evaluate_manual_evidence(base, manual)
        self.assertEqual(result["total_points"], 6)

    def test_ps4_with_manual_pp4_requires_independence_review(self):
        manual = [
            {
                "code": "PP4", "enabled": True,
                "evidence": {
                    "combined_clinical_lr": 350,
                    "source_review_status": "appendix_b",
                    "source_pmid": "31853058",
                    "clinical_data_summary": "Family-history LR; cohort recorded.",
                },
            },
            {
                "code": "PS4", "enabled": True,
                "evidence": {
                    "p_value": 0.01,
                    "odds_ratio": 5,
                    "lower_ci": 2.1,
                    "case_control_country_matched": True,
                    "case_control_ethnicity_matched": True,
                },
            },
        ]
        with self.assertRaisesRegex(ValueError, "PS4 cannot be combined"):
            evaluate_manual_evidence([], manual)

    def test_curated_pvs1_rna_removes_predictive_pp3(self):
        result = evaluate_manual_evidence(
            [{"name": "PP3", "applies": True, "strength": "Supporting", "points": 1}],
            [{
                "code": "PVS1_RNA", "enabled": True,
                "evidence": {
                    "assay_scope": "mrna_only", "rna_conclusion": "damaging",
                    "functional_transcript_remaining": "absent_or_minimal",
                    "curated_strength": "Very Strong", "transcript_accession": "NM_007294.4",
                    "tissue_or_cell_type": "blood", "nmd_assessed": "yes",
                },
            }],
        )
        self.assertEqual(result["total_points"], 8)
        self.assertEqual(
            result["evidence_interactions"][0]["suppressed"],
            ["PP3"],
        )

    def test_pvs1_rna_replaces_all_applicable_bioinformatic_codes(self):
        base = [
            {"name": "PP3", "applies": True, "strength": "Supporting", "points": 1},
            {"name": "BP1", "applies": True, "strength": "Strong", "points": -4},
            {"name": "BP4", "applies": True, "strength": "Supporting", "points": -1},
            {"name": "BP7", "applies": True, "strength": "Supporting", "points": -1},
            {"name": "PS1", "applies": True, "strength": "Strong", "points": 4},
        ]
        result = evaluate_manual_evidence(
            base,
            [{
                "code": "PVS1_RNA", "enabled": True,
                "evidence": {
                    "assay_scope": "mrna_only", "rna_conclusion": "damaging",
                    "functional_transcript_remaining": "absent_or_minimal",
                    "curated_strength": "Very Strong",
                    "transcript_accession": "NM_007294.4",
                    "tissue_or_cell_type": "blood", "nmd_assessed": "yes",
                },
            }],
        )
        self.assertEqual(result["total_points"], 8)
        warning = result["evidence_interactions"][0]
        self.assertEqual(warning["status"], "deduplicated")
        self.assertEqual(
            set(warning["suppressed"]),
            {"PP3", "BP1", "BP4", "BP7", "PS1"},
        )

    def test_bp7_rna_upgrades_bp7_but_retains_bp4(self):
        result = evaluate_manual_evidence(
            [
                {"name": "BP4", "applies": True, "strength": "Supporting", "points": -1},
                {"name": "BP7", "applies": True, "strength": "Supporting", "points": -1},
            ],
            [{
                "code": "BP7_RNA", "enabled": True,
                "evidence": {
                    "assay_scope": "mrna_only",
                    "rna_conclusion": "no_damaging_effect",
                    "transcript_accession": "NM_007294.4",
                    "tissue_or_cell_type": "blood",
                    "nmd_assessed": "not_applicable",
                },
            }],
            variant_context={
                "gene": "BRCA1",
                "c_notation": "c.4185G>A",
                "p_notation": "p.(Gln1395=)",
            },
        )
        self.assertEqual(result["total_points"], -5)
        warning = result["evidence_interactions"][0]
        self.assertEqual(warning["retained"], ["BP7_RNA"])
        self.assertEqual(warning["suppressed"], ["BP7"])

    def test_bp7_rna_and_splice_pp3_are_retained_as_conflict(self):
        result = evaluate_manual_evidence(
            [{"name": "PP3", "applies": True, "strength": "Supporting", "points": 1}],
            [{
                "code": "BP7_RNA", "enabled": True,
                "evidence": {
                    "assay_scope": "mrna_only",
                    "rna_conclusion": "no_damaging_effect",
                    "transcript_accession": "NM_007294.4",
                    "tissue_or_cell_type": "blood",
                    "nmd_assessed": "not_applicable",
                },
            }],
            variant_context={
                "gene": "BRCA1",
                "c_notation": "c.4987-6T>G",
                "p_notation": "p.?",
            },
        )
        self.assertEqual(result["total_points"], -3)
        warning = result["evidence_interactions"][0]
        self.assertEqual(warning["status"], "conflict")
        self.assertTrue(warning["review_required"])
        self.assertEqual(set(warning["retained"]), {"BP7_RNA", "PP3"})

    def test_manual_evidence_creates_separate_amended_result(self):
        base = [
            {
                "name": "PS3",
                "applies": True,
                "strength": "Strong",
                "points": 4,
                "reason": "calibrated functional evidence",
            }
        ]
        manual = [
            {
                "code": "PM3",
                "enabled": True,
                "evidence": {
                    "evidence_points": 1,
                    "cooccurring_variant_classification_basis": "vcep_specifications",
                    "vua_benign_population_review": "does_not_meet",
                },
                "notes": "one PM3 evidence point",
                "references": ["PMID:1"],
            },
            {
                "code": "PP1",
                "enabled": True,
                "evidence": {"likelihood_ratio": 2.08},
                "notes": "quantitative segregation",
                "references": ["PMID:2"],
            },
        ]

        result = evaluate_manual_evidence(base, manual)

        self.assertEqual(result["predicted_class"], 4)
        self.assertEqual(result["total_points"], 6)
        self.assertEqual(base[0]["points"], 4)

    def test_reviewer_cannot_override_a_rule_derived_strength(self):
        with self.assertRaisesRegex(ValueError, "overrides are not permitted"):
            evaluate_manual_evidence(
                [],
                [
                    {
                        "code": "PP1",
                        "enabled": True,
                        "evidence": {"likelihood_ratio": 4.3},
                        "override_strength": "Supporting",
                        "notes": "attempted conservative adjustment",
                        "references": ["PMID:3"],
                    }
                ],
            )

    def test_invalid_strength_for_criterion_is_rejected(self):
        with self.assertRaises(ValueError):
            evaluate_manual_evidence(
                [],
                [
                    {
                        "code": "PS4",
                        "enabled": True,
                        "evidence": {},
                        "override_strength": "Moderate",
                    }
                ],
            )

    def test_pvs1_rna_adds_pathogenic_points_to_amended_result(self):
        result = evaluate_manual_evidence(
            [],
            [
                {
                    "code": "PVS1_RNA",
                    "enabled": True,
                    "evidence": {
                        "assay_scope": "mrna_only",
                        "rna_conclusion": "damaging",
                        "functional_transcript_remaining": "absent_or_minimal",
                        "curated_strength": "Very Strong",
                        "transcript_accession": "NM_007294.4",
                        "tissue_or_cell_type": "fibroblasts",
                        "nmd_assessed": "yes",
                        "transcript_result_summary": "aberrant out-of-frame transcript with minimal normal transcript",
                    },
                    "notes": "curated RNA assay",
                    "references": ["PMID:4"],
                }
            ],
        )
        criterion = result["manual_criteria"][0]
        self.assertTrue(criterion["applies"])
        self.assertEqual(criterion["selected_strength"], "Very Strong")
        self.assertEqual(criterion["points"], 8)

    def test_bp7_rna_adds_benign_strong_points_to_amended_result(self):
        result = evaluate_manual_evidence(
            [],
            [
                {
                    "code": "BP7_RNA",
                    "enabled": True,
                    "evidence": {
                        "assay_scope": "mrna_only",
                        "rna_conclusion": "no_damaging_effect",
                        "transcript_accession": "NM_000059.4",
                        "tissue_or_cell_type": "blood",
                        "nmd_assessed": "not_applicable",
                        "transcript_result_summary": "normal transcript profile within assay sensitivity",
                    },
                    "notes": "curated RNA assay",
                    "references": ["PMID:5"],
                }
            ],
            variant_context={
                "gene": "BRCA2",
                "c_notation": "c.8953+3T>C",
                "p_notation": "p.?",
            },
        )
        criterion = result["manual_criteria"][0]
        self.assertTrue(criterion["applies"])
        self.assertEqual(criterion["selected_strength"], "Strong")
        self.assertEqual(criterion["points"], -4)

    def test_pvs1_init_adds_pathogenic_points_to_amended_result(self):
        result = evaluate_manual_evidence(
            [],
            [
                {
                    "code": "PVS1_INIT",
                    "enabled": True,
                    "evidence": {
                        "met1_loss_confirmed": True,
                        "alternative_start_assessed": "yes",
                        "nearest_alternative_start": "p.Met16",
                        "upstream_pathogenic_evidence": "yes",
                        "functional_domain_impact": "yes",
                        "curated_strength": "Moderate",
                        "initiation_flowchart_rationale": "start-loss flowchart supports PVS1_Moderate",
                    },
                    "notes": "curated initiation-codon flowchart review",
                    "references": ["ENIGMA Table 4"],
                }
            ],
        )
        criterion = result["manual_criteria"][0]
        self.assertTrue(criterion["applies"])
        self.assertEqual(criterion["selected_strength"], "Moderate")
        self.assertEqual(criterion["points"], 2)

    def test_ps1_splice_adds_pathogenic_points_to_amended_result(self):
        result = evaluate_manual_evidence(
            [],
            [
                {
                    "code": "PS1_SPLICE",
                    "enabled": True,
                    "evidence": {
                        "reference_variant": "BRCA1 c.4185G>A",
                        "reference_classification": "Pathogenic",
                        "reference_classification_source": "ENIGMA ST2/ST7",
                        "same_splice_event_confirmed": True,
                        "vua_splice_event": "exon 12 deletion",
                        "reference_splice_event": "exon 12 deletion",
                        "prediction_strength_comparison": "stronger",
                        "curated_strength": "Strong",
                        "ps1_splice_rationale": "same splice event and stronger predicted splice impact",
                    },
                    "notes": "curated PS1(splice) review",
                    "references": ["ENIGMA Supplementary Table 2 row 67"],
                }
            ],
        )
        criterion = result["manual_criteria"][0]
        self.assertTrue(criterion["applies"])
        self.assertEqual(criterion["selected_strength"], "Strong")
        self.assertEqual(criterion["points"], 4)

    def test_ps1_protein_manual_review_adds_ps1_once(self):
        evidence = {
            "reference_variant": "BRCA1 c.122A>G",
            "reference_p_notation": "p.(His41Arg)",
            "reference_classification": "Pathogenic",
            "classification_verification": "external_vcep_assertion",
            "classification_source": "ClinGen ENIGMA expert panel assertion",
            "same_missense_confirmed": True,
            "different_nucleotide_change_confirmed": True,
            "vua_spliceai_score": 0.01,
            "reference_spliceai_score": 0.02,
            "splice_source_check_completed": True,
            "splice_sources_checked": [
                "ENIGMA Specifications Table 9 v1.2",
                "ENIGMA Supplementary Table 2 v1.2",
            ],
            "vua_confirmed_splice_status": "none_identified",
            "reference_confirmed_splice_status": "none_identified",
            "reference_classification_used_ps1": "unknown",
            "ps1_protein_rationale": "Complete ENIGMA protein PS1 review.",
        }
        result = evaluate_manual_evidence(
            [],
            [{
                "code": "PS1_PROTEIN",
                "enabled": True,
                "evidence": evidence,
                "notes": "Reviewed protein PS1 candidate",
                "references": ["ClinGen assertion"],
            }],
        )
        criterion = result["manual_criteria"][0]
        self.assertTrue(criterion["applies"])
        self.assertEqual(criterion["selected_strength"], "Strong")
        self.assertEqual(result["total_points"], 4)

        with self.assertRaisesRegex(ValueError, "cannot be counted twice"):
            evaluate_manual_evidence(
                [{"name": "PS1", "applies": True, "strength": "Strong", "points": 4}],
                [{"code": "PS1_PROTEIN", "enabled": True, "evidence": evidence}],
            )

    def test_ps1_protein_manual_review_requires_every_defined_splice_source(self):
        evidence = {
            "reference_variant": "BRCA1 c.122A>G",
            "reference_p_notation": "p.(His41Arg)",
            "reference_classification": "Pathogenic",
            "classification_verification": "external_vcep_assertion",
            "classification_source": "ClinGen ENIGMA expert panel assertion",
            "same_missense_confirmed": True,
            "different_nucleotide_change_confirmed": True,
            "vua_spliceai_score": 0.01,
            "reference_spliceai_score": 0.02,
            "splice_source_check_completed": True,
            "splice_sources_checked": ["ENIGMA Specifications Table 9 v1.2"],
            "vua_confirmed_splice_status": "none_identified",
            "reference_confirmed_splice_status": "none_identified",
            "reference_classification_used_ps1": "unknown",
            "ps1_protein_rationale": "Incomplete source review.",
        }
        self.assertIsNone(suggest_strength("PS1_PROTEIN", evidence))

    def test_rna_override_without_complete_record_is_rejected(self):
        with self.assertRaises(ValueError):
            evaluate_manual_evidence(
                [],
                [
                    {
                        "code": "PVS1_RNA",
                        "enabled": True,
                        "evidence": {},
                        "override_strength": "Very Strong",
                    }
                ],
            )

    def test_pvs1_init_override_without_complete_record_is_rejected(self):
        with self.assertRaises(ValueError):
            evaluate_manual_evidence(
                [],
                [
                    {
                        "code": "PVS1_INIT",
                        "enabled": True,
                        "evidence": {},
                        "override_strength": "Moderate",
                    }
                ],
            )

    def test_ps1_splice_override_without_complete_record_is_rejected(self):
        with self.assertRaises(ValueError):
            evaluate_manual_evidence(
                [],
                [
                    {
                        "code": "PS1_SPLICE",
                        "enabled": True,
                        "evidence": {},
                        "override_strength": "Strong",
                    }
                ],
            )


class ClinVarReviewStarTests(unittest.TestCase):
    def test_official_review_status_mapping(self):
        self.assertEqual(clinvar_review_stars("practice guideline"), 4)
        self.assertEqual(clinvar_review_stars("reviewed by expert panel"), 3)
        self.assertEqual(
            clinvar_review_stars(
                "criteria provided, multiple submitters, no conflicts"
            ),
            2,
        )
        self.assertEqual(
            clinvar_review_stars("criteria provided, single submitter"), 1
        )
        self.assertEqual(clinvar_review_stars("no assertion criteria provided"), 0)


if __name__ == "__main__":
    unittest.main()
