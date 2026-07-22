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
                {"p_value": 0.05, "odds_ratio": 4, "lower_ci": 2.01},
            ),
            "Strong",
        )
        self.assertIsNone(
            suggest_strength(
                "PS4",
                {"p_value": 0.05, "odds_ratio": 4, "lower_ci": 2.0},
            )
        )

    def test_pm3_and_bs2_point_thresholds(self):
        for code in ("PM3", "BS2"):
            self.assertEqual(
                suggest_strength(code, {"evidence_points": 1}), "Supporting"
            )
            self.assertEqual(
                suggest_strength(code, {"evidence_points": 2}), "Moderate"
            )
            self.assertEqual(
                suggest_strength(code, {"evidence_points": 4}), "Strong"
            )

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
        self.assertEqual(
            suggest_strength("PP1", {"likelihood_ratio": 350}), "Very Strong"
        )

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
        evidence["assay_scope"] = "combined_mrna_protein"
        self.assertIsNone(suggest_strength("PVS1_RNA", evidence))

    def test_bp7_rna_is_strong_only_with_eligibility_stipulation(self):
        evidence = {
            "assay_scope": "mrna_only",
            "rna_conclusion": "no_damaging_effect",
            "bp7_rna_eligible": True,
            "transcript_accession": "NM_000059.4",
            "tissue_or_cell_type": "blood",
            "nmd_assessed": "not_applicable",
        }
        self.assertEqual(suggest_strength("BP7_RNA", evidence), "Strong")
        evidence["bp7_rna_eligible"] = False
        self.assertIsNone(suggest_strength("BP7_RNA", evidence))

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


class ManualEvidenceClassificationTests(unittest.TestCase):
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
                "evidence": {"evidence_points": 1},
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

    def test_reviewer_can_override_with_an_allowed_strength(self):
        result = evaluate_manual_evidence(
            [],
            [
                {
                    "code": "PP1",
                    "enabled": True,
                    "evidence": {"likelihood_ratio": 4.3},
                    "override_strength": "Supporting",
                    "notes": "conservative reviewer adjustment",
                    "references": ["PMID:3"],
                }
            ],
        )
        criterion = result["manual_criteria"][0]
        self.assertEqual(criterion["suggested_strength"], "Moderate")
        self.assertEqual(criterion["selected_strength"], "Supporting")
        self.assertTrue(criterion["overridden"])

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
                        "bp7_rna_eligible": True,
                        "transcript_accession": "NM_000059.4",
                        "tissue_or_cell_type": "blood",
                        "nmd_assessed": "not_applicable",
                        "transcript_result_summary": "normal transcript profile within assay sensitivity",
                    },
                    "notes": "curated RNA assay",
                    "references": ["PMID:5"],
                }
            ],
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
