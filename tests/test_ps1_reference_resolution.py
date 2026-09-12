import unittest

from backend.services.ps1_reference_resolution import (
    Ps1ReferenceDependencies,
    resolve_ps1_reference,
)


class Ps1ReferenceResolutionTests(unittest.IsolatedAsyncioTestCase):
    @staticmethod
    def dependencies(*, clinvar, clingen=None, erepo_registry=None, registry=None):
        scores = {
            "c.5217T>A": 0.03,
            "c.5217T>G": 0.00,
            "c.5216A>T": 0.01,
            "c.123A>G": 0.01,
            "c.122A>G": 0.01,
            "c.131_132delinsCT": 0.01,
            "c.130T>A": 0.01,
        }
        statuses = {
            key: {"status": "ok", "reason": "test score"}
            for key in scores
        }
        return Ps1ReferenceDependencies(
            spliceai_lookup=lambda gene, c_notation: scores.get(c_notation),
            spliceai_status=lambda gene, c_notation: statuses.get(c_notation, {}),
            clinvar_lookup=lambda gene, c_notation: clinvar,
            clingen_lookup=lambda gene, c_notation: clingen or {"status": "not_found"},
            erepo_registry_lookup=lambda gene, c_notation: (
                erepo_registry or {"status": "not_found"}
            ),
            registry_lookup=lambda gene, c_notation: registry,
        )

    async def test_one_star_clinvar_prefills_facts_but_is_not_vcep_verified(self):
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.5217T>A",
            "c.5217T>G",
            dependencies=self.dependencies(
                clinvar={
                    "status": "ok",
                    "variation_id": "55467",
                    "accession": "VCV000055467.8",
                    "aggregate": {
                        "classification": "Likely pathogenic",
                        "review_status": "criteria provided, single submitter",
                    },
                }
            ),
        )

        self.assertEqual(result["reference"]["p_notation"], "p.(Asp1739Glu)")
        self.assertTrue(result["same_missense_substitution"])
        self.assertTrue(result["different_nucleotide_change"])
        self.assertEqual(result["assessed"]["spliceai_score"], 0.03)
        self.assertEqual(result["reference"]["spliceai_score"], 0.00)
        self.assertEqual(result["clinvar_stars"], 1)
        self.assertEqual(result["classification"], "")
        self.assertEqual(result["classification_verification"], "unresolved")
        self.assertEqual(result["classification_source"], "")
        self.assertFalse(result["objective_ps1_checks_pass"])

    async def test_clinvar_enigma_expert_panel_assertion_is_display_only(self):
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.5217T>A",
            "c.5217T>G",
            dependencies=self.dependencies(
                clinvar={
                    "status": "ok",
                    "variation_id": "55467",
                    "aggregate": {
                        "classification": "Likely pathogenic",
                        "review_status": "reviewed by expert panel",
                    },
                    "enigma_submission": {
                        "class": "Likely pathogenic",
                        "scv": "SCV000000001",
                    },
                }
            ),
        )

        self.assertEqual(result["classification"], "")
        self.assertEqual(result["classification_verification"], "unresolved")
        self.assertFalse(result["objective_ps1_checks_pass"])
        self.assertIn("three-star", result["historical_expert_panel_warning"])

    async def test_current_local_erepo_v12_assertion_is_recognized(self):
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.5217T>A",
            "c.5217T>G",
            dependencies=self.dependencies(
                clinvar={"status": "not_found"},
                erepo_registry={
                    "status": "current_vcep_assertion",
                    "record": {
                        "uuid": "11111111-2222-3333-4444-555555555555",
                        "classification": "Likely Pathogenic",
                        "assertion_method_version": "1.2.0",
                        "erepo_url": "https://erepo.example/current",
                    },
                },
            ),
        )

        self.assertEqual(result["classification"], "Likely Pathogenic")
        self.assertEqual(result["classification_verification"], "external_vcep_assertion")
        self.assertEqual(result["erepo_registry_status"], "current_vcep_assertion")
        self.assertTrue(result["objective_ps1_checks_pass"])

    async def test_historical_erepo_assertion_is_not_used_automatically(self):
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.5217T>A",
            "c.5217T>G",
            dependencies=self.dependencies(
                clinvar={"status": "not_found"},
                erepo_registry={
                    "status": "historical_vcep_assertion",
                    "record": {
                        "uuid": "11111111-2222-3333-4444-555555555555",
                        "classification": "Likely Pathogenic",
                        "assertion_method_version": "1.0.0",
                        "erepo_url": "https://erepo.example/historical",
                    },
                },
            ),
        )

        self.assertEqual(result["classification"], "")
        self.assertEqual(result["classification_verification"], "unresolved")
        self.assertFalse(result["objective_ps1_checks_pass"])
        self.assertIn("not the active v1.2", result["historical_expert_panel_warning"])

    async def test_different_protein_consequence_is_reported_not_guessed(self):
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.5217T>A",
            "c.5216A>T",
            dependencies=self.dependencies(clinvar={"status": "not_found"}),
        )

        self.assertEqual(result["reference"]["p_notation"], "p.(Asp1739Val)")
        self.assertFalse(result["same_missense_substitution"])
        self.assertFalse(result["objective_ps1_checks_pass"])

    async def test_external_lookup_failure_is_reported_without_losing_sequence_facts(self):
        def unavailable(gene, c_notation):
            raise TimeoutError("test timeout")

        dependencies = self.dependencies(clinvar={"status": "not_found"})
        dependencies = Ps1ReferenceDependencies(
            spliceai_lookup=dependencies.spliceai_lookup,
            spliceai_status=dependencies.spliceai_status,
            clinvar_lookup=unavailable,
            clingen_lookup=unavailable,
            erepo_registry_lookup=lambda gene, c_notation: {"status": "not_found"},
            registry_lookup=lambda gene, c_notation: None,
        )
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.5217T>A",
            "c.5217T>G",
            dependencies=dependencies,
        )

        self.assertEqual(result["reference"]["p_notation"], "p.(Asp1739Glu)")
        self.assertEqual(result["clinvar_status"], "api_error")
        self.assertEqual(result["clingen_status"], "api_error")
        self.assertIn("TimeoutError", result["clinvar_error"])
        self.assertIn("ClinVar was unavailable", result["review_message"])

    async def test_st7_reference_is_an_accepted_prefilled_classification_basis(self):
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.5217T>A",
            "c.5217T>G",
            dependencies=self.dependencies(
                clinvar={"status": "not_found"},
                registry={
                    "classification": "Pathogenic",
                    "classification_basis": "enigma_st7_v1_2_reference_set",
                    "classification_source": "Parsons et al. 2019",
                    "reference_status": "approved",
                    "reference_splice_evidence_status": "normal",
                    "reference_splice_sources_checked": [
                        "ENIGMA Specifications Table 9 v1.2",
                        "ENIGMA Supplementary Table 2 v1.2",
                    ],
                    "classification_ps1_dependency_used": False,
                },
            ),
        )

        self.assertEqual(result["classification"], "Pathogenic")
        self.assertEqual(
            result["classification_verification"],
            "enigma_st7_v1_2_reference_set",
        )
        self.assertEqual(result["reference_classification_used_ps1"], "no")
        self.assertTrue(result["objective_ps1_checks_pass"])

    async def test_dna_delins_missense_can_use_a_separately_approved_ps1_reference(self):
        result = await resolve_ps1_reference(
            "BRCA1",
            "c.131_132delinsCT",
            "c.130T>A",
            dependencies=self.dependencies(
                clinvar={"status": "not_found"},
                registry={
                    "classification": "Pathogenic",
                    "classification_basis": "enigma_st7_v1_2_reference_set",
                    "classification_source": "ENIGMA Supplementary Table 7 v1.2",
                    "reference_status": "approved",
                    "reference_splice_evidence_status": "normal",
                    "reference_splice_sources_checked": [
                        "ENIGMA Specifications Table 9 v1.2",
                        "ENIGMA Supplementary Table 2 v1.2",
                    ],
                    "classification_ps1_dependency_used": False,
                },
            ),
        )

        self.assertEqual(result["assessed"]["p_notation"], "p.(Cys44Ser)")
        self.assertEqual(result["reference"]["p_notation"], "p.(Cys44Ser)")
        self.assertTrue(result["same_missense_substitution"])
        self.assertTrue(result["different_nucleotide_change"])
        self.assertTrue(result["objective_ps1_checks_pass"])


if __name__ == "__main__":
    unittest.main()
