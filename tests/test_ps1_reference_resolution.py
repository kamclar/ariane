import unittest

from backend.services.ps1_reference_resolution import (
    Ps1ReferenceDependencies,
    resolve_ps1_reference,
)


class Ps1ReferenceResolutionTests(unittest.IsolatedAsyncioTestCase):
    @staticmethod
    def dependencies(*, clinvar, clingen=None):
        scores = {
            "c.5217T>A": 0.03,
            "c.5217T>G": 0.00,
            "c.5216A>T": 0.01,
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
        self.assertEqual(
            result["classification_verification"], "historical_classification_only"
        )
        self.assertFalse(result["objective_ps1_checks_pass"])

    async def test_enigma_expert_panel_assertion_is_recognized(self):
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

        self.assertEqual(result["classification"], "Likely Pathogenic")
        self.assertEqual(result["classification_verification"], "external_vcep_assertion")
        self.assertTrue(result["objective_ps1_checks_pass"])

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


if __name__ == "__main__":
    unittest.main()
