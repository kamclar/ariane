"""Regression matrix for variant-type routing through ENIGMA Module 1."""

import unittest

from tests.dag_test_support import classify_with_dag as evaluate_variant
from backend.modules.pp4_bp5 import evaluate_pp4_bp5
from backend.modules.table9 import table9_lookup_ps3_bs3
from backend.modules.variant_type import infer_variant_type


class VariantTypeInferenceMatrixTests(unittest.TestCase):
    """Every supported notation family must keep its normalized variant type."""

    CASES = (
        ("nonsense", "c.303T>G", "p.(Tyr101Ter)", "nonsense"),
        ("frameshift duplication", "c.5266dup", "p.(Gln1756ProfsTer74)", "frameshift"),
        ("missense", "c.5366C>T", "p.(Ala1789Val)", "missense"),
        ("synonymous", "c.4185G>A", "p.(Gln1395=)", "synonymous"),
        ("intronic SNV", "c.548-9A>G", "p.(?)", "intronic"),
        ("canonical splice site", "c.8953+2T>C", "p.(?)", "splice_site"),
        ("in-frame deletion", "c.3891_3893del", "p.(Ser1298del)", "inframe_deletion"),
        (
            "exon deletion",
            "c.(793+1_794-1)_(1909+1_1910-1)del",
            "p.(?)",
            "exon_deletion",
        ),
        (
            "exon duplication",
            "c.(80+1_81-1)_(134+1_135-1)dup",
            "p.(?)",
            "exon_duplication",
        ),
        ("initiation codon", "c.1A>G", "p.(Met1Val)", "initiation_codon"),
        ("initiation codon unknown protein", "c.1A>G", "p.(Met1?)", "initiation_codon"),
        ("stop lost", "c.5590T>A", "p.(Ter1864ArgextTer39)", "stop_lost"),
        ("5-prime UTR", "c.-10A>G", "p.(?)", "5utr"),
        ("3-prime UTR", "c.*10A>G", "p.(?)", "3utr"),
    )

    def test_variant_type_inference_matrix(self):
        for label, c_notation, p_notation, expected_type in self.CASES:
            with self.subTest(case=label):
                self.assertEqual(
                    infer_variant_type(c_notation, p_notation),
                    expected_type,
                )


class VariantTypeRuleRoutingMatrixTests(unittest.TestCase):
    """Representative variants must enter only the intended rule branches."""

    CASES = (
        {
            "label": "nonsense PTC",
            "gene": "BRCA1",
            "c": "c.303T>G",
            "p": "p.(Tyr101Ter)",
            "spliceai": 0.90,
            "pp4": True,
            "required": {"PVS1", "PM5_PTC", "PP4"},
            "forbidden": {"PP3", "BP1", "BP4", "BP7"},
            "class": 5,
        },
        {
            "label": "frameshift PTC",
            "gene": "BRCA1",
            "c": "c.5266dup",
            "p": "p.(Gln1756ProfsTer74)",
            "pp4": True,
            "required": {"PVS1", "PM5_PTC"},
            "forbidden": {"PP4", "PP3", "BP1", "BP4", "BP7"},
            "class": 5,
        },
        {
            "label": "missense in functional domain",
            "gene": "BRCA1",
            "c": "c.5366C>T",
            "p": "p.(Ala1789Val)",
            "spliceai": 0.05,
            "bayesdel": 0.50,
            "required": {"PP3"},
            "forbidden": {"PVS1", "PM5_PTC", "BP1", "BP4", "BP7"},
            "class": 3,
        },
        {
            "label": "missense outside functional domains with low SpliceAI",
            "gene": "BRCA1",
            "c": "c.3247A>G",
            "p": "p.(Met1083Val)",
            "spliceai": 0.01,
            "required": {"BP1"},
            "forbidden": {"PVS1", "PM5_PTC", "PP3", "BP4", "BP7"},
            "class": 2,
        },
        {
            "label": "synonymous in functional domain",
            "gene": "BRCA1",
            "c": "c.4188G>A",
            "p": "p.(Gln1396=)",
            "spliceai": 0.05,
            "required": {"BP4", "BP7"},
            "forbidden": {"PVS1", "PM5_PTC", "PP3", "BP1"},
            "class": 2,
        },
        {
            "label": "synonymous with ENIGMA patient RNA evidence",
            "gene": "BRCA1",
            "c": "c.4185G>A",
            "p": "p.(Gln1395=)",
            "spliceai": 0.05,
            "required": {"PVS1_RNA"},
            "forbidden": {"PVS1", "PM5_PTC", "PP3", "BP1", "BP4", "BP7"},
            "class": 3,
        },
        {
            "label": "synonymous outside functional domains",
            "gene": "BRCA1",
            "c": "c.306A>C",
            "p": "p.(Ala102=)",
            "spliceai": 0.05,
            "required": {"BP1"},
            "forbidden": {"PVS1", "PM5_PTC", "PP3", "BP4", "BP7"},
            "class": 2,
        },
        {
            "label": "intronic predicted splice effect",
            "gene": "BRCA1",
            "c": "c.548-9A>G",
            "p": "p.(?)",
            "spliceai": 0.86,
            "required": {"PP3"},
            "forbidden": {"PVS1", "PM5_PTC", "BP1", "BP4", "BP7"},
            "class": 3,
        },
        {
            "label": "canonical splice site requiring RNA review",
            "gene": "BRCA1",
            "c": "c.212+1G>T",
            "p": "p.(?)",
            "spliceai": 0.90,
            "required": set(),
            "forbidden": {"PVS1", "PP3", "BP1", "BP4", "BP7"},
            "rna_review": True,
            "class": 3,
        },
        {
            "label": "Table 9 in-frame deletion",
            "gene": "BRCA1",
            "c": "c.3891_3893del",
            "p": "p.(Ser1298del)",
            "spliceai": 0.15,
            "table9": True,
            "required": {"BS3"},
            "forbidden": {"PVS1", "PM5_PTC", "PP3", "BP1", "BP4", "BP7"},
            "class": 3,
        },
        {
            "label": "confirmed tandem exon duplication",
            "gene": "BRCA1",
            "c": "c.(80+1_81-1)_(134+1_135-1)dup",
            "p": "p.(?)",
            "dup_type": "Tandem",
            "required": {"PVS1"},
            "forbidden": {"PM5_PTC", "PP3", "BP1", "BP4", "BP7"},
            "class": 3,
        },
        {
            "label": "initiation codon remains manual",
            "gene": "BRCA1",
            "c": "c.1A>G",
            "p": "p.(Met1?)",
            "spliceai": 0.01,
            "required": set(),
            "forbidden": {"PVS1", "PM5_PTC", "PP3", "BP1", "BP4", "BP7"},
            "initiation_review": True,
            "class": 3,
        },
        {
            "label": "stop-loss extension is not a nonsense PTC",
            "gene": "BRCA1",
            "c": "c.5590T>A",
            "p": "p.(Ter1864ArgextTer39)",
            "required": set(),
            "forbidden": {"PVS1", "PM5_PTC", "PP3", "BP1", "BP4", "BP7"},
            "class": 3,
        },
    )

    @staticmethod
    def _applied_codes(result):
        return {
            code
            for code, criterion in result["criteria"].items()
            if criterion.get("applies", True)
        }

    def test_rule_routing_matrix(self):
        for case in self.CASES:
            with self.subTest(case=case["label"]):
                table9_result = (
                    table9_lookup_ps3_bs3(case["gene"], case["c"])
                    if case.get("table9")
                    else None
                )
                pp4_result = (
                    evaluate_pp4_bp5(case["gene"], case["c"])
                    if case.get("pp4")
                    else None
                )
                result = evaluate_variant(
                    gene=case["gene"],
                    variant_type=infer_variant_type(case["c"], case["p"]),
                    p_notation=case["p"],
                    c_notation=case["c"],
                    spliceai_score=case.get("spliceai"),
                    bayesdel_score=case.get("bayesdel"),
                    table9_result=table9_result,
                    pp4_bp5_result=pp4_result,
                    dup_type=case.get("dup_type", "Unknown"),
                )
                applied = self._applied_codes(result)
                self.assertTrue(case["required"].issubset(applied), applied)
                self.assertTrue(case["forbidden"].isdisjoint(applied), applied)
                self.assertEqual(result["predicted_class"], case["class"])
                if case.get("rna_review"):
                    self.assertTrue(result["rna_review"]["recommended"])
                if case.get("initiation_review"):
                    self.assertTrue(result["initiation_review"]["recommended"])


if __name__ == "__main__":
    unittest.main()
