import unittest

from backend.modules.criterion_order import criterion_sort_key, sorted_criterion_items


class CriterionOrderTests(unittest.TestCase):
    def test_guideline_group_order_and_internal_qualifiers(self):
        codes = [
            "BP7_RNA",
            "BS1_Supporting",
            "BA1",
            "PP4",
            "PM5_PTC",
            "PM2_Supporting",
            "PS3",
            "PVS1",
        ]
        self.assertEqual(
            sorted(codes, key=criterion_sort_key),
            [
                "PVS1",
                "PS3",
                "PM2_Supporting",
                "PM5_PTC",
                "PP4",
                "BA1",
                "BS1_Supporting",
                "BP7_RNA",
            ],
        )

    def test_base_code_precedes_qualified_code(self):
        criteria = {
            "PS1_PROTEIN": {},
            "PS1": {},
            "PVS1_RNA": {},
            "PVS1": {},
        }
        self.assertEqual(
            [code for code, _ in sorted_criterion_items(criteria)],
            ["PVS1", "PVS1_RNA", "PS1", "PS1_PROTEIN"],
        )


if __name__ == "__main__":
    unittest.main()
