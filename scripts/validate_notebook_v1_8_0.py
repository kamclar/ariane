import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "BRCA_ACMG_Criteria_Module1_v1_8_0.ipynb"


def main():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    cells = notebook["cells"]

    for index, cell in enumerate(cells):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        compile(source, f"{NOTEBOOK.name}:cell-{index}", "exec")

    text = "\n".join("".join(cell.get("source", [])) for cell in cells)
    required = [
        "def split_combined_hgvs(",
        "def infer_variant_type(",
        "def evaluate_ps1(",
        "def evaluate_pp4_bp5(",
        "def classify_by_enigma_combination(",
        '\"RING\": (2, 101)',
        "c_notation=c_notation,",
        "## Summary - Version 1.8.0",
        "PVS1 (RNA) is outside the automated Module 1 scope",
        'variant_type.lower() == "intronic"',
    ]
    for marker in required:
        if marker not in text:
            raise AssertionError(f"Missing notebook marker: {marker}")

    parity_cell = next(
        "".join(cell.get("source", []))
        for cell in cells
        if "v1.8.0 application parity helpers" in "".join(cell.get("source", []))
    )
    namespace = {"PROJECT_ROOT": ROOT / "backend", "Path": Path}
    exec(parity_cell, namespace)

    c_notation, p_notation = namespace["split_combined_hgvs"](
        "c.6147_6149del (p.Val2050del)"
    )
    assert c_notation == "c.6147_6149del"
    assert p_notation == "p.(Val2050del)"
    assert namespace["infer_variant_type"](c_notation, p_notation) == "inframe_deletion"

    bs3_only = {"BS3": {"strength": "Strong", "points": -4}}
    assert namespace["classify_by_enigma_combination"](bs3_only, -4)[:2] == (3, "VUS")

    pvs1_only = {"PVS1": {"strength": "Very Strong", "points": 8}}
    pvs1_result = namespace["classify_by_enigma_combination"](pvs1_only, 8)
    assert pvs1_result[:2] == (3, "VUS")
    assert "requires at least one additional Supporting" in pvs1_result[2]

    bp1_only = {"BP1": {"strength": "Strong", "points": -4}}
    assert namespace["classify_by_enigma_combination"](bp1_only, -4)[:2] == (
        2,
        "Likely Benign",
    )

    # Execute the notebook orchestrator with deterministic local stubs.
    domain_cell = next(
        "".join(cell.get("source", []))
        for cell in cells
        if "FUNCTIONAL_DOMAINS =" in "".join(cell.get("source", []))
    )
    exec(domain_cell, namespace)

    for marker in ("def evaluate_bp1(", "BAYESDEL_THRESHOLDS =", "def evaluate_bp7("):
        cell_source = next(
            "".join(cell.get("source", []))
            for cell in cells
            if marker in "".join(cell.get("source", []))
        )
        exec(cell_source, namespace)

    orchestrator_cell = next(
        "".join(cell.get("source", []))
        for cell in cells
        if "def evaluate_variant(variant:" in "".join(cell.get("source", []))
    )
    orchestrator_function = orchestrator_cell.split(
        "# ============================================================\n"
        "# Run evaluation on all test variants",
        1,
    )[0]

    def empty_pvs1(*args, **kwargs):
        return {
            "applies": False,
            "strength": None,
            "points": 0,
            "reason": "",
            "requires_rna": False,
            "pm5_code": None,
            "pm5_strength": None,
            "pm5_points": 0,
        }

    namespace.update(
        {
            "gnomad_status_summary": lambda data: "not queried",
            "evaluate_frequency_criteria": lambda data, variant_type: {},
            "get_bayesdel_score": lambda gene, c_notation: None,
            "evaluate_pvs1": empty_pvs1,
            "normalize_variant_type": lambda variant_type: (variant_type or "").strip().lower(),
            "variant_type_allows_spliceai_pp3": lambda variant_type: variant_type in {
                "synonymous", "silent", "missense",
                "inframe_deletion", "inframe_insertion", "inframe_delins", "delins",
                "intronic",
            },
            "spliceai_predicts_splice_effect": lambda score: score is not None and score >= 0.2,
            "table9_lookup_ps3_bs3": lambda gene, c_notation: (
                {"applies": True, "code": "BS3", "strength": "Strong", "points": -4, "reason": "Table 9"}
                if (gene, c_notation) == ("BRCA1", "c.3891_3893del")
                else {"applies": False, "code": None, "strength": None, "points": 0, "reason": ""}
            ),
        }
    )
    exec(orchestrator_function, namespace)

    namespace["get_spliceai_score"] = lambda gene, c_notation: None
    ser1298 = namespace["evaluate_variant"](
        {
            "gene": "BRCA1",
            "c_notation": "c.3891_3893del",
            "p_notation": "p.(Ser1298del)",
        },
        use_gnomad_cache=False,
    )
    assert ser1298["predicted_class"] == 3
    assert ser1298["total_points"] == -4

    namespace["get_spliceai_score"] = lambda gene, c_notation: 0.05
    val2050 = namespace["evaluate_variant"](
        {
            "gene": "BRCA2",
            "c_notation": "c.6147_6149del (p.Val2050del)",
        },
        use_gnomad_cache=False,
    )
    assert val2050["predicted_class"] == 2
    assert val2050["total_points"] == -4
    assert val2050["criteria"]["BP1"]["strength"] == "Strong"

    print(f"notebook_json_ok cells={len(cells)}")
    print("all_code_cells_compile")
    print("parity_helpers_ok")
    print("notebook_orchestrator_smoke_ok")


if __name__ == "__main__":
    main()
