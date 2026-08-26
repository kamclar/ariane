import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from backend.modules.bp1 import evaluate_bp1
from backend.modules.bp7 import evaluate_bp7
from backend.modules.enigma_rules import (
    clinical_annotations_for_variant,
    get_decision_tree,
    public_catalog,
    search_reference_table,
    search_table9,
    validate_rule_catalog,
)
from backend.modules.classifier import evaluate_variant
from backend.modules.pvs1_rna import evaluate_pvs1_rna
from backend.modules.pvs1 import evaluate_pvs1
from backend.modules.table9 import table9_lookup_ps3_bs3
from backend.modules.pp3_bp4 import evaluate_pp3_bp4


ROOT = Path(__file__).resolve().parents[1]


def test_catalog_is_valid_and_exposes_no_local_paths():
    validate_rule_catalog()
    payload = public_catalog()
    serialized = json.dumps(payload)

    assert payload["criteria_specification_version"] == "1.2.0"
    assert {source["id"] for source in payload["sources"]} >= {
        "enigma-v1.2-specifications",
        "enigma-v1.2-table9",
    }
    assert "local_path" not in serialized
    assert "file_path" not in serialized
    assert "F:\\" not in serialized
    assert "/home/" not in serialized
    assert all(source["official_url"].startswith("https://") for source in payload["sources"])
    assert all(len(source["sha256"]) == 64 for source in payload["sources"])
    assert len(payload["figures"]) == 14
    assert payload["tables"]["table_count"] == 42
    assert {
        group["id"]: sum(
            table["group"] == group["id"]
            for table in payload["tables"]["items"]
        )
        for group in payload["tables"]["groups"]
    } == {"specification": 9, "appendix": 17, "supplementary": 16}
    assert {
        role["id"]: role["table_count"] for role in payload["tables"]["roles"]
    } == {
        "used_by_ariane": 15,
        "expert_review": 4,
        "supporting_reference": 23,
    }
    assert all(table["role"] for table in payload["tables"]["items"])
    used_tables = [
        table for table in payload["tables"]["items"]
        if table["role"] == "used_by_ariane"
    ]
    assert len(used_tables) == 15
    assert all(table["consumer"] for table in used_tables)
    assert all(
        table["consumer"] is None
        for table in payload["tables"]["items"]
        if table["role"] != "used_by_ariane"
    )
    for figure in payload["figures"]:
        assert figure["asset_url"].startswith("/static/enigma/")
        assert (
            ROOT / "frontend" / "static" / figure["asset_url"].removeprefix("/static/")
        ).is_file()


def test_figure_1a_is_a_closed_and_referentially_valid_graph():
    tree = get_decision_tree("figure-1a")
    assert tree is not None
    assert {branch["id"] for branch in tree["branches"]} == {
        "missense-inframe", "synonymous", "intronic"
    }
    condition_ids = {item["id"] for item in tree["condition_definitions"]}
    for branch in tree["branches"]:
        nodes = {node["id"]: node for node in branch["nodes"]}
        assert branch["entry_node"] in nodes
        assert nodes[branch["entry_node"]]["kind"] == "entry"
        assert any(
            edge["from"] == branch["entry_node"] and edge["result"] == "start"
            for edge in branch["edges"]
        )
        for node in nodes.values():
            if node["kind"] == "decision":
                assert node["condition_id"] in condition_ids
        for edge in branch["edges"]:
            assert edge["from"] in nodes
            assert edge["to"] in nodes
        assert set(branch["layout"]["positions"]) == set(nodes)

    missense = next(branch for branch in tree["branches"] if branch["id"] == "missense-inframe")
    assert {node["id"] for node in missense["nodes"]} >= {
        "mi-protein-pp3-unknown",
        "mi-protein-pp3-low",
        "mi-no-code-protein-unknown",
        "mi-no-code-protein-low",
    }
    assert tree["diagram_provenance"] == "official_redraw"
    assert tree["original_figure_ids"]
    assert any(edge.get("points") for edge in missense["edges"])
    assert any(len(position) == 4 for position in missense["layout"]["positions"].values())


def test_all_official_templates_and_derived_paths_are_explicitly_distinguished():
    payload = public_catalog()
    trees = {tree["id"]: tree for tree in payload["decision_trees"]}

    assert {
        "figure-1a",
        "figure-1b",
        "figure-1c",
        "appendix-pvs1-nonsplice",
        "appendix-pvs1-splice",
        "appendix-figure-9",
        "appendix-exon-ptc-maps",
    } <= set(trees)
    assert {
        "population-evidence-path",
        "ps1-reference-eligibility-path",
        "clinical-lr-path",
        "classification-combination-path",
    } <= set(trees)

    for tree in trees.values():
        if tree["diagram_provenance"] == "official_redraw":
            assert tree["original_figure_ids"]
        else:
            assert tree["diagram_provenance"] == "ariane_derived"
            assert tree["original_figure_ids"] == []


def test_table9_search_is_paginated_filtered_and_whitelisted():
    result = search_table9(gene="BRCA1", query="c.509G>A", page=1, page_size=5)
    assert result["total"] == 1
    assert result["page_size"] == 5
    item = result["items"][0]
    assert item["gene"] == "BRCA1"
    assert item["c_notation"] == "c.509G>A"
    assert "result_1" not in item
    assert isinstance(item["results"], list)
    assert not any("path" in key for key in item)


def test_all_enigma_tables_are_paginated_and_searchable_without_source_paths():
    result = search_reference_table(
        "supplementary-table-7",
        query="c.1001C>A",
        page=1,
        page_size=5,
    )
    assert result is not None
    assert result["total"] == 1
    assert result["items"][0]["source_row"] == 5
    assert result["items"][0]["cells"][:4] == [
        "BRCA1", "c.1001C>A", "p.(Pro334His)", "(Likely) Benign",
    ]
    serialized = json.dumps(result)
    assert "local_path" not in serialized
    assert "file_path" not in serialized
    assert "F:\\" not in serialized

    assert search_reference_table("missing-table") is None
    assert search_reference_table(
        "specification-table-4", section_id="missing-sheet"
    ) is None


def test_reduced_penetrance_annotation_comes_from_enigma_appendix_table_11():
    annotations = clinical_annotations_for_variant("BRCA1", "c.5096G>A")
    assert len(annotations) == 1
    annotation = annotations[0]
    assert annotation["label"] == "Variant with reduced penetrance"
    assert "proven reduced penetrance" in annotation["evidence"].lower()
    assert annotation["source"] == "ENIGMA BRCA1/2 VCEP v1.2 Appendix Table 11"
    assert {item["pmid"] for item in annotation["publications"]} == {
        "22889855", "28490613",
    }
    assert annotation["affects_classification"] is False
    assert clinical_annotations_for_variant("BRCA1", "c.509G>A") == []


@pytest.mark.parametrize(
    "result,criterion,branch,outcome",
    [
        (
            evaluate_pp3_bp4("BRCA1", "missense", "p.(Ala1789Val)", 0.35, 0.03)["PP3"],
            "PP3", "missense-inframe", "mi-protein-pp3-low",
        ),
        (
            evaluate_pp3_bp4("BRCA1", "missense", "p.(Ala1789Val)", 0.10, 0.03)["BP4"],
            "BP4", "missense-inframe", "mi-bp4",
        ),
        (
            evaluate_bp1("BRCA1", "missense", "p.(Arg500Gln)", 0.03),
            "BP1", "missense-inframe", "mi-bp1",
        ),
        (
            evaluate_bp7("intronic", 0.03, False, True, "c.548-21A>G", gene="BRCA1"),
            "BP7", "intronic", "int-bp4-bp7",
        ),
    ],
)
def test_applied_figure1a_criteria_carry_the_actual_decision_path(
    result, criterion, branch, outcome
):
    assert result["applies"] is True
    path = result["decision_path"]
    assert path["tree_id"] == "figure-1a"
    assert path["criterion"] == criterion
    assert path["branch_id"] == branch
    assert path["outcome_node"] == outcome
    assert path["steps"]
    assert path["sources"][0]["source_id"] == "enigma-v1.2-specifications"
    assert path["sources"][0]["figure_url"].startswith("/static/enigma/figure-1a-")


def test_frontend_contains_complete_tables_page_and_expandable_decision_path():
    html = (ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
    javascript = (ROOT / "frontend" / "static" / "js" / "app.js").read_text(encoding="utf-8")

    assert "ENIGMA rules" in html
    assert "Show decision path" in html
    assert "ENIGMA tables used by ARIANE are shown first" in html
    assert "Specification Tables" not in html  # supplied by the validated catalogue
    assert "referenceTableRoles()" in html
    assert "referenceTablesForRole()" in html
    assert "Original ENIGMA figures" in html
    assert "decisionPathSvg(c.decision_path)" in html
    assert "Selected variant path" in html
    assert "Official ENIGMA figure redrawn" in html
    assert "ARIANE decision path derived from ENIGMA rules" in html
    assert "selectRuleTree(selectedRuleTreeId)" in html
    assert "/api/rules/trees/${encodeURIComponent(selectedId)}" in javascript
    assert "/api/rules/tables/table9" in javascript
    assert "/api/rules/tables/${encodeURIComponent(this.selectedReferenceTableId)}" in javascript
    assert "decisionTreeSvg(branch, path = null)" in javascript
    assert "const observedAnnotations = Object.fromEntries" in javascript
    assert "this.wrapGraphText(value, maximumCharacters)" in javascript
    assert "this.compactGraphText(observed[node.id]" not in javascript
    assert "this.compactGraphText(node.observed" not in javascript
    assert "observedLinesByNode" in javascript
    assert "tree-edge-label-active" in javascript
    assert "tree-has-active-path" in javascript
    assert "edgeDisplayLabel(value)" in javascript
    assert "edgeLabelLines(value)" in javascript
    assert "graph-node-entry" in javascript


def test_intronic_trace_starts_at_the_figure_1a_entry_node():
    result = evaluate_pp3_bp4(
        "BRCA1",
        "intronic",
        "p.(?)",
        spliceai_score=0.9,
        c_notation="c.548-9A>G",
    )["PP3"]

    assert [item["node_id"] for item in result["decision_path"]["steps"]] == [
        "int-canonical-site",
        "int-splice-impact",
    ]


@pytest.mark.parametrize(
    ("gene", "variant_type", "p_notation", "bayesdel_score", "c_notation"),
    (
        ("BRCA1", "missense", "p.(Ala1789Val)", 0.50, "c.5366C>T"),
        ("BRCA1", "inframe_deletion", "p.(Ser1715del)", 0.50, "c.5143_5145del"),
        ("BRCA1", "synonymous", "p.(Val1653=)", None, "c.4959A>G"),
        ("BRCA1", "intronic", "p.?", None, "c.548-9A>G"),
    ),
)
def test_figure1a_does_not_apply_codes_when_spliceai_is_unavailable(
    gene, variant_type, p_notation, bayesdel_score, c_notation
):
    pp3_bp4 = evaluate_pp3_bp4(
        gene,
        variant_type,
        p_notation,
        bayesdel_score=bayesdel_score,
        spliceai_score=None,
        c_notation=c_notation,
    )
    bp1 = evaluate_bp1(
        gene,
        variant_type,
        p_notation,
        spliceai_score=None,
    )
    bp7 = evaluate_bp7(
        variant_type,
        spliceai_score=None,
        in_domain=True,
        bp4_met=False,
        c_notation=c_notation,
        gene=gene,
    )

    assert not any(item.get("applies") for item in pp3_bp4.values())
    assert bp1["applies"] is False
    assert bp7["applies"] is False


def test_measured_spliceai_not_informative_band_still_uses_official_bayesdel_route():
    result = evaluate_pp3_bp4(
        "BRCA1",
        "missense",
        "p.(Ala1789Val)",
        bayesdel_score=0.50,
        spliceai_score=0.15,
        c_notation="c.5366C>T",
    )["PP3"]

    assert result["applies"] is True
    assert result["decision_path"]["steps"][0]["result"] == "not_informative"
    assert result["decision_path"]["outcome_node"] == "mi-protein-pp3-unknown"


def test_missing_spliceai_produces_explicit_figure1a_unavailable_warning():
    result = evaluate_variant(
        gene="BRCA1",
        variant_type="missense",
        p_notation="p.(Ala1789Val)",
        c_notation="c.5366C>T",
        spliceai_score=None,
        bayesdel_score=0.50,
    )

    assert "PP3" not in result["criteria"]
    warning = next(
        item for item in result["warnings"]
        if "Figure 1A bioinformatic result unavailable" in item
    )
    assert "Missing data was not treated as an ENIGMA prediction band" in warning
    assert "PP3, BP4, BP1 and BP7 were not applied" in warning


def test_rna_and_functional_evidence_link_to_their_official_figures():
    rna = evaluate_pvs1_rna("BRCA1", "c.4185G>A")
    assert rna["applies"] is True
    assert rna["decision_path"]["tree_id"] == "figure-1b"
    assert rna["decision_path"]["outcome_node"] == "rna-other-aberrant"

    functional = table9_lookup_ps3_bs3("BRCA1", "c.509G>A")
    result = evaluate_variant(
        gene="BRCA1",
        variant_type="missense",
        p_notation="p.(Arg170Gln)",
        c_notation="c.509G>A",
        spliceai_score=0.05,
        bayesdel_score=0.10,
        table9_result=functional,
    )
    path = result["criteria"]["BS3"]["decision_path"]
    assert path["tree_id"] == "figure-1c"
    assert path["branch_id"] == "exonic-missense-inframe"
    assert path["outcome_node"] == "func-protein-code"


def test_automated_pvs1_links_to_the_gene_specific_appendix_figure():
    ptc = evaluate_pvs1(
        "BRCA1", "frameshift", "p.(Cys1225SerfsTer10)", "c.3668_3671dup"
    )
    assert ptc["applies"] is True
    assert ptc["decision_path"]["tree_id"] == "appendix-pvs1-nonsplice"
    assert ptc["decision_path"]["branch_id"] == "nonsense-frameshift"
    assert ptc["decision_path"]["outcome_node"] == "pvs-ptc-table4"
    assert ptc["decision_path"]["sources"][0]["location"] == "Figure 3"

    splice = evaluate_pvs1(
        "BRCA2", "splice_site", "p.?", "c.8953+2T>A", spliceai_score=0.9
    )
    assert splice["applies"] is True
    assert splice["decision_path"]["tree_id"] == "appendix-pvs1-splice"
    assert splice["decision_path"]["branch_id"] == "canonical-splice-outcome"
    assert splice["decision_path"]["sources"][0]["location"] == "Figure 6"

    deletion = evaluate_pvs1(
        "BRCA1",
        "exon_deletion",
        "p.?",
        "c.(670+1_671-1)_(4096+1_4097-1)del",
    )
    assert deletion["applies"] is True
    assert deletion["decision_path"]["branch_id"] == "exon-deletion"

    duplication = evaluate_pvs1(
        "BRCA1",
        "exon_duplication",
        "p.?",
        "c.(670+1_671-1)_(4096+1_4097-1)dup",
        dup_type="Tandem",
    )
    assert duplication["applies"] is True
    assert duplication["decision_path"]["branch_id"] == "duplication"


def test_rules_http_api_uses_the_validated_catalog():
    from backend.main import app

    client = TestClient(app)
    catalog = client.get("/api/rules")
    tree = client.get("/api/rules/trees/figure-1a")
    table = client.get(
        "/api/rules/tables/table9",
        params={"gene": "BRCA1", "query": "c.509G>A", "page_size": 5},
    )
    supplementary_table = client.get(
        "/api/rules/tables/supplementary-table-7",
        params={"query": "c.1001C>A", "page_size": 5},
    )

    assert catalog.status_code == 200
    assert tree.status_code == 200
    assert tree.json()["source_location"] == "Figure 1A"
    assert table.status_code == 200
    assert table.json()["items"][0]["c_notation"] == "c.509G>A"
    assert supplementary_table.status_code == 200
    assert supplementary_table.json()["items"][0]["source_row"] == 5
    assert client.get("/api/rules/trees/not-a-tree").status_code == 404
    assert client.get("/api/rules/tables/not-a-table").status_code == 404
    assert client.get("/api/rules/tables/table9", params={"page_size": 101}).status_code == 422

    for tree_summary in catalog.json()["decision_trees"]:
        response = client.get(f"/api/rules/trees/{tree_summary['id']}")
        assert response.status_code == 200
        assert response.json()["diagram_provenance"] in {"official_redraw", "ariane_derived"}
