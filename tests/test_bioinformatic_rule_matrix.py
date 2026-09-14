"""Regression matrix for ENIGMA Figure 1A and its RNA interactions."""

import pytest

from backend.criteria.bp1 import evaluate_bp1
from backend.criteria.bp7 import evaluate_bp7
from backend.criteria.evidence_interactions import (
    apply_rna_interactions,
    automatic_functional_interactions,
    rna_interaction_codes,
)
from backend.review.manual_evidence import suggest_strength
from backend.criteria.pp3_bp4 import evaluate_pp3_bp4
from backend.criteria.pvs1_rna import evaluate_pvs1_rna
from backend.reference_data.table9 import table9_lookup_ps3_bs3
from backend.variant_processing.variant_type import infer_variant_type


def applied_codes(results):
    return {
        code
        for code, result in results.items()
        if isinstance(result, dict) and result.get("applies")
    }


@pytest.mark.parametrize(
    "variant_type",
    (
        "nonsense",
        "frameshift",
        "splice_site",
        "initiation_codon",
        "stop_lost",
        "exon_deletion",
        "exon_duplication",
        "deletion",
        "insertion",
        "duplication",
        "delins",
        "unknown",
    ),
)
def test_unsupported_variant_types_cannot_enter_figure1a(variant_type):
    pp3_bp4 = evaluate_pp3_bp4(
        "BRCA1",
        variant_type,
        "p.(Cys61Gly)",
        bayesdel_score=0.90,
        spliceai_score=0.90,
        c_notation="c.181T>G",
    )
    bp1 = evaluate_bp1(
        "BRCA1",
        variant_type,
        "p.(Ala102Gly)",
        spliceai_score=0.01,
    )
    bp7 = evaluate_bp7(
        variant_type,
        spliceai_score=0.01,
        in_domain=True,
        bp4_met=True,
        c_notation="c.306A>G",
        gene="BRCA1",
    )

    assert not applied_codes(pp3_bp4)
    assert not bp1["applies"]
    assert not bp7["applies"]


@pytest.mark.parametrize("c_notation", ("c.212+1G>T", "c.212+2T>C", "c.213-1G>A", "c.213-2A>G"))
def test_canonical_splice_positions_fail_closed_if_mislabeled_intronic(c_notation):
    damaging = evaluate_pp3_bp4(
        "BRCA1",
        "intronic",
        "p.(?)",
        spliceai_score=0.90,
        c_notation=c_notation,
    )
    benign = evaluate_pp3_bp4(
        "BRCA1",
        "intronic",
        "p.(?)",
        spliceai_score=0.01,
        c_notation=c_notation,
    )

    assert "PP3" not in applied_codes(damaging)
    assert "BP4" not in applied_codes(benign)


@pytest.mark.parametrize(
    ("c_notation", "bp7_expected"),
    (
        ("c.500+6A>G", False),
        ("c.500+7A>G", True),
        ("c.501-20A>G", False),
        ("c.501-21A>G", True),
    ),
)
def test_intronic_bp7_motif_boundaries_match_figure1a(c_notation, bp7_expected):
    bioinformatic = evaluate_pp3_bp4(
        "BRCA1",
        "intronic",
        "p.(?)",
        spliceai_score=0.10,
        c_notation=c_notation,
    )
    assert applied_codes(bioinformatic) == {"BP4"}

    bp7 = evaluate_bp7(
        "intronic",
        spliceai_score=0.10,
        bp4_met=True,
        c_notation=c_notation,
        gene="BRCA1",
    )
    assert bp7["applies"] is bp7_expected


@pytest.mark.parametrize(
    ("spliceai_score", "expected"),
    (
        (0.10, {"BP4"}),
        (0.100001, set()),
        (0.199999, set()),
        (0.20, {"PP3"}),
    ),
)
def test_spliceai_boundaries_for_synonymous_domain_variant(spliceai_score, expected):
    result = evaluate_pp3_bp4(
        "BRCA1",
        "synonymous",
        "p.(Gln1396=)",
        spliceai_score=spliceai_score,
        c_notation="c.4188G>A",
    )
    assert applied_codes(result) == expected


def test_inframe_and_unresolved_indel_types_are_distinguished():
    assert infer_variant_type("c.3891_3893del", "p.(Ser1298del)") == "inframe_deletion"
    assert infer_variant_type("c.3891_3893del", "p.(?)") == "deletion"
    assert infer_variant_type("c.922_923delinsGC", "p.(?)") == "delins"

    unresolved = evaluate_pp3_bp4(
        "BRCA1",
        "delins",
        "p.(?)",
        bayesdel_score=0.90,
        spliceai_score=0.90,
        c_notation="c.922_923delinsGC",
    )
    assert not applied_codes(unresolved)


@pytest.mark.parametrize("variant_type", ("missense", "Missense", " MISSENSE "))
def test_bp1_normalizes_variant_type_at_rule_boundary(variant_type):
    result = evaluate_bp1(
        "BRCA1",
        variant_type,
        "p.(Arg500Gln)",
        spliceai_score=0.01,
    )

    assert result["applies"] is True
    assert result["strength"] == "Strong"
    assert result["points"] == -4
    assert result["decision_path"]["branch_id"] == "missense-inframe"


@pytest.mark.parametrize(
    ("variant_type", "c_notation", "in_domain", "branch_id"),
    (
        ("Intronic", "c.500+7A>G", False, "intronic"),
        (" INTRONIC ", "c.500+7A>G", False, "intronic"),
        ("Synonymous", "c.306A>G", True, "synonymous"),
        (" SYNONYMOUS ", "c.306A>G", True, "synonymous"),
    ),
)
def test_bp7_normalizes_variant_type_at_rule_boundary(
    variant_type,
    c_notation,
    in_domain,
    branch_id,
):
    result = evaluate_bp7(
        variant_type,
        spliceai_score=0.01,
        in_domain=in_domain,
        bp4_met=True,
        c_notation=c_notation,
        gene="BRCA1",
    )

    assert result["applies"] is True
    assert result["points"] == -1
    assert result["decision_path"]["branch_id"] == branch_id


@pytest.mark.parametrize("bioinformatic_code", ("PP3", "BP4", "BP7", "BP1"))
def test_pvs1_rna_replaces_each_figure1a_code(bioinformatic_code):
    combined = {
        "PVS1_RNA": {"points": 8},
        bioinformatic_code: {"points": 1 if bioinformatic_code == "PP3" else -1},
    }
    interactions = apply_rna_interactions(combined, {"PVS1_RNA"})

    assert bioinformatic_code not in combined
    assert interactions[0]["status"] == "deduplicated"
    assert interactions[0]["suppressed"] == [bioinformatic_code]


def test_automatic_rna_interaction_dispatch_discovers_all_applied_rna_codes():
    combined = {
        "PVS1_RNA": {"points": 4},
        "BP7_RNA": {"points": -4},
        "PP4": {"points": 4},
    }

    assert rna_interaction_codes(combined) == {"PVS1_RNA", "BP7_RNA"}


@pytest.mark.parametrize("bioinformatic_code", ("BP1", "BP4"))
def test_bp7_rna_retains_applicable_figure1a_benign_code_with_audit(bioinformatic_code):
    combined = {
        "BP7_RNA": {"points": -4},
        bioinformatic_code: {"points": -4 if bioinformatic_code == "BP1" else -1},
    }
    interactions = apply_rna_interactions(combined, {"BP7_RNA"})

    assert set(combined) == {"BP7_RNA", bioinformatic_code}
    assert len(interactions) == 1
    assert interactions[0]["status"] == "info"
    assert interactions[0]["retained"] == ["BP7_RNA", bioinformatic_code]


@pytest.mark.parametrize(
    ("functional_code", "bioinformatic_code", "expected_status"),
    (
        ("PS3", "PP3", "info"),
        ("PS3", "BP4", "conflict"),
        ("PS3", "BP7", "conflict"),
        ("PS3", "BP1", "conflict"),
        ("BS3", "PP3", "conflict"),
        ("BS3", "BP4", "info"),
        ("BS3", "BP7", "info"),
        ("BS3", "BP1", "info"),
    ),
)
def test_functional_and_figure1a_interaction_matrix(
    functional_code, bioinformatic_code, expected_status
):
    criteria = {
        functional_code: {"points": 4 if functional_code == "PS3" else -4},
        bioinformatic_code: {"points": 1 if bioinformatic_code == "PP3" else -1},
    }
    interactions = automatic_functional_interactions(criteria)

    assert len(interactions) == 1
    assert interactions[0]["status"] == expected_status
    assert interactions[0]["review_required"] is (expected_status == "conflict")
    assert set(interactions[0]["retained"]) == {functional_code, bioinformatic_code}


@pytest.mark.parametrize(
    ("gene", "c_notation"),
    (
        ("BRCA1", "c.5074G>A"),
        ("BRCA1", "c.5074G>C"),
        ("BRCA2", "c.7976G>A"),
        ("BRCA2", "c.7976G>C"),
    ),
)
def test_last_exon_base_examples_keep_table9_ps3_and_route_st2_rna_to_review(
    gene, c_notation
):
    functional = table9_lookup_ps3_bs3(gene, c_notation)
    rna = evaluate_pvs1_rna(gene, c_notation)

    assert functional["applies"]
    assert functional["code"] == "PS3"
    assert functional["strength"] == "Strong"
    assert functional["points"] == 4
    assert "mRNA and protein" in functional["source_record"]["standardised_text"]
    assert not rna["applies"]
    assert rna["points"] == 0
    assert rna["review_required"]
    assert rna["manual_review_prefill"]["assay_scope"] == "mrna_only"
    assert rna["source_record"]["supporting_st3_references"]
    assert rna["manual_review_prefill"]["source_references"]
    assert all(
        url.startswith("https://pubmed.ncbi.nlm.nih.gov/")
        for url in rna["manual_review_prefill"]["source_references"]
    )


def test_mrna_only_assay_cannot_create_manual_ps3():
    evidence = {
        "assay_scope": "mrna_only",
        "functional_conclusion": "abnormal",
        "calibration_status": "reviewed_under_enigma_vcep",
        "pathogenic_and_benign_controls_confirmed": True,
        "curated_strength": "Strong",
        "assay_name": "RNA assay",
        "source_citation": "PMID:test",
        "calibration_summary": "Reviewed controls and assay calibration.",
        "variant_result_summary": "Aberrant transcript detected.",
        "functional_reviewed_by": "Reviewer",
    }

    assert suggest_strength("PS3", evidence) is None
    evidence["assay_scope"] = "combined_mrna_protein"
    assert suggest_strength("PS3", evidence) == "Strong"


def test_pvs1_rna_and_table9_ps3_are_retained_with_pending_vcep_warning():
    combined = {
        "PVS1_RNA": {"points": 8},
        "PS3": {
            "points": 4,
            "table9_audit": {
                "assigned_code": "PS3",
                "standardised_text": "Combined mRNA and protein assay.",
            },
        },
    }

    interactions = apply_rna_interactions(combined, {"PVS1_RNA"})

    assert set(combined) == {"PVS1_RNA", "PS3"}
    assert len(interactions) == 1
    assert interactions[0]["status"] == "review_required"
    assert interactions[0]["retained"] == ["PVS1_RNA", "PS3"]
    assert "beyond merely detecting an aberrant transcript" in interactions[0]["reason"]
    assert "VCEP clarification" in interactions[0]["reason"]
