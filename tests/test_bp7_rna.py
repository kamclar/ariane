"""Regression tests for ENIGMA BP7 Strong (RNA) variant stipulations."""

from backend.modules.bp7_rna import evaluate_bp7_rna_variant_context


def _context(c_notation: str, p_notation: str, gene: str = "BRCA1"):
    return {
        "gene": gene,
        "c_notation": c_notation,
        "p_notation": p_notation,
    }


def _table9_bs3():
    return [{
        "name": "BS3",
        "applies": True,
        "strength": "Strong",
        "points": -4,
        "decision_path": {
            "sources": [{"source_id": "enigma-v1.2-table9"}],
        },
    }]


def test_intronic_variant_is_context_eligible_without_bs3():
    result = evaluate_bp7_rna_variant_context(
        _context("c.4987-6T>G", "p.?"), []
    )
    assert result["eligible"] is True
    assert result["requires_bs3"] is False


def test_missense_outside_domain_is_context_eligible_without_bs3():
    result = evaluate_bp7_rna_variant_context(
        _context("c.509G>A", "p.(Arg170Gln)"), []
    )
    assert result["eligible"] is True
    assert result["functional_domains"] == []


def test_missense_inside_domain_requires_table9_bs3():
    context = _context("c.5123C>T", "p.(Ala1708Val)")
    missing = evaluate_bp7_rna_variant_context(context, [])
    present = evaluate_bp7_rna_variant_context(context, _table9_bs3())
    assert missing["eligible"] is False
    assert missing["requires_bs3"] is True
    assert missing["bs3_met"] is False
    assert present["eligible"] is True
    assert present["bs3_met"] is True


def test_bs3_without_table9_provenance_does_not_satisfy_domain_requirement():
    result = evaluate_bp7_rna_variant_context(
        _context("c.5123C>T", "p.(Ala1708Val)"),
        [{
            "name": "BS3",
            "applies": True,
            "strength": "Strong",
            "points": -4,
        }],
    )
    assert result["eligible"] is False
    assert result["bs3_met"] is False


def test_inframe_interval_overlapping_domain_is_not_treated_as_outside():
    result = evaluate_bp7_rna_variant_context(
        _context("c.4945_4953del", "p.(Arg1649_Ser1651del)"), []
    )
    assert result["eligible"] is False
    assert result["protein_interval"] == [1649, 1651]
    assert result["functional_domains"] == ["BRCT"]


def test_missing_variant_context_fails_closed():
    result = evaluate_bp7_rna_variant_context(None, [])
    assert result["eligible"] is False
    assert result["status"] == "unavailable"
