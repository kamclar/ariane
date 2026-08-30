"""Approved clinical regression corpus for the production DAG.

Expected classes, points and criteria are explicit. No second classifier is
used as an oracle, so a change in clinical output requires a deliberate review
of this corpus.
"""

import pytest

from backend.classification_dag import ClassificationInputs, execute_classification
from backend.modules.exon_cnv_evidence import lookup_exon_cnv_evidence
from backend.population_frequency.policy import classification_policy_for_gene
from backend.modules.pp4_bp5 import evaluate_pp4_bp5
from backend.modules.table9 import table9_lookup_ps3_bs3
from backend.modules.variant_type import infer_variant_type
from backend.population_frequency.indel_size import is_indel_allele


REGRESSION_VARIANTS = (
    ("BRCA1", "c.509G>A", "p.(Arg170Gln)", 0.05, None, "Unknown", 1, -8,
     {"BS3": ("Strong", -4), "BP1": ("Strong", -4)}, False, ()),
    ("BRCA1", "c.1534C>T", "p.(Leu512Phe)", 0.02, None, "Unknown", 1, -8,
     {"BS3": ("Strong", -4), "BP1": ("Strong", -4)}, False, ()),
    ("BRCA1", "c.3668_3671dup", "p.(Cys1225SerfsTer10)", None, None,
     "Unknown", 5, 16, {"PVS1": ("Very Strong", 8),
     "PM5_PTC": ("Strong", 4), "PP4": ("Strong", 4)}, False, ("PM2",)),
    ("BRCA2", "c.9097del", "p.(Thr3033fs)", None, None, "Unknown", 5, 12,
     {"PVS1": ("Very Strong", 8), "PM5_PTC": ("Strong", 4)}, False,
     ("PM2",)),
    ("BRCA1", "c.5551_5552insT", "p.(Asp1851ValfsTer29)", None, None,
     "Unknown", 5, 12, {"PVS1": ("Very Strong", 8),
     "PM5_PTC": ("Strong", 4)}, False, ("PM2",)),
    ("BRCA2", "c.(793+1_794-1)_(1909+1_1910-1)del", "p.(?)", None, None,
     "Unknown", 3, 1, {"PM2_Supporting": ("Supporting", 1)}, False,
     ("PVS1",)),
    ("BRCA2", "c.6147_6149del", "p.(Val2050del)", 0.01, None, "Unknown",
     2, -4, {"BP1": ("Strong", -4)}, False, ("PM2",)),
    ("BRCA1", "c.3891_3893del", "p.(Ser1298del)", 0.15, None, "Unknown",
     3, -4, {"BS3": ("Strong", -4)}, False, ("PM2",)),
    ("BRCA1", "c.4185G>A", "p.(Gln1395=)", 0.95, None, "Unknown", 4, 8,
     {"PVS1_RNA": ("Strong", 4), "PP4": ("Strong", 4)}, False, ()),
    ("BRCA1", "c.628C>T", "p.(Gln210Ter)", None, None, "Unknown", 3, 0,
     {}, False, ("PVS1",)),
    ("BRCA2", "c.8953+2T>C", "p.(?)", 0.90, None, "Unknown", 3, 0, {},
     False, ("PVS1",)),
    ("BRCA1", "c.5366C>T", "p.(Ala1789Val)", 0.05, 0.50, "Unknown", 3,
     5, {"PS3": ("Strong", 4), "PP3": ("Supporting", 1)}, False, ()),
    ("BRCA1", "c.3247A>G", "p.(Met1083Val)", 0.01, None, "Unknown", 2,
     -3, {"PP4": ("Supporting", 1), "BP1": ("Strong", -4)}, True, ()),
    ("BRCA1", "c.5556_5560del", "p.(Gln1853fs)", None, None, "Unknown", 5,
     12, {"PVS1": ("Very Strong", 8), "PM5_PTC": ("Strong", 4)}, False,
     ("PM2",)),
    ("BRCA1", "c.5533_5534insG", "p.(Tyr1845Ter)", None, None, "Unknown",
     5, 12, {"PVS1": ("Very Strong", 8), "PM5_PTC": ("Strong", 4)},
     False, ("PM2",)),
    ("BRCA2", "c.9891_9894dup", "p.(Gln3299IlefsTer29)", None, None,
     "Unknown", 5, 11, {"PVS1": ("Very Strong", 8),
     "PM5_PTC": ("Strong", 4), "BP5": ("Supporting", -1)}, True,
     ("PM2",)),
)


def _inputs(gene, c_notation, p_notation, spliceai, bayesdel, dup_type):
    variant_type = infer_variant_type(c_notation, p_notation)
    return ClassificationInputs(
        gene=gene,
        variant_type=variant_type,
        c_notation=c_notation,
        p_notation=p_notation,
        spliceai_score=spliceai,
        bayesdel_score=bayesdel,
        frequency_policy=classification_policy_for_gene(gene),
        table9_result=table9_lookup_ps3_bs3(gene, c_notation),
        pp4_bp5_result=evaluate_pp4_bp5(gene, c_notation),
        exon_cnv_result=(
            lookup_exon_cnv_evidence(gene, c_notation)
            if is_indel_allele(c_notation)
            else None
        ),
        dup_type=dup_type,
    )


def _criterion_summary(result):
    return {
        code: (criterion["strength"], criterion["points"])
        for code, criterion in result["criteria"].items()
    }


@pytest.mark.parametrize(
    (
        "gene,c_notation,p_notation,spliceai,bayesdel,dup_type,expected_class,"
        "expected_points,expected_criteria,expected_mixed,expected_na"
    ),
    REGRESSION_VARIANTS,
)
def test_approved_variant_corpus_has_stable_dag_results(
    gene,
    c_notation,
    p_notation,
    spliceai,
    bayesdel,
    dup_type,
    expected_class,
    expected_points,
    expected_criteria,
    expected_mixed,
    expected_na,
):
    result = execute_classification(
        _inputs(gene, c_notation, p_notation, spliceai, bayesdel, dup_type)
    ).result

    assert result["predicted_class"] == expected_class
    assert result["total_points"] == expected_points
    assert _criterion_summary(result) == expected_criteria
    assert result["mixed_evidence"] is expected_mixed
    assert tuple(sorted(result["not_applicable_criteria"])) == expected_na


@pytest.mark.parametrize(
    "ps1_result,expected_class,expected_points,expected_criteria,review_required",
    (
        (
            {
                "applies": True,
                "strength": "Strong",
                "points": 4,
                "reason": "Approved PS1 reference used by regression test.",
            },
            4,
            8,
            {"PS3": ("Strong", 4), "PS1": ("Strong", 4)},
            False,
        ),
        (
            {
                "applies": False,
                "review_required": True,
                "application_status": "review_required",
                "reason": "Candidate reference requires review.",
            },
            3,
            4,
            {"PS3": ("Strong", 4)},
            True,
        ),
    ),
)
def test_protein_ps1_applied_and_review_paths_are_explicit_regressions(
    ps1_result,
    expected_class,
    expected_points,
    expected_criteria,
    review_required,
):
    result = execute_classification(
        ClassificationInputs(
            gene="BRCA1",
            variant_type="missense",
            c_notation="c.131G>C",
            p_notation="p.(Cys44Ser)",
            spliceai_score=0.01,
            table9_result=table9_lookup_ps3_bs3("BRCA1", "c.131G>C"),
            ps1_result=ps1_result,
        )
    ).result

    assert result["predicted_class"] == expected_class
    assert result["total_points"] == expected_points
    assert _criterion_summary(result) == expected_criteria
    assert result["protein_ps1_review"]["recommended"] is review_required


def _frequency_input(*, max_af=None, found=False, absent=False):
    status = "found" if found else "absent"
    policy = classification_policy_for_gene("BRCA1")
    return {
        "policy_id": policy["policy_id"],
        "classification_policy": policy,
        "frequency_policy": policy["frequency_criteria"],
        "status": "found" if found else "absent_with_coverage",
        "found": found,
        "max_af": max_af,
        "frequency_metric": "faf95",
        "pm2_absence_established": absent,
        "founder_exception": {
            "status": "not_found",
            "is_pathogenic_founder": False,
            "reason": "regression fixture is not a pathogenic founder variant",
            "snapshot_version": "test",
        },
        "datasets": {
            name: {
                "status": status,
                "max_af": max_af if found else None,
                "coverage": {"mean_depth": 30.0},
                "quality_filter_passed": True if found else None,
                "non_founder_allele_count": 4 if found else 0,
            }
            for name in ("v2_1_non_cancer", "v3_1_non_cancer")
        },
    }


@pytest.mark.parametrize(
    "gnomad_data,expected_class,expected_points,expected_criteria,expected_mixed",
    (
        (_frequency_input(max_af=0.002, found=True), 1, -99,
         {"BA1": ("Stand-alone", -99)}, False),
        (_frequency_input(max_af=0.00005, found=True), 3, 4,
         {"BS1_Supporting": ("Supporting", -1), "PS3": ("Strong", 4),
          "PP3": ("Supporting", 1)}, True),
        (_frequency_input(absent=True), 4, 6,
         {"PM2_Supporting": ("Supporting", 1), "PS3": ("Strong", 4),
          "PP3": ("Supporting", 1)}, False),
    ),
)
def test_population_terminal_mixed_and_absent_paths_are_explicit_regressions(
    gnomad_data,
    expected_class,
    expected_points,
    expected_criteria,
    expected_mixed,
):
    result = execute_classification(
        ClassificationInputs(
            gene="BRCA1",
            variant_type="missense",
            c_notation="c.5366C>T",
            p_notation="p.(Ala1789Val)",
            spliceai_score=0.05,
            bayesdel_score=0.50,
            gnomad_data=gnomad_data,
            frequency_policy=classification_policy_for_gene("BRCA1"),
            table9_result=table9_lookup_ps3_bs3("BRCA1", "c.5366C>T"),
        )
    ).result

    assert result["predicted_class"] == expected_class
    assert result["total_points"] == expected_points
    assert _criterion_summary(result) == expected_criteria
    assert result["mixed_evidence"] is expected_mixed
