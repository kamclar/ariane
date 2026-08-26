"""Clinical parity corpus for incremental migration from legacy to DAG nodes.

The legacy side of these comparisons must remain independent. When a criterion
family moves into native DAG nodes, update only the DAG construction, not the
legacy call or the corpus inputs.
"""

import pytest

from backend.classification_dag import ClassificationInputs, execute_classification
from backend.classification_dag.runtime import compare_classification_results
from backend.modules.classifier import evaluate_variant
from backend.modules.exon_cnv_evidence import lookup_exon_cnv_evidence
from backend.modules.pp4_bp5 import evaluate_pp4_bp5
from backend.modules.table9 import table9_lookup_ps3_bs3
from backend.modules.variant_type import infer_variant_type
from backend.modules.frequency import _classification_policy_for_gene


PARITY_VARIANTS = (
    ("BRCA1", "c.509G>A", "p.(Arg170Gln)", 0.05, None, "Unknown"),
    ("BRCA1", "c.1534C>T", "p.(Leu512Phe)", 0.02, None, "Unknown"),
    ("BRCA1", "c.3668_3671dup", "p.(Cys1225SerfsTer10)", None, None, "Unknown"),
    ("BRCA2", "c.9097del", "p.(Thr3033fs)", None, None, "Unknown"),
    ("BRCA1", "c.5551_5552insT", "p.(Asp1851ValfsTer29)", None, None, "Unknown"),
    (
        "BRCA2",
        "c.(793+1_794-1)_(1909+1_1910-1)del",
        "p.(?)",
        None,
        None,
        "Unknown",
    ),
    ("BRCA2", "c.6147_6149del", "p.(Val2050del)", 0.01, None, "Unknown"),
    ("BRCA1", "c.3891_3893del", "p.(Ser1298del)", 0.15, None, "Unknown"),
    ("BRCA1", "c.4185G>A", "p.(Gln1395=)", 0.95, None, "Unknown"),
    ("BRCA1", "c.628C>T", "p.(Gln210Ter)", None, None, "Unknown"),
    ("BRCA2", "c.8953+2T>C", "p.(?)", 0.90, None, "Unknown"),
    ("BRCA1", "c.5366C>T", "p.(Ala1789Val)", 0.05, 0.50, "Unknown"),
    ("BRCA1", "c.3247A>G", "p.(Met1083Val)", 0.01, None, "Unknown"),
    ("BRCA1", "c.5556_5560del", "p.(Gln1853fs)", None, None, "Unknown"),
    ("BRCA1", "c.5533_5534insG", "p.(Tyr1845Ter)", None, None, "Unknown"),
    ("BRCA2", "c.9891_9894dup", "p.(Gln3299IlefsTer29)", None, None, "Unknown"),
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
        table9_result=table9_lookup_ps3_bs3(gene, c_notation),
        pp4_bp5_result=evaluate_pp4_bp5(gene, c_notation),
        exon_cnv_result=(
            lookup_exon_cnv_evidence(gene, c_notation)
            if variant_type in {"exon_deletion", "exon_duplication"}
            else None
        ),
        dup_type=dup_type,
    )


def _oracle_kwargs(inputs):
    return {
        "gene": inputs.gene,
        "variant_type": inputs.variant_type,
        "p_notation": inputs.p_notation,
        "c_notation": inputs.c_notation,
        "spliceai_score": inputs.spliceai_score,
        "bayesdel_score": inputs.bayesdel_score,
        "gnomad_data": inputs.gnomad_data,
        "table9_result": inputs.table9_result,
        "pp4_bp5_result": inputs.pp4_bp5_result,
        "ps1_result": inputs.ps1_result,
        "exon_cnv_result": inputs.exon_cnv_result,
        "residue_info": inputs.residue_info,
        "dup_type": inputs.dup_type,
    }


def _frequency_input(*, gene="BRCA1", max_af=None, found=False, absent=False):
    status = "found" if found else "absent"
    policy = _classification_policy_for_gene(gene)
    return {
        "policy_id": policy["policy_id"],
        "frequency_policy": policy["frequency_criteria"],
        "status": "found" if found else "absent_with_coverage",
        "found": found,
        "max_af": max_af,
        "frequency_metric": "faf95",
        "pm2_absence_established": absent,
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
    "gene,c_notation,p_notation,spliceai,bayesdel,dup_type",
    PARITY_VARIANTS,
)
def test_tutorial_and_regression_variants_have_exact_legacy_dag_parity(
    gene, c_notation, p_notation, spliceai, bayesdel, dup_type
):
    inputs = _inputs(
        gene, c_notation, p_notation, spliceai, bayesdel, dup_type
    )
    legacy_result = evaluate_variant(**_oracle_kwargs(inputs))
    dag_result = execute_classification(
        inputs,
        mode="dag",
    ).result

    assert compare_classification_results(legacy_result, dag_result) == ()


@pytest.mark.parametrize(
    "criterion",
    (
        {
            "code": "PS1",
            "result": {
                "applies": True,
                "strength": "Strong",
                "points": 4,
                "reason": "Approved PS1 reference used by parity test.",
            },
        },
        {
            "code": "PS1",
            "result": {
                "applies": False,
                "review_required": True,
                "application_status": "review_required",
                "reason": "Candidate reference requires review.",
            },
        },
    ),
)
def test_ps1_applied_and_review_paths_have_exact_legacy_dag_parity(criterion):
    inputs = ClassificationInputs(
        gene="BRCA1",
        variant_type="missense",
        c_notation="c.131G>C",
        p_notation="p.(Cys44Ser)",
        spliceai_score=0.01,
        table9_result=table9_lookup_ps3_bs3("BRCA1", "c.131G>C"),
        ps1_result=criterion["result"],
    )
    legacy_result = evaluate_variant(**_oracle_kwargs(inputs))
    dag_result = execute_classification(
        inputs,
        mode="dag",
    ).result
    assert compare_classification_results(legacy_result, dag_result) == ()


@pytest.mark.parametrize(
    "gnomad_data",
    (
        _frequency_input(max_af=0.002, found=True),
        _frequency_input(max_af=0.00005, found=True),
        _frequency_input(absent=True),
    ),
)
def test_population_terminal_mixed_and_absent_paths_have_exact_parity(gnomad_data):
    inputs = ClassificationInputs(
        gene="BRCA1",
        variant_type="missense",
        c_notation="c.5366C>T",
        p_notation="p.(Ala1789Val)",
        spliceai_score=0.05,
        bayesdel_score=0.50,
        gnomad_data=gnomad_data,
        table9_result=table9_lookup_ps3_bs3("BRCA1", "c.5366C>T"),
    )
    legacy_result = evaluate_variant(**_oracle_kwargs(inputs))
    dag_result = execute_classification(inputs, mode="dag").result
    assert compare_classification_results(legacy_result, dag_result) == ()
