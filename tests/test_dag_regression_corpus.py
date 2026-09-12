"""Approved clinical regression corpus for the production DAG.

Expected classes, points and criteria are explicit. No second classifier is
used as an oracle, so a change in clinical output requires a deliberate review
of this corpus.
"""

import json

import pytest

from backend.classification_dag import ClassificationInputs, execute_classification
from backend.config import ST2_SPLICE_EVIDENCE_PATH
from backend.lookups.founder_variants import (
    FOUNDER_VARIANT_SNAPSHOT,
    lookup_pathogenic_founder_variant,
)
from backend.modules.exon_cnv_evidence import lookup_exon_cnv_evidence
from backend.population_frequency.policy import classification_policy_for_gene
from backend.modules.pp4_bp5 import evaluate_pp4_bp5
from backend.modules.pvs1 import evaluate_pvs1
from backend.modules.pvs1_rna import (
    evaluate_pvs1_rna,
    lookup_st2_pvs1_rna_evidence,
)
from backend.modules.table9 import table9_lookup_ps3_bs3
from backend.modules.variant_type import infer_variant_type
from backend.population_frequency.indel_size import is_indel_allele


REGRESSION_VARIANTS = (
    ("BRCA1", "c.509G>A", "p.(Arg170Gln)", 0.05, None, "Unknown", 1, -12,
     {"BS3": ("Strong", -4), "BP5": ("Strong", -4),
      "BP1": ("Strong", -4)}, False, ()),
    ("BRCA1", "c.1534C>T", "p.(Leu512Phe)", 0.02, None, "Unknown", 1, -16,
     {"BS3": ("Strong", -4), "BP5": ("Very Strong", -8),
      "BP1": ("Strong", -4)}, False, ()),
    ("BRCA1", "c.3668_3671dup", "p.(Cys1225SerfsTer10)", None, None,
     "Unknown", 5, 12, {"PVS1": ("Very Strong", 8),
     "PM5_PTC": ("Strong", 4)}, False, ("PM2",)),
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
     1, -8, {"BS3": ("Strong", -4), "BP5": ("Strong", -4)}, False,
     ("PM2",)),
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


FOUNDER_REGRESSION_VARIANTS = tuple(
    (
        record["gene"],
        record["canonical_c_notation"],
        record["protein_notation"],
    )
    for record in json.loads(
        FOUNDER_VARIANT_SNAPSHOT.read_text(encoding="utf-8")
    )["variants"]
)


UNQUANTIFIED_PATIENT_RNA_VARIANTS = tuple(
    (record["gene"], record["c_notation"])
    for record in json.loads(
        ST2_SPLICE_EVIDENCE_PATH.read_text(encoding="utf-8")
    )["variants"]
    if " ".join(
        str(record.get("splicing_assay_result_category") or "").split()
    ).lower()
    == (
        "patient not allele-specific; aberrant transcripts consistent with "
        "loss of function"
    )
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


def _frequency_input(*, gene="BRCA1", max_af=None, found=False, absent=False):
    status = "found" if found else "absent"
    policy = classification_policy_for_gene(gene)
    return {
        "policy_id": policy["policy_id"],
        "classification_policy": policy,
        "frequency_policy": policy["frequency_criteria"],
        "status": "found" if found else "absent_with_coverage",
        "found": found,
        "max_af": max_af,
        "frequency_metric": "faf95",
        "pm2_absence_established": absent,
        "pm2_coverage_method": {
            "status": "approved_test_fixture",
            "automatic_assignment_allowed": True,
            "reason": "explicit approved method fixture",
        },
        "founder_exception": {
            "status": "reviewed_not_found",
            "is_pathogenic_founder": False,
            "reason": "regression fixture records an authoritative negative review",
            "snapshot_version": "test",
        },
        "datasets": {
            name: {
                "status": status,
                "max_af": max_af if found else None,
                "coverage": {
                    "mean_depth": 30.0,
                    "classification_compatible": True,
                },
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


def test_c4185_uses_curated_st2_rna_and_replaces_pp3():
    automatic = execute_classification(
        _inputs(
            "BRCA1",
            "c.4185G>A",
            "p.(Gln1395=)",
            0.95,
            None,
            "Unknown",
        )
    ).result

    assert automatic["predicted_class"] == 4
    assert automatic["total_points"] == 8
    assert _criterion_summary(automatic) == {
        "PVS1_RNA": ("Strong", 4),
        "PP4": ("Strong", 4),
    }
    assert automatic["rna_review"]["recommended"] is False
    assert automatic["evidence_interactions"] == []


@pytest.mark.parametrize("gene,c_notation", UNQUANTIFIED_PATIENT_RNA_VARIANTS)
def test_unquantified_patient_rna_is_scored_only_for_unambiguous_st2_table4_path(
    gene,
    c_notation,
):
    result = evaluate_pvs1_rna(gene, c_notation)

    if result["application_status"] == "applied_from_curated_st2":
        assert result["applies"] is True
        assert result["strength"] == "Strong"
        assert result["points"] == 4
        assert result["review_required"] is False
        assert result["table4_exon"]
        assert result["source_record"]["supporting_st3_references"]
    else:
        assert result["applies"] is False
        assert result["points"] == 0
        assert result["review_required"] is True
        assert result["application_status"] == "review_required"
        assert result["appendix_branch"] == (
            "unquantified_patient_mrna_requires_consensus_review"
        )
        assert "curated_strength" not in result["manual_review_prefill"]


def test_unquantified_patient_rna_scope_is_explicit_and_complete():
    outcomes = [
        evaluate_pvs1_rna(gene, c_notation)["application_status"]
        for gene, c_notation in UNQUANTIFIED_PATIENT_RNA_VARIANTS
    ]

    assert len(outcomes) == 26
    assert outcomes.count("applied_from_curated_st2") == 16
    assert outcomes.count("review_required") == 10


@pytest.mark.parametrize(
    "baseline_strength,baseline_code,expected_strength,expected_points",
    (
        ("Very Strong", "PVS1", "Strong", 4),
        ("Strong", "PVS1_Strong", "Moderate", 2),
        ("Moderate", "PVS1_Moderate", "Moderate", 2),
        ("Supporting", "PVS1_Supporting", "Supporting", 1),
    ),
)
def test_unquantified_patient_rna_uses_appendix_e_weight_matrix(
    monkeypatch,
    baseline_strength,
    baseline_code,
    expected_strength,
    expected_points,
):
    curated = lookup_st2_pvs1_rna_evidence("BRCA1", "c.4185G>A")
    baseline = {
        "found": True,
        "pvs1_strength": baseline_strength,
        "pvs1_code": baseline_code,
    }
    curated = {
        **curated,
        "table4_baseline": baseline,
    }
    monkeypatch.setattr(
        "backend.modules.pvs1_rna.table4_lookup_deletion",
        lambda gene, exon: baseline,
    )

    result = evaluate_pvs1_rna(
        "BRCA1",
        "c.4185G>A",
        st2_evidence_result=curated,
    )

    assert result["applies"] is True
    assert result["strength"] == expected_strength
    assert result["points"] == expected_points


@pytest.mark.parametrize(
    "gene,c_notation,p_notation",
    FOUNDER_REGRESSION_VARIANTS,
)
def test_all_registered_pathogenic_founders_suppress_ba1_and_bs1(
    gene,
    c_notation,
    p_notation,
):
    frequency = _frequency_input(gene=gene, max_af=0.002, found=True)
    frequency["founder_exception"] = lookup_pathogenic_founder_variant(
        gene, c_notation
    )

    result = execute_classification(
        ClassificationInputs(
            gene=gene,
            variant_type=infer_variant_type(c_notation, p_notation),
            c_notation=c_notation,
            p_notation=p_notation,
            gnomad_data=frequency,
            frequency_policy=classification_policy_for_gene(gene),
        )
    ).result

    assert frequency["founder_exception"]["status"] == "pathogenic_founder"
    assert "BA1" not in result["criteria"]
    assert not any(code.startswith("BS1") for code in result["criteria"])
    assert result["excluded_criteria"]["BA1"]["applies"] is False
    assert "pathogenic founder" in result["excluded_criteria"]["BA1"][
        "reason"
    ]


@pytest.mark.parametrize(
    (
        "gene,c_notation,p_notation,expected_applies,"
        "expected_pvs1_code,expected_pm5_code"
    ),
    (
        (
            "BRCA1", "c.5560C>T", "p.(Leu1854Ter)", True,
            "PVS1", "PM5_Strong (PTC)",
        ),
        (
            "BRCA1", "c.5563C>T", "p.(Ile1855Ter)", False,
            "PVS1_N/A", None,
        ),
        (
            "BRCA2", "c.9925G>T", "p.(Glu3309Ter)", True,
            "PVS1", "PM5_Strong (PTC)",
        ),
        (
            "BRCA2", "c.9928A>T", "p.(Lys3310Ter)", False,
            "PVS1_N/A", None,
        ),
    ),
)
def test_last_exon_ptc_boundaries_are_pinned_for_both_genes(
    gene,
    c_notation,
    p_notation,
    expected_applies,
    expected_pvs1_code,
    expected_pm5_code,
):
    result = evaluate_pvs1(
        gene,
        "nonsense",
        p_notation,
        c_notation,
    )

    assert result["applies"] is expected_applies
    assert result["pvs1_code"] == expected_pvs1_code
    assert result["pm5_code"] == expected_pm5_code
