from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
from unittest.mock import patch

import pytest

from backend.classification_dag.runtime import ClassificationInputs, execute_classification
from backend.gene_policy import (
    GENE_POLICY_MANIFEST_PATH,
    GenePolicyError,
    active_genes,
    bayesdel_thresholds,
    load_gene_policy_manifest,
    reference_transcript,
    runtime_policy_id,
    validate_gene_policy_payload,
    validate_policy_source_bindings,
)
from backend.modules.pp3_bp4 import evaluate_pp3_bp4
from backend.modules.manual_evidence import (
    manual_criteria_for_gene,
    resource_links_for_gene,
)
from backend.modules.pvs1 import evaluate_pvs1


def _modified_manifest():
    return deepcopy(load_gene_policy_manifest())


def test_checked_manifest_and_source_bindings_are_valid():
    manifest = load_gene_policy_manifest()
    assert active_genes() == ("BRCA1", "BRCA2")
    assert manifest["manifest_version"] == "2026.08.25"
    assert reference_transcript("BRCA1") == "NM_007294.4"
    assert runtime_policy_id("BRCA2") == "ENIGMA_BRCA_VCEP_1.2"
    assert manifest["genes"]["BRCA1"]["decision_assets"]["PVS1"][
        "splice"
    ]["figure_number"] == "4"
    assert manifest["genes"]["BRCA2"]["functional_domains"]["DBD"][
        "description"
    ].startswith("DNA-binding domain")
    validate_policy_source_bindings()


def test_unknown_gene_has_no_threshold_fallback():
    with pytest.raises(GenePolicyError, match="No active gene policy"):
        bayesdel_thresholds("NOT_A_CONFIGURED_GENE")


def test_incomplete_gene_threshold_set_is_rejected():
    manifest = _modified_manifest()
    del manifest["genes"]["BRCA1"]["thresholds"]["bayesdel_noaf"]["pp3_min_inclusive"]
    with pytest.raises(GenePolicyError, match="incomplete"):
        validate_gene_policy_payload(manifest)


def test_checksum_drift_is_rejected():
    manifest = _modified_manifest()
    raw = json.dumps(manifest, sort_keys=True).encode("utf-8")
    metadata = {
        "schema_version": 1,
        "manifest_id": manifest["manifest_id"],
        "manifest_version": manifest["manifest_version"],
        "manifest_sha256": "0" * 64,
        "active_gene_count": 2,
        "active_policy_count": 1,
        "validation_status": "approved",
    }
    assert hashlib.sha256(raw).hexdigest() != metadata["manifest_sha256"]
    with pytest.raises(GenePolicyError, match="checksum mismatch"):
        validate_gene_policy_payload(manifest, metadata, manifest_bytes=raw)


def test_bayesdel_decision_changes_from_manifest_without_python_change():
    baseline = evaluate_pp3_bp4(
        "BRCA1", "missense", "p.(Ala1789Val)",
        bayesdel_score=0.438, spliceai_score=0.03,
    )
    assert baseline["PP3"]["applies"] is True

    manifest = _modified_manifest()
    manifest["genes"]["BRCA1"]["thresholds"]["bayesdel_noaf"]["pp3_min_inclusive"] = 0.50
    with patch("backend.gene_policy.load_gene_policy_manifest", return_value=manifest):
        changed = evaluate_pp3_bp4(
            "BRCA1", "missense", "p.(Ala1789Val)",
            bayesdel_score=0.438, spliceai_score=0.03,
        )
    assert "PP3" not in changed or changed["PP3"]["applies"] is False


def test_non_applicable_rule_is_excluded_before_classification():
    manifest = _modified_manifest()
    manifest["genes"]["BRCA1"]["applicable_rules"].remove("PP3")
    inputs = ClassificationInputs(
        gene="BRCA1",
        variant_type="missense",
        c_notation="c.5366C>T",
        p_notation="p.(Ala1789Val)",
        reference_transcript="NM_007294.4",
        spliceai_score=0.03,
        bayesdel_score=0.438,
    )
    with patch("backend.gene_policy.load_gene_policy_manifest", return_value=manifest):
        result = execute_classification(inputs).result
    assert "PP3" not in result["criteria"]
    assert result["excluded_criteria"]["PP3"]["points"] == 0
    assert "not applicable under the configured VCEP policy" in result["excluded_criteria"]["PP3"]["reason"]


def test_unimplemented_policy_profile_fails_closed_before_classification():
    manifest = _modified_manifest()
    manifest["policies"]["enigma_brca_vcep_1_2"][
        "implementation_profile"
    ] = "future_vcep_profile"
    inputs = ClassificationInputs(
        gene="BRCA1",
        variant_type="missense",
        c_notation="c.5366C>T",
        p_notation="p.(Ala1789Val)",
        reference_transcript="NM_007294.4",
    )
    with patch("backend.gene_policy.load_gene_policy_manifest", return_value=manifest):
        with pytest.raises(RuntimeError, match="No classification DAG is implemented"):
            execute_classification(inputs)


def test_gene_specific_vcep_sources_come_from_the_manifest():
    brca1_manual = manual_criteria_for_gene("BRCA1")
    brca2_manual = manual_criteria_for_gene("BRCA2")
    assert "GN092" in brca1_manual["PS4"]["source_url"]
    assert "GN097" in brca2_manual["PS4"]["source_url"]
    assert "GN097" in resource_links_for_gene("BRCA2")[0]["url"]
    assert "GN097" in evaluate_pvs1(
        "BRCA2", "missense", "p.(Ala1Val)", "c.2C>T"
    )["source"]


def test_generic_runtime_files_do_not_branch_on_brca_gene_symbols():
    project_root = Path(__file__).resolve().parents[1]
    files = (
        "backend/main.py",
        "backend/modules/variant_input.py",
        "backend/modules/hgvs_engine.py",
        "backend/modules/pvs1.py",
        "backend/modules/narrative.py",
        "backend/modules/table4.py",
        "backend/modules/enigma_rules.py",
        "backend/lookups/founder_variants.py",
        "backend/modules/exon_cnv_evidence.py",
    )
    forbidden = (
        re.compile(r'gene\s*==\s*["\']BRCA[12]["\']'),
        re.compile(r'gene\s+not\s+in\s+\{[^}]*BRCA1[^}]*BRCA2'),
        re.compile(r'BRCA\\\[12\\\]'),
    )
    for relative in files:
        source = (project_root / relative).read_text(encoding="utf-8")
        assert not any(pattern.search(source) for pattern in forbidden), relative
