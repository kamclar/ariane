"""Runner contracts: aligned arrays, persisted results and honest denominators."""

from types import SimpleNamespace

import numpy as np
import pytest

from scripts import spliceai_distance_scan as scan
from scripts.compare_spliceai_distances import align_predictions
from scripts.spliceai_distance_store import ResultStore


def test_simple_indel_alignment_matches_reference_coordinates():
    ref = np.zeros((7, 3))
    deleted = np.ones((5, 3))
    _, aligned = align_predictions(ref, deleted, 3, 1, 3)
    np.testing.assert_array_equal(aligned[4:6], np.zeros((2, 3)))
    assert aligned.shape == ref.shape
    inserted = np.arange(27).reshape(9, 3)
    _, aligned = align_predictions(ref, inserted, 1, 3, 3)
    np.testing.assert_array_equal(aligned[3], inserted[3:6].max(axis=0))
    assert aligned.shape == ref.shape
    with pytest.raises(ValueError, match="outside"):
        align_predictions(ref, ref, 2, 2, 3)


def case(identifier):
    return {"id": identifier, "variant": "chr17-100-A-C", "gene": "BRCA1",
            "variant_kind": "snv", "strata": ["snv", "exon"], "locus_id": "tx:100"}


def test_completed_evaluation_can_be_checkpointed_and_resumed(tmp_path, monkeypatch):
    def fake_infer(_annotator, _case, distance):
        ref = np.zeros((2 * distance + 1, 3), dtype=np.float32)
        return (ref, ref.copy()), {}, 0.01

    monkeypatch.setattr(scan.core, "infer", fake_infer)
    monkeypatch.setattr(scan.core, "validate_capture", lambda *_: None)
    args = SimpleNamespace(distances=[50, 4999, 5000, 10000], repeats=2,
                           timing_repeat_cases=1, tolerance=1e-6,
                           save_arrays="none", output=tmp_path)
    result = scan.evaluate_case(None, case("a"), 0, args)
    assert result["id"] == "a" and result["status"] == "ok"
    assert len(result["runs"]) == 6  # D50 is fixture-only; repeats are not new cases.
    path = tmp_path / "results.sqlite3"
    with ResultStore(path, {"test": True}) as store:
        store.put("a", result)
    with ResultStore(path, {"test": True}, resume=True) as store:
        assert store.get("a") == result
        assert store.counts() == {"total": 1, "ok": 1, "error": 0}


def test_summary_never_counts_error_or_unfinished_as_agreement():
    plan = {"cohort": {"cases": [case("a"), case("b"), case("c")], "coverage": {}}}
    rows = [("a", {"status": "ok", "flags": {"delta_changed": False}}),
            ("b", {"status": "error", "error": "model failed"})]
    result = scan.summarize(plan, rows, "paused")
    assert result["counts"] == {"planned": 3, "ok": 1, "error": 1, "remaining": 1}
    assert result["groups"]["snv"] == {
        "planned": 3, "ok": 1, "error": 1, "delta_changed": 0, "band_changed": 0}
    assert result["distinct_loci"] == {"planned": 1, "with_successful_results": 1}
    assert result["first_anomalies"][0]["status"] == "error"


def test_retry_requires_explicit_resume():
    args = SimpleNamespace(positions_per_region=3, repeats=2, threads=4,
                           progress_every=100, max_errors=5, timing_repeat_cases=12,
                           tolerance=1e-6, max_cases=None, max_seconds=None,
                           distances=[4999, 5000, 10000], retry_errors=True, resume=False)
    with pytest.raises(ValueError, match="--resume"):
        scan.validate_options(args)
