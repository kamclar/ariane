"""Contracts for the D=50/500 extension; no model or network access."""

import json
import sqlite3

from scripts.spliceai_default_distance_scan import (
    baseline_projection,
    baseline_summary,
    differing_paths,
    flags_for,
)


def test_tolerant_nested_baseline_comparison_reports_exact_paths():
    expected = {"maximum_delta": 0.2, "components": {"AG": {"position": 4}}}
    observed = {"maximum_delta": 0.2000004, "components": {"AG": {"position": 5}}}
    assert differing_paths(expected, observed, 1e-6) == ["components.AG.position"]
    assert differing_paths(expected, {"maximum_delta": 0.3}, 1e-6) == [
        "components", "maximum_delta"
    ]


def test_base_checkpoint_requires_one_successful_first_repeat(tmp_path):
    database = sqlite3.connect(tmp_path / "base.sqlite3")
    database.execute("CREATE TABLE cases(case_id TEXT, status TEXT, result_json TEXT)")
    result = {
        "runs": [
            {"distance": 4999, "repeat": 1, "seconds": 2.0, "maximum_delta": 0.25},
            {"distance": 4999, "repeat": 2, "seconds": 1.9, "maximum_delta": 0.25},
        ]
    }
    database.execute("INSERT INTO cases VALUES (?, ?, ?)", ("case-a", "ok", json.dumps(result)))
    database.commit()
    assert baseline_summary(database, "case-a")["maximum_delta"] == 0.25


def test_baseline_projection_ignores_argmax_tie_position_and_probabilities():
    first = {
        "distance": 4999, "maximum_delta": 0.03, "band_raw": "le_0.1",
        "band_rounded_3dp": "le_0.1", "components": {
            "AL": {"delta": 0.0002, "position": 2696, "ref": 0.03, "alt": 0.0298}
        }, "edge_deltas": {}, "outside_model_radius": {"positions_examined": 0},
    }
    second = json.loads(json.dumps(first))
    second["components"]["AL"].update(position=2796, ref=0.22, alt=0.2198)
    assert differing_paths(baseline_projection(first), baseline_projection(second), 1e-6) == []
    second["components"]["AL"]["delta"] = 0.002
    assert differing_paths(baseline_projection(first), baseline_projection(second), 1e-6) == [
        "component_deltas.AL"
    ]


def comparison(*, delta_equal=True, rounded_equal=True, position=False,
               negligible=False, band_equal=True, overlap=0.0):
    channel = {"position_changed": position, "both_deltas_negligible": negligible}
    return {
        "all_delta_scores_within_tolerance": delta_equal,
        "all_rounded_delta_scores_equal": rounded_equal,
        "raw_band_equal": band_equal,
        "rounded_band_equal": band_equal,
        "overlapping_ref_max_abs_difference": overlap,
        "overlapping_alt_max_abs_difference": 0.0,
        "channels": {name: channel for name in ("AG", "AL", "DG", "DL")},
    }


def test_flags_ignore_arbitrary_argmax_for_negligible_zero_scores():
    result = flags_for(
        [comparison(position=True, negligible=True)],
        [{"outside_model_radius": {"maximum_absolute_delta": None}}],
        "snv",
        1e-6,
    )
    assert result == {
        "delta_changed": False,
        "rounded_delta_changed": False,
        "position_changed": False,
        "band_changed": False,
        "overlap_changed": False,
        "snv_outside_radius": False,
    }


def test_flags_separate_numeric_rounding_band_and_overlap_changes():
    result = flags_for(
        [comparison(delta_equal=False, rounded_equal=False, position=True,
                    band_equal=False, overlap=2e-6)],
        [{"outside_model_radius": {"maximum_absolute_delta": 3e-6}}],
        "snv",
        1e-6,
    )
    assert all(result.values())
