"""Check the experiment's analysis using artificial arrays, without model I/O."""

import numpy as np
import pytest

from scripts.compare_spliceai_distances import (
    compare_arrays,
    score_band,
    summarize_arrays,
    validate_capture,
)


def arrays(distance):
    ref = np.zeros((2 * distance + 1, 3), dtype=np.float32)
    return ref, ref.copy()


def test_detects_new_signal_exactly_at_5000_not_in_4999_window():
    small = arrays(4999)
    large = arrays(5000)
    large[1][10000, 1] = 0.3
    result = compare_arrays(small, large, 4999, 5000, 1e-6)
    assert not result["all_delta_scores_within_tolerance"]
    assert result["channels"]["AG"]["large_position"] == 5000
    assert result["overlapping_alt_max_abs_difference"] == 0
    assert not result["raw_band_equal"]


def test_zero_delta_argmax_moves_without_changing_score():
    result = compare_arrays(arrays(4999), arrays(10000), 4999, 10000, 1e-6)
    assert result["all_delta_scores_within_tolerance"]
    assert result["channels"]["AG"]["position_changed"]
    assert result["channels"]["AG"]["both_deltas_negligible"]


def test_unrounded_changes_are_not_hidden_by_equal_rounded_scores():
    small = arrays(4999)
    large = arrays(10000)
    small[1][4999, 2] = 0.2061
    large[1][10000, 2] = 0.2062
    result = compare_arrays(small, large, 4999, 10000, 1e-6)
    assert result["all_rounded_delta_scores_equal"]
    assert not result["all_delta_scores_within_tolerance"]
    assert result["overlapping_alt_max_abs_difference"] > 1e-6


def test_loss_channel_and_model_boundary_are_checked_separately():
    ref, alt = arrays(10000)
    ref[4999, 2] = 0.5  # donor loss at -5001, outside receptive field
    alt[15000, 1] = 0.4  # acceptor gain at +5000, on the boundary
    result = summarize_arrays(ref, alt, 10000, 1e-6)
    assert result["components"]["DL"]["position"] == -5001
    assert result["outside_model_radius"]["positions_above_tolerance"] == 1
    assert result["outside_model_radius"]["maximum_absolute_delta"] == 0.5


def test_missing_outside_window_is_not_reported_as_observed_zero():
    result = summarize_arrays(*arrays(4999), 4999, 1e-6)
    assert result["outside_model_radius"]["positions_examined"] == 0
    assert result["outside_model_radius"]["maximum_absolute_delta"] is None


def test_nonfinite_predictions_are_rejected():
    ref, alt = arrays(50)
    alt[0, 0] = np.nan
    with pytest.raises(ValueError, match="Non-finite"):
        summarize_arrays(ref, alt, 50, 1e-6)


def test_threshold_boundaries_and_capture_orientation_guard():
    assert score_band(0.1) == "le_0.1"
    assert score_band(0.100001) == "between_0.1_and_0.2"
    assert score_band(0.2) == "ge_0.2"
    summary = summarize_arrays(*arrays(50), 50, 1e-6)
    row = {}
    for channel, item in summary["components"].items():
        row.update({f"DS_{channel}": 0, f"DS_{channel}_REF": 0,
                    f"DS_{channel}_ALT": 0, f"DP_{channel}": item["position"]})
    validate_capture(summary, row)
    row["DP_AG"] = 50
    with pytest.raises(ValueError, match="Captured position"):
        validate_capture(summary, row)
