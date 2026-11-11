import pytest

from backend.policy.criterion_points import points_for_strength


@pytest.mark.parametrize(
    ("strength", "points"),
    (("Supporting", 1), ("Moderate", 2), ("Strong", 4), ("Very Strong", 8)),
)
def test_points_for_strength(strength: str, points: int) -> None:
    assert points_for_strength(strength) == points
    assert points_for_strength(strength, benign=True) == -points


def test_missing_strength_has_no_points() -> None:
    assert points_for_strength(None) == 0


def test_points_for_strength_rejects_unknown_strength() -> None:
    with pytest.raises(KeyError):
        points_for_strength("Stand-alone")
