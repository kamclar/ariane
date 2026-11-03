from __future__ import annotations

import subprocess
from pathlib import Path

from scripts.prepare_source_update import _clinical_lr_diff, _increment_patch_version
from scripts.sync_source_watch_issue import MARKER_PREFIX, build_issue_body


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_issue_body_links_report_and_manual_candidate_workflow() -> None:
    report = {
        "review_fingerprint": "abc123",
        "sources": [
            {
                "review_required": True,
                "candidate_workflow": "draft_pull_request",
            }
        ],
    }

    body = build_issue_body(
        report,
        "# Report\n\nA source changed.\n",
        repository="owner/repository",
        run_url="https://github.test/run/1",
    )

    assert "# Report" in body
    assert "https://github.test/run/1" in body
    assert "actions/workflows/prepare-source-update.yml" in body
    assert f"{MARKER_PREFIX}abc123 -->" in body
    assert "does not deploy" in body


def test_candidate_version_increment_changes_only_patch(tmp_path: Path) -> None:
    version_path = tmp_path / "version.py"
    version_path.write_text('ARIANE_VERSION = "1.9.36"\n', encoding="utf-8")

    previous, candidate = _increment_patch_version(version_path)

    assert previous == "1.9.36"
    assert candidate == "1.9.37"
    assert version_path.read_text(encoding="utf-8") == 'ARIANE_VERSION = "1.9.37"\n'


def test_generated_source_update_review_is_not_gitignored() -> None:
    result = subprocess.run(
        [
            "git",
            "check-ignore",
            "--quiet",
            "docs/source_update_reviews/example.md",
        ],
        cwd=REPOSITORY_ROOT,
        check=False,
    )

    assert result.returncode == 1


def test_clinical_lr_diff_highlights_threshold_and_direction_changes() -> None:
    pinned = {
        "BRCA1:c.1A>G": {
            "combined_lr": 1.0,
            "criterion": None,
            "strength": None,
            "points": 0,
            "automatic_application_status": "eligible",
            "overlap_status": "source_curated_combination",
        },
        "BRCA2:c.2T>C": {"combined_lr": 0.1, "criterion": "BP5", "strength": "Moderate", "points": -2},
        "BRCA2:c.3G>A": {"combined_lr": 10.0, "criterion": "PP4", "strength": "Moderate", "points": 2},
    }
    candidate = {
        "BRCA1:c.1A>G": {
            "combined_lr": 20.0,
            "criterion": "PP4",
            "strength": "Strong",
            "points": 4,
            "automatic_application_status": "eligible",
            "overlap_status": "source_curated_combination",
        },
        "BRCA2:c.2T>C": {"combined_lr": 3.0, "criterion": "PP4", "strength": "Supporting", "points": 1},
        "BRCA1:c.4C>T": {"combined_lr": 1.0, "criterion": None, "strength": None, "points": 0},
    }

    diff = _clinical_lr_diff(pinned, candidate)

    assert diff["added"] == ["BRCA1:c.4C>T"]
    assert diff["removed"] == ["BRCA2:c.3G>A"]
    assert diff["threshold_crossings"] == ["BRCA1:c.1A>G", "BRCA2:c.2T>C"]
    assert diff["direction_changes"] == ["BRCA1:c.1A>G", "BRCA2:c.2T>C"]
