from __future__ import annotations

from pathlib import Path

from scripts.prepare_source_update import _increment_patch_version
from scripts.sync_source_watch_issue import MARKER_PREFIX, build_issue_body


def test_issue_body_links_report_and_manual_candidate_workflow() -> None:
    report = {"review_fingerprint": "abc123"}

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
