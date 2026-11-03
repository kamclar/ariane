from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from scripts import source_watch


class FakeResponse(io.BytesIO):
    def __init__(self, payload: bytes):
        super().__init__(payload)
        self.headers = {"Content-Length": str(len(payload))}

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()


def _json_bytes(value: object) -> bytes:
    return json.dumps(value).encode("utf-8")


def _cspec_payload(*, version: str = "1.2", revision: str = "rev-1") -> dict:
    return {
        "data": {
            "SequenceVariantInterpretation": [
                {
                    "entId": "GN092",
                    "rev": revision,
                    "modified": "2026-07-20T17:44:33Z",
                    "entContent": {
                        "version": version,
                        "approvedOn": "2025-01-09T01:11:34Z",
                        "title": "ENIGMA BRCA1 specification",
                        "states": [{"name": "Released", "current": True}],
                    },
                }
            ]
        }
    }


def _cspec_source() -> dict:
    return {
        "records": [
            {
                "id": "GN092",
                "gene": "BRCA1",
                "api_url": "https://example.test/GN092",
                "pinned_version": "1.2.0",
                "observed_revision": "rev-1",
                "observed_modified": "2026-07-20T17:44:33Z",
            }
        ]
    }


def test_repository_source_watch_registry_is_valid_and_fail_closed() -> None:
    registry = source_watch.load_registry()

    assert registry["policy"]["automatic_activation"] is False
    assert registry["policy"]["runtime_network_dependency"] is False
    assert {source["checker"] for source in registry["sources"]} >= {
        "cspec_records",
        "erepo_export",
        "clinical_lr",
        "gnomad_releases",
        "manual",
    }


def test_cspec_checker_accepts_equivalent_version_format() -> None:
    payload = _json_bytes(_cspec_payload())

    result = source_watch.check_cspec_records(
        _cspec_source(),
        1,
        opener=lambda *_args, **_kwargs: FakeResponse(payload),
    )

    assert result["status"] == "current"
    assert result["review_required"] is False
    assert result["records"][0]["remote_version"] == "1.2"


@pytest.mark.parametrize(
    ("payload", "expected_status"),
    [
        (_cspec_payload(version="1.3"), "update_available"),
        (_cspec_payload(revision="rev-2"), "upstream_changed"),
    ],
)
def test_cspec_checker_requires_review_for_version_or_identity_change(
    payload: dict, expected_status: str
) -> None:
    result = source_watch.check_cspec_records(
        _cspec_source(),
        1,
        opener=lambda *_args, **_kwargs: FakeResponse(_json_bytes(payload)),
    )

    assert result["status"] == expected_status
    assert result["review_required"] is True


def test_erepo_checker_compares_export_and_service_version(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    export = (
        b"Uuid\tExpert Panel\tHGNC Gene Symbol\tVariation\n"
        b"1\tPanel\tBRCA1\ttest\n"
    )
    metadata_path = tmp_path / "erepo.metadata.json"
    metadata_path.write_text(
        json.dumps(
            {
                "source_export_url": "https://example.test/export.csv",
                "source_export_sha256": hashlib.sha256(export).hexdigest(),
                "source_service_url": "https://example.test/service",
                "source_api_version": "1.6.1",
                "retrieved_at": "2026-09-12T00:00:00Z",
                "record_count": 1,
                "expert_panel": "Panel",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(source_watch, "ROOT", tmp_path)

    def opener(request: Any, **_kwargs: object) -> FakeResponse:
        if request.full_url.endswith("export.csv"):
            return FakeResponse(export)
        return FakeResponse(_json_bytes({"data": {"version": "1.6.1"}}))

    selection = source_watch.erepo_selection_identity(
        export, expert_panel="Panel", genes={"BRCA1"}
    )
    result = source_watch.check_erepo_export(
        {
            "metadata_path": "erepo.metadata.json",
            "genes": ["BRCA1"],
            "selected_rows_sha256": selection["sha256"],
            "selected_row_count": selection["row_count"],
        },
        1,
        opener=opener,
    )

    assert result["status"] == "current"
    assert result["remote_export_sha256"] == result["pinned_export_sha256"]


def test_erepo_checker_detects_changed_export(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    metadata_path = tmp_path / "erepo.metadata.json"
    metadata_path.write_text(
        json.dumps(
            {
                "source_export_url": "https://example.test/export.csv",
                "source_export_sha256": hashlib.sha256(b"old").hexdigest(),
                "source_service_url": "https://example.test/service",
                "source_api_version": "1.6.1",
                "expert_panel": "Panel",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(source_watch, "ROOT", tmp_path)

    def opener(request: Any, **_kwargs: object) -> FakeResponse:
        payload = (
            b"Uuid\tExpert Panel\tHGNC Gene Symbol\tVariation\n"
            b"1\tPanel\tBRCA1\tnew\n"
            if request.full_url.endswith("export.csv")
            else _json_bytes({"data": {"version": "1.6.1"}})
        )
        return FakeResponse(payload)

    old_selection = b"Uuid\tExpert Panel\tHGNC Gene Symbol\tVariation\n1\tPanel\tBRCA1\told\n"
    old_identity = source_watch.erepo_selection_identity(
        old_selection, expert_panel="Panel", genes={"BRCA1"}
    )
    result = source_watch.check_erepo_export(
        {
            "metadata_path": "erepo.metadata.json",
            "genes": ["BRCA1"],
            "selected_rows_sha256": old_identity["sha256"],
            "selected_row_count": old_identity["row_count"],
        },
        1,
        opener=opener,
    )

    assert result["status"] == "update_available"
    assert result["review_required"] is True


def test_erepo_checker_ignores_changes_outside_enigma_scope(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    selected = (
        b"Uuid\tExpert Panel\tHGNC Gene Symbol\tVariation\n"
        b"1\tPanel\tBRCA1\tselected\n"
    )
    remote = selected + b"2\tOther panel\tOTHER\tunrelated change\n"
    identity = source_watch.erepo_selection_identity(
        selected, expert_panel="Panel", genes={"BRCA1"}
    )
    metadata_path = tmp_path / "erepo.metadata.json"
    metadata_path.write_text(
        json.dumps(
            {
                "source_export_url": "https://example.test/export.csv",
                "source_export_sha256": hashlib.sha256(b"old full export").hexdigest(),
                "source_service_url": "https://example.test/service",
                "source_api_version": "1.6.1",
                "expert_panel": "Panel",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(source_watch, "ROOT", tmp_path)

    def opener(request: Any, **_kwargs: object) -> FakeResponse:
        payload = (
            remote
            if request.full_url.endswith("export.csv")
            else _json_bytes({"data": {"version": "1.6.1"}})
        )
        return FakeResponse(payload)

    result = source_watch.check_erepo_export(
        {
            "metadata_path": "erepo.metadata.json",
            "genes": ["BRCA1"],
            "selected_rows_sha256": identity["sha256"],
            "selected_row_count": identity["row_count"],
        },
        1,
        opener=opener,
    )

    assert result["status"] == "current"
    assert result["full_export_changed"] is True
    assert result["review_required"] is False


def test_erepo_record_diff_reports_added_modified_and_removed_records() -> None:
    export = (
        "Uuid\tExpert Panel\tHGNC Gene Symbol\tVariation\tAssertion\t"
        "Applied Evidence Codes (Met)\tApplied Evidence Codes (Not Met)\t"
        "Summary of interpretation\tApproval Date\tPublished Date\tRetracted\t"
        "Evidence Repo Link\n"
        "kept\tPanel\tBRCA1\tNM_007294.4(BRCA1):c.1A>G\tLikely Pathogenic\t"
        "PS3\t\tnew summary\t2026-10-30\t2026-10-29\tfalse\t"
        "https://example.test/kept\n"
        "added\tPanel\tBRCA2\tNM_000059.4(BRCA2):c.2T>C\tUncertain Significance\t"
        "\tPP3\tadded summary\t2026-10-30\t2026-10-29\tfalse\t"
        "https://example.test/added\n"
    ).encode()
    pinned = {
        "records": [
            {
                "uuid": "kept",
                "gene": "BRCA1",
                "c_notation": "c.1A>G",
                "p_notation": "",
                "classification": "Uncertain Significance",
                "criteria_met": [],
                "criteria_not_met": ["PS3"],
                "interpretation_summary": "old summary",
                "approval_date": "2026-01-01",
                "published_date": "2026-01-01",
                "retracted": False,
            },
            {
                "uuid": "removed",
                "gene": "BRCA2",
                "c_notation": "c.3G>A",
                "p_notation": "",
                "classification": "Benign",
                "erepo_url": "https://example.test/removed",
            },
        ]
    }

    diff = source_watch.erepo_record_diff(
        export,
        pinned_registry=pinned,
        expert_panel="Panel",
        genes={"BRCA1", "BRCA2"},
    )

    assert diff["added_count"] == 1
    assert diff["modified_count"] == 1
    assert diff["removed_count"] == 1
    assert diff["added"][0]["uuid"] == "added"
    assert diff["removed"][0]["uuid"] == "removed"
    assert diff["modified"][0]["uuid"] == "kept"
    assert set(diff["modified"][0]["changes"]) >= {
        "classification",
        "criteria_met",
        "criteria_not_met",
    }


def test_run_checks_distinguishes_review_from_operational_failure() -> None:
    registry = {
        "policy": {
            "default_timeout_seconds": 1,
            "required_before_activation": ["review"],
        },
        "sources": [
            {"id": "current", "name": "Current", "checker": "first"},
            {"id": "changed", "name": "Changed", "checker": "second"},
        ],
    }
    checkers = {
        "first": lambda _source, _timeout: {
            "status": "current",
            "review_required": False,
        },
        "second": lambda _source, _timeout: {
            "status": "update_available",
            "review_required": True,
        },
    }

    review_report = source_watch.run_checks(registry, checkers=checkers)

    assert review_report["overall_status"] == "review_required"
    assert review_report["exit_code"] == 2
    assert len(review_report["review_fingerprint"]) == 64

    def failing(_source: dict, _timeout: float) -> dict:
        raise TimeoutError("upstream timed out")

    failed_report = source_watch.run_checks(
        {**registry, "sources": [registry["sources"][0]]},
        checkers={"first": failing},
    )
    assert failed_report["overall_status"] == "check_failed"
    assert failed_report["exit_code"] == 3
    assert failed_report["sources"][0]["error_type"] == "TimeoutError"


def test_gnomad_identity_change_is_review_notice_not_check_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class IdentityChanged(RuntimeError):
        detail = {
            "dataset": "gnomad_v3",
            "source_type": "sites",
            "field": "etag",
            "expected": "old",
            "observed": "new",
        }

    def changed(_path: Path) -> None:
        raise IdentityChanged("identity changed")

    module = SimpleNamespace(
        GnomadSourceIdentityChangedError=IdentityChanged,
        check_updates=changed,
    )
    monkeypatch.setattr(source_watch, "_load_script_module", lambda *_args: module)

    result = source_watch.check_gnomad_releases(
        {"manifest_path": "backend/data/gnomad/gnomad_panel_manifest.json"},
        1,
    )

    assert result["status"] == "upstream_changed"
    assert result["review_required"] is True
    assert result["identity_change"]["field"] == "etag"


def test_markdown_explains_review_notice_and_candidate_diff() -> None:
    report = {
        "overall_status": "review_required",
        "generated_utc": "2026-11-03T00:00:00Z",
        "sources": [
            {
                "name": "ClinGen ERepo ENIGMA VCEP assertions",
                "checker": "erepo_export",
                "status": "update_available",
                "review_required": True,
                "candidate_workflow": "draft_pull_request",
                "official_url": "https://example.test",
                "activation_command": "python build.py",
                "candidate_diff": {
                    "pinned_count": 1,
                    "remote_count": 2,
                    "added_count": 1,
                    "modified_count": 0,
                    "removed_count": 0,
                    "added": [
                        {
                            "gene": "BRCA1",
                            "variant": "c.1A>G",
                            "classification": "Likely Pathogenic",
                            "published_date": "2026-10-29",
                            "erepo_url": "https://example.test/record",
                        }
                    ],
                    "modified": [],
                    "removed": [],
                },
            }
        ],
    }

    markdown = source_watch.render_markdown(report)

    assert "Candidate record diff" in markdown
    assert "[c.1A>G](https://example.test/record)" in markdown
    assert "review notice, not a CI failure" in markdown
    assert "Prepare source update" in markdown


def test_manual_source_is_visible_but_does_not_claim_freshness() -> None:
    result = source_watch.check_manual_source(
        {
            "review_trigger": "policy change",
            "tracked_files": ["data/example.json"],
        },
        1,
    )

    assert result["status"] == "manual_only"
    assert result["review_required"] is False
    assert result["review_trigger"] == "policy change"
