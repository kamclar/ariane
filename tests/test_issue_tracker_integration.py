from fastapi.testclient import TestClient

from backend.main import app
from backend.version import ARIANE_VERSION


def test_resources_expose_configured_issue_tracker_context(monkeypatch):
    monkeypatch.setenv("ARIANE_ISSUE_TRACKER_URL", "https://bugs.example.org/bug_report_page.php")
    monkeypatch.setenv("ARIANE_BUILD_REVISION", "abc1234")

    response = TestClient(app).get("/api/resources")

    assert response.status_code == 200
    payload = response.json()
    assert payload["version"] == ARIANE_VERSION
    assert payload["build_revision"] == "abc1234"
    assert payload["issue_tracker_url"] == "https://bugs.example.org/bug_report_page.php"


def test_resources_hide_unconfigured_issue_tracker(monkeypatch):
    monkeypatch.delenv("ARIANE_ISSUE_TRACKER_URL", raising=False)
    monkeypatch.delenv("ARIANE_BUILD_REVISION", raising=False)

    payload = TestClient(app).get("/api/resources").json()

    assert payload["build_revision"] == ""
    assert payload["issue_tracker_url"] == ""
