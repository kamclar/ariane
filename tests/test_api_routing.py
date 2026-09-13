"""Contract tests for paired browser and public API route registration."""

from __future__ import annotations

from fastapi.testclient import TestClient


PAIRED_PATHS = (
    "/classify",
    "/resources",
    "/rules",
    "/rules/trees/{tree_id}",
    "/rules/tables/table9",
    "/rules/tables/{table_id}",
    "/audit/client-validation",
    "/manual-evidence/evaluate",
    "/manual-evidence/status",
    "/manual-evidence/resolve-ps1-reference",
    "/normalize",
)


def test_paired_routes_register_both_paths_and_hide_ui_from_openapi() -> None:
    from backend.main import app

    paths = {route.path for route in app.routes}
    for path in PAIRED_PATHS:
        assert f"/api{path}" in paths
        assert f"/ui-api{path}" in paths

    openapi_paths = app.openapi()["paths"]
    assert not any(path.startswith("/ui-api/") for path in openapi_paths)


def test_paired_routes_keep_the_two_authentication_gates_separate() -> None:
    from backend.main import app

    client = TestClient(app)
    assert client.get("/api/rules").status_code == 401
    assert client.get("/ui-api/rules").status_code == 401

    browser = TestClient(app)
    assert browser.get("/").status_code == 200
    assert browser.get("/ui-api/rules").status_code == 200
    assert browser.get("/api/rules").status_code == 401
