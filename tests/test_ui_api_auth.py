from fastapi.testclient import TestClient
import pytest

from backend.models import ClassificationResult
from backend.ui_session import (
    UI_SESSION_AUTHENTICATOR,
    UI_SESSION_COOKIE,
    UI_SESSION_MAX_AGE_SECONDS,
    UI_SESSION_SECRET_ENVIRONMENT,
    UiSessionConfigurationError,
    UiSessionAuthenticator,
)


API_HEADERS = {"X-ARIANE-API-Key": "test-api-key"}


def _result() -> ClassificationResult:
    return ClassificationResult(
        variant="BRCA1 c.4185G>A p.(Gln1395=)",
        gene="BRCA1",
        c_notation="c.4185G>A",
        p_notation="p.(Gln1395=)",
        predicted_class=3,
        predicted_label="VUS",
        total_points=0,
    )


def _mock_classification(monkeypatch) -> None:
    from backend import main

    async def classify_cached(*args, **kwargs):
        return _result(), "hit", "test-fingerprint"

    monkeypatch.setattr(main, "_classify_one_cached", classify_cached)
    monkeypatch.setattr(main, "CLASSIFICATION_USAGE", None)
    monkeypatch.setattr(main, "_audit", lambda *args, **kwargs: None)


def _allow_test_api_key(monkeypatch) -> None:
    from backend.api_auth import PUBLIC_API_KEY_AUTHENTICATOR

    monkeypatch.setattr(
        PUBLIC_API_KEY_AUTHENTICATOR,
        "authenticate",
        lambda value: "test-client" if value == "test-api-key" else None,
    )


def test_ui_session_token_expires():
    authenticator = UiSessionAuthenticator()
    token = authenticator.issue(now=1_000)

    assert authenticator.validate(token, now=1_000)
    assert authenticator.validate(
        token,
        now=1_000 + UI_SESSION_MAX_AGE_SECONDS,
    )
    assert not authenticator.validate(
        token,
        now=1_001 + UI_SESSION_MAX_AGE_SECONDS,
    )
    assert not authenticator.validate(f"{token}changed", now=1_000)


def test_ui_session_secret_is_required(monkeypatch):
    monkeypatch.delenv(UI_SESSION_SECRET_ENVIRONMENT, raising=False)

    with pytest.raises(UiSessionConfigurationError):
        UiSessionAuthenticator()


def test_root_issues_httponly_same_site_ui_session():
    from backend import main

    response = TestClient(main.app).get("/")

    assert response.status_code == 200
    cookie = response.headers["set-cookie"]
    assert f"{UI_SESSION_COOKIE}=" in cookie
    assert "HttpOnly" in cookie
    assert "SameSite=strict" in cookie


def test_ui_classification_requires_session_even_with_api_key(monkeypatch):
    from backend import main

    _mock_classification(monkeypatch)
    response = TestClient(main.app).post(
        "/ui-api/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers=API_HEADERS,
    )

    assert response.status_code == 401
    assert response.json()["detail"]["code"] == "ui_session_required"


def test_ui_classification_works_after_loading_page(monkeypatch):
    from backend import main

    _mock_classification(monkeypatch)
    client = TestClient(main.app)
    assert client.get("/").status_code == 200

    response = client.post(
        "/ui-api/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers={
            "Origin": "http://testserver",
            "Sec-Fetch-Site": "same-origin",
        },
    )

    assert response.status_code == 200
    assert response.json()["predicted_class"] == 3


def test_ui_classification_rejects_cross_site_browser_request(monkeypatch):
    from backend import main

    _mock_classification(monkeypatch)
    client = TestClient(main.app)
    client.cookies.set(UI_SESSION_COOKIE, UI_SESSION_AUTHENTICATOR.issue())

    response = client.post(
        "/ui-api/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers={
            "Origin": "https://external.example.org",
            "Sec-Fetch-Site": "cross-site",
        },
    )

    assert response.status_code == 403
    assert response.json()["detail"]["code"] == "cross_site_ui_request_rejected"


def test_legacy_classification_requires_api_key(monkeypatch):
    from backend import main

    _mock_classification(monkeypatch)
    client = TestClient(main.app)
    assert client.post(
        "/api/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
    ).status_code == 401


def test_legacy_classification_accepts_valid_api_key(monkeypatch):
    from backend import main

    _mock_classification(monkeypatch)
    _allow_test_api_key(monkeypatch)
    response = TestClient(main.app).post(
        "/api/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers=API_HEADERS,
    )

    assert response.status_code == 200


def test_non_health_api_data_routes_require_api_key():
    from backend import main

    client = TestClient(main.app)
    assert client.get("/api/health").status_code == 200
    assert client.get("/api/v1/capabilities").status_code == 200
    assert client.get("/api/resources").status_code == 401
    assert client.get("/api/rules").status_code == 401
    assert client.post(
        "/api/normalize",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
    ).status_code == 401


def test_openapi_marks_unversioned_api_routes_as_key_protected():
    from backend import main

    schema = main.app.openapi()
    for path, method in (
        ("/api/classify", "post"),
        ("/api/classify/batch", "post"),
        ("/api/normalize", "post"),
        ("/api/resources", "get"),
        ("/api/rules", "get"),
    ):
        assert schema["paths"][path][method]["security"] == [
            {"ArianeApiKey": []}
        ]
    assert not any(path.startswith("/ui-api/") for path in schema["paths"])
