from fastapi.testclient import TestClient

from backend.models import ClassificationResult


API_HEADERS = {"X-ARIANE-API-Key": "test-api-key"}


def _allow_test_api_key(monkeypatch):
    from backend.api_auth import PUBLIC_API_KEY_AUTHENTICATOR

    monkeypatch.setattr(
        PUBLIC_API_KEY_AUTHENTICATOR,
        "authenticate",
        lambda value: "test-client" if value == "test-api-key" else None,
    )


def _classification_result() -> ClassificationResult:
    return ClassificationResult(
        variant="BRCA1 c.4185G>A p.(Gln1395=)",
        gene="BRCA1",
        c_notation="c.4185G>A",
        p_notation="p.(Gln1395=)",
        predicted_class=3,
        predicted_label="VUS",
        total_points=0,
    )


def test_batch_returns_item_error_without_discarding_valid_variants(monkeypatch):
    from backend import main

    _allow_test_api_key(monkeypatch)

    calls: list[str] = []

    async def classify_cached(gene, c_notation, *args, **kwargs):
        calls.append(f"{gene}:{c_notation}")
        return _classification_result(), "miss", "test-fingerprint"

    monkeypatch.setattr(main, "_classify_one_cached", classify_cached)
    monkeypatch.setattr(main, "CLASSIFICATION_USAGE", None)
    monkeypatch.setattr(main, "_audit", lambda *args, **kwargs: None)

    response = TestClient(main.app).post(
        "/api/classify/batch",
        json={
            "variants": [
                {"gene": "BRCA1", "c_notation": "c.4185G>A"},
                {"gene": "BRCA1", "c_notation": "not-hgvs"},
            ]
        },
        headers=API_HEADERS,
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["total"] == 2
    assert payload["success_count"] == 1
    assert payload["error_count"] == 1
    assert [item["status"] for item in payload["results"]] == ["ok", "error"]
    assert payload["results"][1]["variant"] == "BRCA1 not-hgvs"
    assert "variant description" in payload["results"][1]["error"]
    assert "input_value" not in payload["results"][1]["error"]
    assert calls == ["BRCA1:c.4185G>A"]


def test_batch_keeps_request_level_item_limit(monkeypatch):
    from backend import main

    _allow_test_api_key(monkeypatch)

    response = TestClient(main.app).post(
        "/api/classify/batch",
        json={
            "variants": [
                {"gene": "BRCA1", "c_notation": "c.4185G>A"}
                for _ in range(201)
            ]
        },
        headers=API_HEADERS,
    )

    assert response.status_code == 422
    assert "Maximum 200 variants per batch" in response.text


def test_legacy_batch_requires_api_key():
    from backend import main

    response = TestClient(main.app).post(
        "/api/classify/batch",
        json={"variants": []},
    )

    assert response.status_code == 401
