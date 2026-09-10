from fastapi import HTTPException
from fastapi.testclient import TestClient

from backend.models import ClassificationResult, RnaReviewRecommendation


API_HEADERS = {"X-ARIANE-API-Key": "test-api-key"}


def _result_with_review() -> ClassificationResult:
    return ClassificationResult(
        variant="BRCA1 c.4185G>A p.(Gln1395=)",
        gene="BRCA1",
        c_notation="c.4185G>A",
        p_notation="p.(Gln1395=)",
        predicted_class=3,
        predicted_label="VUS",
        total_points=5,
        rna_review=RnaReviewRecommendation(
            recommended=True,
            priority="high",
            title="RNA evidence review recommended",
            summary="RNA evidence may clarify the transcript consequence.",
            potential_branches=["PVS1 (RNA)", "BP7 (RNA)"],
            is_evidence_criterion=False,
            manual_review_prefill={"assay_scope": "mrna_only"},
        ),
    )


def _client(monkeypatch) -> TestClient:
    from backend import main
    from backend.api_auth import PUBLIC_API_KEY_AUTHENTICATOR

    async def classify_cached(*args, **kwargs):
        return _result_with_review(), "hit", "test-fingerprint"

    monkeypatch.setattr(main, "_classify_one_cached", classify_cached)
    monkeypatch.setattr(main, "CLASSIFICATION_USAGE", None)
    monkeypatch.setattr(main, "_audit", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        PUBLIC_API_KEY_AUTHENTICATOR,
        "authenticate",
        lambda value: "test-client" if value == "test-api-key" else None,
    )
    return TestClient(main.app)


def test_public_api_classification_has_versioned_metadata_and_separate_review(
    monkeypatch,
):
    response = _client(monkeypatch).post(
        "/api/v1/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers={**API_HEADERS, "X-Request-ID": "public-api-test"},
    )

    assert response.status_code == 200
    assert response.headers["X-ARIANE-API-Version"] == "1.0"
    assert response.headers["X-Request-ID"] == "public-api-test"
    payload = response.json()
    assert payload["metadata"]["api_version"] == "1.0"
    assert payload["metadata"]["classifier_fingerprint"] == "test-fingerprint"
    assert payload["metadata"]["cache_status"] == "hit"
    assert payload["metadata"]["policy"]["policy_id"] == "ENIGMA_BRCA_VCEP_1.2"
    assert "rna_review" not in payload["classification"]
    assert payload["classification"]["predicted_class"] == 3
    assert payload["manual_review"]["recommended"] is True
    assert payload["manual_review"]["affects_automatic_classification"] is False
    assert payload["manual_review"]["items"][0]["kind"] == "rna"
    assert payload["manual_review"]["items"][0]["details"][
        "is_evidence_criterion"
    ] is False


def test_public_api_validation_error_is_machine_readable(monkeypatch):
    response = _client(monkeypatch).post(
        "/api/v1/classify",
        json={"gene": "BRCA1", "c_notation": "not-hgvs"},
        headers=API_HEADERS,
    )

    assert response.status_code == 422
    payload = response.json()
    assert payload["error"]["code"] == "invalid_request"
    assert payload["error"]["retryable"] is False
    assert payload["metadata"]["api_version"] == "1.0"


def test_public_api_replaces_unsafe_request_id(monkeypatch):
    response = _client(monkeypatch).post(
        "/api/v1/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers={**API_HEADERS, "X-Request-ID": "unsafe request id for logs"},
    )

    assert response.status_code == 200
    request_id = response.headers["X-Request-ID"]
    assert request_id != "unsafe request id for logs"
    assert len(request_id) == 32
    assert request_id == response.json()["metadata"]["request_id"]


def test_public_api_batch_keeps_valid_items_and_codes_invalid_items(monkeypatch):
    response = _client(monkeypatch).post(
        "/api/v1/classify/batch",
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
    assert payload["success_count"] == 1
    assert payload["error_count"] == 1
    assert payload["results"][0]["status"] == "ok"
    assert payload["results"][0]["manual_review"]["recommended"] is True
    assert payload["results"][1]["status"] == "error"
    assert payload["results"][1]["error"]["code"] == "invalid_variant"
    assert payload["results"][1]["error"]["retryable"] is False


def test_public_api_marks_required_spliceai_timeout_as_retryable_item_error(
    monkeypatch,
):
    from backend import main
    from backend.api_auth import PUBLIC_API_KEY_AUTHENTICATOR

    async def unavailable(*args, **kwargs):
        raise HTTPException(
            status_code=503,
            detail="SpliceAI is required but temporarily unavailable",
            headers={
                "X-ARIANE-Error-Code": "spliceai_temporarily_unavailable",
                "X-ARIANE-Retryable": "true",
                "Retry-After": "5",
            },
        )

    monkeypatch.setattr(main, "_classify_one_cached", unavailable)
    monkeypatch.setattr(main, "CLASSIFICATION_USAGE", None)
    monkeypatch.setattr(main, "_audit", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        PUBLIC_API_KEY_AUTHENTICATOR,
        "authenticate",
        lambda value: "test-client" if value == "test-api-key" else None,
    )
    response = TestClient(main.app).post(
        "/api/v1/classify/batch",
        json={
            "variants": [{"gene": "BRCA1", "c_notation": "c.5366C>T"}],
        },
        headers=API_HEADERS,
    )

    assert response.status_code == 200
    item = response.json()["results"][0]
    assert item["status"] == "error"
    assert item["error"]["code"] == "spliceai_temporarily_unavailable"
    assert item["error"]["retryable"] is True
    assert item["error"]["message"] == (
        "SpliceAI is required but temporarily unavailable"
    )
    assert "classification" not in item


def test_public_api_capabilities_publish_limits_and_policy():
    from backend import main

    response = TestClient(main.app).get("/api/v1/capabilities")

    assert response.status_code == 200
    payload = response.json()
    assert payload["status"] == "beta"
    assert payload["limits"]["maximum_batch_items"] == 10
    assert payload["limits"]["recommended_uncached_batch_items"] == 5
    assert payload["limits"]["per_key_requests_per_minute"] == 30
    assert payload["limits"]["per_key_request_burst"] == 3
    assert payload["authentication_required"] is True
    assert payload["authentication_header"] == "X-ARIANE-API-Key"
    assert {gene["symbol"] for gene in payload["supported_genes"]} == {
        "BRCA1",
        "BRCA2",
    }


def test_openapi_documents_api_key_on_classification_routes():
    from backend import main

    schema = main.app.openapi()
    assert schema["paths"]["/api/v1/classify"]["post"]["security"] == [
        {"ArianeApiKey": []}
    ]
    assert "security" not in schema["paths"]["/api/v1/capabilities"]["get"]
    security_scheme = schema["components"]["securitySchemes"]["ArianeApiKey"]
    assert security_scheme["in"] == "header"
    assert security_scheme["name"] == "X-ARIANE-API-Key"


def test_public_api_rejects_more_than_synchronous_batch_limit(monkeypatch):
    from backend import main
    from backend.api_auth import PUBLIC_API_KEY_AUTHENTICATOR

    monkeypatch.setattr(
        PUBLIC_API_KEY_AUTHENTICATOR,
        "authenticate",
        lambda value: "test-client" if value == "test-api-key" else None,
    )

    response = TestClient(main.app).post(
        "/api/v1/classify/batch",
        json={
            "variants": [
                {"gene": "BRCA1", "c_notation": "c.4185G>A"}
                for _ in range(11)
            ]
        },
        headers=API_HEADERS,
    )

    assert response.status_code == 422
    assert response.json()["error"]["code"] == "invalid_request"
    assert "at most 10 variants" in response.json()["error"]["message"]


def test_public_api_requires_key_for_classification(monkeypatch):
    response = _client(monkeypatch).post(
        "/api/v1/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
    )

    assert response.status_code == 401
    assert response.json()["error"]["code"] == "api_key_required"


def test_public_api_rejects_invalid_key(monkeypatch):
    response = _client(monkeypatch).post(
        "/api/v1/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers={"X-ARIANE-API-Key": "wrong-key"},
    )

    assert response.status_code == 401
    assert response.json()["error"]["code"] == "invalid_api_key"


def test_public_api_fails_closed_when_key_registry_is_missing(
    tmp_path,
    monkeypatch,
):
    from backend import main

    monkeypatch.setenv(
        "ARIANE_API_KEYS_FILE",
        str(tmp_path / "missing-api-keys.json"),
    )
    response = TestClient(main.app).post(
        "/api/v1/classify",
        json={"gene": "BRCA1", "c_notation": "c.4185G>A"},
        headers={"X-ARIANE-API-Key": "unverifiable-key"},
    )

    assert response.status_code == 503
    assert response.json()["error"]["code"] == "api_authentication_unavailable"
    assert response.json()["error"]["retryable"] is False
