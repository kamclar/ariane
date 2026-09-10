import json

import pytest

from scripts.ariane_api_client import (
    classify_batch_with_retries,
    read_variants,
    request_json,
)


class _Response:
    def __init__(self, status_code, payload, headers=None):
        self.status_code = status_code
        self._payload = payload
        self.headers = headers or {}
        self.text = json.dumps(payload)
        self.ok = 200 <= status_code < 300

    def json(self):
        return self._payload


class _Session:
    def __init__(self, responses):
        self.responses = iter(responses)
        self.calls = 0

    def request(self, *args, **kwargs):
        self.calls += 1
        return next(self.responses)


def test_reference_client_reads_required_and_optional_tsv_fields(tmp_path):
    source = tmp_path / "variants.tsv"
    source.write_text(
        "gene\tc_notation\tp_notation\n"
        "BRCA1\tc.4185G>A\tp.(Gln1395=)\n",
        encoding="utf-8",
    )

    assert read_variants(source) == [{
        "gene": "BRCA1",
        "c_notation": "c.4185G>A",
        "p_notation": "p.(Gln1395=)",
    }]


def test_reference_client_retries_only_temporary_http_status(monkeypatch):
    session = _Session([
        _Response(503, {"error": {"message": "Evidence unavailable"}}),
        _Response(200, {"status": "ok"}),
    ])
    monkeypatch.setattr("scripts.ariane_api_client.time.sleep", lambda value: None)

    result = request_json(
        session,
        "GET",
        "https://example.test/capabilities",
        request_id="client-test",
    )

    assert result == {"status": "ok"}
    assert session.calls == 2


def test_reference_client_does_not_retry_invalid_request():
    session = _Session([
        _Response(422, {"error": {"message": "Invalid request"}}),
        _Response(200, {"status": "unexpected"}),
    ])

    with pytest.raises(RuntimeError, match="Invalid request"):
        request_json(
            session,
            "POST",
            "https://example.test/classify",
            payload={},
            request_id="client-test",
        )

    assert session.calls == 1


def test_reference_client_retries_only_retryable_batch_items(monkeypatch):
    session = _Session([
        _Response(200, {
            "results": [
                {"index": 0, "status": "ok", "variant": "BRCA1 c.1A>G"},
                {
                    "index": 1,
                    "status": "error",
                    "variant": "BRCA1 c.2A>G",
                    "error": {
                        "code": "spliceai_temporarily_unavailable",
                        "message": "temporary",
                        "retryable": True,
                    },
                },
                {
                    "index": 2,
                    "status": "error",
                    "variant": "BRCA1 c.3A>G",
                    "error": {
                        "code": "spliceai_coordinates_unavailable",
                        "message": "coordinates missing",
                        "retryable": False,
                    },
                },
            ],
        }),
        _Response(200, {
            "results": [
                {"index": 0, "status": "ok", "variant": "BRCA1 c.2A>G"},
            ],
        }),
    ])
    monkeypatch.setattr("scripts.ariane_api_client.time.sleep", lambda value: None)

    results = classify_batch_with_retries(
        session,
        "https://example.test/api/v1",
        [
            {"gene": "BRCA1", "c_notation": "c.1A>G"},
            {"gene": "BRCA1", "c_notation": "c.2A>G"},
            {"gene": "BRCA1", "c_notation": "c.3A>G"},
        ],
        offset=10,
    )

    assert session.calls == 2
    assert [item["input_index"] for item in results] == [10, 11, 12]
    assert [item["status"] for item in results] == ["ok", "ok", "error"]
