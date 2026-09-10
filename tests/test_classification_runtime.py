from datetime import datetime, timedelta, timezone
from pathlib import Path
import sqlite3

from fastapi import FastAPI, Request, Response
from fastapi.testclient import TestClient

from backend.classification_runtime.cache import ClassificationCacheRepository
from backend.classification_runtime.identity import resolve_usage_identity
from backend.classification_runtime.usage import ClassificationUsageRepository
from backend.models import ClassificationResult


def _result() -> ClassificationResult:
    return ClassificationResult(
        variant="BRCA1 c.4185G>A p.(Gln1395=)",
        gene="BRCA1",
        c_notation="c.4185G>A",
        p_notation="p.(Gln1395=)",
        predicted_class=3,
        predicted_label="Uncertain significance",
    )


def _cache_arguments(fingerprint: str = "policy-a") -> dict:
    return {
        "gene": "BRCA1",
        "transcript": "NM_007294.4",
        "c_notation": "c.4185G>A",
        "p_notation": "p.(Gln1395=)",
        "dup_type": "Unknown",
        "fingerprint": fingerprint,
    }


def test_classification_cache_requires_exact_fingerprint(tmp_path: Path):
    repository = ClassificationCacheRepository(
        tmp_path / "cache.sqlite3", max_age_seconds=3600
    )
    repository.put(
        _result(),
        transcript="NM_007294.4",
        dup_type="Unknown",
        fingerprint="policy-a",
    )

    hit = repository.get(**_cache_arguments("policy-a"))
    assert hit.status == "hit"
    assert hit.result is not None
    assert hit.result.predicted_class == 3
    assert repository.get(**_cache_arguments("policy-b")).status == "miss"


def test_classification_cache_rejects_modified_payload(tmp_path: Path):
    database = tmp_path / "cache.sqlite3"
    repository = ClassificationCacheRepository(database, max_age_seconds=3600)
    repository.put(
        _result(),
        transcript="NM_007294.4",
        dup_type="Unknown",
        fingerprint="policy-a",
    )
    with sqlite3.connect(database) as connection:
        connection.execute(
            "UPDATE classification_results SET result_json = ?",
            ('{"predicted_class": 5}',),
        )
        connection.commit()
    assert repository.get(**_cache_arguments()).status == "invalid"


def test_classification_cache_has_no_default_time_expiration(tmp_path: Path, monkeypatch):
    monkeypatch.delenv("ARIANE_CLASSIFICATION_CACHE_MAX_AGE_SECONDS", raising=False)
    database = tmp_path / "cache.sqlite3"
    repository = ClassificationCacheRepository(database)
    repository.put(
        _result(),
        transcript="NM_007294.4",
        dup_type="Unknown",
        fingerprint="policy-a",
    )
    with sqlite3.connect(database) as connection:
        expires_at = connection.execute(
            "SELECT expires_at FROM classification_results"
        ).fetchone()[0]
    assert expires_at == ""
    assert repository.get(**_cache_arguments()).status == "hit"


def test_usage_events_count_every_search_separately(tmp_path: Path):
    repository = ClassificationUsageRepository(tmp_path / "usage.sqlite3")
    common = {
        "actor_id": "jane",
        "actor_type": "account",
        "request_mode": "single",
        "gene": "BRCA1",
        "transcript": "NM_007294.4",
        "c_notation": "c.4185G>A",
        "status": "completed",
        "predicted_class": 3,
        "total_points": 0,
        "duration_ms": 12.0,
        "policy_id": "ENIGMA_BRCA_VCEP_1.2",
        "policy_version": "1.2.0",
        "classifier_fingerprint": "fingerprint",
    }
    repository.record(request_id="one", cache_status="miss", **common)
    repository.record(request_id="two", cache_status="hit", **common)

    now = datetime.now(timezone.utc)
    summary = repository.summary(
        start=now - timedelta(minutes=5),
        end=now + timedelta(minutes=5),
    )
    assert summary["searches"] == 2
    assert summary["actors"] == 1
    assert summary["cache_hits"] == 1
    assert summary["top_variants"] == [
        {"variant_key": "BRCA1:c.4185G>A", "searches": 2}
    ]


def test_identity_uses_cookie_without_storing_ip_and_supports_trusted_header(monkeypatch):
    app = FastAPI()

    @app.get("/identity")
    async def identity(request: Request, response: Response):
        value = resolve_usage_identity(request)
        value.set_cookie(response)
        return {"actor_id": value.actor_id, "actor_type": value.actor_type}

    client = TestClient(app)
    first = client.get("/identity")
    second = client.get("/identity")
    assert first.json()["actor_type"] == "visitor"
    assert second.json()["actor_id"] == first.json()["actor_id"]

    monkeypatch.setenv("ARIANE_TRUSTED_USER_HEADER", "X-ARIANE-User")
    account = client.get("/identity", headers={"X-ARIANE-User": "jana@example.org"})
    assert account.json() == {"actor_id": "jana@example.org", "actor_type": "account"}


def test_classification_endpoint_caches_result_but_counts_both_searches(tmp_path, monkeypatch):
    from backend import main

    cache = ClassificationCacheRepository(
        tmp_path / "cache.sqlite3", max_age_seconds=3600
    )
    usage = ClassificationUsageRepository(tmp_path / "usage.sqlite3")
    calls = 0

    async def classify_once(*args, **kwargs):
        nonlocal calls
        calls += 1
        return _result()

    monkeypatch.setattr(main, "CLASSIFICATION_CACHE", cache)
    monkeypatch.setattr(main, "CLASSIFICATION_USAGE", usage)
    monkeypatch.setattr(main, "_classify_one", classify_once)
    monkeypatch.setattr(main, "classification_fingerprint", lambda gene, mode: "test-fingerprint")

    client = TestClient(main.app)
    payload = {
        "gene": "BRCA1",
        "c_notation": "c.4185G>A",
        "p_notation": "p.(Gln1395=)",
    }
    assert client.post("/api/classify", json=payload).status_code == 200
    assert client.post("/api/classify", json=payload).status_code == 200
    assert calls == 1

    now = datetime.now(timezone.utc)
    summary = usage.summary(
        start=now - timedelta(minutes=5),
        end=now + timedelta(minutes=5),
    )
    assert summary["searches"] == 2
    assert summary["actors"] == 1
    assert summary["cache_hits"] == 1
