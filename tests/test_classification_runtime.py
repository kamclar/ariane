from datetime import datetime, timedelta, timezone
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import sqlite3

from fastapi import FastAPI, Request, Response
from fastapi.testclient import TestClient

from backend.classification_runtime.cache import ClassificationCacheRepository
from backend.classification_runtime.identity import resolve_usage_identity
from backend.classification_runtime.usage import ClassificationUsageRepository
from backend.contracts import ClassificationResult, SpliceAIAudit


def _result() -> ClassificationResult:
    return ClassificationResult(
        variant="BRCA1 c.4185G>A p.(Gln1395=)",
        gene="BRCA1",
        c_notation="c.4185G>A",
        p_notation="p.(Gln1395=)",
        predicted_class=3,
        predicted_label="Uncertain significance",
    )


def _incomplete_figure1a_result() -> ClassificationResult:
    return ClassificationResult(
        variant="BRCA1 c.5366C>T p.(Ala1789Val)",
        gene="BRCA1",
        c_notation="c.5366C>T",
        p_notation="p.(Ala1789Val)",
        variant_type="missense",
        predicted_class=3,
        predicted_label="Uncertain significance",
        spliceai_audit=SpliceAIAudit(
            status="api_error",
            score=None,
            required_for_classification=True,
            retryable=True,
            reason="TimeoutError: The read operation timed out",
        ),
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


def test_classification_cache_does_not_store_incomplete_required_spliceai(tmp_path: Path):
    repository = ClassificationCacheRepository(
        tmp_path / "cache.sqlite3", max_age_seconds=3600
    )
    result = _incomplete_figure1a_result()
    repository.put(
        result,
        transcript="NM_007294.4",
        dup_type="Unknown",
        fingerprint="policy-a",
    )

    hit = repository.get(
        gene="BRCA1",
        transcript="NM_007294.4",
        c_notation="c.5366C>T",
        p_notation="p.(Ala1789Val)",
        dup_type="Unknown",
        fingerprint="policy-a",
    )
    assert hit.status == "miss"


def test_classification_cache_does_not_store_incomplete_ps1_reference_lookup(
    tmp_path: Path,
):
    repository = ClassificationCacheRepository(
        tmp_path / "cache.sqlite3", max_age_seconds=3600
    )
    result = _result().model_copy(
        update={
            "variant_type": "synonymous",
            "spliceai_audit": SpliceAIAudit(
                status="ok",
                score=0.01,
                required_for_classification=True,
                reference_lookup_required=True,
                reference_lookup_complete=False,
                reference_variant_statuses={
                    "c.4185G>C": {
                        "status": "api_error",
                        "score": None,
                        "retryable": True,
                    }
                },
            ),
        }
    )
    repository.put(
        result,
        transcript="NM_007294.4",
        dup_type="Unknown",
        fingerprint="policy-a",
    )

    assert repository.get(**_cache_arguments()).status == "miss"


def test_classification_cache_rejects_preexisting_incomplete_result(tmp_path: Path):
    database = tmp_path / "cache.sqlite3"
    repository = ClassificationCacheRepository(database, max_age_seconds=3600)
    repository.put(
        _result(),
        transcript="NM_007294.4",
        dup_type="Unknown",
        fingerprint="policy-a",
    )
    incomplete = _incomplete_figure1a_result()
    serialized = json.dumps(
        incomplete.model_dump(mode="json"),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    checksum = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
    with sqlite3.connect(database) as connection:
        connection.execute(
            """
            UPDATE classification_results
            SET gene = ?, c_notation = ?, p_notation = ?,
                result_json = ?, result_sha256 = ?
            """,
            (
                incomplete.gene,
                incomplete.c_notation,
                incomplete.p_notation,
                serialized,
                checksum,
            ),
        )
        connection.commit()

    assert repository.get(**_cache_arguments()).status == "incomplete"


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


def test_usage_events_accept_authenticated_api_key_identity(tmp_path: Path):
    repository = ClassificationUsageRepository(tmp_path / "usage.sqlite3")
    repository.record(
        actor_id="external-client-01",
        actor_type="api_key",
        request_id="api-request",
        request_mode="single",
        gene="BRCA1",
        transcript="NM_007294.4",
        c_notation="c.4185G>A",
        status="completed",
        cache_status="hit",
        predicted_class=3,
        total_points=0,
        duration_ms=4.0,
        policy_id="ENIGMA_BRCA_VCEP_1.2",
        policy_version="1.2.0",
        classifier_fingerprint="fingerprint",
    )

    now = datetime.now(timezone.utc)
    summary = repository.summary(
        start=now - timedelta(minutes=5),
        end=now + timedelta(minutes=5),
    )
    assert summary["top_actors"] == [
        {
            "actor_id": "external-client-01",
            "actor_type": "api_key",
            "searches": 1,
        }
    ]


def test_public_api_daily_quota_counts_variants_and_resets_by_utc_day(tmp_path: Path):
    repository = ClassificationUsageRepository(tmp_path / "usage.sqlite3")
    first_day = datetime(2026, 9, 11, 23, 59, tzinfo=timezone.utc)

    first = repository.reserve_public_api_classifications(
        api_key_id="external-client-01",
        units=3,
        limit=5,
        now=first_day,
    )
    rejected = repository.reserve_public_api_classifications(
        api_key_id="external-client-01",
        units=3,
        limit=5,
        now=first_day,
    )
    next_day = repository.reserve_public_api_classifications(
        api_key_id="external-client-01",
        units=3,
        limit=5,
        now=first_day + timedelta(minutes=2),
    )

    assert first.allowed is True
    assert first.used == 3
    assert first.remaining == 2
    assert rejected.allowed is False
    assert rejected.used == 3
    assert rejected.remaining == 2
    assert next_day.allowed is True
    assert next_day.used == 3


def test_public_api_daily_quota_reservation_is_atomic(tmp_path: Path):
    database = tmp_path / "usage.sqlite3"
    repository = ClassificationUsageRepository(database)
    now = datetime(2026, 9, 11, 12, 0, tzinfo=timezone.utc)

    def reserve_once(_: int) -> bool:
        return repository.reserve_public_api_classifications(
            api_key_id="parallel-client",
            units=1,
            limit=10,
            now=now,
        ).allowed

    with ThreadPoolExecutor(max_workers=8) as executor:
        decisions = list(executor.map(reserve_once, range(24)))

    assert sum(decisions) == 10


def test_usage_database_migrates_actor_type_constraint_for_api_keys(tmp_path: Path):
    database = tmp_path / "usage.sqlite3"
    with sqlite3.connect(database) as connection:
        connection.executescript(
            """
            CREATE TABLE classification_usage_events (
                event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                schema_version INTEGER NOT NULL,
                occurred_at TEXT NOT NULL,
                actor_id TEXT NOT NULL,
                actor_type TEXT NOT NULL CHECK (actor_type IN ('account', 'visitor')),
                request_id TEXT NOT NULL,
                request_mode TEXT NOT NULL CHECK (request_mode IN ('single', 'batch')),
                variant_key TEXT NOT NULL,
                gene TEXT NOT NULL,
                transcript TEXT NOT NULL,
                c_notation TEXT NOT NULL,
                status TEXT NOT NULL CHECK (status IN ('completed', 'error')),
                cache_status TEXT NOT NULL,
                predicted_class INTEGER,
                total_points INTEGER,
                duration_ms REAL NOT NULL,
                app_version TEXT NOT NULL,
                policy_id TEXT NOT NULL,
                policy_version TEXT NOT NULL,
                classifier_fingerprint TEXT NOT NULL
            );
            INSERT INTO classification_usage_events (
                schema_version, occurred_at, actor_id, actor_type, request_id,
                request_mode, variant_key, gene, transcript, c_notation,
                status, cache_status, duration_ms, app_version, policy_id,
                policy_version, classifier_fingerprint
            ) VALUES (
                1, '2026-09-11T00:00:00+00:00', 'visitor-1', 'visitor', 'old',
                'single', 'BRCA1:c.1A>G', 'BRCA1', 'NM_007294.4', 'c.1A>G',
                'completed', 'miss', 1.0, '1.9.4', 'policy', '1.2.0', 'fp'
            );
            """
        )

    repository = ClassificationUsageRepository(database)
    with sqlite3.connect(database) as connection:
        rows = connection.execute(
            "SELECT schema_version, actor_id, actor_type "
            "FROM classification_usage_events ORDER BY event_id"
        ).fetchall()
        table_sql = connection.execute(
            "SELECT sql FROM sqlite_master "
            "WHERE type = 'table' AND name = 'classification_usage_events'"
        ).fetchone()[0]

    assert rows == [(2, "visitor-1", "visitor")]
    assert "'api_key'" in table_sql
    assert repository.database_path == database


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


def test_usage_identity_prefers_authenticated_api_key_id():
    app = FastAPI()

    @app.get("/identity")
    async def identity(request: Request):
        request.state.api_key_id = "batch-client-01"
        value = resolve_usage_identity(request)
        return {"actor_id": value.actor_id, "actor_type": value.actor_type}

    response = TestClient(app).get("/identity")
    assert response.json() == {
        "actor_id": "batch-client-01",
        "actor_type": "api_key",
    }


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
    monkeypatch.setattr(main.CLASSIFICATION_API, "classify_uncached", classify_once)
    monkeypatch.setattr(
        "backend.services.variant_classification_service.classification_fingerprint",
        lambda gene: "test-fingerprint",
    )

    client = TestClient(main.app)
    assert client.get("/").status_code == 200
    payload = {
        "gene": "BRCA1",
        "c_notation": "c.4185G>A",
        "p_notation": "p.(Gln1395=)",
    }
    headers = {"Origin": "http://testserver"}
    assert client.post(
        "/ui-api/classify", json=payload, headers=headers
    ).status_code == 200
    assert client.post(
        "/ui-api/classify", json=payload, headers=headers
    ).status_code == 200
    assert calls == 1

    now = datetime.now(timezone.utc)
    summary = usage.summary(
        start=now - timedelta(minutes=5),
        end=now + timedelta(minutes=5),
    )
    assert summary["searches"] == 2
    assert summary["actors"] == 1
    assert summary["cache_hits"] == 1
