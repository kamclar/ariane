"""Version-bound cache of completed automatic classification results."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import hashlib
import json
import os
from pathlib import Path
import sqlite3
from typing import Iterator

from backend.contracts import ClassificationResult
from backend.infrastructure.cache_registry import register_runtime_cache
from backend.policy.classification_scope import automatic_classification_is_supported
from backend.policy.spliceai import spliceai_result_is_complete
from backend.infrastructure.runtime_cache import runtime_cache_path


SCHEMA_VERSION = 2
DEFAULT_DATABASE_NAME = "classification_results.sqlite3"


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _is_complete(result: ClassificationResult) -> bool:
    if not automatic_classification_is_supported(result.variant_type):
        return False
    audit = result.spliceai_audit
    assessed_complete = spliceai_result_is_complete(
        result.variant_type,
        status=audit.status if audit is not None else "",
        score=audit.score if audit is not None else None,
    )
    if not assessed_complete:
        return False
    if audit is None:
        return True
    return not audit.reference_lookup_required or audit.reference_lookup_complete


def _cache_key(
    *,
    gene: str,
    transcript: str,
    c_notation: str,
    p_notation: str,
    dup_type: str,
    fingerprint: str,
) -> str:
    payload = _canonical_json({
        "gene": gene,
        "transcript": transcript,
        "c_notation": c_notation,
        "p_notation": p_notation,
        "dup_type": dup_type,
        "fingerprint": fingerprint,
    })
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _configured_max_age() -> int | None:
    raw = os.getenv("ARIANE_CLASSIFICATION_CACHE_MAX_AGE_SECONDS", "").strip()
    if not raw:
        return None
    try:
        value = int(raw)
    except ValueError:
        return None
    return max(1, value)


@dataclass(frozen=True)
class ClassificationCacheResult:
    status: str
    result: ClassificationResult | None = None


class ClassificationCacheRepository:
    """Store only completed Module 1 results and reject stale cache entries."""

    def __init__(
        self,
        database_path: Path | None = None,
        max_age_seconds: int | None = None,
    ):
        self.database_path = database_path or runtime_cache_path(DEFAULT_DATABASE_NAME)
        self.max_age_seconds = (
            _configured_max_age()
            if max_age_seconds is None
            else max(1, max_age_seconds)
        )
        self.enabled = os.getenv("ARIANE_CLASSIFICATION_CACHE_ENABLED", "1").strip() != "0"
        if self.enabled:
            self.database_path.parent.mkdir(parents=True, exist_ok=True)
            self._initialize()

    @contextmanager
    def _connection(self) -> Iterator[sqlite3.Connection]:
        connection = sqlite3.connect(self.database_path, timeout=10)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA busy_timeout = 10000")
        try:
            yield connection
        finally:
            connection.close()

    def _initialize(self) -> None:
        with self._connection() as connection:
            connection.executescript(
                """
                PRAGMA journal_mode = WAL;
                CREATE TABLE IF NOT EXISTS classification_results (
                    cache_key TEXT PRIMARY KEY,
                    schema_version INTEGER NOT NULL,
                    variant_key TEXT NOT NULL,
                    gene TEXT NOT NULL,
                    transcript TEXT NOT NULL,
                    c_notation TEXT NOT NULL,
                    p_notation TEXT NOT NULL,
                    dup_type TEXT NOT NULL,
                    fingerprint TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    expires_at TEXT NOT NULL,
                    last_accessed_at TEXT NOT NULL,
                    hit_count INTEGER NOT NULL DEFAULT 0,
                    result_json TEXT NOT NULL,
                    result_sha256 TEXT NOT NULL
                );
                CREATE INDEX IF NOT EXISTS classification_results_variant_idx
                    ON classification_results(variant_key, created_at DESC);
                CREATE INDEX IF NOT EXISTS classification_results_expiry_idx
                    ON classification_results(expires_at);
                """
            )
            connection.commit()

    def get(
        self,
        *,
        gene: str,
        transcript: str,
        c_notation: str,
        p_notation: str,
        dup_type: str,
        fingerprint: str,
    ) -> ClassificationCacheResult:
        if not self.enabled:
            return ClassificationCacheResult("disabled")
        key = _cache_key(
            gene=gene,
            transcript=transcript,
            c_notation=c_notation,
            p_notation=p_notation,
            dup_type=dup_type,
            fingerprint=fingerprint,
        )
        now = _utc_now()
        with self._connection() as connection:
            row = connection.execute(
                "SELECT * FROM classification_results WHERE cache_key = ?",
                (key,),
            ).fetchone()
            if row is None:
                return ClassificationCacheResult("miss")
            if row["schema_version"] != SCHEMA_VERSION:
                return ClassificationCacheResult("invalid")
            raw_expiry = str(row["expires_at"] or "").strip()
            if raw_expiry:
                try:
                    expires_at = datetime.fromisoformat(raw_expiry)
                except ValueError:
                    return ClassificationCacheResult("invalid")
                if expires_at <= now:
                    return ClassificationCacheResult("expired")
            result_json = str(row["result_json"])
            if hashlib.sha256(result_json.encode("utf-8")).hexdigest() != row["result_sha256"]:
                return ClassificationCacheResult("invalid")
            try:
                result = ClassificationResult.model_validate_json(result_json)
            except Exception:
                return ClassificationCacheResult("invalid")
            if not _is_complete(result):
                return ClassificationCacheResult("incomplete")
            connection.execute(
                """
                UPDATE classification_results
                SET last_accessed_at = ?, hit_count = hit_count + 1
                WHERE cache_key = ?
                """,
                (now.isoformat(), key),
            )
            connection.commit()
        return ClassificationCacheResult("hit", result)

    def put(
        self,
        result: ClassificationResult,
        *,
        transcript: str,
        dup_type: str,
        fingerprint: str,
    ) -> None:
        if not self.enabled:
            return
        if not _is_complete(result):
            return
        now = _utc_now()
        expires_at = (
            now + timedelta(seconds=self.max_age_seconds)
            if self.max_age_seconds is not None
            else None
        )
        # ClinVar and ClinGen are volatile, read-only comparisons. They are
        # refreshed independently on a cache hit and must not become part of a
        # long-lived Module 1 assertion snapshot.
        cacheable_result = result.model_copy(update={"external": None})
        result_json = _canonical_json(cacheable_result.model_dump(mode="json"))
        result_sha256 = hashlib.sha256(result_json.encode("utf-8")).hexdigest()
        key = _cache_key(
            gene=result.gene,
            transcript=transcript,
            c_notation=result.c_notation,
            p_notation=result.p_notation,
            dup_type=dup_type,
            fingerprint=fingerprint,
        )
        with self._connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            connection.execute(
                """
                INSERT INTO classification_results (
                    cache_key, schema_version, variant_key, gene, transcript,
                    c_notation, p_notation, dup_type, fingerprint, created_at,
                    expires_at, last_accessed_at, hit_count, result_json,
                    result_sha256
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0, ?, ?)
                ON CONFLICT(cache_key) DO UPDATE SET
                    schema_version = excluded.schema_version,
                    created_at = excluded.created_at,
                    expires_at = excluded.expires_at,
                    last_accessed_at = excluded.last_accessed_at,
                    hit_count = 0,
                    result_json = excluded.result_json,
                    result_sha256 = excluded.result_sha256
                """,
                (
                    key,
                    SCHEMA_VERSION,
                    f"{result.gene}:{result.c_notation}",
                    result.gene,
                    transcript,
                    result.c_notation,
                    result.p_notation,
                    dup_type,
                    fingerprint,
                    now.isoformat(),
                    expires_at.isoformat() if expires_at is not None else "",
                    now.isoformat(),
                    result_json,
                    result_sha256,
                ),
            )
            if self.max_age_seconds is not None:
                connection.execute(
                    """
                    DELETE FROM classification_results
                    WHERE expires_at != '' AND expires_at < ?
                    """,
                    ((now - timedelta(days=7)).isoformat(),),
                )
            connection.commit()

    def clear(self) -> None:
        if not self.enabled:
            return
        with self._connection() as connection:
            connection.execute("DELETE FROM classification_results")
            connection.commit()


def _clear_classification_cache() -> None:
    ClassificationCacheRepository().clear()


register_runtime_cache("classification_results", _clear_classification_cache)
