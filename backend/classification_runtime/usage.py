"""Persistent, privacy-limited events for classification usage statistics."""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import os
from pathlib import Path
import sqlite3
from typing import Any, Iterator

from backend.runtime_data import runtime_data_path
from backend.version import ARIANE_VERSION


SCHEMA_VERSION = 2
DEFAULT_DATABASE_NAME = "classification_usage.sqlite3"


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _retention_days() -> int:
    try:
        return max(1, int(os.getenv("ARIANE_USAGE_RETENTION_DAYS", "365")))
    except ValueError:
        return 365


class ClassificationUsageRepository:
    """Append one event per attempted classification request."""

    def __init__(self, database_path: Path | None = None):
        self.database_path = database_path or runtime_data_path(DEFAULT_DATABASE_NAME)
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
            connection.execute("PRAGMA journal_mode = WAL")
            existing = connection.execute(
                "SELECT sql FROM sqlite_master "
                "WHERE type = 'table' AND name = 'classification_usage_events'"
            ).fetchone()
            if existing is not None and "'api_key'" not in str(existing["sql"]):
                connection.executescript(
                    """
                    BEGIN IMMEDIATE;
                    CREATE TABLE classification_usage_events_v2 (
                        event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                        schema_version INTEGER NOT NULL,
                        occurred_at TEXT NOT NULL,
                        actor_id TEXT NOT NULL,
                        actor_type TEXT NOT NULL CHECK (
                            actor_type IN ('account', 'visitor', 'api_key')
                        ),
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
                    INSERT INTO classification_usage_events_v2 (
                        event_id, schema_version, occurred_at, actor_id, actor_type,
                        request_id, request_mode, variant_key, gene, transcript,
                        c_notation, status, cache_status, predicted_class,
                        total_points, duration_ms, app_version, policy_id,
                        policy_version, classifier_fingerprint
                    )
                    SELECT event_id, 2, occurred_at, actor_id, actor_type,
                        request_id, request_mode, variant_key, gene, transcript,
                        c_notation, status, cache_status, predicted_class,
                        total_points, duration_ms, app_version, policy_id,
                        policy_version, classifier_fingerprint
                    FROM classification_usage_events;
                    DROP TABLE classification_usage_events;
                    ALTER TABLE classification_usage_events_v2
                        RENAME TO classification_usage_events;
                    COMMIT;
                    """
                )
            connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS classification_usage_events (
                    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    schema_version INTEGER NOT NULL,
                    occurred_at TEXT NOT NULL,
                    actor_id TEXT NOT NULL,
                    actor_type TEXT NOT NULL CHECK (
                        actor_type IN ('account', 'visitor', 'api_key')
                    ),
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
                CREATE INDEX IF NOT EXISTS classification_usage_time_idx
                    ON classification_usage_events(occurred_at DESC);
                CREATE INDEX IF NOT EXISTS classification_usage_variant_idx
                    ON classification_usage_events(variant_key, occurred_at DESC);
                CREATE INDEX IF NOT EXISTS classification_usage_actor_idx
                    ON classification_usage_events(actor_id, occurred_at DESC);
                """
            )
            connection.commit()

    def record(
        self,
        *,
        actor_id: str,
        actor_type: str,
        request_id: str,
        request_mode: str,
        gene: str,
        transcript: str,
        c_notation: str,
        status: str,
        cache_status: str,
        predicted_class: int | None,
        total_points: int | None,
        duration_ms: float,
        policy_id: str,
        policy_version: str,
        classifier_fingerprint: str,
    ) -> None:
        now = _utc_now()
        cutoff = now - timedelta(days=_retention_days())
        with self._connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            connection.execute(
                """
                INSERT INTO classification_usage_events (
                    schema_version, occurred_at, actor_id, actor_type, request_id,
                    request_mode, variant_key, gene, transcript, c_notation,
                    status, cache_status, predicted_class, total_points,
                    duration_ms, app_version, policy_id, policy_version,
                    classifier_fingerprint
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    SCHEMA_VERSION,
                    now.isoformat(),
                    actor_id,
                    actor_type,
                    request_id,
                    request_mode,
                    f"{gene}:{c_notation}",
                    gene,
                    transcript,
                    c_notation,
                    status,
                    cache_status,
                    predicted_class,
                    total_points,
                    round(max(0.0, duration_ms), 1),
                    ARIANE_VERSION,
                    policy_id,
                    policy_version,
                    classifier_fingerprint,
                ),
            )
            connection.execute(
                "DELETE FROM classification_usage_events WHERE occurred_at < ?",
                (cutoff.isoformat(),),
            )
            connection.commit()

    def summary(self, *, start: datetime, end: datetime) -> dict[str, Any]:
        parameters = (start.isoformat(), end.isoformat())
        with self._connection() as connection:
            totals = connection.execute(
                """
                SELECT COUNT(*) AS searches,
                       COUNT(DISTINCT actor_type || ':' || actor_id) AS actors,
                       SUM(CASE WHEN cache_status = 'hit' THEN 1 ELSE 0 END) AS cache_hits,
                       SUM(CASE WHEN status = 'error' THEN 1 ELSE 0 END) AS errors
                FROM classification_usage_events
                WHERE occurred_at >= ? AND occurred_at <= ?
                """,
                parameters,
            ).fetchone()
            top_variants = connection.execute(
                """
                SELECT variant_key, COUNT(*) AS searches
                FROM classification_usage_events
                WHERE occurred_at >= ? AND occurred_at <= ?
                GROUP BY variant_key ORDER BY searches DESC, variant_key LIMIT 10
                """,
                parameters,
            ).fetchall()
            top_actors = connection.execute(
                """
                SELECT actor_id, actor_type, COUNT(*) AS searches
                FROM classification_usage_events
                WHERE occurred_at >= ? AND occurred_at <= ?
                GROUP BY actor_type, actor_id ORDER BY searches DESC, actor_id LIMIT 10
                """,
                parameters,
            ).fetchall()
        return {
            "searches": int(totals["searches"] or 0),
            "actors": int(totals["actors"] or 0),
            "cache_hits": int(totals["cache_hits"] or 0),
            "errors": int(totals["errors"] or 0),
            "top_variants": [dict(row) for row in top_variants],
            "top_actors": [dict(row) for row in top_actors],
        }
