"""Crash-safe, bounded-memory checkpoints for the SpliceAI distance experiment.

Only completed cases are written. Both successful and failed cases are retained;
the experiment runner decides whether to explicitly retry a failed case. A
successful result is immutable, except for an identical, idempotent write.

The caller supplies all inference and cohort provenance in ``identity``. Resume
requires the exact same canonical JSON identity, including any pinned manifest,
input hashes, code hashes, and run configuration supplied by the caller.
"""

from __future__ import annotations

from collections.abc import Iterator
import hashlib
import json
from pathlib import Path
import sqlite3
from typing import Any


SCHEMA_VERSION = 1


def _canonical_json(value: dict[str, Any]) -> str:
    """Reject lossy or non-standard JSON rather than weakening identity checks."""
    def validate(item: Any) -> None:
        if isinstance(item, dict):
            for key, child in item.items():
                if not isinstance(key, str):
                    raise ValueError("JSON object keys must be strings")
                validate(child)
        elif isinstance(item, list):
            for child in item:
                validate(child)
        elif item is not None and not isinstance(item, (str, int, float, bool)):
            raise ValueError(f"Unsupported JSON value: {type(item).__name__}")

    if not isinstance(value, dict):
        raise ValueError("Expected a JSON object")
    try:
        validate(value)
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, RecursionError) as exc:
        raise ValueError("Expected finite, non-recursive JSON data") from exc


def _validate_case_id(case_id: str) -> None:
    if not isinstance(case_id, str) or not case_id.strip():
        raise ValueError("case_id must be a non-empty string")


def _validate_result(case_id: str, result: dict[str, Any]) -> None:
    _validate_case_id(case_id)
    if not isinstance(result, dict) or result.get("id") != case_id:
        raise ValueError("Result id must match case_id")
    if result.get("status") not in ("ok", "error"):
        raise ValueError("Only completed results with status 'ok' or 'error' may be stored")


class ResultStore:
    """SQLite-backed per-case results, committed independently with full durability.

    ``resume=False`` refuses an existing path; ``resume=True`` refuses a missing,
    incompatible, or differently identified database. Failed results may be
    replaced by an explicit ``put``; the store never silently filters them out.
    """

    def __init__(self, path: Path, identity: dict[str, Any], *, resume: bool = False):
        self.path = Path(path).resolve()
        self._connection: sqlite3.Connection | None = None
        identity_json = _canonical_json(identity)
        identity_digest = hashlib.sha256(identity_json.encode("utf-8")).hexdigest()

        if resume:
            if not self.path.is_file():
                raise FileNotFoundError(f"Cannot resume: checkpoint does not exist: {self.path}")
            # Validate before opening writable or setting journal pragmas. A bad
            # identity must not initialize or otherwise modify an existing DB.
            readonly = sqlite3.connect(f"{self.path.as_uri()}?mode=ro", uri=True)
            try:
                self._check_identity(readonly, identity_json, identity_digest)
            finally:
                readonly.close()
        else:
            try:
                # Exclusive creation prevents accidentally reopening an existing
                # run even if another process creates it at the same instant.
                with self.path.open("xb"):
                    pass
            except FileExistsError as exc:
                raise FileExistsError(
                    f"Checkpoint already exists; use resume=True: {self.path}"
                ) from exc

        connection = sqlite3.connect(f"{self.path.as_uri()}?mode=rw", uri=True)
        try:
            if resume:
                self._check_identity(connection, identity_json, identity_digest)
            journal_mode = connection.execute("PRAGMA journal_mode=WAL").fetchone()[0]
            if journal_mode != "wal":
                raise ValueError("Checkpoint filesystem does not support SQLite WAL durability")
            connection.execute("PRAGMA synchronous=FULL")
            if not resume:
                with connection:
                    connection.execute("BEGIN IMMEDIATE")
                    connection.execute(
                        "CREATE TABLE metadata ("
                        "singleton INTEGER PRIMARY KEY CHECK (singleton = 1), "
                        "schema_version INTEGER NOT NULL, identity_json TEXT NOT NULL, "
                        "identity_sha256 TEXT NOT NULL)"
                    )
                    connection.execute(
                        "INSERT INTO metadata VALUES (1, ?, ?, ?)",
                        (SCHEMA_VERSION, identity_json, identity_digest),
                    )
                    connection.execute(
                        "CREATE TABLE cases (case_id TEXT PRIMARY KEY NOT NULL, "
                        "status TEXT NOT NULL CHECK (status IN ('ok', 'error')), "
                        "result_json TEXT NOT NULL)"
                    )
                    # Counts scan this compact covering index instead of all
                    # potentially large per-case JSON pages.
                    connection.execute("CREATE INDEX cases_status ON cases (status)")
            self._connection = connection
        except BaseException:
            connection.close()
            raise

    @staticmethod
    def _check_identity(
        connection: sqlite3.Connection, identity_json: str, identity_digest: str
    ) -> None:
        try:
            rows = connection.execute(
                "SELECT singleton, schema_version, identity_json, identity_sha256 FROM metadata"
            ).fetchall()
            connection.execute("SELECT case_id, status, result_json FROM cases LIMIT 0")
        except sqlite3.DatabaseError as exc:
            raise ValueError("Not a valid SpliceAI distance checkpoint database") from exc
        if len(rows) != 1 or rows[0][0] != 1 or rows[0][1] != SCHEMA_VERSION:
            raise ValueError("Unsupported or invalid checkpoint schema version")
        if rows[0][2:] != (identity_json, identity_digest):
            raise ValueError("Checkpoint identity mismatch; use the original inputs/configuration")

    def _db(self) -> sqlite3.Connection:
        if self._connection is None:
            raise RuntimeError("Result store is closed")
        return self._connection

    @staticmethod
    def _decode(case_id: str, status: str, payload: str) -> dict[str, Any]:
        result = json.loads(payload)
        _validate_result(case_id, result)
        if result["status"] != status:
            raise ValueError(f"Checkpoint status disagrees with stored result for {case_id}")
        return result

    def get(self, case_id: str) -> dict[str, Any] | None:
        """Return a completed result, including a failure, or None if unfinished."""
        _validate_case_id(case_id)
        row = self._db().execute(
            "SELECT status, result_json FROM cases WHERE case_id = ?", (case_id,)
        ).fetchone()
        return None if row is None else self._decode(case_id, *row)

    def put(self, case_id: str, result: dict[str, Any]) -> None:
        """Commit a whole case atomically; caller must explicitly authorize retries."""
        _validate_result(case_id, result)
        payload = _canonical_json(result)
        connection = self._db()
        with connection:
            # Serialize the read/update pair so a second writer cannot replace
            # a successful result created after our initial read.
            connection.execute("BEGIN IMMEDIATE")
            previous = connection.execute(
                "SELECT status, result_json FROM cases WHERE case_id = ?", (case_id,)
            ).fetchone()
            if previous is not None:
                self._decode(case_id, *previous)
                if previous[0] == "ok":
                    if previous[1] != payload:
                        raise ValueError(f"Cannot overwrite successful result for {case_id}")
                    return
            connection.execute(
                "INSERT INTO cases (case_id, status, result_json) VALUES (?, ?, ?) "
                "ON CONFLICT(case_id) DO UPDATE SET "
                "status = excluded.status, result_json = excluded.result_json",
                (case_id, result["status"], payload),
            )

    def results(self) -> Iterator[tuple[str, dict[str, Any]]]:
        """Stream all results in insertion order without loading the cohort."""
        cursor = self._db().execute(
            "SELECT case_id, status, result_json FROM cases ORDER BY rowid"
        )
        try:
            for case_id, status, payload in cursor:
                yield case_id, self._decode(case_id, status, payload)
        finally:
            cursor.close()

    def counts(self) -> dict[str, int]:
        """Count both successes and failures; total never means successful only."""
        counts = {"total": 0, "ok": 0, "error": 0}
        for status, count in self._db().execute("SELECT status, COUNT(*) FROM cases GROUP BY status"):
            if status not in ("ok", "error"):
                raise ValueError(f"Invalid stored result status: {status}")
            counts[status] = count
            counts["total"] += count
        return counts

    def close(self) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def __enter__(self) -> ResultStore:
        self._db()
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()
