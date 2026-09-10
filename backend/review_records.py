"""Immutable, versioned persistence for manually reviewed classifications."""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import uuid
from typing import Any, Iterator, Mapping

from backend.gene_policy import get_gene_policy
from backend.runtime_data import runtime_data_path
from backend.version import ARIANE_VERSION


SCHEMA_VERSION = "1.0"
DEFAULT_DATABASE_NAME = "review_records.sqlite3"


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class ReviewRecordRepository:
    """Store each draft and approval as a new immutable database row."""

    def __init__(self, database_path: Path | None = None):
        self.database_path = database_path or runtime_data_path(DEFAULT_DATABASE_NAME)
        self.database_path.parent.mkdir(parents=True, exist_ok=True)
        self._initialize()

    @contextmanager
    def _connection(self) -> Iterator[sqlite3.Connection]:
        connection = sqlite3.connect(self.database_path, timeout=10)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
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
                CREATE TABLE IF NOT EXISTS review_records (
                    record_id TEXT PRIMARY KEY,
                    schema_version TEXT NOT NULL,
                    variant_key TEXT NOT NULL,
                    version INTEGER NOT NULL,
                    parent_record_id TEXT,
                    status TEXT NOT NULL CHECK (status IN ('draft', 'approved')),
                    created_at TEXT NOT NULL,
                    created_by_account TEXT NOT NULL,
                    reviewer_name TEXT NOT NULL,
                    reviewer_role TEXT NOT NULL,
                    review_date TEXT NOT NULL,
                    approval_comment TEXT NOT NULL DEFAULT '',
                    approved_at TEXT,
                    app_version TEXT NOT NULL,
                    policy_id TEXT NOT NULL,
                    policy_name TEXT NOT NULL,
                    policy_version TEXT NOT NULL,
                    policy_source_url TEXT NOT NULL,
                    module1_json TEXT NOT NULL,
                    manual_evidence_json TEXT NOT NULL,
                    amended_result_json TEXT NOT NULL,
                    module1_sha256 TEXT NOT NULL,
                    manual_evidence_sha256 TEXT NOT NULL,
                    amended_result_sha256 TEXT NOT NULL,
                    content_sha256 TEXT NOT NULL,
                    record_sha256 TEXT NOT NULL,
                    FOREIGN KEY(parent_record_id) REFERENCES review_records(record_id),
                    UNIQUE(variant_key, version)
                );
                CREATE INDEX IF NOT EXISTS review_records_variant_idx
                    ON review_records(variant_key, version DESC);
                CREATE INDEX IF NOT EXISTS review_records_status_idx
                    ON review_records(status, created_at DESC);
                """
            )
            connection.commit()

    @staticmethod
    def _verify_row(row: sqlite3.Row) -> None:
        section_values = {
            "module1_sha256": row["module1_json"],
            "manual_evidence_sha256": row["manual_evidence_json"],
            "amended_result_sha256": row["amended_result_json"],
        }
        for hash_field, value in section_values.items():
            if _sha256_text(value) != row[hash_field]:
                raise RuntimeError(
                    f"Review record {row['record_id']} failed {hash_field} verification"
                )
        expected_content = _sha256_text(_canonical_json({
            "variant_key": row["variant_key"],
            "module1_sha256": row["module1_sha256"],
            "manual_evidence_sha256": row["manual_evidence_sha256"],
            "amended_result_sha256": row["amended_result_sha256"],
            "policy_id": row["policy_id"],
            "policy_version": row["policy_version"],
            "app_version": row["app_version"],
        }))
        if expected_content != row["content_sha256"]:
            raise RuntimeError(
                f"Review record {row['record_id']} failed content verification"
            )
        record_fields = {
            "record_id": row["record_id"],
            "variant_key": row["variant_key"],
            "version": row["version"],
            "parent_record_id": row["parent_record_id"],
            "status": row["status"],
            "created_at": row["created_at"],
            "created_by_account": row["created_by_account"],
            "reviewer_name": row["reviewer_name"],
            "reviewer_role": row["reviewer_role"],
            "content_sha256": row["content_sha256"],
        }
        if row["status"] == "draft":
            record_fields["review_date"] = row["review_date"]
        else:
            record_fields["approval_comment"] = row["approval_comment"]
        if _sha256_text(_canonical_json(record_fields)) != row["record_sha256"]:
            raise RuntimeError(
                f"Review record {row['record_id']} failed record verification"
            )

    @classmethod
    def _public_record(cls, row: sqlite3.Row) -> dict[str, Any]:
        cls._verify_row(row)
        return {
            "record_id": row["record_id"],
            "schema_version": row["schema_version"],
            "variant_key": row["variant_key"],
            "version": row["version"],
            "parent_record_id": row["parent_record_id"],
            "status": row["status"],
            "created_at": row["created_at"],
            "created_by_account": row["created_by_account"],
            "reviewer_name": row["reviewer_name"],
            "reviewer_role": row["reviewer_role"],
            "review_date": row["review_date"],
            "approval_comment": row["approval_comment"],
            "approved_at": row["approved_at"],
            "app_version": row["app_version"],
            "policy": {
                "id": row["policy_id"],
                "name": row["policy_name"],
                "version": row["policy_version"],
                "source_url": row["policy_source_url"],
            },
            "module1_result": json.loads(row["module1_json"]),
            "manual_evidence": json.loads(row["manual_evidence_json"]),
            "amended_result": json.loads(row["amended_result_json"]),
            "hashes": {
                "module1_sha256": row["module1_sha256"],
                "manual_evidence_sha256": row["manual_evidence_sha256"],
                "amended_result_sha256": row["amended_result_sha256"],
                "content_sha256": row["content_sha256"],
                "record_sha256": row["record_sha256"],
            },
            "integrity_status": "verified",
        }

    def create_draft(
        self,
        *,
        variant_context: Mapping[str, Any],
        module1_result: Mapping[str, Any],
        manual_evidence: Mapping[str, Any],
        amended_result: Mapping[str, Any],
        reviewer_name: str,
        reviewer_role: str,
        review_date: str,
        authenticated_account: str,
    ) -> dict[str, Any]:
        gene = str(variant_context["gene"]).strip().upper()
        c_notation = str(variant_context["c_notation"]).strip()
        variant_key = f"{gene}:{c_notation}"
        policy = get_gene_policy(gene)["policy"]
        module1_json = _canonical_json(module1_result)
        manual_json = _canonical_json(manual_evidence)
        amended_json = _canonical_json(amended_result)
        module1_sha = _sha256_text(module1_json)
        manual_sha = _sha256_text(manual_json)
        amended_sha = _sha256_text(amended_json)
        content_sha = _sha256_text(_canonical_json({
            "variant_key": variant_key,
            "module1_sha256": module1_sha,
            "manual_evidence_sha256": manual_sha,
            "amended_result_sha256": amended_sha,
            "policy_id": policy["runtime_policy_id"],
            "policy_version": policy["version"],
            "app_version": ARIANE_VERSION,
        }))
        record_id = str(uuid.uuid4())
        created_at = _utc_now()

        with self._connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            latest = connection.execute(
                """
                SELECT record_id, version FROM review_records
                WHERE variant_key = ? ORDER BY version DESC LIMIT 1
                """,
                (variant_key,),
            ).fetchone()
            version = (int(latest["version"]) + 1) if latest else 1
            parent_record_id = latest["record_id"] if latest else None
            record_sha = _sha256_text(_canonical_json({
                "record_id": record_id,
                "variant_key": variant_key,
                "version": version,
                "parent_record_id": parent_record_id,
                "status": "draft",
                "created_at": created_at,
                "created_by_account": authenticated_account,
                "reviewer_name": reviewer_name,
                "reviewer_role": reviewer_role,
                "review_date": review_date,
                "content_sha256": content_sha,
            }))
            connection.execute(
                """
                INSERT INTO review_records (
                    record_id, schema_version, variant_key, version,
                    parent_record_id, status, created_at, created_by_account,
                    reviewer_name, reviewer_role, review_date, approval_comment,
                    approved_at, app_version, policy_id, policy_name,
                    policy_version, policy_source_url, module1_json,
                    manual_evidence_json, amended_result_json, module1_sha256,
                    manual_evidence_sha256, amended_result_sha256,
                    content_sha256, record_sha256
                ) VALUES (?, ?, ?, ?, ?, 'draft', ?, ?, ?, ?, ?, '', NULL,
                          ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record_id, SCHEMA_VERSION, variant_key, version,
                    parent_record_id, created_at, authenticated_account,
                    reviewer_name, reviewer_role, review_date,
                    ARIANE_VERSION, policy["runtime_policy_id"], policy["name"],
                    policy["version"], policy["source_url"], module1_json,
                    manual_json, amended_json, module1_sha, manual_sha,
                    amended_sha, content_sha, record_sha,
                ),
            )
            connection.commit()
        return self.get(record_id)

    def approve(
        self,
        record_id: str,
        *,
        approver_name: str,
        approver_role: str,
        approval_comment: str,
        authenticated_account: str,
    ) -> dict[str, Any]:
        approved_at = _utc_now()
        approved_id = str(uuid.uuid4())
        with self._connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            source = connection.execute(
                "SELECT * FROM review_records WHERE record_id = ?",
                (record_id,),
            ).fetchone()
            if source is None:
                raise KeyError("Review record not found")
            if source["status"] != "draft":
                raise ValueError("Only a draft review record can be approved")
            latest_version = connection.execute(
                "SELECT MAX(version) FROM review_records WHERE variant_key = ?",
                (source["variant_key"],),
            ).fetchone()[0]
            version = int(latest_version) + 1
            record_sha = _sha256_text(_canonical_json({
                "record_id": approved_id,
                "variant_key": source["variant_key"],
                "version": version,
                "parent_record_id": record_id,
                "status": "approved",
                "created_at": approved_at,
                "created_by_account": authenticated_account,
                "reviewer_name": approver_name,
                "reviewer_role": approver_role,
                "approval_comment": approval_comment,
                "content_sha256": source["content_sha256"],
            }))
            connection.execute(
                """
                INSERT INTO review_records (
                    record_id, schema_version, variant_key, version,
                    parent_record_id, status, created_at, created_by_account,
                    reviewer_name, reviewer_role, review_date, approval_comment,
                    approved_at, app_version, policy_id, policy_name,
                    policy_version, policy_source_url, module1_json,
                    manual_evidence_json, amended_result_json, module1_sha256,
                    manual_evidence_sha256, amended_result_sha256,
                    content_sha256, record_sha256
                ) VALUES (?, ?, ?, ?, ?, 'approved', ?, ?, ?, ?, ?, ?, ?,
                          ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    approved_id, source["schema_version"], source["variant_key"],
                    version, record_id, approved_at, authenticated_account,
                    approver_name, approver_role, source["review_date"],
                    approval_comment, approved_at, source["app_version"],
                    source["policy_id"], source["policy_name"],
                    source["policy_version"], source["policy_source_url"],
                    source["module1_json"], source["manual_evidence_json"],
                    source["amended_result_json"], source["module1_sha256"],
                    source["manual_evidence_sha256"],
                    source["amended_result_sha256"], source["content_sha256"],
                    record_sha,
                ),
            )
            connection.commit()
        return self.get(approved_id)

    def get(self, record_id: str) -> dict[str, Any]:
        with self._connection() as connection:
            row = connection.execute(
                "SELECT * FROM review_records WHERE record_id = ?",
                (record_id,),
            ).fetchone()
        if row is None:
            raise KeyError("Review record not found")
        return self._public_record(row)

    def list(self, *, limit: int = 100) -> list[dict[str, Any]]:
        with self._connection() as connection:
            rows = connection.execute(
                """
                SELECT * FROM review_records
                ORDER BY created_at DESC, version DESC LIMIT ?
                """,
                (max(1, min(limit, 500)),),
            ).fetchall()
        return [self._public_record(row) for row in rows]
