"""Checkpoint durability and provenance tests, without inference/model imports."""

import json
from pathlib import Path
import sqlite3
import subprocess
import sys

import pytest

from scripts.spliceai_distance_store import ResultStore


@pytest.fixture
def identity():
    return {
        "manifest_sha256": "fixed-manifest-digest",
        "case_sha256": "fixed-cohort-digest",
        "config": {"distances": [4999, 5000, 10000], "tolerance": 1e-6},
    }


def test_completed_results_and_failures_survive_resume(tmp_path, identity):
    path = tmp_path / "results.sqlite"
    success = {"id": "case-a", "status": "ok", "measurements": {"delta": 0.123456789}}
    failure = {"id": "case-b", "status": "error", "error": "fixture mismatch"}
    with ResultStore(path, identity) as store:
        assert store.counts() == {"total": 0, "ok": 0, "error": 0}
        assert store.get("unfinished") is None
        store.put("case-a", success)
        store.put("case-b", failure)
        assert store.counts() == {"total": 2, "ok": 1, "error": 1}

    with ResultStore(path, identity, resume=True) as store:
        assert store.get("case-a") == success
        assert store.get("case-b") == failure
        assert list(store.results()) == [("case-a", success), ("case-b", failure)]
        assert store.get("unfinished") is None


def test_new_store_refuses_existing_database_without_modification(tmp_path, identity):
    path = tmp_path / "results.sqlite"
    with ResultStore(path, identity) as store:
        store.put("a", {"id": "a", "status": "ok"})
    before = path.read_bytes()
    with pytest.raises(FileExistsError, match="already exists"):
        ResultStore(path, identity)
    assert path.read_bytes() == before


def test_resume_requires_existing_database(tmp_path, identity):
    path = tmp_path / "missing.sqlite"
    with pytest.raises(FileNotFoundError, match="does not exist"):
        ResultStore(path, identity, resume=True)
    assert not path.exists()


@pytest.mark.parametrize("changed_key", ["manifest_sha256", "case_sha256", "config"])
def test_identity_mismatch_is_rejected_without_modification(tmp_path, identity, changed_key):
    path = tmp_path / "results.sqlite"
    with ResultStore(path, identity) as store:
        store.put("a", {"id": "a", "status": "error", "error": "original"})
    before = path.read_bytes()
    changed = dict(identity, **{changed_key: "different"})
    with pytest.raises(ValueError, match="identity mismatch"):
        ResultStore(path, changed, resume=True)
    assert path.read_bytes() == before
    with ResultStore(path, identity, resume=True) as store:
        assert store.get("a")["error"] == "original"


def test_identity_is_canonical_and_copied(tmp_path, identity):
    path = tmp_path / "results.sqlite"
    reordered = dict(reversed(list(identity.items())))
    with ResultStore(path, identity):
        identity["config"]["tolerance"] = 0.5
    original = json.loads(json.dumps(reordered))
    original["config"]["tolerance"] = 1e-6
    with ResultStore(path, original, resume=True) as store:
        assert store.counts()["total"] == 0


def test_success_is_immutable_but_identical_write_is_idempotent(tmp_path, identity):
    with ResultStore(tmp_path / "results.sqlite", identity) as store:
        result = {"id": "a", "status": "ok", "value": 0.1234}
        store.put("a", result)
        store.put("a", dict(reversed(list(result.items()))))
        with pytest.raises(ValueError, match="overwrite successful"):
            store.put("a", dict(result, value=1))
        with pytest.raises(ValueError, match="overwrite successful"):
            store.put("a", {"id": "a", "status": "error", "error": "later"})
        assert store.get("a") == result
        assert store.counts() == {"total": 1, "ok": 1, "error": 0}


def test_explicit_retry_can_replace_a_failure(tmp_path, identity):
    with ResultStore(tmp_path / "results.sqlite", identity) as store:
        store.put("a", {"id": "a", "status": "error", "error": "first attempt"})
        store.put("a", {"id": "a", "status": "error", "error": "second attempt"})
        assert store.get("a")["error"] == "second attempt"
        store.put("a", {"id": "a", "status": "ok", "value": 0})
        assert store.counts() == {"total": 1, "ok": 1, "error": 0}


@pytest.mark.parametrize("result", [
    {},
    {"id": "other", "status": "ok"},
    {"id": "a"},
    {"id": "a", "status": "running"},
    {"id": "a", "status": "pending"},
    {"id": "a", "status": "ok", "value": float("nan")},
    {"id": "a", "status": "ok", "value": float("inf")},
    {"id": "a", "status": "ok", "nested": {1: "lossy key"}},
])
def test_invalid_or_incomplete_results_do_not_clobber_failures(tmp_path, identity, result):
    with ResultStore(tmp_path / "results.sqlite", identity) as store:
        failure = {"id": "a", "status": "error", "error": "retained"}
        store.put("a", failure)
        with pytest.raises(ValueError):
            store.put("a", result)
        assert store.get("a") == failure
        assert store.counts() == {"total": 1, "ok": 0, "error": 1}


@pytest.mark.parametrize("case_id", ["", " ", None, 123])
def test_invalid_case_id_is_rejected(tmp_path, identity, case_id):
    with ResultStore(tmp_path / "results.sqlite", identity) as store:
        with pytest.raises(ValueError, match="case_id"):
            store.put(case_id, {"id": case_id, "status": "ok"})
        assert store.counts()["total"] == 0


@pytest.mark.parametrize("invalid_identity", [[], {1: "not a string"}, {"value": float("nan")}])
def test_invalid_identity_creates_no_file(tmp_path, invalid_identity):
    path = tmp_path / "results.sqlite"
    with pytest.raises(ValueError):
        ResultStore(path, invalid_identity)
    assert not path.exists()


def test_unrelated_database_is_not_initialized_on_resume(tmp_path, identity):
    path = tmp_path / "unrelated.sqlite"
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE unrelated (value TEXT)")
        connection.execute("INSERT INTO unrelated VALUES ('retained')")
    before = path.read_bytes()
    with pytest.raises(ValueError, match="Not a valid"):
        ResultStore(path, identity, resume=True)
    assert path.read_bytes() == before


@pytest.mark.parametrize("column, value, message", [
    ("schema_version", 999, "schema version"),
    ("identity_sha256", "corrupt-digest", "identity mismatch"),
])
def test_unsupported_or_corrupt_metadata_is_not_modified(tmp_path, identity, column, value, message):
    path = tmp_path / "results.sqlite"
    with ResultStore(path, identity):
        pass
    with sqlite3.connect(path) as connection:
        connection.execute(f"UPDATE metadata SET {column} = ?", (value,))
    before = path.read_bytes()
    with pytest.raises(ValueError, match=message):
        ResultStore(path, identity, resume=True)
    assert path.read_bytes() == before


def test_committed_put_survives_abrupt_process_exit(tmp_path, identity):
    path = tmp_path / "results.sqlite"
    script = """
import json, os, sys
from pathlib import Path
from scripts.spliceai_distance_store import ResultStore
store = ResultStore(Path(sys.argv[1]), json.loads(sys.argv[2]))
store.put('done', {'id': 'done', 'status': 'ok'})
os._exit(17)
"""
    process = subprocess.run(
        [sys.executable, "-c", script, str(path), json.dumps(identity)],
        check=False, cwd=Path(__file__).resolve().parents[1],
    )
    assert process.returncode == 17
    with ResultStore(path, identity, resume=True) as store:
        assert store.get("done") == {"id": "done", "status": "ok"}


def test_commit_visible_without_closing_and_crashed_transaction_is_rolled_back(tmp_path, identity):
    path = tmp_path / "results.sqlite"
    with ResultStore(path, identity) as store:
        store.put("committed", {"id": "committed", "status": "ok"})
        # A separate process sees the committed case, then dies after an
        # uncommitted update. No Python cleanup or connection.close is run.
        script = """
import os, sqlite3, sys
db = sqlite3.connect(sys.argv[1])
assert db.execute('SELECT COUNT(*) FROM cases').fetchone()[0] == 1
db.execute('BEGIN IMMEDIATE')
db.execute('DELETE FROM cases')
db.execute("INSERT INTO cases VALUES ('incomplete', 'ok', '{}')")
os._exit(17)
"""
        process = subprocess.run([sys.executable, "-c", script, str(path)], check=False)
        assert process.returncode == 17
        assert store.get("committed") == {"id": "committed", "status": "ok"}
        assert store.get("incomplete") is None
    with ResultStore(path, identity, resume=True) as store:
        assert store.counts() == {"total": 1, "ok": 1, "error": 0}


def test_sql_failure_rolls_back_case_transaction(tmp_path, identity):
    path = tmp_path / "results.sqlite"
    with ResultStore(path, identity) as store:
        failure = {"id": "a", "status": "error", "error": "retained"}
        store.put("a", failure)
        with sqlite3.connect(path) as connection:
            connection.execute(
                "CREATE TRIGGER reject_update BEFORE UPDATE ON cases "
                "BEGIN SELECT RAISE(ABORT, 'simulated write failure'); END"
            )
        with pytest.raises(sqlite3.IntegrityError, match="simulated write failure"):
            store.put("a", {"id": "a", "status": "ok"})
        assert store.get("a") == failure
        store.put("b", {"id": "b", "status": "ok"})
        assert store.counts() == {"total": 2, "ok": 1, "error": 1}


def test_wal_full_durability_and_compact_json(tmp_path, identity):
    path = tmp_path / "results.sqlite"
    with ResultStore(path, identity) as store:
        store.put("a", {"id": "a", "status": "ok", "list": [1, 2, 3]})
        assert store._db().execute("PRAGMA journal_mode").fetchone()[0] == "wal"
        assert store._db().execute("PRAGMA synchronous").fetchone()[0] == 2
        payload = store._db().execute("SELECT result_json FROM cases").fetchone()[0]
        assert payload == '{"id":"a","list":[1,2,3],"status":"ok"}'


def test_close_is_idempotent_and_use_after_close_is_clear(tmp_path, identity):
    store = ResultStore(tmp_path / "results.sqlite", identity)
    store.close()
    store.close()
    with pytest.raises(RuntimeError, match="closed"):
        store.counts()
