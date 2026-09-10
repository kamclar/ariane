#!/usr/bin/env python3
"""Create a transactionally consistent SQLite database backup."""

from pathlib import Path
import sqlite3
import sys


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: backup-review-records.py <source.sqlite3> <destination.sqlite3>", file=sys.stderr)
        return 2
    source = Path(sys.argv[1]).resolve()
    destination = Path(sys.argv[2]).resolve()
    if not source.is_file():
        print(f"Database not found: {source}", file=sys.stderr)
        return 3
    if destination.exists():
        print(f"Destination already exists: {destination}", file=sys.stderr)
        return 4
    destination.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(source) as source_db, sqlite3.connect(destination) as target_db:
        source_db.backup(target_db)
        result = target_db.execute("PRAGMA integrity_check").fetchone()[0]
        if result != "ok":
            raise RuntimeError(f"SQLite integrity check failed: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
