#!/usr/bin/env python3
"""Create, list and disable ARIANE public API keys."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import secrets
import sys
import tempfile

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from backend.api.auth import API_KEY_SCHEMA_VERSION


KEY_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")


def _empty_registry() -> dict:
    return {"schema_version": API_KEY_SCHEMA_VERSION, "keys": []}


def _read(path: Path) -> dict:
    if not path.exists():
        return _empty_registry()
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(payload, dict)
        or payload.get("schema_version") != API_KEY_SCHEMA_VERSION
        or not isinstance(payload.get("keys"), list)
    ):
        raise ValueError("Unsupported or malformed API-key registry")
    return payload


def _write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    previous = path.stat() if path.exists() else None
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    )
    temporary = Path(handle.name)
    try:
        with handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, previous.st_mode & 0o777 if previous else 0o600)
        if previous and hasattr(os, "chown"):
            try:
                os.chown(temporary, previous.st_uid, previous.st_gid)
            except PermissionError:
                pass
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _find(payload: dict, key_id: str) -> dict | None:
    return next((item for item in payload["keys"] if item.get("id") == key_id), None)


def create_key(path: Path, key_id: str, description: str) -> str:
    if not KEY_ID_RE.fullmatch(key_id):
        raise ValueError(
            "Key id must contain 1 to 64 letters, digits, dots, underscores or hyphens"
        )
    payload = _read(path)
    if _find(payload, key_id) is not None:
        raise ValueError(f"API-key id already exists: {key_id}")
    secret = f"ariane_v1_{secrets.token_urlsafe(32)}"
    payload["keys"].append({
        "id": key_id,
        "secret_sha256": hashlib.sha256(secret.encode("utf-8")).hexdigest(),
        "enabled": True,
        "description": description.strip(),
        "created_at": datetime.now(timezone.utc).isoformat(),
    })
    _write(path, payload)
    return secret


def disable_key(path: Path, key_id: str) -> None:
    payload = _read(path)
    record = _find(payload, key_id)
    if record is None:
        raise ValueError(f"Unknown API-key id: {key_id}")
    record["enabled"] = False
    _write(path, payload)


def list_keys(path: Path) -> list[dict]:
    return [
        {
            "id": str(item.get("id") or ""),
            "enabled": bool(item.get("enabled")),
            "description": str(item.get("description") or ""),
            "created_at": str(item.get("created_at") or ""),
        }
        for item in _read(path)["keys"]
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_path = Path(
        os.getenv(
            "ARIANE_API_KEYS_FILE",
            str(PROJECT_ROOT / ".runtime-data" / "api_keys.json"),
        )
    )
    parser.add_argument("--file", type=Path, default=default_path)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("create", help="Create a new key")
    create.add_argument("--id", required=True, dest="key_id")
    create.add_argument("--description", default="")
    disable = commands.add_parser("disable", help="Disable an existing key")
    disable.add_argument("--id", required=True, dest="key_id")
    commands.add_parser("list", help="List key IDs and states")
    args = parser.parse_args()

    if args.command == "create":
        secret = create_key(args.file, args.key_id, args.description)
        print(f"Created API key: {args.key_id}")
        print("Copy this value now. It is not stored in recoverable form:")
        print(secret)
    elif args.command == "disable":
        disable_key(args.file, args.key_id)
        print(f"Disabled API key: {args.key_id}")
    else:
        for item in list_keys(args.file):
            state = "enabled" if item["enabled"] else "disabled"
            print(
                f"{item['id']}\t{state}\t{item['created_at']}\t"
                f"{item['description']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
