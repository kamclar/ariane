"""API-key authentication for the versioned public API."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import hmac
import json
import os
from pathlib import Path
import re
from typing import Annotated

from fastapi import HTTPException, Request, Security
from fastapi.security import APIKeyHeader

from backend.infrastructure.runtime_data import runtime_data_path


API_KEY_HEADER_NAME = "X-ARIANE-API-Key"
API_KEY_FILE_ENVIRONMENT = "ARIANE_API_KEYS_FILE"
API_KEY_SCHEMA_VERSION = 1
_KEY_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_API_KEY_HEADER = APIKeyHeader(
    name=API_KEY_HEADER_NAME,
    scheme_name="ArianeApiKey",
    description="ARIANE public API key",
    auto_error=False,
)


class ApiKeyConfigurationError(RuntimeError):
    """The API-key registry is missing or invalid."""


@dataclass(frozen=True)
class ApiKeyRecord:
    id: str
    secret_sha256: str
    enabled: bool
    description: str = ""
    created_at: str = ""


def configured_api_key_path() -> Path:
    configured = os.getenv(API_KEY_FILE_ENVIRONMENT, "").strip()
    if configured:
        return Path(configured)
    return runtime_data_path("api_keys.json")


def load_api_key_records(path: Path) -> tuple[ApiKeyRecord, ...]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ApiKeyConfigurationError(
            f"API-key registry is missing: {path}"
        ) from exc
    except (OSError, json.JSONDecodeError) as exc:
        raise ApiKeyConfigurationError(
            f"API-key registry cannot be read: {path}"
        ) from exc
    if not isinstance(payload, dict):
        raise ApiKeyConfigurationError("API-key registry root must be an object")
    if payload.get("schema_version") != API_KEY_SCHEMA_VERSION:
        raise ApiKeyConfigurationError("Unsupported API-key registry schema")
    raw_keys = payload.get("keys")
    if not isinstance(raw_keys, list):
        raise ApiKeyConfigurationError("API-key registry must contain a keys list")

    records: list[ApiKeyRecord] = []
    seen_ids: set[str] = set()
    for index, raw in enumerate(raw_keys):
        if not isinstance(raw, dict):
            raise ApiKeyConfigurationError(
                f"API-key record {index} must be an object"
            )
        key_id = str(raw.get("id") or "")
        digest = str(raw.get("secret_sha256") or "").lower()
        enabled = raw.get("enabled")
        if not _KEY_ID_RE.fullmatch(key_id):
            raise ApiKeyConfigurationError(
                f"API-key record {index} has an invalid id"
            )
        if key_id in seen_ids:
            raise ApiKeyConfigurationError(f"Duplicate API-key id: {key_id}")
        if not _SHA256_RE.fullmatch(digest):
            raise ApiKeyConfigurationError(
                f"API-key record {key_id} has an invalid SHA-256 digest"
            )
        if not isinstance(enabled, bool):
            raise ApiKeyConfigurationError(
                f"API-key record {key_id} has no boolean enabled state"
            )
        seen_ids.add(key_id)
        records.append(ApiKeyRecord(
            id=key_id,
            secret_sha256=digest,
            enabled=enabled,
            description=str(raw.get("description") or ""),
            created_at=str(raw.get("created_at") or ""),
        ))
    return tuple(records)


class ApiKeyAuthenticator:
    """Validate high-entropy API keys against a checksum-only registry."""

    def authenticate(self, presented_key: str) -> str | None:
        value = (presented_key or "").strip()
        if not value or len(value) > 256:
            return None
        presented_digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
        matched_id: str | None = None
        for record in load_api_key_records(configured_api_key_path()):
            digest_matches = hmac.compare_digest(
                presented_digest,
                record.secret_sha256,
            )
            if digest_matches and record.enabled:
                matched_id = record.id
        return matched_id


PUBLIC_API_KEY_AUTHENTICATOR = ApiKeyAuthenticator()


async def require_public_api_key(
    request: Request,
    presented_key: Annotated[str | None, Security(_API_KEY_HEADER)] = None,
) -> str:
    if not presented_key:
        raise HTTPException(
            status_code=401,
            detail={
                "code": "api_key_required",
                "message": f"Provide an API key in the {API_KEY_HEADER_NAME} header",
                "retryable": False,
            },
            headers={"WWW-Authenticate": "ApiKey"},
        )
    try:
        key_id = PUBLIC_API_KEY_AUTHENTICATOR.authenticate(presented_key)
    except ApiKeyConfigurationError as exc:
        raise HTTPException(
            status_code=503,
            detail={
                "code": "api_authentication_unavailable",
                "message": "Public API authentication is not configured",
                "retryable": False,
            },
        ) from exc
    if key_id is None:
        raise HTTPException(
            status_code=401,
            detail={
                "code": "invalid_api_key",
                "message": "The supplied API key is invalid or disabled",
                "retryable": False,
            },
            headers={"WWW-Authenticate": "ApiKey"},
        )
    request.state.api_key_id = key_id
    return key_id
