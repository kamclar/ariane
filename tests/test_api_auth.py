import hashlib
import json

import pytest

from backend.api_auth import (
    API_KEY_FILE_ENVIRONMENT,
    ApiKeyAuthenticator,
    ApiKeyConfigurationError,
    load_api_key_records,
)
from scripts.manage_api_keys import create_key, disable_key, list_keys


def test_authenticator_accepts_enabled_key_and_returns_nonsecret_id(
    tmp_path,
    monkeypatch,
):
    secret = "ariane_v1_test-secret"
    path = tmp_path / "api-keys.json"
    path.write_text(json.dumps({
        "schema_version": 1,
        "keys": [{
            "id": "integration-a",
            "secret_sha256": hashlib.sha256(secret.encode()).hexdigest(),
            "enabled": True,
        }],
    }), encoding="utf-8")
    monkeypatch.setenv(API_KEY_FILE_ENVIRONMENT, str(path))

    authenticator = ApiKeyAuthenticator()
    assert authenticator.authenticate(secret) == "integration-a"
    assert authenticator.authenticate("wrong") is None


def test_disabled_key_is_rejected(tmp_path, monkeypatch):
    secret = "ariane_v1_disabled"
    path = tmp_path / "api-keys.json"
    path.write_text(json.dumps({
        "schema_version": 1,
        "keys": [{
            "id": "disabled-client",
            "secret_sha256": hashlib.sha256(secret.encode()).hexdigest(),
            "enabled": False,
        }],
    }), encoding="utf-8")
    monkeypatch.setenv(API_KEY_FILE_ENVIRONMENT, str(path))

    assert ApiKeyAuthenticator().authenticate(secret) is None


def test_missing_or_malformed_registry_fails_closed(tmp_path, monkeypatch):
    missing = tmp_path / "missing.json"
    monkeypatch.setenv(API_KEY_FILE_ENVIRONMENT, str(missing))
    with pytest.raises(ApiKeyConfigurationError):
        ApiKeyAuthenticator().authenticate("some-key")

    malformed = tmp_path / "malformed.json"
    malformed.write_text('{"schema_version":1,"keys":"wrong"}', encoding="utf-8")
    monkeypatch.setenv(API_KEY_FILE_ENVIRONMENT, str(malformed))
    with pytest.raises(ApiKeyConfigurationError):
        ApiKeyAuthenticator().authenticate("some-key")


def test_management_creates_checksum_only_key_and_disables_it(tmp_path):
    path = tmp_path / "api-keys.json"
    secret = create_key(path, "batch-client", "Batch testing")
    raw = path.read_text(encoding="utf-8")

    assert secret.startswith("ariane_v1_")
    assert secret not in raw
    assert load_api_key_records(path)[0].enabled is True
    assert list_keys(path)[0]["description"] == "Batch testing"

    disable_key(path, "batch-client")
    assert load_api_key_records(path)[0].enabled is False
