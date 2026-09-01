import hashlib
import json
from pathlib import Path

import pytest

from scripts.check_enigma_clinical_lr_update import check_update


TRACK_DB = b"track BRCAmla\nshortLabel ENIGMA\ndataVersion 2026-08-18\n"


class FakeResponse:
    def __init__(self, payload: bytes, headers: dict[str, str] | None = None):
        self.payload = payload
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def read(self, _limit: int = -1) -> bytes:
        return self.payload


def _manifest(path: Path, payload: bytes, *, automatic_activation: bool = False) -> Path:
    content = {
        "schema_version": 3,
        "clinical_lr_data_release": {"release_date": "2026-08-18"},
        "datasets": {
            "active": {
                "url": "https://example.test/BRCAmfa.bb",
                "derived_from_bigbed_sha256": hashlib.sha256(payload).hexdigest(),
            }
        },
        "update_policy": {
            "automatic_release_activation": automatic_activation,
            "required_before_activation": ["expert review"],
        },
    }
    path.write_text(json.dumps(content), encoding="utf-8")
    return path


def _opener_for(bigbed_payload: bytes):
    def opener(request, timeout):
        assert timeout > 0
        if request.full_url.endswith("trackDb.txt"):
            return FakeResponse(TRACK_DB)
        return FakeResponse(
            bigbed_payload,
            {
                "Content-Length": str(len(bigbed_payload)),
                "Last-Modified": "Tue, 18 Aug 2026 15:07:00 GMT",
                "ETag": '"test"',
            },
        )

    return opener


def test_current_release_is_detected_without_writing(tmp_path):
    payload = b"pinned bigbed"
    result = check_update(
        manifest_path=_manifest(tmp_path / "manifest.json", payload),
        opener=_opener_for(payload),
    )

    assert result["status"] == "current"
    assert result["review_required"] is False
    assert result["automatic_activation"] is False
    assert result["remote_data_version"] == "2026-08-18"
    assert list(tmp_path.iterdir()) == [tmp_path / "manifest.json"]


def test_changed_release_is_saved_only_as_inactive_candidate(tmp_path):
    pinned = b"pinned bigbed"
    remote = b"new candidate bigbed"
    candidate_dir = tmp_path / "candidates"
    result = check_update(
        manifest_path=_manifest(tmp_path / "manifest.json", pinned),
        candidate_dir=candidate_dir,
        opener=_opener_for(remote),
    )

    assert result["status"] == "update_available"
    assert result["review_required"] is True
    assert result["automatic_activation"] is False
    candidate = Path(result["candidate_file"])
    audit = Path(result["candidate_audit_file"])
    assert candidate.read_bytes() == remote
    assert audit.exists()
    assert json.loads(audit.read_text(encoding="utf-8"))["status"] == "update_available"


def test_manifest_cannot_enable_automatic_release_activation(tmp_path):
    payload = b"pinned bigbed"
    manifest = _manifest(
        tmp_path / "manifest.json", payload, automatic_activation=True
    )

    with pytest.raises(RuntimeError, match="prohibit automatic activation"):
        check_update(manifest_path=manifest, opener=_opener_for(payload))
