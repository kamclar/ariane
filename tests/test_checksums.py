from __future__ import annotations

import hashlib

from backend.integrity.checksums import (
    canonical_json_sha256,
    file_sha256,
    normalized_text_sha256,
)


def test_file_sha256_hashes_exact_bytes(tmp_path) -> None:
    path = tmp_path / "payload.bin"
    payload = b"\x00ARIANE\r\n\xff"
    path.write_bytes(payload)

    assert file_sha256(path) == hashlib.sha256(payload).hexdigest()


def test_normalized_text_sha256_is_independent_of_checkout_line_endings(
    tmp_path,
) -> None:
    lf_path = tmp_path / "lf.txt"
    crlf_path = tmp_path / "crlf.txt"
    lf_path.write_bytes(b"alpha\nbeta\n")
    crlf_path.write_bytes(b"alpha\r\nbeta\r\n")

    assert normalized_text_sha256(lf_path) == normalized_text_sha256(crlf_path)


def test_canonical_json_sha256_is_order_independent_and_unicode_preserving() -> None:
    first = {"gene": "BRCA1", "label": "Příliš", "values": [2, 1]}
    second = {"values": [2, 1], "label": "Příliš", "gene": "BRCA1"}
    expected = hashlib.sha256(
        '{"gene":"BRCA1","label":"Příliš","values":[2,1]}'.encode("utf-8")
    ).hexdigest()

    assert canonical_json_sha256(first) == canonical_json_sha256(second) == expected
