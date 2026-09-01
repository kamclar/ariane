"""Detect a new UCSC ENIGMA PP4/BP5 data release without activating it."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "data/sources/enigma/clinical_lr_sources.manifest.json"
TRACK_DB_URL = "https://hgdownload.soe.ucsc.edu/hubs/enigma/hg38/trackDb.txt"
MAX_DOWNLOAD_BYTES = 20 * 1024 * 1024


def _load_manifest(path: Path = MANIFEST_PATH) -> tuple[dict, dict]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 3:
        raise RuntimeError("Clinical LR source manifest schema is unsupported")
    datasets = manifest.get("datasets") or {}
    if len(datasets) != 1:
        raise RuntimeError("Clinical LR manifest must pin exactly one active source dataset")
    if manifest.get("update_policy", {}).get("automatic_release_activation") is not False:
        raise RuntimeError("Clinical LR update policy must prohibit automatic activation")
    return manifest, next(iter(datasets.values()))


def _download(
    url: str,
    *,
    timeout: float,
    opener: Callable = urllib.request.urlopen,
) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "ARIANE-clinical-LR-update-check/1.0"},
    )
    with opener(request, timeout=timeout) as response:
        content_length = response.headers.get("Content-Length")
        if content_length and int(content_length) > MAX_DOWNLOAD_BYTES:
            raise RuntimeError("Remote clinical LR source exceeds the safety size limit")
        payload = response.read(MAX_DOWNLOAD_BYTES + 1)
        if len(payload) > MAX_DOWNLOAD_BYTES:
            raise RuntimeError("Remote clinical LR source exceeds the safety size limit")
        headers = {
            "content_length": response.headers.get("Content-Length", ""),
            "last_modified": response.headers.get("Last-Modified", ""),
            "etag": response.headers.get("ETag", ""),
        }
    return payload, headers


def _track_data_version(
    *,
    timeout: float,
    opener: Callable = urllib.request.urlopen,
) -> str:
    try:
        payload, _headers = _download(TRACK_DB_URL, timeout=timeout, opener=opener)
    except Exception:
        return "unavailable"
    text = payload.decode("utf-8", errors="replace")
    match = re.search(
        r"(?ms)^track BRCAmla\s+.*?^dataVersion\s+(.+?)\s*$",
        text,
    )
    return match.group(1).strip() if match else "not declared"


def check_update(
    *,
    manifest_path: Path = MANIFEST_PATH,
    timeout: float = 30.0,
    candidate_dir: Path | None = None,
    opener: Callable = urllib.request.urlopen,
) -> dict:
    manifest, dataset = _load_manifest(manifest_path)
    payload, headers = _download(dataset["url"], timeout=timeout, opener=opener)
    remote_sha256 = hashlib.sha256(payload).hexdigest()
    pinned_sha256 = dataset["derived_from_bigbed_sha256"]
    update_available = remote_sha256 != pinned_sha256
    result = {
        "status": "update_available" if update_available else "current",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": dataset["url"],
        "pinned_release_date": manifest["clinical_lr_data_release"]["release_date"],
        "pinned_bigbed_sha256": pinned_sha256,
        "remote_bigbed_sha256": remote_sha256,
        "remote_content_length": len(payload),
        "remote_headers": headers,
        "remote_data_version": _track_data_version(timeout=timeout, opener=opener),
        "automatic_activation": False,
        "review_required": update_available,
        "activation_requirements": manifest["update_policy"]["required_before_activation"],
    }
    if candidate_dir is not None:
        candidate_dir.mkdir(parents=True, exist_ok=True)
        candidate_name = f"BRCAmfa.{remote_sha256[:12]}.candidate.bb"
        candidate_path = candidate_dir / candidate_name
        candidate_path.write_bytes(payload)
        audit_path = candidate_dir / f"{candidate_name}.json"
        audit_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        result["candidate_file"] = str(candidate_path)
        result["candidate_audit_file"] = str(audit_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument(
        "--candidate-dir",
        type=Path,
        help="Optionally retain the downloaded file as an inactive candidate.",
    )
    args = parser.parse_args()
    result = check_update(
        manifest_path=args.manifest,
        timeout=args.timeout,
        candidate_dir=args.candidate_dir,
    )
    print(json.dumps(result, indent=2))
    return 2 if result["status"] == "update_available" else 0


if __name__ == "__main__":
    raise SystemExit(main())
