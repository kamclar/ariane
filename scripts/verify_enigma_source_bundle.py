"""Verify the pinned ClinGen ENIGMA v1.2 source bundle."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]
BUNDLE_ROOT = ROOT / "docs" / "enigma" / "v1.2"
MANIFEST_PATH = BUNDLE_ROOT / "manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_bundle() -> list[str]:
    errors: list[str] = []
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if manifest.get("criteria_specification_id") != "GN092":
        errors.append("unexpected criteria specification id")
    if manifest.get("criteria_specification_version") != "1.2.0":
        errors.append("unexpected criteria specification version")

    for item in manifest.get("files", []):
        relative_path = Path(item["path"])
        path = BUNDLE_ROOT / relative_path
        if not path.is_file():
            errors.append(f"missing file: {relative_path.as_posix()}")
            continue
        if path.stat().st_size != item.get("size_bytes"):
            errors.append(f"size mismatch: {relative_path.as_posix()}")
        if sha256(path) != str(item.get("sha256", "")).lower():
            errors.append(f"checksum mismatch: {relative_path.as_posix()}")
        try:
            with zipfile.ZipFile(path) as archive:
                damaged_member = archive.testzip()
            if damaged_member:
                errors.append(
                    f"damaged OOXML member in {relative_path.as_posix()}: "
                    f"{damaged_member}"
                )
        except zipfile.BadZipFile:
            errors.append(f"not a valid OOXML ZIP file: {relative_path.as_posix()}")
    return errors


def main() -> None:
    errors = verify_bundle()
    if errors:
        raise SystemExit("ENIGMA source bundle validation failed:\n- " + "\n- ".join(errors))
    print("ENIGMA GN092 v1.2.0 source bundle verified")


if __name__ == "__main__":
    main()
