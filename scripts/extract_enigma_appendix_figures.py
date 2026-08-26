"""Extract the embedded ENIGMA Appendix v1.2 figure assets reproducibly.

The source DOCX is checksum-pinned in ``backend/data/enigma_rule_catalog.json``.
This script does not transform the source images. It gives the embedded media
stable public names so the rule explorer can show the official originals next
to ARIANE's accessible SVG redraws.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs" / "enigma" / "v1.2" / "source" / "Appendix_V1.2.docx"
DESTINATION = ROOT / "frontend" / "static" / "enigma"
EXPECTED_SOURCE_SHA256 = "9a9cc9210178e8d947597089ba8edd5ff0633ff7e4fea268b8f202d46e6c8ebb"

ASSETS = {
    "word/media/image1.jpg": "appendix-figure-1-brca1-gene-map.jpg",
    "word/media/image2.png": "appendix-figure-2-brca2-gene-map.png",
    "word/media/image3.jpg": "appendix-figure-3-brca1-pvs1.jpg",
    "word/media/image4.jpg": "appendix-figure-4-brca1-splice-pvs1.jpg",
    "word/media/image5.png": "appendix-figure-4-brca1-splice-pvs1-continued.png",
    "word/media/image6.jpg": "appendix-figure-5-brca2-pvs1.jpg",
    "word/media/image7.jpg": "appendix-figure-6-brca2-splice-pvs1.jpg",
    "word/media/image8.png": "appendix-figure-6-brca2-splice-pvs1-continued.png",
    "word/media/image9.jpg": "appendix-figure-9-evidence-interactions.jpg",
}


def main() -> None:
    if not SOURCE.is_file():
        raise SystemExit(f"Missing pinned ENIGMA Appendix source: {SOURCE}")
    actual_sha256 = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if actual_sha256 != EXPECTED_SOURCE_SHA256:
        raise SystemExit(
            "ENIGMA Appendix checksum mismatch: "
            f"expected {EXPECTED_SOURCE_SHA256}, got {actual_sha256}"
        )
    DESTINATION.mkdir(parents=True, exist_ok=True)
    with ZipFile(SOURCE) as archive:
        available = set(archive.namelist())
        missing = set(ASSETS) - available
        if missing:
            raise SystemExit(f"Appendix media layout changed; missing: {sorted(missing)}")
        for member, filename in ASSETS.items():
            target = DESTINATION / filename
            target.write_bytes(archive.read(member))
            print(f"{member} -> {target.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
