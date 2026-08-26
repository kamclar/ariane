import re
from pathlib import Path

from backend.version import ARIANE_VERSION


MAIN_SOURCE = (Path(__file__).resolve().parents[1] / "backend" / "main.py").read_text(encoding="utf-8")


def test_application_version_uses_semantic_version_format():
    assert re.fullmatch(r"\d+\.\d+\.\d+", ARIANE_VERSION)


def test_fastapi_and_public_endpoints_use_single_version_constant():
    assert "from backend.version import ARIANE_VERSION" in MAIN_SOURCE
    assert "version=ARIANE_VERSION" in MAIN_SOURCE
    assert MAIN_SOURCE.count('"version": ARIANE_VERSION') == 2
    assert '"version": "1.8.' not in MAIN_SOURCE
