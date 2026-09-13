import re
from pathlib import Path

from backend.version import ARIANE_VERSION


MAIN_SOURCE = (Path(__file__).resolve().parents[1] / "backend" / "main.py").read_text(encoding="utf-8")
SYSTEM_API_SOURCE = (
    Path(__file__).resolve().parents[1] / "backend" / "api" / "system.py"
).read_text(encoding="utf-8")
PUBLIC_API_SOURCE = (
    Path(__file__).resolve().parents[1] / "backend" / "api" / "public.py"
).read_text(encoding="utf-8")


def test_application_version_uses_semantic_version_format():
    assert re.fullmatch(r"\d+\.\d+\.\d+", ARIANE_VERSION)


def test_fastapi_and_public_endpoints_use_single_version_constant():
    assert "from backend.version import ARIANE_VERSION" in MAIN_SOURCE
    assert "version=ARIANE_VERSION" in MAIN_SOURCE
    assert SYSTEM_API_SOURCE.count('"version": ARIANE_VERSION') == 2
    assert "application_version: str = ARIANE_VERSION" in PUBLIC_API_SOURCE
    assert '"version": "1.8.' not in "\n".join(
        (MAIN_SOURCE, SYSTEM_API_SOURCE, PUBLIC_API_SOURCE)
    )
