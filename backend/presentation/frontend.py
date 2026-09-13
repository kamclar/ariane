"""Load the browser shell, feature templates, and versioned static assets."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path


PARTIAL_NAMES = (
    "search",
    "classification-result",
    "manual-review",
    "external-comparison",
    "rules-explorer",
    "batch",
)

NESTED_PARTIAL_NAMES = (
    "manual-review-form",
    "manual-review-clinical-fields",
    "manual-review-functional-fields",
    "manual-review-rna-splice-fields",
    "manual-review-protein-ps1-fields",
    "manual-review-result",
)


@dataclass(frozen=True)
class FrontendAssets:
    html: str
    version: str


def load_frontend_template(frontend_dir: Path) -> str:
    """Compose the page from explicit feature templates."""
    shell_path = frontend_dir / "index.html"
    try:
        html = shell_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise RuntimeError(f"Frontend shell cannot be loaded: {shell_path}: {exc}") from exc

    partial_dir = frontend_dir / "templates"
    for name in (*PARTIAL_NAMES, *NESTED_PARTIAL_NAMES):
        marker = f"<!-- ARIANE_INCLUDE:{name} -->"
        if html.count(marker) != 1:
            raise RuntimeError(f"Frontend templates must contain one {marker} marker")
        path = partial_dir / f"{name}.html"
        try:
            content = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise RuntimeError(f"Frontend template cannot be loaded: {path}: {exc}") from exc
        html = html.replace(marker, content.rstrip())
    if "<!-- ARIANE_INCLUDE:" in html:
        raise RuntimeError("Frontend shell contains an unknown include marker")
    return html


def frontend_asset_version(frontend_dir: Path) -> str:
    """Bind browser asset URLs to the release-independent file contents."""
    static_dir = frontend_dir / "static"
    asset_paths = sorted(
        path
        for path in static_dir.rglob("*")
        if path.is_file()
    )
    digest = hashlib.sha256()
    for path in asset_paths:
        digest.update(path.relative_to(static_dir).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()[:12]


def load_frontend_assets(frontend_dir: Path, application_version: str) -> FrontendAssets:
    digest = frontend_asset_version(frontend_dir)
    return FrontendAssets(
        html=load_frontend_template(frontend_dir),
        version=f"{application_version}-{digest}",
    )
