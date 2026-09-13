"""Browser shell route and static asset mounting."""

from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

from backend.api.session import issue_ui_session
from backend.presentation.frontend import FrontendAssets


def install_frontend(
    app: FastAPI,
    *,
    frontend_dir: Path,
    assets: FrontendAssets,
) -> None:
    app.mount("/static", StaticFiles(directory=frontend_dir / "static"), name="static")
    router = APIRouter()

    @router.get("/", response_class=HTMLResponse)
    async def index(request: Request):
        response = HTMLResponse(
            assets.html.replace("__ARIANE_ASSET_VERSION__", assets.version),
            headers={
                "Cache-Control": "no-store, max-age=0",
                "Pragma": "no-cache",
            },
        )
        issue_ui_session(request, response)
        return response

    app.include_router(router)
