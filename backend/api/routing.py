"""Registration helpers for equivalent public-API and browser-UI routes."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from fastapi import APIRouter, Depends

from backend.api.auth import require_public_api_key
from backend.api.session import require_ui_session


Endpoint = Callable[..., Any]


class PairedApiRoutes:
    """Register one handler behind the two intentionally different auth gates."""

    def __init__(self, router: APIRouter) -> None:
        self._router = router

    def _route(
        self,
        path: str,
        *,
        methods: Sequence[str],
        **route_options: Any,
    ) -> Callable[[Endpoint], Endpoint]:
        if not path.startswith("/"):
            raise ValueError("Paired API paths must start with '/'")
        if "dependencies" in route_options or "include_in_schema" in route_options:
            raise ValueError(
                "Paired route authentication and schema visibility are fixed by policy"
            )

        def register(endpoint: Endpoint) -> Endpoint:
            self._router.add_api_route(
                f"/api{path}",
                endpoint,
                methods=list(methods),
                dependencies=[Depends(require_public_api_key)],
                **route_options,
            )
            self._router.add_api_route(
                f"/ui-api{path}",
                endpoint,
                methods=list(methods),
                dependencies=[Depends(require_ui_session)],
                include_in_schema=False,
                **route_options,
            )
            return endpoint

        return register

    def get(
        self,
        path: str,
        **route_options: Any,
    ) -> Callable[[Endpoint], Endpoint]:
        return self._route(path, methods=("GET",), **route_options)

    def post(
        self,
        path: str,
        **route_options: Any,
    ) -> Callable[[Endpoint], Endpoint]:
        return self._route(path, methods=("POST",), **route_options)
