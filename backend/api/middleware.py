"""Application-wide HTTP middleware registration."""

from __future__ import annotations

from collections.abc import Callable
import time
import uuid

from fastapi import FastAPI, Request

from backend.api.public import PUBLIC_API_VERSION
from backend.version import ARIANE_VERSION


AuditWriter = Callable[..., None]

CONTENT_SECURITY_POLICY = "; ".join((
    "default-src 'self'",
    "script-src 'self' 'unsafe-eval' https://cdn.jsdelivr.net",
    "style-src 'self' 'unsafe-inline'",
    "img-src 'self' data:",
    "connect-src 'self'",
    "font-src 'self'",
    "object-src 'none'",
    "base-uri 'self'",
    "form-action 'self'",
    "frame-ancestors 'none'",
))


def request_id(value: str | None) -> str:
    if value:
        candidate = value.strip()
        if (
            1 <= len(candidate) <= 128
            and all(char.isalnum() or char in "-._:" for char in candidate)
        ):
            return candidate
    return uuid.uuid4().hex


def install_http_middleware(app: FastAPI, *, audit: AuditWriter) -> None:
    """Install security and audit middleware in their deliberate order."""

    @app.middleware("http")
    async def security_headers(request: Request, call_next):
        response = await call_next(request)
        response.headers.setdefault("Content-Security-Policy", CONTENT_SECURITY_POLICY)
        response.headers.setdefault("X-Frame-Options", "DENY")
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("Referrer-Policy", "no-referrer")
        response.headers.setdefault(
            "Permissions-Policy",
            "camera=(), microphone=(), geolocation=(), payment=(), usb=()",
        )
        return response

    @app.middleware("http")
    async def audit_request(request: Request, call_next):
        request.state.request_id = request_id(request.headers.get("x-request-id"))
        started = time.monotonic()
        try:
            response = await call_next(request)
        except Exception as exc:
            audit(
                request,
                "request_exception",
                level="exception",
                duration_ms=round((time.monotonic() - started) * 1000, 1),
                error_type=type(exc).__name__,
                error=str(exc)[:2000],
            )
            raise
        response.headers["X-Request-ID"] = request.state.request_id
        if request.url.path.startswith("/api/v1/"):
            response.headers["X-ARIANE-API-Version"] = PUBLIC_API_VERSION
            response.headers["X-ARIANE-Version"] = ARIANE_VERSION
        log_completion = (
            request.url.path.startswith("/admin/")
            or request.url.path.startswith("/ui-api/")
            or request.url.path.startswith("/api/")
            and request.url.path != "/api/health"
        )
        if log_completion:
            audit(
                request,
                "request_completed",
                level="warning" if response.status_code >= 400 else "info",
                status_code=response.status_code,
                duration_ms=round((time.monotonic() - started) * 1000, 1),
            )
        return response
