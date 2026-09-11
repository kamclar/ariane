"""Short-lived signed sessions for the public browser interface.

The browser session is deliberately separate from public API authentication.
It lets the same-origin web interface call its internal transport routes without
placing a reusable API key in downloadable JavaScript.
"""

from __future__ import annotations

import base64
import hashlib
import hmac
import os
import secrets
import time
from urllib.parse import urlsplit

from fastapi import HTTPException, Request, Response


UI_SESSION_COOKIE = "ariane_ui_session"
UI_SESSION_SECRET_ENVIRONMENT = "ARIANE_UI_SESSION_SECRET"
UI_SESSION_MAX_AGE_SECONDS = 12 * 60 * 60
_TOKEN_VERSION = "v1"


class UiSessionConfigurationError(RuntimeError):
    """The configured browser-session signing secret is invalid."""


def _urlsafe_encode(value: bytes) -> str:
    return base64.urlsafe_b64encode(value).decode("ascii").rstrip("=")


def _urlsafe_decode(value: str) -> bytes:
    padding = "=" * (-len(value) % 4)
    return base64.urlsafe_b64decode(value + padding)


class UiSessionAuthenticator:
    """Issue and validate stateless HMAC-signed browser sessions."""

    def __init__(self) -> None:
        configured = os.getenv(UI_SESSION_SECRET_ENVIRONMENT, "").strip()
        if len(configured.encode("utf-8")) < 32:
            raise UiSessionConfigurationError(
                f"{UI_SESSION_SECRET_ENVIRONMENT} must contain at least 32 bytes"
            )
        self._secret = configured.encode("utf-8")

    def issue(self, *, now: int | None = None) -> str:
        issued_at = int(time.time()) if now is None else int(now)
        payload = f"{_TOKEN_VERSION}.{issued_at}.{secrets.token_urlsafe(18)}"
        signature = hmac.new(
            self._secret,
            payload.encode("ascii"),
            hashlib.sha256,
        ).digest()
        return f"{payload}.{_urlsafe_encode(signature)}"

    def validate(self, token: str, *, now: int | None = None) -> bool:
        if not token or len(token) > 512:
            return False
        try:
            version, issued_text, nonce, signature_text = token.split(".", 3)
            issued_at = int(issued_text)
            signature = _urlsafe_decode(signature_text)
        except (AttributeError, ValueError, TypeError):
            return False
        if version != _TOKEN_VERSION or not nonce:
            return False
        current_time = int(time.time()) if now is None else int(now)
        if issued_at > current_time + 60:
            return False
        if current_time - issued_at > UI_SESSION_MAX_AGE_SECONDS:
            return False
        payload = f"{version}.{issued_at}.{nonce}"
        expected = hmac.new(
            self._secret,
            payload.encode("ascii"),
            hashlib.sha256,
        ).digest()
        return hmac.compare_digest(signature, expected)


UI_SESSION_AUTHENTICATOR = UiSessionAuthenticator()


def _secure_cookie(request: Request) -> bool:
    forwarded = request.headers.get("x-forwarded-proto", "").split(",", 1)[0]
    scheme = forwarded.strip().lower() or request.url.scheme.lower()
    return scheme == "https"


def issue_ui_session(request: Request, response: Response) -> None:
    response.set_cookie(
        UI_SESSION_COOKIE,
        UI_SESSION_AUTHENTICATOR.issue(),
        max_age=UI_SESSION_MAX_AGE_SECONDS,
        httponly=True,
        secure=_secure_cookie(request),
        samesite="strict",
        path="/",
    )


async def require_ui_session(request: Request) -> str:
    token = request.cookies.get(UI_SESSION_COOKIE, "")
    if not UI_SESSION_AUTHENTICATOR.validate(token):
        raise HTTPException(
            status_code=401,
            detail={
                "code": "ui_session_required",
                "message": "Reload the ARIANE page to establish a browser session",
                "retryable": False,
            },
        )

    fetch_site = request.headers.get("sec-fetch-site", "").strip().lower()
    if fetch_site and fetch_site not in {"same-origin", "none"}:
        raise HTTPException(
            status_code=403,
            detail={
                "code": "cross_site_ui_request_rejected",
                "message": "The browser request must originate from ARIANE",
                "retryable": False,
            },
        )
    if request.method not in {"GET", "HEAD", "OPTIONS"}:
        origin = request.headers.get("origin", "").strip()
        forwarded_scheme = request.headers.get(
            "x-forwarded-proto", ""
        ).split(",", 1)[0].strip().lower()
        scheme = forwarded_scheme or request.url.scheme.lower()
        forwarded_host = request.headers.get(
            "x-forwarded-host", ""
        ).split(",", 1)[0].strip().lower()
        host = forwarded_host or request.headers.get("host", "").strip().lower()
        parsed_origin = urlsplit(origin)
        if (
            not origin
            or parsed_origin.scheme.lower() != scheme
            or parsed_origin.netloc.lower() != host
        ):
            raise HTTPException(
                status_code=403,
                detail={
                    "code": "ui_origin_required",
                    "message": "The browser request must originate from ARIANE",
                    "retryable": False,
                },
            )
    request.state.client_surface = "web_ui"
    return "web_ui"
