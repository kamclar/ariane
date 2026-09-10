"""Resolve an analytics identity without treating network addresses as users."""

from __future__ import annotations

from dataclasses import dataclass
import os
import re
import uuid

from fastapi import Request, Response


VISITOR_COOKIE = "ariane_visitor_id"
_ACCOUNT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@+:-]{0,99}$")


@dataclass(frozen=True)
class UsageIdentity:
    actor_id: str
    actor_type: str
    new_visitor_id: str = ""
    secure_cookie: bool = True

    def set_cookie(self, response: Response) -> None:
        if not self.new_visitor_id:
            return
        response.set_cookie(
            VISITOR_COOKIE,
            self.new_visitor_id,
            max_age=365 * 24 * 60 * 60,
            httponly=True,
            secure=self.secure_cookie,
            samesite="lax",
        )


def resolve_usage_identity(request: Request) -> UsageIdentity:
    """Prefer a proxy-authenticated account, otherwise use a browser identifier."""
    trusted_header = os.getenv("ARIANE_TRUSTED_USER_HEADER", "").strip().lower()
    if trusted_header:
        value = request.headers.get(trusted_header, "").strip()
        if _ACCOUNT_RE.fullmatch(value):
            return UsageIdentity(actor_id=value, actor_type="account")

    visitor_id = request.cookies.get(VISITOR_COOKIE, "").strip()
    try:
        visitor_id = str(uuid.UUID(visitor_id))
    except (ValueError, AttributeError):
        visitor_id = str(uuid.uuid4())
        return UsageIdentity(
            actor_id=visitor_id,
            actor_type="visitor",
            new_visitor_id=visitor_id,
            secure_cookie=request.url.scheme == "https",
        )
    return UsageIdentity(
        actor_id=visitor_id,
        actor_type="visitor",
        secure_cookie=request.url.scheme == "https",
    )
