"""Process-wide registry of degraded data sources visible to API users."""
from __future__ import annotations

import threading
import re


_LOCK = threading.Lock()
_ISSUES: dict[str, str] = {}


_PROJECT_PATH_RE = re.compile(
    r"(?:[A-Za-z]:)?[\\/](?:[^\\/\s]+[\\/])*ariane[\\/][^\s,;]+",
    re.IGNORECASE,
)


def _user_safe_reason(reason: str) -> str:
    """Hide deployment-specific prefixes while retaining useful project paths."""
    def shorten(match: re.Match[str]) -> str:
        normalized = match.group(0).replace("\\", "/")
        marker = normalized.lower().rfind("/ariane/")
        return f"…ariane/{normalized[marker + len('/ariane/'):]}"

    return _PROJECT_PATH_RE.sub(shorten, reason)


def register_issue(component: str, reason: str) -> None:
    with _LOCK:
        _ISSUES[component] = reason


def clear_issue(component: str) -> None:
    with _LOCK:
        _ISSUES.pop(component, None)


def get_data_issues() -> list[dict[str, str]]:
    with _LOCK:
        return [
            {"component": component, "reason": _user_safe_reason(reason)}
            for component, reason in sorted(_ISSUES.items())
        ]


def get_user_warnings() -> list[str]:
    return [
        f"Data source degraded: {issue['component']}: {issue['reason']}"
        for issue in get_data_issues()
    ]
