"""Process-wide registry of degraded data sources visible to API users."""
from __future__ import annotations

import threading


_LOCK = threading.Lock()
_ISSUES: dict[str, str] = {}


def register_issue(component: str, reason: str) -> None:
    with _LOCK:
        _ISSUES[component] = reason


def clear_issue(component: str) -> None:
    with _LOCK:
        _ISSUES.pop(component, None)


def get_data_issues() -> list[dict[str, str]]:
    with _LOCK:
        return [
            {"component": component, "reason": reason}
            for component, reason in sorted(_ISSUES.items())
        ]


def get_user_warnings() -> list[str]:
    return [
        f"Data source degraded: {issue['component']}: {issue['reason']}"
        for issue in get_data_issues()
    ]
