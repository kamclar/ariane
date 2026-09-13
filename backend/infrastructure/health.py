"""Thread-safe runtime health state owned by the application composition root."""
from __future__ import annotations

import threading
import re


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


class DataHealthRegistry:
    """Collect degraded-source diagnostics for one application runtime."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._issues: dict[str, str] = {}

    def register(self, component: str, reason: str) -> None:
        with self._lock:
            self._issues[component] = reason

    def clear(self, component: str) -> None:
        with self._lock:
            self._issues.pop(component, None)

    def issues(self) -> list[dict[str, str]]:
        with self._lock:
            return [
                {"component": component, "reason": _user_safe_reason(reason)}
                for component, reason in sorted(self._issues.items())
            ]

    def user_warnings(self) -> list[str]:
        warnings = []
        for issue in self.issues():
            if issue["reason"].startswith((
                "score was obtained and used",
                "coordinates were obtained and used",
            )):
                prefix = "Runtime cache persistence warning"
            else:
                prefix = "Data source degraded"
            warnings.append(f"{prefix}: {issue['component']}: {issue['reason']}")
        return warnings
