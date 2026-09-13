"""Structured request audit logging for HTTP adapters."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import sys
from typing import Any

from fastapi import Request


class RequestAuditLogger:
    """Write transport audit events without coupling the app root to logging I/O."""

    def __init__(self, logger: logging.Logger) -> None:
        self._logger = logger

    @classmethod
    def configured(cls) -> "RequestAuditLogger":
        logger = logging.getLogger("ariane.audit")
        logger.setLevel(logging.INFO)
        logger.propagate = False
        if not logger.handlers:
            stream_handler = logging.StreamHandler(sys.stdout)
            stream_handler.setFormatter(logging.Formatter("%(message)s"))
            logger.addHandler(stream_handler)

        path = Path(os.getenv("ARIANE_AUDIT_LOG", "/var/log/ariane/audit.jsonl"))
        try:
            if path.parent.is_dir() and not any(
                isinstance(handler, logging.FileHandler)
                and Path(handler.baseFilename) == path
                for handler in logger.handlers
            ):
                file_handler = logging.FileHandler(path, encoding="utf-8")
                file_handler.setFormatter(logging.Formatter("%(message)s"))
                logger.addHandler(file_handler)
        except OSError:
            logger.exception("Failed to open the audit log file")
        return cls(logger)

    def __call__(
        self,
        request: Request,
        event: str,
        level: str = "info",
        **fields: Any,
    ) -> None:
        record = {
            "log_type": "ariane_audit",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "event": event,
            "request_id": getattr(request.state, "request_id", ""),
            "source_ip": request.client.host if request.client else "unknown",
            "method": request.method,
            "path": request.url.path,
            "user_agent": request.headers.get("user-agent", "")[:300],
            **fields,
        }
        message = json.dumps(record, ensure_ascii=True, separators=(",", ":"))
        getattr(self._logger, level)(message)
