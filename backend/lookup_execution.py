"""Observable execution wrapper for non-critical external evidence lookups."""
from __future__ import annotations

import asyncio
import logging


LOGGER = logging.getLogger(__name__)
EXTERNAL_LOOKUP_TIMEOUT = 12


async def lookup_or_unavailable(func, default, service, diagnostics, *args):
    try:
        return await asyncio.wait_for(
            asyncio.to_thread(func, *args),
            timeout=EXTERNAL_LOOKUP_TIMEOUT,
        )
    except asyncio.TimeoutError:
        message = f"{service} lookup timed out after {EXTERNAL_LOOKUP_TIMEOUT} seconds"
        LOGGER.warning(message, extra={"lookup_service": service})
        diagnostics.append(message)
        return default
    except Exception as exc:
        message = f"{service} lookup failed: {type(exc).__name__}: {exc}"
        LOGGER.exception(message, extra={"lookup_service": service})
        diagnostics.append(message)
        return default
