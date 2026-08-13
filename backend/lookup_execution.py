"""Observable execution wrapper for non-critical external evidence lookups."""
from __future__ import annotations

import asyncio
import logging
import os


LOGGER = logging.getLogger(__name__)
EXTERNAL_LOOKUP_TIMEOUT = 12
SERVICE_LOOKUP_TIMEOUTS = {
    # The lower-level API call defaults to 25 seconds. This margin covers rate
    # limiting and still keeps the complete request below nginx's 60-second
    # proxy timeout. Offline cache builders do not use this wrapper.
    "SpliceAI": int(os.environ.get("SPLICEAI_LOOKUP_TIMEOUT", "30")),
}


async def lookup_or_unavailable(func, default, service, diagnostics, *args):
    timeout = SERVICE_LOOKUP_TIMEOUTS.get(service, EXTERNAL_LOOKUP_TIMEOUT)
    try:
        return await asyncio.wait_for(
            asyncio.to_thread(func, *args),
            timeout=timeout,
        )
    except asyncio.TimeoutError:
        message = f"{service} lookup timed out after {timeout} seconds"
        LOGGER.warning(message, extra={"lookup_service": service})
        diagnostics.append(message)
        return default
    except Exception as exc:
        message = f"{service} lookup failed: {type(exc).__name__}: {exc}"
        LOGGER.exception(message, extra={"lookup_service": service})
        diagnostics.append(message)
        return default
