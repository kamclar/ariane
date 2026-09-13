"""Shared HTTP exception-to-response mapping."""

from __future__ import annotations

from collections.abc import Callable

from fastapi import FastAPI, HTTPException, Request
from fastapi.encoders import jsonable_encoder
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from backend.api.public import PublicApiError, PublicApiErrorResponse, build_metadata
from backend.domain.classification import CLASSIFICATION_ENGINE_ID


AuditWriter = Callable[..., None]


def install_exception_handlers(app: FastAPI, *, audit: AuditWriter) -> None:
    @app.exception_handler(RequestValidationError)
    async def validation_error_handler(request: Request, exc: RequestValidationError):
        audit(
            request,
            "validation_error",
            level="warning",
            input=jsonable_encoder(exc.body),
            errors=jsonable_encoder(exc.errors()),
        )
        if request.url.path.startswith("/api/v1/"):
            messages = []
            for error in exc.errors():
                location = ".".join(str(value) for value in error.get("loc", ()))
                message = str(error.get("msg", "Invalid request"))
                if message.startswith("Value error, "):
                    message = message[len("Value error, "):]
                rendered = f"{location}: {message}" if location else message
                if rendered not in messages:
                    messages.append(rendered)
            response = PublicApiErrorResponse(
                metadata=build_metadata(
                    request_id=request.state.request_id,
                    engine=CLASSIFICATION_ENGINE_ID,
                ),
                error=PublicApiError(
                    code="invalid_request",
                    message="; ".join(messages) or "Invalid request",
                    retryable=False,
                ),
            )
            return JSONResponse(
                status_code=422,
                content=response.model_dump(mode="json", exclude_none=True),
            )
        return JSONResponse(
            status_code=422,
            content={"detail": jsonable_encoder(exc.errors())},
        )

    @app.exception_handler(HTTPException)
    async def http_error_handler(request: Request, exc: HTTPException):
        if not request.url.path.startswith("/api/v1/"):
            return JSONResponse(
                status_code=exc.status_code,
                content={"detail": jsonable_encoder(exc.detail)},
                headers=exc.headers,
            )
        code_by_status = {
            401: "authentication_required",
            404: "not_found",
            422: "variant_not_classifiable",
            429: "rate_limited",
            503: "evidence_unavailable",
        }
        error_headers = exc.headers or {}
        if isinstance(exc.detail, dict):
            error_code = str(
                exc.detail.get("code")
                or code_by_status.get(exc.status_code, "request_failed")
            )
            error_message = str(exc.detail.get("message") or "Request failed")
            retryable = bool(exc.detail.get("retryable", False))
        else:
            error_code = error_headers.get(
                "X-ARIANE-Error-Code",
                code_by_status.get(exc.status_code, "request_failed"),
            )
            error_message = str(exc.detail)
            retryable_header = error_headers.get("X-ARIANE-Retryable")
            retryable = (
                retryable_header.lower() == "true"
                if retryable_header is not None
                else exc.status_code in {429, 502, 503, 504}
            )
        response = PublicApiErrorResponse(
            metadata=build_metadata(
                request_id=request.state.request_id,
                engine=CLASSIFICATION_ENGINE_ID,
            ),
            error=PublicApiError(
                code=error_code,
                message=error_message,
                retryable=retryable,
            ),
        )
        return JSONResponse(
            status_code=exc.status_code,
            content=response.model_dump(mode="json", exclude_none=True),
            headers=exc.headers,
        )
