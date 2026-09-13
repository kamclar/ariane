"""Browser-side validation reporting schema."""

from __future__ import annotations

import json
from typing import Any, Literal

from pydantic import BaseModel, Field, model_validator


class ClientValidationRequest(BaseModel):
    form: Literal["single", "batch"]
    input: dict[str, Any]
    error: str = Field(min_length=1, max_length=500)

    @model_validator(mode="after")
    def limit_input_size(self):
        if len(json.dumps(self.input, ensure_ascii=False)) > 5000:
            raise ValueError("Client validation input is too large")
        return self
