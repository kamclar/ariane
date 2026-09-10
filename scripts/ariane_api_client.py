#!/usr/bin/env python3
"""Small reference client for the versioned ARIANE public API."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import time
import uuid

import requests


DEFAULT_BASE_URL = "https://ariane-app.duckdns.org/api/v1"
RETRYABLE_HTTP_STATUS = {429, 502, 503, 504}


def read_variants(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene", "c_notation"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Input TSV must contain gene and c_notation columns")
        variants = []
        for row_number, row in enumerate(reader, start=2):
            gene = str(row.get("gene", "")).strip()
            c_notation = str(row.get("c_notation", "")).strip()
            if not gene or not c_notation:
                raise ValueError(f"Input row {row_number} has no gene or c_notation")
            item = {"gene": gene, "c_notation": c_notation}
            for field in ("p_notation", "assembly", "dup_type"):
                value = str(row.get(field, "")).strip()
                if value:
                    item[field] = value
            variants.append(item)
    return variants


def request_json(
    session: requests.Session,
    method: str,
    url: str,
    *,
    payload: dict | None = None,
    request_id: str,
    attempts: int = 3,
) -> dict:
    headers = {"X-Request-ID": request_id}
    for attempt in range(attempts):
        response = session.request(
            method,
            url,
            json=payload,
            headers=headers,
            timeout=180,
        )
        if response.ok:
            return response.json()
        try:
            error_payload = response.json()
        except ValueError:
            error_payload = {"error": {"message": response.text[:500]}}
        message = (
            error_payload.get("error", {}).get("message")
            or f"HTTP {response.status_code}"
        )
        if response.status_code not in RETRYABLE_HTTP_STATUS:
            raise RuntimeError(f"ARIANE request failed: {message}")
        if attempt + 1 == attempts:
            raise RuntimeError(f"ARIANE request failed after retries: {message}")
        retry_after = response.headers.get("Retry-After", "")
        try:
            delay = max(1.0, min(float(retry_after), 30.0))
        except ValueError:
            delay = min(2 ** attempt, 30)
        time.sleep(delay)
    raise RuntimeError("ARIANE request did not complete")


def chunks(values: list[dict[str, str]], size: int):
    for offset in range(0, len(values), size):
        yield offset, values[offset:offset + size]


def classify_batch_with_retries(
    session: requests.Session,
    base_url: str,
    batch: list[dict[str, str]],
    *,
    offset: int,
    attempts: int = 3,
) -> list[dict]:
    pending = list(enumerate(batch))
    completed: dict[int, dict] = {}
    for attempt in range(attempts):
        response = request_json(
            session,
            "POST",
            f"{base_url}/classify/batch",
            payload={"variants": [variant for _, variant in pending]},
            request_id=(
                f"batch-{offset:08d}-attempt-{attempt + 1}-"
                f"{uuid.uuid4().hex[:8]}"
            ),
        )
        retry_items: list[tuple[int, dict[str, str]]] = []
        for response_item, (local_index, variant) in zip(
            response["results"],
            pending,
            strict=True,
        ):
            response_item["index"] = local_index
            response_item["input_index"] = offset + local_index
            error = response_item.get("error") or {}
            if (
                response_item.get("status") == "error"
                and error.get("retryable") is True
                and attempt + 1 < attempts
            ):
                retry_items.append((local_index, variant))
            else:
                completed[local_index] = response_item
        if not retry_items:
            break
        pending = retry_items
        time.sleep(min(2 ** attempt, 30))
    return [completed[index] for index in sorted(completed)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input TSV file")
    parser.add_argument("--output", type=Path, required=True, help="Output JSONL file")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--batch-size", type=int, default=5)
    parser.add_argument(
        "--api-key",
        default=os.getenv("ARIANE_API_KEY", ""),
        help="API key, preferably supplied through ARIANE_API_KEY",
    )
    args = parser.parse_args()

    if not args.api_key:
        raise ValueError(
            "Provide an API key with ARIANE_API_KEY or --api-key"
        )

    base_url = args.base_url.rstrip("/")
    session = requests.Session()
    session.headers.update({
        "User-Agent": "ariane-public-api-client/1.0",
        "X-ARIANE-API-Key": args.api_key,
    })
    capabilities = request_json(
        session,
        "GET",
        f"{base_url}/capabilities",
        request_id=f"capabilities-{uuid.uuid4().hex[:12]}",
    )
    maximum = int(capabilities["limits"]["maximum_batch_items"])
    if args.batch_size < 1 or args.batch_size > maximum:
        raise ValueError(f"Batch size must be between 1 and {maximum}")

    variants = read_variants(args.input)
    with args.output.open("w", encoding="utf-8", newline="\n") as output:
        for offset, batch in chunks(variants, args.batch_size):
            results = classify_batch_with_retries(
                session,
                base_url,
                batch,
                offset=offset,
            )
            for item in results:
                output.write(json.dumps(item, ensure_ascii=False) + "\n")
            output.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
