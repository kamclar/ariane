#!/usr/bin/env python3
"""Validate a running SpliceAI service against the active ARIANE profile."""

from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROFILE_PATH = (
    PROJECT_ROOT / "data" / "spliceai" / "enigma_v1_2_spliceai_profile.json"
)
CASES_PATH = (
    PROJECT_ROOT / "data" / "spliceai" / "local_service_validation_cases.json"
)


class ValidationError(RuntimeError):
    """Raised when the service does not match the versioned scoring profile."""


def load_json_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValidationError(f"Expected a JSON object in {path}")
    return value


def build_url(base_url: str, profile: dict[str, Any], variant: str) -> str:
    base_url = base_url.strip()
    if not base_url.endswith("/"):
        base_url += "/"
    query = urllib.parse.urlencode(
        {
            "hg": 38,
            "variant": variant,
            "distance": profile["max_distance"],
            "mask": profile["mask"],
            "bc": profile["annotation_subset"],
            "show-ref-alt": 1,
        }
    )
    return f"{base_url}?{query}"


def fetch_json(url: str, timeout: float) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "ARIANE-SpliceAI-validator"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        value = json.loads(response.read())
    if not isinstance(value, dict):
        raise ValidationError("SpliceAI response is not a JSON object")
    return value


def _as_float_map(row: dict[str, Any], expected: dict[str, Any]) -> dict[str, float]:
    result: dict[str, float] = {}
    for field in expected:
        if field not in row:
            raise ValidationError(f"Reference transcript row is missing {field}")
        try:
            result[field] = float(row[field])
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"Reference transcript field {field} is not numeric") from exc
    return result


def validate_payload(
    payload: dict[str, Any],
    profile: dict[str, Any],
    case: dict[str, Any],
) -> dict[str, Any]:
    try:
        response_distance = int(payload.get("distance"))
        response_mask = int(payload.get("mask"))
    except (TypeError, ValueError) as exc:
        raise ValidationError("Response does not identify distance and mask") from exc

    assembly = str(payload.get("genomeVersion") or payload.get("hg") or "")
    if assembly != "38":
        raise ValidationError(f"Expected GRCh38, received {assembly!r}")
    if response_distance != int(profile["max_distance"]):
        raise ValidationError(
            f"Expected distance {profile['max_distance']}, received {response_distance}"
        )
    if response_mask != int(profile["mask"]):
        raise ValidationError(f"Expected mask {profile['mask']}, received {response_mask}")
    if str(payload.get("bc") or "") != str(profile["annotation_subset"]):
        raise ValidationError(
            f"Expected annotation {profile['annotation_subset']!r}, "
            f"received {payload.get('bc')!r}"
        )

    rows = payload.get("scores")
    if not isinstance(rows, list) or not rows:
        raise ValidationError("Response has no transcript score rows")
    matching_rows = [
        row
        for row in rows
        if isinstance(row, dict)
        and row.get("t_id") == case["transcript"]
        and case["refseq"] in (row.get("t_refseq_ids") or [])
    ]
    if len(matching_rows) != 1:
        raise ValidationError(
            f"Expected one row for {case['transcript']} / {case['refseq']}, "
            f"received {len(matching_rows)}"
        )
    row = matching_rows[0]

    expected_groups = (
        case["delta_scores"],
        case["reference_scores"],
        case["alternate_scores"],
    )
    for expected in expected_groups:
        observed = _as_float_map(row, expected)
        mismatches = {
            field: {"expected": float(expected[field]), "observed": observed[field]}
            for field in expected
            if observed[field] != float(expected[field])
        }
        if mismatches:
            raise ValidationError(
                f"Score mismatch for {case['id']}: "
                + json.dumps(mismatches, sort_keys=True)
            )

    delta_scores = _as_float_map(row, case["delta_scores"])
    return {
        "id": case["id"],
        "variant": case["variant"],
        "transcript": case["transcript"],
        "score": max(delta_scores.values()),
    }


def validate_service(
    base_url: str,
    profile_path: Path = PROFILE_PATH,
    cases_path: Path = CASES_PATH,
    timeout: float = 180.0,
) -> dict[str, Any]:
    profile = load_json_object(profile_path)
    fixture = load_json_object(cases_path)
    if fixture.get("profile_id") != profile.get("profile_id"):
        raise ValidationError("Validation cases do not belong to the active profile")
    cases = fixture.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValidationError("No SpliceAI validation cases are configured")

    results = []
    for case in cases:
        if not isinstance(case, dict):
            raise ValidationError("Invalid SpliceAI validation case")
        payload = fetch_json(build_url(base_url, profile, case["variant"]), timeout)
        results.append(validate_payload(payload, profile, case))
    return {
        "status": "ok",
        "profile_id": profile["profile_id"],
        "docker_image": profile["approved_engine"]["docker_image"],
        "cases": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://127.0.0.1:8081/spliceai/")
    parser.add_argument("--profile", type=Path, default=PROFILE_PATH)
    parser.add_argument("--cases", type=Path, default=CASES_PATH)
    parser.add_argument("--timeout", type=float, default=180.0)
    args = parser.parse_args()
    try:
        result = validate_service(args.url, args.profile, args.cases, args.timeout)
    except Exception as exc:
        print(f"SpliceAI validation failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
