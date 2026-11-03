"""Build a reviewable source-update candidate without deploying it.

This command is intended for the manually dispatched GitHub workflow. It
updates repository files in the workflow checkout so they can be reviewed in a
draft pull request. It never connects to or deploys the production server.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import source_watch  # noqa: E402
from scripts.build_enigma_erepo_vcep_registry import build as build_erepo  # noqa: E402


SUPPORTED_SOURCE = "enigma_erepo_vcep_assertions"
VERSION_PATH = ROOT / "backend" / "version.py"


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def _increment_patch_version(path: Path = VERSION_PATH) -> tuple[str, str]:
    text = path.read_text(encoding="utf-8")
    match = re.search(r'ARIANE_VERSION = "(\d+)\.(\d+)\.(\d+)"', text)
    if match is None:
        raise RuntimeError(f"Could not find ARIANE_VERSION in {path}")
    old = match.group(0)
    previous = ".".join(match.groups())
    current = f"{match.group(1)}.{match.group(2)}.{int(match.group(3)) + 1}"
    path.write_text(text.replace(old, f'ARIANE_VERSION = "{current}"'), encoding="utf-8")
    return previous, current


def _render_review(
    *,
    watch_result: dict[str, Any],
    candidate: dict[str, Any],
    previous_version: str,
    candidate_version: str,
) -> str:
    diff = watch_result.get("candidate_diff") or {}
    lines = [
        "# ClinGen ERepo ENIGMA candidate update",
        "",
        f"Prepared: `{datetime.now(timezone.utc).isoformat()}`  ",
        f"ARIANE version: `{previous_version}` to `{candidate_version}`  ",
        f"Pinned records: **{diff.get('pinned_count', 0)}**  ",
        f"Candidate records: **{diff.get('remote_count', candidate['record_count'])}**  ",
        f"Added: **{diff.get('added_count', 0)}**; modified: "
        f"**{diff.get('modified_count', 0)}**; removed: "
        f"**{diff.get('removed_count', 0)}**.",
        "",
        "## Candidate method status",
        "",
        f"Assertion-method versions: `{candidate['assertion_method_version_counts']}`  ",
        f"Runtime source statuses: `{candidate['source_status_counts']}`",
        "",
        "Unversioned and historical records remain display-only. A record is not "
        "made authoritative merely because it is present in ERepo.",
        "",
    ]
    source_watch.render_erepo_diff(lines, diff)
    lines.extend(
        [
            "",
            "## Required approval checks",
            "",
            "- [ ] Confirm the official ERepo records and assertion-method versions.",
            "- [ ] Review every changed classification and evidence-code assignment.",
            "- [ ] Confirm that historical and unversioned records remain display-only.",
            "- [ ] Review the targeted and complete CI results.",
            "- [ ] Approve and merge the pull request only after expert review.",
            "",
            "This candidate was built in GitHub Actions. It has not been deployed.",
            "",
        ]
    )
    return "\n".join(lines)


def prepare_erepo(*, review_output: Path) -> None:
    registry_path = source_watch.DEFAULT_REGISTRY
    watch_registry = source_watch.load_registry(registry_path)
    source = next(
        item for item in watch_registry["sources"] if item["id"] == SUPPORTED_SOURCE
    )
    watch_result = source_watch.check_erepo_export(
        source,
        float(watch_registry["policy"]["default_timeout_seconds"]),
    )
    if watch_result["status"] == "current":
        raise RuntimeError("ERepo already matches the approved repository baseline")

    candidate, metadata = build_erepo()
    registry_file = ROOT / str(source["registry_path"])
    metadata_file = ROOT / str(source["metadata_path"])
    _write_json(registry_file, candidate)
    _write_json(metadata_file, metadata)

    source["selected_rows_sha256"] = watch_result["remote_selected_rows_sha256"]
    source["selected_row_count"] = watch_result["remote_selected_row_count"]
    source["selected_rows_observed_utc"] = datetime.now(timezone.utc).isoformat()
    _write_json(registry_path, watch_registry)
    previous_version, candidate_version = _increment_patch_version()

    review_output.parent.mkdir(parents=True, exist_ok=True)
    review_output.write_text(
        _render_review(
            watch_result=watch_result,
            candidate=candidate,
            previous_version=previous_version,
            candidate_version=candidate_version,
        ),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=[SUPPORTED_SOURCE], required=True)
    parser.add_argument("--review-output", type=Path, required=True)
    parser.add_argument(
        "--write",
        action="store_true",
        help="Required acknowledgement that repository candidate files will change.",
    )
    args = parser.parse_args()
    if not args.write:
        parser.error("--write is required; candidate generation changes checkout files")
    prepare_erepo(review_output=args.review_output)
    print(f"Prepared {args.source} candidate in {args.review_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
