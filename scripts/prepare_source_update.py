"""Build a reviewable source-update candidate without deploying it.

This command is intended for the manually dispatched GitHub workflow. It
updates repository files in the workflow checkout so they can be reviewed in a
draft pull request. It never connects to or deploys the production server.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import stat
import subprocess
import sys
import tempfile
import urllib.request
from email.utils import parsedate_to_datetime
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import source_watch  # noqa: E402
from scripts.build_pp4_clinical_lr_snapshot import (  # noqa: E402
    TRANSCRIPTS,
    build as build_clinical_lr,
)
from scripts.build_enigma_erepo_vcep_registry import build as build_erepo  # noqa: E402
from scripts.check_enigma_clinical_lr_update import check_update as check_clinical_lr  # noqa: E402


EREPO_SOURCE = "enigma_erepo_vcep_assertions"
CLINICAL_LR_SOURCE = "enigma_clinical_lr"
SUPPORTED_SOURCES = (EREPO_SOURCE, CLINICAL_LR_SOURCE)
VERSION_PATH = ROOT / "backend" / "version.py"
CLINICAL_LR_SNAPSHOT = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.index.json"
CLINICAL_LR_METADATA = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json"


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clinical_lr_diff(
    pinned: dict[str, Any], candidate: dict[str, Any]
) -> dict[str, Any]:
    added_keys = sorted(candidate.keys() - pinned.keys())
    removed_keys = sorted(pinned.keys() - candidate.keys())
    comparison_fields = (
        "combined_lr",
        "criterion",
        "strength",
        "points",
        "automatic_application_status",
        "overlap_status",
    )
    modified = []
    for key in sorted(pinned.keys() & candidate.keys()):
        changes = {
            field: {"pinned": pinned[key].get(field), "candidate": candidate[key].get(field)}
            for field in comparison_fields
            if pinned[key].get(field) != candidate[key].get(field)
        }
        if changes:
            modified.append({"variant": key, "changes": changes})

    def direction(record: dict[str, Any]) -> str:
        points = int(record.get("points") or 0)
        return "pathogenic" if points > 0 else "benign" if points < 0 else "neutral"

    shared = pinned.keys() & candidate.keys()
    threshold_crossings = [
        key
        for key in sorted(shared)
        if (
            pinned[key].get("criterion"),
            pinned[key].get("strength"),
            pinned[key].get("points"),
        )
        != (
            candidate[key].get("criterion"),
            candidate[key].get("strength"),
            candidate[key].get("points"),
        )
    ]
    direction_changes = [
        key
        for key in sorted(shared)
        if direction(pinned[key]) != direction(candidate[key])
    ]
    return {
        "pinned_count": len(pinned),
        "candidate_count": len(candidate),
        "added": added_keys,
        "removed": removed_keys,
        "modified": modified,
        "threshold_crossings": threshold_crossings,
        "direction_changes": direction_changes,
    }


def _render_clinical_lr_review(
    *,
    check: dict[str, Any],
    diff: dict[str, Any],
    metadata: dict[str, Any],
    previous_version: str,
    candidate_version: str,
) -> str:
    lines = [
        "# UCSC ENIGMA clinical LR candidate update",
        "",
        f"Prepared: `{datetime.now(timezone.utc).isoformat()}`  ",
        f"ARIANE version: `{previous_version}` to `{candidate_version}`  ",
        f"Pinned release: `{check['pinned_release_date']}`  ",
        f"Candidate data version: `{check['remote_data_version']}`  ",
        f"Pinned BigBed SHA-256: `{check['pinned_bigbed_sha256']}`  ",
        f"Candidate BigBed SHA-256: `{check['remote_bigbed_sha256']}`",
        "",
        "## Record impact",
        "",
        f"Pinned records: **{diff['pinned_count']}**; candidate records: "
        f"**{diff['candidate_count']}**.  ",
        f"Added: **{len(diff['added'])}**; modified: **{len(diff['modified'])}**; "
        f"removed: **{len(diff['removed'])}**.  ",
        f"PP4/BP5 threshold crossings: **{len(diff['threshold_crossings'])}**; "
        f"evidence-direction changes: **{len(diff['direction_changes'])}**.",
        "",
        f"Candidate criterion counts: `{metadata.get('criteria', {})}`",
        "",
    ]
    for heading, records in (
        ("Added variants", diff["added"]),
        ("Removed variants", diff["removed"]),
        ("Threshold crossings", diff["threshold_crossings"]),
        ("Evidence-direction changes", diff["direction_changes"]),
    ):
        if records:
            lines.extend([f"### {heading}", "", *[f"- `{item}`" for item in records[:200]], ""])
            if len(records) > 200:
                lines.extend([f"Only the first 200 of {len(records)} are shown.", ""])
    if diff["modified"]:
        lines.extend(["### Modified LR records", ""])
        for item in diff["modified"][:200]:
            fields = ", ".join(
                f"{field}: {values['pinned']} -> {values['candidate']}"
                for field, values in item["changes"].items()
            )
            lines.append(f"- `{item['variant']}`: {fields}")
        lines.append("")
    lines.extend(
        [
            "## Required approval checks",
            "",
            "- [ ] Review UCSC dataVersion, source history and build record.",
            "- [ ] Review all PP4/BP5 threshold and evidence-direction changes.",
            "- [ ] Review added, removed and normalization-conflict records.",
            "- [ ] Confirm complete CI results.",
            "- [ ] Approve and merge only after expert review.",
            "",
            "The publisher-combined LR is retained as one evidence item. Publication "
            "components are not multiplied again. This candidate has not been deployed.",
            "",
        ]
    )
    return "\n".join(lines)


def _download_conversion_tool(manifest: dict[str, Any], destination: Path) -> None:
    conversion = manifest["conversion"]
    request = urllib.request.Request(
        conversion["tool_download_url"],
        headers={"User-Agent": "ARIANE-clinical-LR-candidate/1.0"},
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        destination.write_bytes(response.read())
    if _sha256(destination) != conversion["tool_sha256_used_for_pinned_bed"]:
        raise RuntimeError("Downloaded bigBedToBed checksum does not match the pinned tool")
    destination.chmod(destination.stat().st_mode | stat.S_IXUSR)


def _clinical_lr_row_counts(path: Path) -> tuple[int, int]:
    rows = 0
    without_hgvsc = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            rows += 1
            name = values[3] if len(values) > 3 else ""
            transcript, separator, c_notation = name.partition(":")
            if not separator or transcript not in TRANSCRIPTS or not c_notation.startswith("c."):
                without_hgvsc += 1
    return rows, without_hgvsc


def prepare_clinical_lr(*, review_output: Path) -> None:
    source_manifest_path = ROOT / "data/sources/enigma/clinical_lr_sources.manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    pinned_snapshot = json.loads(CLINICAL_LR_SNAPSHOT.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="ariane-clinical-lr-") as directory:
        workspace = Path(directory)
        check = check_clinical_lr(
            manifest_path=source_manifest_path,
            timeout=90,
            candidate_dir=workspace,
        )
        if check["status"] == "current":
            raise RuntimeError("Clinical LR already matches the approved baseline")
        release = str(check["remote_data_version"])
        if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", release):
            raise RuntimeError(f"Clinical LR dataVersion is not an ISO date: {release!r}")
        candidate_bb = Path(check["candidate_file"])
        source_dir = source_manifest_path.parent
        final_bb = source_dir / f"BRCAmfa.hg38.{release}.bb"
        final_bed = source_dir / f"BRCAmfa.hg38.{release}.bed"
        final_bb.write_bytes(candidate_bb.read_bytes())
        tool = workspace / "bigBedToBed"
        _download_conversion_tool(source_manifest, tool)
        subprocess.run([str(tool), str(final_bb), str(final_bed)], check=True)

    rows, without_hgvsc = _clinical_lr_row_counts(final_bed)
    old_dataset = next(iter(source_manifest["datasets"].values()))
    dataset = dict(old_dataset)
    dataset.update(
        {
            "file": final_bed.name,
            "sha256": _sha256(final_bed),
            "derived_from_bigbed_file": final_bb.name,
            "derived_from_bigbed_sha256": _sha256(final_bb),
            "remote_last_modified": parsedate_to_datetime(
                check["remote_headers"]["last_modified"]
            ).astimezone(timezone.utc).isoformat(),
            "remote_content_length": final_bb.stat().st_size,
            "expected_rows": rows,
            "expected_rows_without_transcript_hgvsc": without_hgvsc,
        }
    )
    dataset_id = f"ucsc_enigma_brcamfa_{release.replace('-', '_')}"
    source_manifest["datasets"] = {dataset_id: dataset}
    source_manifest["clinical_lr_data_release"]["release_date"] = release
    source_manifest["conversion"]["output_sha256_required"] = _sha256(final_bed)
    _write_json(source_manifest_path, source_manifest)

    metadata = build_clinical_lr(
        source=final_bed,
        output=CLINICAL_LR_SNAPSHOT,
        metadata_path=CLINICAL_LR_METADATA,
        source_manifest_path=source_manifest_path,
    )
    candidate_snapshot = json.loads(CLINICAL_LR_SNAPSHOT.read_text(encoding="utf-8"))
    diff = _clinical_lr_diff(pinned_snapshot, candidate_snapshot)
    previous_version, candidate_version = _increment_patch_version()
    review_output.parent.mkdir(parents=True, exist_ok=True)
    review_output.write_text(
        _render_clinical_lr_review(
            check=check,
            diff=diff,
            metadata=metadata,
            previous_version=previous_version,
            candidate_version=candidate_version,
        ),
        encoding="utf-8",
    )


def prepare_erepo(*, review_output: Path) -> None:
    registry_path = source_watch.DEFAULT_REGISTRY
    watch_registry = source_watch.load_registry(registry_path)
    source = next(
        item for item in watch_registry["sources"] if item["id"] == EREPO_SOURCE
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
    parser.add_argument("--source", choices=SUPPORTED_SOURCES, required=True)
    parser.add_argument("--review-output", type=Path, required=True)
    parser.add_argument(
        "--write",
        action="store_true",
        help="Required acknowledgement that repository candidate files will change.",
    )
    args = parser.parse_args()
    if not args.write:
        parser.error("--write is required; candidate generation changes checkout files")
    if args.source == EREPO_SOURCE:
        prepare_erepo(review_output=args.review_output)
    else:
        prepare_clinical_lr(review_output=args.review_output)
    print(f"Prepared {args.source} candidate in {args.review_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
