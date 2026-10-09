"""Read-only monitoring of upstream data sources used by ARIANE.

This module detects upstream changes.  It never downloads into production data
paths, rebuilds snapshots, changes manifests, or activates a new release.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import sys
import urllib.request
from collections.abc import Callable, Iterable
from datetime import date, datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = ROOT / "data" / "sources" / "source_watch_registry.json"
MAX_EREPO_EXPORT_BYTES = 50 * 1024 * 1024
AUTOMATED_CHECKERS = {
    "cspec_records",
    "erepo_export",
    "clinical_lr",
    "gnomad_releases",
}
SUPPORTED_CHECKERS = AUTOMATED_CHECKERS | {"manual"}
REVIEW_STATUSES = {"update_available", "upstream_changed", "manual_review_due"}

UrlOpener = Callable[..., Any]
SourceChecker = Callable[[dict[str, Any], float], dict[str, Any]]


class SourceWatchConfigurationError(RuntimeError):
    """Raised when the source-watch registry violates its safety contract."""


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise SourceWatchConfigurationError(f"Expected a JSON object in {path}")
    return payload


def _project_path(relative_path: str) -> Path:
    path = (ROOT / relative_path).resolve()
    if not path.is_relative_to(ROOT):
        raise SourceWatchConfigurationError(
            f"Source-watch path escapes the repository: {relative_path}"
        )
    return path


def load_registry(path: Path = DEFAULT_REGISTRY) -> dict[str, Any]:
    registry = _load_json(path)
    if registry.get("schema_version") != 1:
        raise SourceWatchConfigurationError("Unsupported source-watch registry schema")

    policy = registry.get("policy")
    if not isinstance(policy, dict):
        raise SourceWatchConfigurationError("Source-watch policy is missing")
    if policy.get("automatic_activation") is not False:
        raise SourceWatchConfigurationError(
            "Source watch must explicitly prohibit automatic activation"
        )
    if policy.get("runtime_network_dependency") is not False:
        raise SourceWatchConfigurationError(
            "Source watch must not become a runtime network dependency"
        )

    sources = registry.get("sources")
    if not isinstance(sources, list) or not sources:
        raise SourceWatchConfigurationError("Source-watch registry has no sources")
    identifiers = [str(source.get("id") or "") for source in sources]
    if any(not identifier for identifier in identifiers) or len(set(identifiers)) != len(
        identifiers
    ):
        raise SourceWatchConfigurationError("Source IDs must be non-empty and unique")

    for source in sources:
        checker = source.get("checker")
        if not str(source.get("name") or "").strip():
            raise SourceWatchConfigurationError(
                f"Source {source.get('id')!r} has no display name"
            )
        if checker not in SUPPORTED_CHECKERS:
            raise SourceWatchConfigurationError(
                f"Unsupported checker {checker!r} for {source.get('id')!r}"
            )
        cadence_days = source.get("cadence_days")
        if checker != "manual" and (
            not isinstance(cadence_days, int) or cadence_days < 1
        ):
            raise SourceWatchConfigurationError(
                f"Automated source {source['id']!r} needs a positive cadence_days"
            )
        for field in ("manifest_path", "metadata_path"):
            if source.get(field):
                local_path = _project_path(str(source[field]))
                if not local_path.is_file():
                    raise SourceWatchConfigurationError(
                        f"Referenced {field} does not exist for {source['id']}: "
                        f"{source[field]}"
                    )
        for relative_path in source.get("tracked_files", []):
            local_path = _project_path(str(relative_path))
            if not local_path.is_file():
                raise SourceWatchConfigurationError(
                    f"Tracked file does not exist for {source['id']}: {relative_path}"
                )
        if checker == "cspec_records":
            records = source.get("records")
            if not isinstance(records, list) or not records:
                raise SourceWatchConfigurationError(
                    f"CSpec source {source['id']!r} has no records"
                )
            for record in records:
                missing = {
                    "gene",
                    "id",
                    "api_url",
                    "pinned_version",
                    "observed_revision",
                    "observed_modified",
                } - record.keys()
                if missing:
                    raise SourceWatchConfigurationError(
                        f"CSpec record in {source['id']!r} is missing "
                        + ", ".join(sorted(missing))
                    )
            manifest = _load_json(_project_path(str(source["manifest_path"])))
            genes_payload = manifest.get("genes", manifest.get("targets", []))
            if isinstance(genes_payload, dict):
                targets = [
                    {"gene": gene, **configuration}
                    for gene, configuration in genes_payload.items()
                ]
            else:
                targets = genes_payload
            active_targets = {
                str(target["vcep_specification"]["id"]): {
                    "gene": str(target["gene"]),
                    "version": str(target["vcep_specification"]["version"]),
                }
                for target in targets
                if target.get("activation_status", "active") == "active"
                and target.get("vcep_specification")
            }
            watched_targets = {
                str(record["id"]): {
                    "gene": str(record["gene"]),
                    "version": str(record["pinned_version"]),
                }
                for record in records
            }
            if set(watched_targets) != set(active_targets) or any(
                watched_targets[source_id]["gene"] != target["gene"]
                or not _versions_equal(
                    watched_targets[source_id]["version"], target["version"]
                )
                for source_id, target in active_targets.items()
            ):
                raise SourceWatchConfigurationError(
                    "CSpec source watch does not match every active gene-policy target"
                )
        if checker == "erepo_export":
            fingerprint = str(source.get("selected_rows_sha256") or "")
            if len(fingerprint) != 64 or set(fingerprint) - set("0123456789abcdef"):
                raise SourceWatchConfigurationError(
                    f"ERepo source {source['id']!r} has no valid selected-row SHA-256"
                )
            if not source.get("genes") or int(source.get("selected_row_count") or 0) < 1:
                raise SourceWatchConfigurationError(
                    f"ERepo source {source['id']!r} has no selected genes or row count"
                )
            metadata = _load_json(_project_path(str(source["metadata_path"])))
            if int(source["selected_row_count"]) != int(metadata.get("record_count") or 0):
                raise SourceWatchConfigurationError(
                    f"ERepo source {source['id']!r} row count differs from active metadata"
                )
    return registry


def _request_json(url: str, *, timeout: float, opener: UrlOpener) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "ARIANE-source-watch/1.0"},
    )
    with opener(request, timeout=timeout) as response:
        payload = json.load(response)
    if not isinstance(payload, dict):
        raise RuntimeError(f"Expected a JSON object from {url}")
    return payload


def _download_bytes(
    url: str,
    *,
    timeout: float,
    opener: UrlOpener,
    maximum_bytes: int,
) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "ARIANE-source-watch/1.0"},
    )
    with opener(request, timeout=timeout) as response:
        content_length = response.headers.get("Content-Length")
        if content_length and int(content_length) > maximum_bytes:
            raise RuntimeError(f"Remote source exceeds {maximum_bytes} bytes: {url}")
        payload = response.read(maximum_bytes + 1)
    if len(payload) > maximum_bytes:
        raise RuntimeError(f"Remote source exceeds {maximum_bytes} bytes: {url}")
    return payload


def _version_parts(value: object) -> tuple[int, ...] | None:
    text = str(value or "").strip()
    parts = text.split(".")
    if not parts or any(not part.isdigit() for part in parts):
        return None
    numbers = tuple(int(part) for part in parts)
    return numbers + (0,) * (3 - len(numbers))


def _versions_equal(first: object, second: object) -> bool:
    first_parts = _version_parts(first)
    second_parts = _version_parts(second)
    if first_parts is not None and second_parts is not None:
        return first_parts == second_parts
    return str(first or "").strip() == str(second or "").strip()


def check_cspec_records(
    source: dict[str, Any],
    timeout: float,
    *,
    opener: UrlOpener = urllib.request.urlopen,
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    changed = False
    version_changed = False
    for expected in source.get("records", []):
        payload = _request_json(expected["api_url"], timeout=timeout, opener=opener)
        candidates = (
            payload.get("data", {}).get("SequenceVariantInterpretation", [])
        )
        matches = [item for item in candidates if item.get("entId") == expected["id"]]
        if len(matches) != 1:
            raise RuntimeError(
                f"CSpec returned {len(matches)} matching records for {expected['id']}"
            )
        record = matches[0]
        content = record.get("entContent") or {}
        current_states = sorted(
            state.get("name")
            for state in content.get("states", [])
            if state.get("current") and state.get("name")
        )
        remote_version = str(content.get("version") or "")
        version_differs = not _versions_equal(remote_version, expected["pinned_version"])
        identity_differs = (
            str(record.get("rev") or "") != expected["observed_revision"]
            or str(record.get("modified") or "") != expected["observed_modified"]
            or current_states != ["Released"]
        )
        changed = changed or version_differs or identity_differs
        version_changed = version_changed or version_differs
        records.append(
            {
                "id": expected["id"],
                "gene": expected["gene"],
                "pinned_version": expected["pinned_version"],
                "remote_version": remote_version,
                "pinned_revision": expected["observed_revision"],
                "remote_revision": record.get("rev"),
                "pinned_modified": expected["observed_modified"],
                "remote_modified": record.get("modified"),
                "approved_on": content.get("approvedOn"),
                "current_states": current_states,
                "title": content.get("title"),
            }
        )
    if not records:
        raise SourceWatchConfigurationError("CSpec checker has no records")
    status = "update_available" if version_changed else "upstream_changed" if changed else "current"
    return {
        "status": status,
        "review_required": changed,
        "records": records,
    }


def check_erepo_export(
    source: dict[str, Any],
    timeout: float,
    *,
    opener: UrlOpener = urllib.request.urlopen,
) -> dict[str, Any]:
    metadata = _load_json(_project_path(source["metadata_path"]))
    export_url = str(metadata["source_export_url"])
    export = _download_bytes(
        export_url,
        timeout=timeout,
        opener=opener,
        maximum_bytes=MAX_EREPO_EXPORT_BYTES,
    )
    remote_sha256 = hashlib.sha256(export).hexdigest()
    pinned_sha256 = str(metadata["source_export_sha256"])
    selected_identity = erepo_selection_identity(
        export,
        expert_panel=str(metadata["expert_panel"]),
        genes=set(source["genes"]),
    )
    pinned_selected_sha256 = str(source["selected_rows_sha256"])
    service = _request_json(
        str(metadata["source_service_url"]), timeout=timeout, opener=opener
    )
    remote_api_version = str(service.get("data", {}).get("version") or "")
    pinned_api_version = str(metadata.get("source_api_version") or "")
    selected_rows_changed = (
        selected_identity["sha256"] != pinned_selected_sha256
        or selected_identity["row_count"] != int(source["selected_row_count"])
    )
    api_changed = not _versions_equal(remote_api_version, pinned_api_version)
    status = (
        "update_available"
        if selected_rows_changed
        else "upstream_changed"
        if api_changed
        else "current"
    )
    return {
        "status": status,
        "review_required": selected_rows_changed or api_changed,
        "source_url": export_url,
        "pinned_export_sha256": pinned_sha256,
        "remote_export_sha256": remote_sha256,
        "full_export_changed": remote_sha256 != pinned_sha256,
        "remote_export_bytes": len(export),
        "pinned_selected_rows_sha256": pinned_selected_sha256,
        "remote_selected_rows_sha256": selected_identity["sha256"],
        "pinned_selected_row_count": source["selected_row_count"],
        "remote_selected_row_count": selected_identity["row_count"],
        "pinned_api_version": pinned_api_version,
        "remote_api_version": remote_api_version,
        "pinned_retrieved_at": metadata.get("retrieved_at"),
        "pinned_record_count": metadata.get("record_count"),
    }


def erepo_selection_identity(
    export: bytes,
    *,
    expert_panel: str,
    genes: set[str],
) -> dict[str, Any]:
    """Fingerprint only the ERepo rows in ARIANE scope.

    The official download contains every ClinGen panel.  Hashing the complete
    file alone would alert whenever an unrelated panel changed.  This stable
    projection still records the full-export checksum for provenance but uses
    only ENIGMA BRCA1/2 rows to decide whether ARIANE review is required.
    """
    text = export.decode("utf-8-sig")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    selected = [
        {str(key): str(value or "").strip() for key, value in row.items()}
        for row in rows
        if str(row.get("Expert Panel") or "").strip() == expert_panel
        and str(row.get("HGNC Gene Symbol") or "").strip() in genes
    ]
    if not selected:
        raise RuntimeError("ERepo export contains no selected ENIGMA BRCA1/2 rows")
    selected.sort(key=lambda row: (row.get("Uuid", ""), row.get("Variation", "")))
    canonical = json.dumps(
        selected,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return {"row_count": len(selected), "sha256": hashlib.sha256(canonical).hexdigest()}


def _load_script_module(name: str, filename: str) -> ModuleType:
    path = ROOT / "scripts" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load checker module {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_clinical_lr(source: dict[str, Any], timeout: float) -> dict[str, Any]:
    module = _load_script_module(
        "ariane_clinical_lr_update_check", "check_enigma_clinical_lr_update.py"
    )
    result = module.check_update(
        manifest_path=_project_path(source["manifest_path"]),
        timeout=timeout,
    )
    result["review_required"] = result.get("status") == "update_available"
    return result


def check_gnomad_releases(source: dict[str, Any], timeout: float) -> dict[str, Any]:
    del timeout  # The specialised checker owns safe per-request timeouts.
    module = _load_script_module(
        "ariane_gnomad_update_check", "refresh_gnomad_panel_snapshot.py"
    )
    specialised = module.check_updates(_project_path(source["manifest_path"]))
    review_required = specialised.get("activation_status") == "review_required"
    return {
        "status": "update_available" if review_required else "current",
        "review_required": review_required,
        "active_releases": specialised.get("active_releases", []),
        "available_releases": specialised.get("available_releases", []),
        "detected_other_releases": specialised.get("detected_other_releases", []),
        "newer_equivalent_product_releases": specialised.get(
            "newer_equivalent_product_releases", {}
        ),
        "equivalent_product_checks": specialised.get("equivalent_product_checks", []),
        "source_identities": specialised.get("source_identities", {}),
        "automatic_activation": specialised.get("automatic_activation"),
    }


def check_manual_source(source: dict[str, Any], timeout: float) -> dict[str, Any]:
    del timeout
    due_on = source.get("next_review_due")
    overdue = bool(due_on and date.fromisoformat(str(due_on)) < _utc_now().date())
    return {
        "status": "manual_review_due" if overdue else "manual_only",
        "review_required": overdue,
        "next_review_due": due_on,
        "review_trigger": source.get("review_trigger"),
        "tracked_files": source.get("tracked_files", []),
        "note": source.get("note"),
    }


DEFAULT_CHECKERS: dict[str, SourceChecker] = {
    "cspec_records": check_cspec_records,
    "erepo_export": check_erepo_export,
    "clinical_lr": check_clinical_lr,
    "gnomad_releases": check_gnomad_releases,
    "manual": check_manual_source,
}


def run_checks(
    registry: dict[str, Any],
    *,
    selected_sources: Iterable[str] | None = None,
    timeout: float | None = None,
    checkers: dict[str, SourceChecker] | None = None,
) -> dict[str, Any]:
    sources = registry["sources"]
    selected = set(selected_sources or [])
    known = {source["id"] for source in sources}
    unknown = selected - known
    if unknown:
        raise SourceWatchConfigurationError(
            "Unknown source IDs: " + ", ".join(sorted(unknown))
        )
    effective_timeout = float(
        timeout or registry["policy"].get("default_timeout_seconds", 60)
    )
    handlers = checkers or DEFAULT_CHECKERS
    results: list[dict[str, Any]] = []
    for source in sources:
        if selected and source["id"] not in selected:
            continue
        checked_at = _utc_now().isoformat()
        base = {
            "id": source["id"],
            "name": source["name"],
            "checker": source["checker"],
            "cadence_days": source.get("cadence_days"),
            "official_url": source.get("official_url"),
            "checked_utc": checked_at,
        }
        try:
            detail = handlers[source["checker"]](source, effective_timeout)
            base.update(detail)
        except Exception as exc:
            base.update(
                {
                    "status": "check_failed",
                    "review_required": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )
        results.append(base)

    statuses = {result["status"] for result in results}
    if "check_failed" in statuses:
        overall = "check_failed"
        exit_code = 3
    elif statuses & REVIEW_STATUSES:
        overall = "review_required"
        exit_code = 2
    else:
        overall = "current"
        exit_code = 0
    return {
        "schema_version": 1,
        "generated_utc": _utc_now().isoformat(),
        "overall_status": overall,
        "exit_code": exit_code,
        "automatic_activation": False,
        "runtime_network_dependency": False,
        "summary": {
            "checked": len(results),
            "current": sum(result["status"] == "current" for result in results),
            "manual_only": sum(
                result["status"] == "manual_only" for result in results
            ),
            "review_required": sum(
                result["status"] in REVIEW_STATUSES for result in results
            ),
            "failed": sum(result["status"] == "check_failed" for result in results),
        },
        "sources": results,
        "activation_policy": registry["policy"]["required_before_activation"],
    }


def render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# ARIANE source-watch report",
        "",
        f"Overall status: **{report['overall_status']}**  ",
        f"Generated: `{report['generated_utc']}`",
        "",
        "| Source | Checker | Status | Action |",
        "|---|---|---|---|",
    ]
    for source in report["sources"]:
        action = "review required" if source.get("review_required") else "none"
        if source["status"] == "manual_only":
            action = f"manual trigger: {source.get('review_trigger', '')}"
        if source["status"] == "check_failed":
            action = f"check failed: {source.get('error', '')}"
        action = str(action).replace("|", "\\|").replace("\n", " ")
        lines.append(
            f"| {source['name']} | `{source['checker']}` | "
            f"`{source['status']}` | {action} |"
        )
    lines.extend(
        [
            "",
            "The check is read-only. No production dataset or manifest was changed.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument(
        "--source",
        action="append",
        dest="sources",
        help="Check only this source ID; may be repeated.",
    )
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--output", type=Path, help="Write the JSON report here.")
    parser.add_argument(
        "--markdown-output", type=Path, help="Write a short Markdown summary here."
    )
    parser.add_argument(
        "--validate-registry",
        action="store_true",
        help="Validate local configuration without contacting upstream services.",
    )
    args = parser.parse_args(argv)

    try:
        registry = load_registry(args.registry)
        if args.validate_registry:
            print(f"Source-watch registry is valid: {args.registry}")
            return 0
        report = run_checks(
            registry,
            selected_sources=args.sources,
            timeout=args.timeout,
        )
    except Exception as exc:
        print(f"Source-watch configuration failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 4

    rendered = json.dumps(report, ensure_ascii=False, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(render_markdown(report), encoding="utf-8")
    print(rendered)
    return int(report["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
