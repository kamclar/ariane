#!/usr/bin/env python3
"""Extend a completed distance scan with the original and current UI defaults.

This is a separate, resumable experiment over the exact saved cohort. It scores
D=50, D=500 and the shared D=4999 anchor without changing the completed high-
distance result. The D=4999 summary for every case must reproduce the base run.
"""

from __future__ import annotations

import argparse
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import signal
import sqlite3
import sys
import time

from scripts import compare_spliceai_distances as core
from scripts import spliceai_distance_scan as scan
from scripts.spliceai_distance_store import ResultStore


DISTANCES = (50, 500, 4999)
PAIRS = ((50, 500), (500, 4999), (50, 4999))


def differing_paths(expected, observed, tolerance, path=""):
    """Compare persisted numeric summaries without requiring bitwise equality."""
    if isinstance(expected, dict) and isinstance(observed, dict):
        differences = []
        for key in sorted(set(expected) | set(observed)):
            child = f"{path}.{key}" if path else key
            if key not in expected or key not in observed:
                differences.append(child)
            else:
                differences.extend(differing_paths(expected[key], observed[key], tolerance, child))
        return differences
    if isinstance(expected, float) or isinstance(observed, float):
        try:
            return [] if abs(float(expected) - float(observed)) <= tolerance else [path]
        except (TypeError, ValueError):
            return [path]
    return [] if expected == observed else [path]


def baseline_projection(run):
    """Keep reproducibility-critical scores, not unstable argmax tie details."""
    return {
        "distance": run["distance"],
        "maximum_delta": run["maximum_delta"],
        "band_raw": run["band_raw"],
        "band_rounded_3dp": run["band_rounded_3dp"],
        "component_deltas": {
            name: component["delta"] for name, component in run["components"].items()
        },
        "edge_deltas": run["edge_deltas"],
        "outside_model_radius": run["outside_model_radius"],
    }


def baseline_summary(connection, case_id):
    row = connection.execute(
        "SELECT status, result_json FROM cases WHERE case_id = ?", (case_id,)
    ).fetchone()
    if row is None or row[0] != "ok":
        raise ValueError(f"Base experiment lacks a successful result for {case_id}")
    result = json.loads(row[1])
    candidates = [run for run in result.get("runs", [])
                  if run.get("distance") == 4999 and run.get("repeat") == 1]
    if len(candidates) != 1:
        raise ValueError(f"Base experiment lacks one D=4999 baseline for {case_id}")
    return candidates[0]


def flags_for(comparisons, runs, kind, tolerance):
    return {
        "delta_changed": any(not item["all_delta_scores_within_tolerance"]
                             for item in comparisons),
        "rounded_delta_changed": any(not item["all_rounded_delta_scores_equal"]
                                     for item in comparisons),
        "position_changed": any(
            channel["position_changed"] and not channel["both_deltas_negligible"]
            for item in comparisons for channel in item["channels"].values()
        ),
        "band_changed": any(not item["raw_band_equal"] or not item["rounded_band_equal"]
                            for item in comparisons),
        "overlap_changed": any(
            max(item["overlapping_ref_max_abs_difference"],
                item["overlapping_alt_max_abs_difference"]) > tolerance
            for item in comparisons
        ),
        "snv_outside_radius": kind == "snv" and any(
            (run["outside_model_radius"]["maximum_absolute_delta"] or 0) > tolerance
            for run in runs
        ),
    }


def evaluate_case(annotator, base_db, case, index, args):
    import numpy as np

    repeats = args.repeats if index < args.timing_repeat_cases else 1
    first, runs, comparisons = {}, [], []
    baseline_mismatches = []
    expected_baseline = baseline_summary(base_db, case["id"])
    baseline_position_changes = []
    for repeat in range(repeats):
        shift = index % len(DISTANCES)
        order = list(DISTANCES[shift:] + DISTANCES[:shift])
        if repeat % 2:
            order.reverse()
        for distance in order:
            arrays, row, elapsed = core.infer(annotator, case, distance)
            summary = core.summarize_arrays(*arrays, distance, args.tolerance)
            core.validate_capture(summary, row)
            run = {"distance": distance, "repeat": repeat + 1, "seconds": elapsed, **summary}
            runs.append(run)
            if distance in first:
                comparisons.append(core.compare_arrays(
                    first[distance], arrays, distance, distance, args.tolerance
                ))
            else:
                first[distance] = arrays
            if distance == 4999 and repeat == 0:
                baseline_mismatches = differing_paths(
                    baseline_projection(expected_baseline), baseline_projection(run),
                    args.tolerance,
                )
                baseline_position_changes = sorted(
                    name for name in core.CHANNELS
                    if expected_baseline["components"][name]["position"]
                    != run["components"][name]["position"]
                )
    for smaller, larger in PAIRS:
        comparisons.append(core.compare_arrays(
            first[smaller], first[larger], smaller, larger, args.tolerance
        ))
    flags = flags_for(comparisons, runs, case["variant_kind"], args.tolerance)
    artifacts = []
    retain = args.save_arrays == "all" or (
        args.save_arrays == "changed" and (any(flags.values()) or baseline_mismatches)
    )
    if retain:
        directory = args.output / "arrays"
        directory.mkdir(exist_ok=True)
        for distance, arrays in first.items():
            path = directory / f"{case['id']}-d{distance}.npz"
            np.savez_compressed(path, ref=arrays[0], alt=arrays[1])
            artifacts.append({"path": str(path.relative_to(args.output)),
                              "sha256": core.sha256(path)})
    return {
        "id": case["id"],
        "status": "error" if baseline_mismatches else "ok",
        "case": case,
        "error": "d4999_baseline_mismatch" if baseline_mismatches else None,
        "baseline_mismatch_fields": baseline_mismatches,
        "baseline_argmax_position_changes": baseline_position_changes,
        "runs": runs,
        "comparisons": comparisons,
        "flags": flags,
        "raw_arrays": artifacts,
    }


def publish(output, summary, base_summary):
    summary["base_high_distance_result"] = {
        "status": base_summary["status"],
        "counts": base_summary["counts"],
        "flags": base_summary["flags"],
    }
    scan.write_json(output / "summary.json", summary)
    counts = summary["counts"]
    lines = [
        "# SpliceAI default-distance extension", "",
        f"Status: {summary['status']}", "",
        f"Planned: {counts['planned']}; successful: {counts['ok']}; "
        f"errors: {counts['error']}; remaining: {counts['remaining']}.", "",
        "Distances: D=50, D=500 and shared anchor D=4999.", "",
        summary["statistical_scope"], "", "## Observed differences", "",
    ]
    lines.extend(f"- {key}: {value}" for key, value in summary["flags"].items())
    lines.extend([
        "", "Every successful case reproduced its saved D=4999 summary from the base run.",
        "Complete per-case data and failures: `results.sqlite3`.",
        "Coverage and group counts: `summary.json`. Selected raw arrays: `arrays/`.", "",
    ])
    temporary = output / "summary.md.tmp"
    temporary.write_text("\n".join(lines), encoding="utf-8")
    os.replace(temporary, output / "summary.md")


def validate_args(args):
    for name in ("threads", "repeats", "progress_every", "max_errors"):
        if getattr(args, name) < 1:
            raise ValueError(f"{name} must be positive")
    if args.timing_repeat_cases < 0:
        raise ValueError("timing_repeat_cases must be nonnegative")
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        raise ValueError("tolerance must be positive and finite")
    if args.max_cases is not None and args.max_cases < 1:
        raise ValueError("max_cases must be positive")
    if args.max_seconds is not None and (
        not math.isfinite(args.max_seconds) or args.max_seconds <= 0
    ):
        raise ValueError("max_seconds must be positive and finite")
    if args.retry_errors and not args.resume:
        raise ValueError("Use --resume together with --retry-errors")


def prepare(args):
    base_plan_path = args.base_output / "plan.json"
    base_summary_path = args.base_output / "summary.json"
    base_database = args.base_output / "results.sqlite3"
    if not all(path.is_file() for path in (base_plan_path, base_summary_path, base_database)):
        raise ValueError("Completed base experiment files are missing")
    base_plan = json.loads(base_plan_path.read_text(encoding="utf-8"))
    base_summary = json.loads(base_summary_path.read_text(encoding="utf-8"))
    cases = base_plan["cohort"]["cases"]
    if (base_summary.get("status") != "completed"
            or base_summary.get("counts", {}).get("ok") != len(cases)
            or base_summary.get("counts", {}).get("error") != 0):
        raise ValueError("Base experiment is not complete and error-free")
    with sqlite3.connect(f"file:{base_database.resolve().as_posix()}?mode=ro", uri=True) as db:
        rows = dict(db.execute("SELECT status, COUNT(*) FROM cases GROUP BY status"))
    if rows != {"ok": len(cases)}:
        raise ValueError("Base checkpoint does not match its completed summary")

    profile_path = core.ROOT / "data/spliceai/enigma_v1_2_spliceai_profile.json"
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    engine = profile["approved_engine"]
    annotation = core.ROOT / engine["reference_transcript_annotation"]
    if (args.image_digest != engine["docker_image"]
            or os.environ.get("MODEL_COMMIT") != engine["model_commit"]
            or version("spliceai") != engine["python_package_version"]):
        raise ValueError("Container image/package/commit does not match the pinned profile")
    if core.sha256(annotation) != engine["reference_transcript_annotation_sha256"]:
        raise ValueError("Annotation checksum does not match the profile")
    base_identity = base_plan["identity"]
    if (base_identity["image_digest"] != args.image_digest
            or base_identity["fasta_sha256"] != core.sha256(args.fasta)
            or base_identity["annotation_sha256"] != core.sha256(annotation)):
        raise ValueError("Base experiment used different model/reference inputs")
    for name, expected in base_identity["sources_sha256"].items():
        if core.sha256(core.ROOT / "scripts" / name) != expected:
            raise ValueError(f"Base experiment source changed: {name}")

    import_source = None
    if args.import_output is not None:
        import_plan_path = args.import_output / "plan.json"
        import_database = args.import_output / "results.sqlite3"
        if not import_plan_path.is_file() or not import_database.is_file():
            raise ValueError("Import checkpoint plan/database is missing")
        import_plan = json.loads(import_plan_path.read_text(encoding="utf-8"))
        if (import_plan.get("identity", {}).get("base_plan_sha256")
                != core.sha256(base_plan_path)
                or import_plan.get("identity", {}).get("settings", {}).get("distances")
                != list(DISTANCES)):
            raise ValueError("Import checkpoint belongs to a different cohort or distance set")
        import_source = {
            "directory": args.import_output.name,
            "plan_sha256": core.sha256(import_plan_path),
            "results_sha256": core.sha256(import_database),
            "policy": "copy successful cases only; recompute failures and unfinished cases",
        }

    identity = {
        "schema_version": 1,
        "experiment": "spliceai-default-distance-extension",
        "base_plan_sha256": core.sha256(base_plan_path),
        "base_results_sha256": core.sha256(base_database),
        "base_identity_sha256": scan.canonical_hash(base_identity),
        "extension_source_sha256": core.sha256(Path(__file__)),
        "import_source": import_source,
        "image_digest": args.image_digest,
        "model_commit": os.environ["MODEL_COMMIT"],
        "settings": {
            "distances": list(DISTANCES), "mask": 0, "threads": args.threads,
            "tolerance": args.tolerance, "repeats": args.repeats,
            "timing_repeat_cases": args.timing_repeat_cases,
            "save_arrays": args.save_arrays, "device": "CPU",
        },
    }
    plan = {
        "identity": identity,
        "cohort": base_plan["cohort"],
        "base_experiment": {
            "directory": args.base_output.name,
            "completed_at": base_summary.get("updated_at"),
            "distances": [4999, 5000, 10000],
        },
        "import_source": import_source,
    }
    plan_path = args.output / "plan.json"
    if plan_path.exists():
        if json.loads(plan_path.read_text(encoding="utf-8")) != plan:
            raise ValueError("Plan/configuration/source changed; use a new output directory")
    elif args.resume:
        raise ValueError("Cannot resume without a saved extension plan")
    else:
        scan.write_json(plan_path, plan)
    return plan, annotation, base_summary, base_database


def execute(args, plan, annotation, base_summary, base_database):
    database = args.output / "results.sqlite3"
    base_uri = f"file:{base_database.resolve().as_posix()}?mode=ro"
    with ResultStore(database, plan["identity"], resume=args.resume) as store, \
            sqlite3.connect(base_uri, uri=True) as base_db:
        if not args.resume and args.import_output is not None:
            import_database = args.import_output / "results.sqlite3"
            import_uri = f"file:{import_database.resolve().as_posix()}?mode=ro"
            imported = 0
            with sqlite3.connect(import_uri, uri=True) as source:
                for case_id, payload in source.execute(
                    "SELECT case_id, result_json FROM cases WHERE status = 'ok' ORDER BY rowid"
                ):
                    store.put(case_id, json.loads(payload))
                    imported += 1
            print(json.dumps({"phase": "import", "successful_cases": imported}), flush=True)
        cases = plan["cohort"]["cases"]
        pending = []
        for index, case in enumerate(cases):
            existing = store.get(case["id"])
            if existing is None or (args.retry_errors and existing["status"] == "error"):
                pending.append((index, case))
        if not pending:
            state = "completed_with_errors" if store.counts()["error"] else "completed"
            summary = scan.summarize(plan, store.results(), state)
            publish(args.output, summary, base_summary)
            print(json.dumps({"status": state, **summary["counts"]}), flush=True)
            return 1 if state == "completed_with_errors" else 0

        import tensorflow as tf
        from spliceai.utils import Annotator

        tf.config.set_visible_devices([], "GPU")
        tf.config.threading.set_intra_op_parallelism_threads(args.threads)
        tf.config.threading.set_inter_op_parallelism_threads(1)
        annotator = Annotator(str(args.fasta), str(annotation))
        annotator.models = [core.RecordingModel(model) for model in annotator.models]
        for distance in DISTANCES:
            arrays, row, _ = core.infer(annotator, cases[0], distance)
            core.validate_capture(
                core.summarize_arrays(*arrays, distance, args.tolerance), row
            )

        stopped = False

        def request_stop(_signum, _frame):
            nonlocal stopped
            stopped = True

        previous = {item: signal.signal(item, request_stop)
                    for item in (signal.SIGINT, signal.SIGTERM)}
        started = time.monotonic()
        processed = errors = 0
        state = "running"
        try:
            publish(args.output, scan.summarize(plan, store.results(), state), base_summary)
            for index, case in pending:
                if (stopped
                        or (args.max_cases is not None and processed >= args.max_cases)
                        or (args.max_seconds is not None
                            and time.monotonic() - started >= args.max_seconds)):
                    state = "paused"
                    break
                try:
                    result = evaluate_case(annotator, base_db, case, index, args)
                except Exception as exc:
                    result = {"id": case["id"], "status": "error", "case": case,
                              "error": f"{type(exc).__name__}: {exc}"}
                store.put(case["id"], result)
                processed += 1
                errors += result["status"] == "error"
                if processed % args.progress_every == 0 or result["status"] == "error":
                    summary = scan.summarize(plan, store.results(), state)
                    publish(args.output, summary, base_summary)
                    print(json.dumps({"status": state, **summary["counts"],
                                      "flags": summary["flags"]}), flush=True)
                if result.get("baseline_mismatch_fields") or errors >= args.max_errors:
                    state = "failed"
                    break
            else:
                state = "completed_with_errors" if store.counts()["error"] else "completed"
        except BaseException:
            state = "failed"
            raise
        finally:
            for item, handler in previous.items():
                signal.signal(item, handler)
            summary = scan.summarize(plan, store.results(), state)
            publish(args.output, summary, base_summary)
        print(json.dumps({"status": state, **summary["counts"],
                          "flags": summary["flags"]}), flush=True)
        return 1 if state in ("failed", "completed_with_errors") else 0


def run(args):
    try:
        validate_args(args)
        if args.resume and not args.output.exists():
            raise ValueError("Resume directory does not exist")
        if args.output.exists() and not (args.output / "plan.json").exists():
            raise ValueError("Output must be new or an existing extension plan")
        args.output.mkdir(parents=True, exist_ok=True)
        with scan.exclusive_run(args.output):
            plan, annotation, base_summary, base_database = prepare(args)
            evaluations = len(plan["cohort"]["cases"]) * len(DISTANCES)
            evaluations += min(args.timing_repeat_cases, len(plan["cohort"]["cases"])) \
                * len(DISTANCES) * (args.repeats - 1)
            print(json.dumps({
                "phase": "plan", "cases": len(plan["cohort"]["cases"]),
                "distances": list(DISTANCES), "model_evaluations": evaluations,
                "rough_hours_at_1.5_seconds_per_evaluation": round(evaluations * 1.5 / 3600, 2),
            }), flush=True)
            if args.prepare_only:
                return 0
            return execute(args, plan, annotation, base_summary, base_database)
    except Exception as exc:
        print(f"Default-distance scan failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--base-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--import-output", type=Path)
    parser.add_argument("--image-digest", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--timing-repeat-cases", type=int, default=12)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--max-errors", type=int, default=5)
    parser.add_argument("--max-cases", type=int)
    parser.add_argument("--max-seconds", type=float)
    parser.add_argument("--save-arrays", choices=("all", "changed", "none"), default="changed")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--retry-errors", action="store_true")
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
