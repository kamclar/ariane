"""Resumable, bounded cohort execution for the local distance experiment."""

from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import statistics
import sys
import time
from types import SimpleNamespace

from scripts import compare_spliceai_distances as core
from scripts.spliceai_distance_cases import build_cohort
from scripts.spliceai_distance_store import ResultStore


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8")
    os.replace(temporary, path)


@contextmanager
def exclusive_run(output):
    # Execution is inside the Linux container. OS locks are released even if
    # the process dies; a stale lock file cannot prevent a later resume.
    import fcntl

    with (output / ".run.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("Another process is using this experiment directory") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def annotation_view(annotation, fasta_path):
    from pyfastx import Fasta

    with annotation.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return SimpleNamespace(
        genes=[r["#NAME"] for r in rows], chroms=[r["CHROM"] for r in rows],
        strands=[r["STRAND"] for r in rows],
        tx_starts=[int(r["TX_START"]) + 1 for r in rows],
        tx_ends=[int(r["TX_END"]) for r in rows],
        exon_starts=[[int(v) + 1 for v in r["EXON_START"].split(",") if v] for r in rows],
        exon_ends=[[int(v) for v in r["EXON_END"].split(",") if v] for r in rows],
        ref_fasta=Fasta(str(fasta_path)),
    )


def validate_options(args):
    for name in ("positions_per_region", "repeats", "threads", "progress_every", "max_errors"):
        if getattr(args, name) < 1:
            raise ValueError(f"{name} must be positive")
    if args.timing_repeat_cases < 0:
        raise ValueError("timing_repeat_cases must be nonnegative")
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        raise ValueError("Tolerance must be positive and finite")
    if args.max_cases is not None and args.max_cases < 1:
        raise ValueError("max_cases must be positive")
    if args.max_seconds is not None and (not math.isfinite(args.max_seconds) or args.max_seconds <= 0):
        raise ValueError("max_seconds must be positive and finite")
    if not {4999, 5000, 10000}.issubset(args.distances) or any(
        d < 1 or d > 10000 for d in args.distances
    ):
        raise ValueError("Distances must include 4999, 5000 and 10000, within 1..10000")
    if args.retry_errors and not args.resume:
        raise ValueError("Use --resume together with --retry-errors")


def prepare(args):
    profile_path = core.ROOT / "data/spliceai/enigma_v1_2_spliceai_profile.json"
    fixture_path = core.ROOT / "data/spliceai/local_service_validation_cases.json"
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
    engine = profile["approved_engine"]
    annotation = core.ROOT / engine["reference_transcript_annotation"]
    if (args.image_digest != engine["docker_image"]
            or os.environ.get("MODEL_COMMIT") != engine["model_commit"]
            or version("spliceai") != engine["python_package_version"]):
        raise ValueError("Container image/package/commit does not match the pinned profile")
    if core.sha256(annotation) != engine["reference_transcript_annotation_sha256"]:
        raise ValueError("Annotation checksum does not match the profile")
    if fixture["profile_id"] != profile["profile_id"]:
        raise ValueError("Control fixture/profile mismatch")
    view = annotation_view(annotation, args.fasta)
    cohort = build_cohort(view, profile, fixture, seed=args.seed,
                          positions_per_region=args.positions_per_region,
                          include_indels=not args.no_indels)
    spec = importlib.util.find_spec("spliceai")
    model_dir = Path(spec.origin).parent / "models"
    identity = {
        "schema_version": 2,
        "cohort_sha256": canonical_hash(cohort), "profile_sha256": core.sha256(profile_path),
        "fixture_sha256": core.sha256(fixture_path),
        "annotation_sha256": core.sha256(annotation),
        "fasta_sha256": core.sha256(args.fasta),
        "image_digest": args.image_digest, "model_commit": os.environ["MODEL_COMMIT"],
        "model_weights_sha256": {f"spliceai{i}.h5": core.sha256(model_dir / f"spliceai{i}.h5")
                                 for i in range(1, 6)},
        "sources_sha256": {name: core.sha256(core.ROOT / "scripts" / name) for name in (
            "compare_spliceai_distances.py", "spliceai_distance_scan.py",
            "spliceai_distance_cases.py", "spliceai_distance_store.py")},
        "packages": {name: version(name) for name in ("spliceai", "tensorflow-cpu", "numpy")},
        "settings": {"distances": sorted(set(args.distances)), "mask": 0,
                     "tolerance": args.tolerance, "threads": args.threads,
                     "repeats": args.repeats, "timing_repeat_cases": args.timing_repeat_cases,
                     "save_arrays": args.save_arrays, "device": "CPU",
                     "small_windows": "fixture controls only"},
    }
    plan = {"identity": identity, "cohort": cohort}
    path = args.output / "plan.json"
    if path.exists():
        if json.loads(path.read_text(encoding="utf-8")) != plan:
            raise ValueError("Plan/configuration/source changed; use a new experiment directory")
    else:
        if args.resume:
            raise ValueError("Cannot resume without a saved plan")
        write_json(path, plan)
    return plan, annotation


def case_distances(case, distances):
    # D50 proves the comparison can detect a real window effect, but need not
    # be recomputed for thousands of other variants.
    return sorted({d for d in distances if d >= 4999 or "delta_scores" in case})


def comparison_flags(comparisons, runs, kind, tolerance):
    relevant = [c for c in comparisons if c["smaller_distance"] >= 4999]
    return {
        "delta_changed": any(not c["all_delta_scores_within_tolerance"] for c in relevant),
        "rounded_delta_changed": any(not c["all_rounded_delta_scores_equal"] for c in relevant),
        "position_changed": any(v["position_changed"] for c in relevant for v in c["channels"].values()),
        "band_changed": any(not c["raw_band_equal"] or not c["rounded_band_equal"] for c in relevant),
        "overlap_changed": any(max(c["overlapping_ref_max_abs_difference"],
                                    c["overlapping_alt_max_abs_difference"]) > tolerance
                               for c in relevant),
        "snv_outside_radius": kind == "snv" and any(
            (r["outside_model_radius"]["maximum_absolute_delta"] or 0) > tolerance for r in runs
        ),
    }


def evaluate_case(annotator, case, index, args):
    import numpy as np

    distances = case_distances(case, args.distances)
    repeats = args.repeats if index < args.timing_repeat_cases else 1
    first = {}
    runs, comparisons = [], []
    baseline_errors = []
    for repeat in range(repeats):
        shift = index % len(distances)
        order = distances[shift:] + distances[:shift]
        if repeat % 2:
            order.reverse()
        for distance in order:
            arrays, row, elapsed = core.infer(annotator, case, distance)
            summary = core.summarize_arrays(*arrays, distance, args.tolerance)
            core.validate_capture(summary, row)
            if distance == 10000 and "delta_scores" in case:
                for field, expected in {**case["delta_scores"], **case["reference_scores"],
                                        **case["alternate_scores"]}.items():
                    if float(row.get(field, "nan")) != expected:
                        baseline_errors.append(field)
            runs.append({"distance": distance, "repeat": repeat + 1, "seconds": elapsed,
                         **summary})
            if distance in first:
                comparisons.append(core.compare_arrays(first[distance], arrays, distance,
                                                        distance, args.tolerance))
            else:
                first[distance] = arrays
    for small, large in ((4999, 5000), (4999, 10000), (5000, 10000)):
        comparisons.append(core.compare_arrays(first[small], first[large], small, large,
                                               args.tolerance))
    if 50 in first:
        comparisons.append(core.compare_arrays(first[50], first[10000], 50, 10000, args.tolerance))
    flags = comparison_flags(comparisons, runs, case["variant_kind"], args.tolerance)
    artifacts = []
    retain = args.save_arrays == "all" or (
        args.save_arrays == "changed" and (any(flags.values()) or baseline_errors)
    )
    if retain:
        directory = args.output / "arrays"
        directory.mkdir(exist_ok=True)
        for distance, arrays in first.items():
            path = directory / f"{case['id']}-d{distance}.npz"
            np.savez_compressed(path, ref=arrays[0], alt=arrays[1])
            artifacts.append({"path": str(path.relative_to(args.output)), "sha256": core.sha256(path)})
    return {"id": case["id"], "status": "error" if baseline_errors else "ok", "case": case,
            "error": "baseline_mismatch" if baseline_errors else None,
            "baseline_mismatch_fields": baseline_errors, "runs": runs,
            "comparisons": comparisons, "flags": flags, "raw_arrays": artifacts}


def summarize(plan, rows, state):
    cases = plan["cohort"]["cases"]
    groups = defaultdict(lambda: {"planned": 0, "ok": 0, "error": 0,
                                 "delta_changed": 0, "band_changed": 0})
    planned_loci, finished_loci = set(), set()
    by_id = {c["id"]: c for c in cases}
    for case in cases:
        planned_loci.add(case["locus_id"])
        for group in sorted({case["gene"], case["variant_kind"], *case["strata"]}):
            groups[group]["planned"] += 1
    counts = {"planned": len(cases), "ok": 0, "error": 0}
    flags = defaultdict(int)
    elapsed = defaultdict(list)
    anomalies = []
    for case_id, result in rows:
        case = by_id[case_id]
        status = result["status"]
        counts[status] += 1
        if status == "ok":
            finished_loci.add(case["locus_id"])
        for key, value in result.get("flags", {}).items():
            flags[key] += bool(value)
        for group in sorted({case["gene"], case["variant_kind"], *case["strata"]}):
            groups[group][status] += 1
            if status == "ok":
                for key in ("delta_changed", "band_changed"):
                    groups[group][key] += bool(result.get("flags", {}).get(key))
        for run in result.get("runs", []):
            elapsed[str(run["distance"])].append(run["seconds"])
        if status == "error" or any(result.get("flags", {}).values()):
            if len(anomalies) < 20:
                anomalies.append({"id": case_id, "variant": case["variant"],
                                  "status": status, "error": result.get("error"),
                                  "flags": result.get("flags", {})})
    counts["remaining"] = counts["planned"] - counts["ok"] - counts["error"]
    return {"status": state, "updated_at": datetime.now(timezone.utc).isoformat(),
            "counts": counts, "flags": dict(flags), "groups": dict(groups),
            "distinct_loci": {"planned": len(planned_loci), "with_successful_results": len(finished_loci)},
            "coverage": plan["cohort"]["coverage"],
            "statistical_scope": "Stratified engineering validation, with targeted controls. "
            "Three ALT alleles at one position are correlated; repeats are not new samples. "
            "Counts/fractions are empirical, not population-wide confidence limits. "
            "Strata overlap. Errors are shown separately and are not agreements.",
            "timing_seconds": {d: {"n": len(values), "median": statistics.median(values),
                                    "min": min(values), "max": max(values)}
                               for d, values in elapsed.items()},
            "first_anomalies": anomalies}


def publish_summary(output, summary):
    write_json(output / "summary.json", summary)
    counts = summary["counts"]
    lines = ["# SpliceAI distance scan", "", f"Status: {summary['status']}", "",
             f"Planned: {counts['planned']}; successful: {counts['ok']}; "
             f"errors: {counts['error']}; remaining: {counts['remaining']}.", "",
             summary["statistical_scope"], "", "## Differences among D=4999/5000/10000", ""]
    lines += [f"- {key}: {value}" for key, value in summary["flags"].items()]
    lines += ["", "Complete per-case data and failures: `results.sqlite3`.",
              "Coverage and group counts: `summary.json`. Selected raw arrays: `arrays/`.", ""]
    path = output / "summary.md"
    temporary = path.with_suffix(".md.tmp")
    temporary.write_text("\n".join(lines), encoding="utf-8")
    os.replace(temporary, path)


def execute(args, plan, annotation):
    database = args.output / "results.sqlite3"
    with ResultStore(database, plan["identity"], resume=args.resume) as store:
        cases = plan["cohort"]["cases"]
        pending = [(i, c) for i, c in enumerate(cases)
                   if store.get(c["id"]) is None or
                   (args.retry_errors and store.get(c["id"])["status"] == "error")]
        if not pending:
            state = "completed_with_errors" if store.counts()["error"] else "completed"
            summary = summarize(plan, store.results(), state)
            publish_summary(args.output, summary)
            print(json.dumps({"status": state, **summary["counts"]}), flush=True)
            return 1 if state == "completed_with_errors" else 0
        import tensorflow as tf
        from spliceai.utils import Annotator

        tf.config.set_visible_devices([], "GPU")
        tf.config.threading.set_intra_op_parallelism_threads(args.threads)
        tf.config.threading.set_inter_op_parallelism_threads(1)
        annotator = Annotator(str(args.fasta), str(annotation))
        annotator.models = [core.RecordingModel(m) for m in annotator.models]
        for distance in sorted(set(args.distances)):
            arrays, row, _ = core.infer(annotator, cases[0], distance)
            core.validate_capture(core.summarize_arrays(*arrays, distance, args.tolerance), row)
        stopped = False

        def request_stop(_signum, _frame):
            nonlocal stopped
            stopped = True

        previous = {s: signal.signal(s, request_stop) for s in (signal.SIGINT, signal.SIGTERM)}
        start = time.monotonic()
        processed = errors = 0
        state = "running"
        try:
            publish_summary(args.output, summarize(plan, store.results(), state))
            for index, case in pending:
                if stopped or (args.max_cases is not None and processed >= args.max_cases) or (
                    args.max_seconds is not None and time.monotonic() - start >= args.max_seconds
                ):
                    state = "paused"
                    break
                try:
                    result = evaluate_case(annotator, case, index, args)
                except Exception as exc:
                    result = {"id": case["id"], "status": "error", "case": case,
                              "error": f"{type(exc).__name__}: {exc}"}
                store.put(case["id"], result)
                processed += 1
                errors += result["status"] == "error"
                if processed % args.progress_every == 0 or result["status"] == "error":
                    summary = summarize(plan, store.results(), state)
                    publish_summary(args.output, summary)
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
            for sig, handler in previous.items():
                signal.signal(sig, handler)
            summary = summarize(plan, store.results(), state)
            publish_summary(args.output, summary)
        print(json.dumps({"status": state, **summary["counts"], "flags": summary["flags"]}), flush=True)
        return 1 if state in ("failed", "completed_with_errors") else 0


def run_scan(args):
    try:
        validate_options(args)
        existing = args.output.exists()
        if args.resume and not existing:
            raise ValueError("Resume directory does not exist")
        if existing and not (args.output / "plan.json").exists():
            raise ValueError("Output must be a new directory or an existing scan plan")
        args.output.mkdir(parents=True, exist_ok=True)
        with exclusive_run(args.output):
            plan, annotation = prepare(args)
            cases = plan["cohort"]["cases"]
            calls = sum(len(case_distances(c, args.distances)) * (
                args.repeats if i < args.timing_repeat_cases else 1
            ) for i, c in enumerate(cases))
            print(json.dumps({"phase": "plan", "cases": len(cases), "model_evaluations": calls,
                              "rough_hours_at_2_seconds_per_evaluation": round(calls * 2 / 3600, 2),
                              "note": "Planning estimate only; no statistical extrapolation."}), flush=True)
            if args.prepare_only:
                return 0
            return execute(args, plan, annotation)
    except Exception as exc:
        print(f"Distance scan failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
