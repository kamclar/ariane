#!/usr/bin/env python3
"""Compare SpliceAI search distances using fresh, unrounded local inference.

Run inside the profile-pinned SpliceAI container. No HTTP requests, application
caches, or classification settings are used or modified. NumPy/TensorFlow and
SpliceAI are imported only when needed so analysis helpers remain testable.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import statistics
import sys
import time
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
CHANNELS = {"AG": (1, 1), "AL": (1, -1), "DG": (2, 1), "DL": (2, -1)}
MODEL_RADIUS = 5000


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def score_band(score):
    """Only the SpliceAI threshold band, not an assigned evidence criterion."""
    if score <= 0.1:
        return "le_0.1"
    if score >= 0.2:
        return "ge_0.2"
    return "between_0.1_and_0.2"


def summarize_arrays(ref, alt, distance, tolerance):
    import numpy as np

    expected_shape = (2 * distance + 1, 3)
    if ref.shape != expected_shape or alt.shape != expected_shape:
        raise ValueError(f"Unexpected prediction shapes: {ref.shape}, {alt.shape}")
    if not np.isfinite(ref).all() or not np.isfinite(alt).all():
        raise ValueError("Non-finite raw model predictions")
    delta = alt - ref
    offsets = np.arange(-distance, distance + 1)
    components = {}
    for channel, (column, sign) in CHANNELS.items():
        changes = sign * delta[:, column]
        index = int(changes.argmax())
        components[channel] = {
            "delta": float(changes[index]),
            "position": int(offsets[index]),
            "ref": float(ref[index, column]),
            "alt": float(alt[index, column]),
        }
    outside = abs(offsets) > MODEL_RADIUS
    outside_values = np.abs(delta[outside, 1:])
    maximum = max(item["delta"] for item in components.values())
    return {
        "components": components,
        "maximum_delta": maximum,
        "band_raw": score_band(maximum),
        "band_rounded_3dp": score_band(float(f"{maximum:.3f}")),
        "outside_model_radius": {
            "positions_examined": int(outside.sum()),
            "maximum_absolute_delta": float(outside_values.max()) if outside.any() else None,
            "positions_with_any_nonzero_delta": int(
                np.any(outside_values != 0, axis=1).sum()
            ) if outside.any() else None,
            "positions_above_tolerance": int(
                np.any(outside_values > tolerance, axis=1).sum()
            ) if outside.any() else None,
        },
        "edge_deltas": {
            str(offset): {
                channel: float(sign * delta[offset + distance, column])
                for channel, (column, sign) in CHANNELS.items()
            }
            for offset in (-5001, -5000, -4999, 4999, 5000, 5001)
            if abs(offset) <= distance
        },
    }


def compare_arrays(smaller, larger, d_small, d_large, tolerance):
    import numpy as np

    small_summary = summarize_arrays(*smaller, d_small, tolerance)
    large_summary = summarize_arrays(*larger, d_large, tolerance)
    start = d_large - d_small
    stop = start + 2 * d_small + 1
    per_channel = {}
    for channel in CHANNELS:
        small = small_summary["components"][channel]
        large = large_summary["components"][channel]
        difference = large["delta"] - small["delta"]
        per_channel[channel] = {
            "delta_difference": difference,
            "within_tolerance": abs(difference) <= tolerance,
            "rounded_3dp_equal": f"{small['delta']:.3f}" == f"{large['delta']:.3f}",
            "position_changed": small["position"] != large["position"],
            "small_position": small["position"],
            "large_position": large["position"],
            "both_deltas_negligible": max(abs(small["delta"]), abs(large["delta"]))
            <= tolerance,
        }
    return {
        "smaller_distance": d_small,
        "larger_distance": d_large,
        "channels": per_channel,
        "all_delta_scores_within_tolerance": all(
            item["within_tolerance"] for item in per_channel.values()
        ),
        "all_rounded_delta_scores_equal": all(
            item["rounded_3dp_equal"] for item in per_channel.values()
        ),
        "raw_band_equal": small_summary["band_raw"] == large_summary["band_raw"],
        "rounded_band_equal": small_summary["band_rounded_3dp"]
        == large_summary["band_rounded_3dp"],
        "overlapping_ref_max_abs_difference": float(
            np.abs(smaller[0] - larger[0][start:stop]).max()
        ),
        "overlapping_alt_max_abs_difference": float(
            np.abs(smaller[1] - larger[1][start:stop]).max()
        ),
    }


class RecordingModel:
    """Capture the two actual model calls made by the upstream scorer."""

    def __init__(self, model):
        self.model = model
        self.outputs = []

    def predict(self, inputs, **kwargs):
        kwargs.setdefault("verbose", 0)
        outputs = self.model.predict(inputs, **kwargs)
        self.outputs.append(outputs)
        return outputs


def infer(annotator, case, distance):
    import numpy as np
    from spliceai.utils import get_delta_scores

    for model in annotator.models:
        model.outputs.clear()
    chrom, pos, ref, alt = case["variant"].split("-")
    record = SimpleNamespace(chrom=chrom, pos=int(pos), ref=ref, alts=(alt,))
    started = time.perf_counter()
    rows = get_delta_scores(record, annotator, distance, 0)
    elapsed = time.perf_counter() - started
    if len(rows) != 1 or rows[0].get("NAME") != case["transcript"]:
        raise ValueError(f"Expected exactly one score row for {case['transcript']}")
    if any(len(model.outputs) != 2 for model in annotator.models):
        raise ValueError("Expected exactly one REF and one ALT call per ensemble model")
    ref_scores = np.mean([model.outputs[0] for model in annotator.models], axis=0)[0]
    alt_scores = np.mean([model.outputs[1] for model in annotator.models], axis=0)[0]
    if rows[0]["STRAND"] == "-":
        ref_scores, alt_scores = ref_scores[::-1], alt_scores[::-1]
    ref_scores, alt_scores = align_predictions(
        ref_scores, alt_scores, len(ref), len(alt), distance
    )
    return (ref_scores.copy(), alt_scores.copy()), rows[0], elapsed


def align_predictions(ref_scores, alt_scores, ref_length, alt_length, distance):
    """Reproduce the pinned upstream alignment for SNVs and simple indels."""
    import numpy as np

    if ref_length > 1 and alt_length > 1:
        raise ValueError("Complex delins/MNV alignment is outside this experiment")
    if ref_length > 1:
        cut = distance + alt_length
        alt_scores = np.concatenate((
            alt_scores[:cut], np.zeros((ref_length - alt_length, 3)), alt_scores[cut:]
        ))
    elif alt_length > 1:
        alt_scores = np.concatenate((
            alt_scores[:distance],
            np.max(alt_scores[distance:distance + alt_length], axis=0)[None, :],
            alt_scores[distance + alt_length:],
        ))
    return ref_scores, alt_scores


def make_cases(annotator, fixture, profile):
    cases = [dict(case, provenance="local_service_validation_cases.json")
             for case in fixture["cases"]]
    for gene, notation, allele in (
        ("BRCA1", "c.134+10A>C", "chr17-43115716-T-G"),
        ("BRCA1", "c.135-13T>G", "chr17-43106546-A-C"),
        ("BRCA2", "c.1909+1G>A", "chr13-32333388-G-A"),
    ):
        cases.append({
            "id": f"{gene}_{notation}", "variant": allele, "gene": gene,
            "transcript": profile["reference_transcripts"][gene]["ensembl"],
            "provenance": "data/coordinates/brca_intronic_snv_coordinates.json; "
            "historical D50 scores used only for case selection",
        })
    # Deliberate synthetic SNV halfway through the longest intron of BRCA2.
    # Derive the position and REF from this run's annotation and FASTA.
    transcript = profile["reference_transcripts"]["BRCA2"]["ensembl"]
    index = list(annotator.genes).index(transcript)
    intervals = [(int(end) + 1, int(start) - 1) for end, start in zip(
        annotator.exon_ends[index][:-1], annotator.exon_starts[index][1:]
    )]
    start, end = max(intervals, key=lambda pair: pair[1] - pair[0] + 1)
    pos = (start + end) // 2
    chrom = str(annotator.chroms[index])
    ref = annotator.ref_fasta[chrom][pos - 1:pos].seq.upper()
    alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
    cases.append({
        "id": "BRCA2_synthetic_long_intron_midpoint", "gene": "BRCA2",
        "variant": f"{chrom}-{pos}-{ref}-{alt}", "transcript": transcript,
        "provenance": "synthetic transition, position and REF derived from pinned annotation/FASTA",
        "intron_interval_1based": [start, end],
    })
    for case in cases:
        chrom, pos, ref, alt = case["variant"].split("-")
        observed = annotator.ref_fasta[chrom][int(pos) - 1:int(pos)].seq.upper()
        if observed != ref or ref == alt or len(ref) != 1 or len(alt) != 1:
            raise ValueError(f"Invalid SNV/reference for {case['id']}: {observed}")
    return cases


def save_report(output, report):
    (output / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )


def validate_capture(summary, service_row):
    """Ensure raw capture, genomic orientation and upstream aggregation agree."""
    for channel, component in summary["components"].items():
        for suffix, key in (("", "delta"), ("_REF", "ref"), ("_ALT", "alt")):
            field = f"DS_{channel}{suffix}"
            if float(f"{component[key]:.3f}") != float(service_row[field]):
                raise ValueError(f"Captured prediction disagrees with upstream {field}")
        if component["position"] != int(service_row[f"DP_{channel}"]):
            raise ValueError(f"Captured position disagrees with upstream DP_{channel}")


def run(args):
    import numpy as np
    import tensorflow as tf
    from pkg_resources import resource_filename
    from spliceai.utils import Annotator
    import spliceai.utils

    if args.repeats < 1 or not math.isfinite(args.tolerance) or args.tolerance <= 0:
        raise ValueError("Repeats and finite tolerance must be positive")
    if args.threads < 1:
        raise ValueError("Threads must be positive")
    distances = sorted(set(args.distances))
    if not {4999, 5000, 10000}.issubset(distances) or any(d < 1 or d > 10000 for d in distances):
        raise ValueError("Distances must include 4999, 5000 and 10000, within 1..10000")
    args.output.mkdir(parents=True, exist_ok=False)
    args.output_created = True
    profile_path = ROOT / "data/spliceai/enigma_v1_2_spliceai_profile.json"
    fixture_path = ROOT / "data/spliceai/local_service_validation_cases.json"
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
    engine = profile["approved_engine"]
    annotation = ROOT / engine["reference_transcript_annotation"]
    if sha256(annotation) != engine["reference_transcript_annotation_sha256"]:
        raise ValueError("Annotation does not match pinned profile")
    if fixture["profile_id"] != profile["profile_id"]:
        raise ValueError("Validation fixture belongs to a different profile")
    if args.image_digest != engine["docker_image"]:
        raise ValueError("Launch with the exact image recorded in the active profile")
    if version("spliceai") != engine["python_package_version"]:
        raise ValueError("Installed SpliceAI package does not match active profile")
    if os.environ.get("MODEL_COMMIT") != engine["model_commit"]:
        raise ValueError("Container model commit does not match active profile")
    tf.config.threading.set_intra_op_parallelism_threads(args.threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    started = time.perf_counter()
    annotator = Annotator(str(args.fasta), str(annotation))
    initialization_seconds = time.perf_counter() - started
    annotator.models = [RecordingModel(model) for model in annotator.models]
    cases = make_cases(annotator, fixture, profile)
    report = {
        "schema_version": 1, "status": "running",
        "started_at": datetime.now(timezone.utc).isoformat(),
        "scope": "small SNV pilot; cannot establish equivalence for all variants or indels",
        "profile": profile, "profile_sha256": sha256(profile_path),
        "fixture_sha256": sha256(fixture_path), "script_sha256": sha256(__file__),
        "image_digest": args.image_digest,
        "python": sys.version, "platform": platform.platform(),
        "packages": {"spliceai": version("spliceai"), "tensorflow": tf.__version__,
                     "numpy": np.__version__},
        "model_commit": os.environ["MODEL_COMMIT"],
        "model_weights_sha256": {
            f"spliceai{i}.h5": sha256(resource_filename("spliceai", f"models/spliceai{i}.h5"))
            for i in range(1, 6)
        },
        "spliceai_utils_sha256": sha256(spliceai.utils.__file__),
        "reference_fasta": str(args.fasta),
        "reference_fasta_bytes": args.fasta.stat().st_size,
        "reference_fasta_sha256": sha256(args.fasta),
        "annotation_sha256": sha256(annotation),
        "tensorflow_devices": [device.name for device in tf.config.list_physical_devices()],
        "threads": args.threads, "mask": 0, "distances": distances,
        "absolute_tolerance": args.tolerance, "repeats": args.repeats,
        "cache": "none; direct fresh REF/ALT inference for every run",
        "model_initialization_seconds": initialization_seconds,
        "timing_note": "One untimed warmup per shape; cyclic/reversed distance order. "
        "Wall time includes upstream scoring and formatting, not application/network latency.",
        "cases": [], "warmups": [],
    }
    save_report(args.output, report)
    for distance in distances:
        _, _, elapsed = infer(annotator, cases[0], distance)
        report["warmups"].append({"distance": distance, "seconds": elapsed})
        print(f"Warmup D={distance}: {elapsed:.2f}s", flush=True)
    failed_baselines = []
    for case_index, case in enumerate(cases):
        entry = {"case": case, "runs": [], "comparisons": [], "repeat_comparisons": []}
        first_arrays = {}
        samples = {distance: [] for distance in distances}
        for repeat in range(args.repeats):
            shift = case_index % len(distances)
            order = distances[shift:] + distances[:shift]
            if repeat % 2:
                order = list(reversed(order))
            for distance in order:
                arrays, service_row, elapsed = infer(annotator, case, distance)
                summary = summarize_arrays(*arrays, distance, args.tolerance)
                validate_capture(summary, service_row)
                artifact = f"case-{case_index + 1}-d{distance}-repeat{repeat + 1}.npz"
                np.savez_compressed(args.output / artifact, ref=arrays[0], alt=arrays[1])
                baseline = None
                if distance == profile["max_distance"] and "delta_scores" in case:
                    expected = {**case["delta_scores"], **case["reference_scores"],
                                **case["alternate_scores"]}
                    baseline = {field: {"expected": value, "observed": service_row.get(field)}
                                for field, value in expected.items()
                                if float(service_row.get(field, "nan")) != value}
                    if baseline:
                        failed_baselines.append(case["id"])
                entry["runs"].append({
                    "distance": distance, "repeat": repeat + 1, "seconds": elapsed,
                    **summary,
                    "service_fields": {key: value for key, value in service_row.items()
                                       if key.startswith(("DS_", "DP_"))},
                    "baseline_mismatches": baseline,
                    "raw_arrays": artifact, "raw_arrays_sha256": sha256(args.output / artifact),
                })
                if distance in first_arrays:
                    entry["repeat_comparisons"].append(compare_arrays(
                        first_arrays[distance], arrays, distance, distance, args.tolerance
                    ))
                else:
                    first_arrays[distance] = arrays
                samples[distance].append(elapsed)
                print(f"{case['id']} D={distance} repeat={repeat + 1}: "
                      f"max={summary['maximum_delta']:.9f}, {elapsed:.2f}s", flush=True)
        for smaller in distances:
            if smaller != max(distances):
                entry["comparisons"].append(compare_arrays(
                    first_arrays[smaller], first_arrays[max(distances)],
                    smaller, max(distances), args.tolerance,
                ))
        entry["comparisons"].append(compare_arrays(
            first_arrays[4999], first_arrays[5000], 4999, 5000, args.tolerance,
        ))
        entry["timing_seconds"] = {
            str(distance): {"median": statistics.median(values), "min": min(values),
                            "max": max(values)}
            for distance, values in samples.items()
        }
        report["cases"].append(entry)
        save_report(args.output, report)
    report["status"] = "baseline_mismatch" if failed_baselines else "completed"
    report["baseline_mismatch_cases"] = sorted(set(failed_baselines))
    report["completed_at"] = datetime.now(timezone.utc).isoformat()
    save_report(args.output, report)
    return 1 if failed_baselines else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--image-digest", required=True, help="Exact docker image used to launch")
    parser.add_argument("--distances", nargs="+", type=int, default=[50, 4999, 5000, 10000])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--cohort", choices=("pilot", "stratified"), default="pilot")
    parser.add_argument("--positions-per-region", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--no-indels", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--retry-errors", action="store_true")
    parser.add_argument("--max-cases", type=int, help="New/retried cases in this invocation")
    parser.add_argument("--max-seconds", type=float, help="Pause between cases after this duration")
    parser.add_argument("--max-errors", type=int, default=5)
    parser.add_argument("--timing-repeat-cases", type=int, default=12)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--save-arrays", choices=("all", "changed", "none"), default="changed")
    args = parser.parse_args()
    if args.cohort == "stratified":
        if str(ROOT) not in sys.path:
            sys.path.insert(0, str(ROOT))
        from scripts.spliceai_distance_scan import run_scan
        return run_scan(args)
    if args.resume or args.prepare_only or args.retry_errors or args.max_cases or args.max_seconds:
        parser.error("Prepare, resume and run limits require --cohort stratified")
    args.output_created = False
    try:
        return run(args)
    except Exception as exc:
        report_path = args.output / "report.json"
        if args.output_created and report_path.exists():
            report = json.loads(report_path.read_text(encoding="utf-8"))
            if report.get("status") == "running":
                report.update(status="failed", error=f"{type(exc).__name__}: {exc}")
                save_report(args.output, report)
        print(f"Experiment failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
