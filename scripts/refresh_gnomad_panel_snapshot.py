"""Refresh the panel gnomAD frequency and coverage snapshots from official data.

The active releases and panel intervals live in
``backend/data/gnomad/gnomad_panel_manifest.json``.  The build reads the
official gnomAD Hail Tables, restricts them to the configured intervals and
publishes JSON only after every source and output validation succeeds.

This is a build-time tool.  Hail and Java are not application runtime
dependencies.  Install ``requirements-data.txt`` in a separate environment or
run the documented Hail container command.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
GNOMAD_DIR = ROOT / "backend" / "data" / "gnomad"
DEFAULT_MANIFEST = GNOMAD_DIR / "gnomad_panel_manifest.json"
DEFAULT_VARIANTS = GNOMAD_DIR / "gnomad_brca_frequency_snapshot.json"
DEFAULT_COVERAGE = GNOMAD_DIR / "gnomad_brca_coverage_snapshot.json"
FAF_CONFIDENCE = 0.95


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _active_policies(manifest: dict) -> dict[str, dict]:
    policies = manifest.get("classification_policies") or {}
    return {
        policy_id: policy
        for policy_id, policy in policies.items()
        if policy.get("status") == "active"
    }


def _scored_ancestries(manifest: dict) -> list[str]:
    ancestries = {
        ancestry
        for policy in _active_policies(manifest).values()
        for ancestry in (
            policy.get("frequency_criteria", {}).get(
                "scored_non_founder_ancestries", []
            )
        )
    }
    return sorted(ancestries)


def _validate_enigma_brca_policy(policy: dict) -> None:
    frequency = policy.get("frequency_criteria") or {}
    ba1 = frequency.get("ba1") or {}
    bs1 = frequency.get("bs1") or {}
    pm2 = frequency.get("pm2") or {}
    observations = frequency.get("outbred_observation_count") or {}
    if set(frequency.get("scored_non_founder_ancestries") or []) != {
        "afr", "amr", "eas", "nfe", "sas"
    }:
        raise RuntimeError("ENIGMA v1.2 non-founder ancestry set changed")
    strong = bs1.get("strong") or {}
    supporting = bs1.get("supporting") or {}
    automation = policy.get("automation_scope") or {}
    if (
        set(policy.get("applicable_genes") or []) != {"BRCA1", "BRCA2"}
        or policy.get("policy_scope") != "gnomad_frequency_criteria_only"
        or automation.get("criteria_defined_here") != ["BA1", "BS1", "PM2"]
        or automation.get("inherit_other_vcep_criteria") is not False
        or frequency.get("metric") != "faf95"
        or frequency.get("subset") != "non_cancer"
        or frequency.get("dataset_combination")
        != "assess_exome_and_genome_separately_then_use_highest_qualifying_faf95"
        or observations.get("single_observation") != "not_informative"
        or observations.get("minimum_observations_for_ba1_bs1") != 2
        or observations.get("pm2_requires_zero_observations") is not True
        or ba1.get("operator") != ">"
        or ba1.get("threshold") != 0.001
        or ba1.get("strength") != "Stand-alone"
        or ba1.get("points") != -99
        or ba1.get("minimum_mean_depth") != 20.0
        or ba1.get("requires_pass_record") is not True
        or ba1.get("exclude_well_established_pathogenic_founder_variants") is not True
        or bs1.get("minimum_mean_depth") != 20.0
        or bs1.get("requires_pass_record") is not True
        or bs1.get("exclude_well_established_pathogenic_founder_variants") is not True
        or strong.get("operator") != ">"
        or strong.get("threshold") != 0.0001
        or strong.get("strength") != "Strong"
        or strong.get("points") != -4
        or supporting.get("lower_operator") != ">"
        or supporting.get("lower_threshold") != 0.00002
        or supporting.get("upper_operator") != "<="
        or supporting.get("upper_threshold") != 0.0001
        or supporting.get("strength") != "Supporting"
        or supporting.get("points") != -1
        or pm2.get("strength") != "Supporting"
        or pm2.get("points") != 1
        or pm2.get("minimum_mean_depth") != 25.0
        or pm2.get("single_outbred_observation") != "not_informative"
        or pm2.get("coverage_scope") != "variant_reference_span"
    ):
        raise RuntimeError("ENIGMA BRCA v1.2 frequency policy is incomplete or changed")


def _load_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 2:
        raise RuntimeError("unsupported gnomAD panel manifest schema")
    if not manifest.get("datasets") or not manifest.get("targets"):
        raise RuntimeError("gnomAD panel manifest has no datasets or targets")
    policies = manifest.get("classification_policies") or {}
    active_policies = _active_policies(manifest)
    if not active_policies:
        raise RuntimeError("gnomAD panel manifest has no active classification policy")
    if any(
        policy.get("automatic_release_activation") is not False
        for policy in policies.values()
    ):
        raise RuntimeError("automatic gnomAD release activation must be disabled")
    if "enigma_brca_v1_2" in active_policies:
        _validate_enigma_brca_policy(active_policies["enigma_brca_v1_2"])
    scored_ancestries = _scored_ancestries(manifest)
    if not scored_ancestries:
        raise RuntimeError("active policies define no scored population groups")

    dataset_keys = [item.get("dataset_key") for item in manifest["datasets"]]
    if None in dataset_keys or len(dataset_keys) != len(set(dataset_keys)):
        raise RuntimeError("gnomAD dataset keys are missing or duplicated")
    runtime_keys = [item.get("runtime_key") for item in manifest["datasets"]]
    if None in runtime_keys or len(runtime_keys) != len(set(runtime_keys)):
        raise RuntimeError("gnomAD runtime dataset keys are missing or duplicated")
    available_runtime_keys = set(runtime_keys)
    for policy_id, policy in active_policies.items():
        frequency = policy.get("frequency_criteria") or {}
        required = set(frequency.get("required_dataset_runtime_keys") or [])
        pm2_required = set(
            (frequency.get("pm2") or {}).get(
                "required_absence_dataset_runtime_keys", []
            )
        )
        if not required or not required.issubset(available_runtime_keys):
            raise RuntimeError(f"policy {policy_id} refers to unavailable datasets")
        if pm2_required != required:
            raise RuntimeError(
                f"policy {policy_id} PM2 dataset requirements are inconsistent"
            )
    genes = [item.get("gene") for item in manifest["targets"]]
    if None in genes or len(genes) != len(set(genes)):
        raise RuntimeError("gnomAD target genes are missing or duplicated")
    for target in manifest["targets"]:
        policy_id = target.get("classification_policy_id")
        if policy_id not in active_policies:
            raise RuntimeError(
                f"target {target.get('gene')} has no active gene-specific policy"
            )
        if target.get("gene") not in (
            active_policies[policy_id].get("applicable_genes") or []
        ):
            raise RuntimeError(
                f"policy {policy_id} is not approved for target {target.get('gene')}"
            )
        if target.get("activation_status") != "active":
            raise RuntimeError(
                f"target {target.get('gene')} is not explicitly active"
            )
        if not target.get("reference_transcript") or not target.get(
            "vcep_specification"
        ):
            raise RuntimeError(
                f"target {target.get('gene')} lacks transcript or VCEP provenance"
            )

    for dataset in manifest["datasets"]:
        assembly = dataset.get("assembly")
        if dataset.get("faf95_mode") not in {
            "published_hail_faf",
            "hail_native_from_non_cancer_ac_an",
        }:
            raise RuntimeError(f"unsupported FAF95 mode for {dataset['dataset_key']}")
        context_codes = [
            item.get("code") for item in dataset.get("excluded_population_context", [])
        ]
        if (
            None in context_codes
            or len(context_codes) != len(set(context_codes))
            or set(context_codes).intersection(scored_ancestries)
        ):
            raise RuntimeError(
                f"invalid excluded population context for {dataset['dataset_key']}"
            )
        if any(
            item.get("category") not in {
                "founder_population",
                "not_enigma_non_founder_population",
            }
            for item in dataset.get("excluded_population_context", [])
        ):
            raise RuntimeError(
                f"population context category is missing for {dataset['dataset_key']}"
            )
        for target in manifest["targets"]:
            interval = (target.get("intervals") or {}).get(assembly)
            if not interval or int(interval["start"]) > int(interval["end"]):
                raise RuntimeError(
                    f"target {target['gene']} has no valid {assembly} interval"
                )
    return manifest


def _http_source_identity(url: str) -> dict:
    request = urllib.request.Request(url, method="HEAD")
    request.add_header("User-Agent", "ARIANE gnomAD snapshot builder/2.0")
    with urllib.request.urlopen(request, timeout=90) as response:
        headers = response.headers
        hash_headers = headers.get_all("x-goog-hash") or []
        return {
            "url": url,
            "etag": headers.get("ETag"),
            "x_goog_hash": ", ".join(hash_headers) or None,
            "content_length": headers.get("Content-Length"),
            "last_modified": headers.get("Last-Modified"),
            "generation": headers.get("x-goog-generation"),
        }


def _verify_sources(manifest: dict) -> dict[str, dict[str, dict]]:
    verified: dict[str, dict[str, dict]] = {}
    for dataset in manifest["datasets"]:
        dataset_verified = {}
        for source_type in ("sites", "coverage"):
            url = dataset[f"{source_type}_metadata_url"]
            observed = _http_source_identity(url)
            expected = dataset[f"{source_type}_metadata_identity"]
            for field in ("etag", "x_goog_hash"):
                observed_value = observed.get(field)
                expected_value = expected.get(field)
                if field == "x_goog_hash":
                    observed_value = _hash_header_parts(observed_value)
                    expected_value = _hash_header_parts(expected_value)
                if observed_value != expected_value:
                    raise RuntimeError(
                        f"official {source_type} source identity changed for "
                        f"{dataset['dataset_key']}: {field} expected "
                        f"{expected.get(field)!r}, observed {observed.get(field)!r}"
                    )
            dataset_verified[source_type] = observed
        verified[dataset["dataset_key"]] = dataset_verified
    return verified


def _hash_header_parts(value: str | None) -> frozenset[str]:
    return frozenset(part.strip() for part in (value or "").split(",") if part.strip())


def _discover_releases(manifest: dict) -> list[str]:
    discovery = manifest["release_discovery"]
    query = urllib.parse.urlencode(
        {
            "prefix": discovery["prefix"],
            "delimiter": "/",
            "maxResults": 1000,
        }
    )
    request = urllib.request.Request(
        discovery["gcs_bucket_api"] + "?" + query,
        headers={"User-Agent": "ARIANE gnomAD release checker/2.0"},
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        payload = json.load(response)
    releases = []
    for prefix in payload.get("prefixes", []):
        match = re.fullmatch(r"release/([^/]+)/", prefix)
        if match:
            releases.append(match.group(1))
    return sorted(set(releases), key=_version_key)


def _candidate_sites_metadata_url(dataset: dict, candidate_release: str) -> str:
    """Return the same sites product path for a candidate release.

    A release directory alone is not evidence that a newer release contains
    the small-variant Hail Table used by this pipeline.  For example, gnomAD
    3.1.3 is an STR release and has no v3.1.3 sites Table.  Replacing the pinned
    release in the complete metadata URL checks the equivalent product rather
    than merely comparing release labels.
    """
    active_release = str(dataset["release"])
    active_url = str(dataset["sites_metadata_url"])
    if active_release not in active_url:
        raise RuntimeError(
            f"sites metadata URL for {dataset['dataset_key']} does not contain "
            f"its pinned release {active_release}"
        )
    return active_url.replace(active_release, candidate_release)


def _source_exists(url: str) -> tuple[bool, int | None]:
    request = urllib.request.Request(
        url,
        method="HEAD",
        headers={"User-Agent": "ARIANE gnomAD release checker/2.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            return 200 <= response.status < 300, response.status
    except urllib.error.HTTPError as exc:
        return False, exc.code


def _version_key(value: str) -> tuple:
    parts = re.split(r"([0-9]+)", value)
    return tuple(int(part) if part.isdigit() else part for part in parts)


def check_updates(manifest_path: Path) -> dict:
    manifest = _load_manifest(manifest_path)
    verified = _verify_sources(manifest)
    available = _discover_releases(manifest)
    active = sorted({item["release"] for item in manifest["datasets"]}, key=_version_key)
    candidates = [release for release in available if release not in active]
    same_major_newer_directories = {
        release: [
            candidate
            for candidate in candidates
            if _numeric_version(candidate)[:1] == _numeric_version(release)[:1]
            and _numeric_version(candidate) > _numeric_version(release)
        ]
        for release in active
    }
    product_checks = []
    compatible_updates: dict[str, list[str]] = {release: [] for release in active}
    for dataset in manifest["datasets"]:
        release = dataset["release"]
        for candidate in same_major_newer_directories[release]:
            url = _candidate_sites_metadata_url(dataset, candidate)
            exists, status = _source_exists(url)
            product_checks.append(
                {
                    "dataset": dataset["dataset_key"],
                    "active_release": release,
                    "candidate_release": candidate,
                    "equivalent_sites_product_url": url,
                    "http_status": status,
                    "equivalent_sites_product_available": exists,
                }
            )
            if exists:
                compatible_updates[release].append(candidate)
    return {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_sha256": _canonical_sha256(manifest),
        "active_releases": active,
        "available_releases": available,
        "detected_other_releases": candidates,
        "newer_same_major_release_directories": same_major_newer_directories,
        "equivalent_product_checks": product_checks,
        "newer_equivalent_product_releases": compatible_updates,
        "automatic_activation": False,
        "activation_status": (
            "review_required"
            if any(compatible_updates.values())
            else "up_to_date"
        ),
        "source_identities": verified,
    }


def _numeric_version(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in re.findall(r"[0-9]+", value))


def _require_hail():
    try:
        import hail as hl
    except ImportError as exc:  # pragma: no cover - build environment only
        raise SystemExit(
            "Hail is required for refresh. Install requirements-data.txt in a "
            "separate build environment or use the documented Docker command."
        ) from exc
    return hl


def _hail_intervals(hl, manifest: dict, assembly: str) -> list:
    padding = int(manifest["region_padding_bp"])
    intervals = []
    contig_prefix = "chr" if assembly == "GRCh38" else ""
    for target in manifest["targets"]:
        item = target["intervals"][assembly]
        start = int(item["start"]) - padding
        end = int(item["end"]) + padding
        intervals.append(
            hl.parse_locus_interval(
                f"{contig_prefix}{item['chrom']}:{start}-{end}",
                reference_genome=assembly,
            )
        )
    return intervals


def _frequency_indices(globals_value, dataset: dict, ancestries: Iterable[str]) -> tuple:
    index = dict(globals_value.freq_index_dict)
    style = dataset["frequency_index_style"]
    subset = dataset["subset"]
    if style == "v2_underscore":
        overall_key = subset
        population_keys = {pop: f"{subset}_{pop}" for pop in ancestries}
    elif style == "v3_hyphen_adj":
        overall_key = f"{subset}-adj"
        population_keys = {pop: f"{subset}-{pop}-adj" for pop in ancestries}
    else:
        raise RuntimeError(f"unsupported frequency index style {style!r}")
    missing = [key for key in [overall_key, *population_keys.values()] if key not in index]
    if missing:
        raise RuntimeError(
            f"official frequency indices missing for {dataset['dataset_key']}: {missing}"
        )
    return index[overall_key], {pop: index[key] for pop, key in population_keys.items()}


def _population_frequency_key(dataset: dict, population: str) -> str:
    subset = dataset["subset"]
    style = dataset["frequency_index_style"]
    if style == "v2_underscore":
        return f"{subset}_{population}"
    if style == "v3_hyphen_adj":
        return f"{subset}-{population}-adj"
    raise RuntimeError(f"unsupported frequency index style {style!r}")


def _local_uri(path: Path) -> str:
    return path.resolve().as_uri()


def _export_frequencies(hl, manifest: dict, dataset: dict, destination: Path) -> None:
    table = hl.read_table(dataset["sites_hail_uri"])
    table = hl.filter_intervals(
        table, _hail_intervals(hl, manifest, dataset["assembly"])
    )
    globals_value = hl.eval(table.index_globals())
    ancestries = tuple(_scored_ancestries(manifest))
    overall_index, pop_indices = _frequency_indices(
        globals_value, dataset, ancestries
    )
    frequency_index = dict(globals_value.freq_index_dict)
    context_config = tuple(dataset.get("excluded_population_context", []))
    context_indices = {}
    for item in context_config:
        key = _population_frequency_key(dataset, item["code"])
        if key not in frequency_index:
            raise RuntimeError(
                f"official context population index {key!r} is missing from "
                f"{dataset['dataset_key']}"
            )
        context_indices[item["code"]] = frequency_index[key]
    overall = table.freq[overall_index]
    table = table.filter(hl.is_defined(overall.AC) & (overall.AC > 0) & (overall.AN > 0))
    # Hail expressions are tied to their source Table.  Rebind after filter so
    # no pre-filter expression is mixed into the exported Table.
    overall = table.freq[overall_index]

    faf_expressions = {}
    faf_index = {}
    if dataset["faf95_mode"] == "published_hail_faf":
        faf_index = dict(globals_value.faf_index_dict)
        for pop in ancestries:
            key = f"{dataset['subset']}_{pop}"
            if key not in faf_index:
                raise RuntimeError(
                    f"official published FAF95 index {key!r} is missing from "
                    f"{dataset['dataset_key']}"
                )
            faf_expressions[pop] = table.faf[faf_index[key]].faf95
    else:
        for pop, index in pop_indices.items():
            frequency = table.freq[index]
            faf_expressions[pop] = hl.if_else(
                hl.is_defined(frequency.AC) & hl.is_defined(frequency.AN) & (frequency.AN > 0),
                hl.experimental.filtering_allele_frequency(
                    frequency.AC, frequency.AN, FAF_CONFIDENCE
                ),
                hl.missing(hl.tfloat64),
            )

    populations = {pop: table.freq[index] for pop, index in pop_indices.items()}
    context_populations = {
        pop: table.freq[index] for pop, index in context_indices.items()
    }
    context_faf_expressions = {}
    for pop, frequency in context_populations.items():
        published_key = _population_frequency_key(dataset, pop)
        if dataset["faf95_mode"] == "published_hail_faf" and published_key in faf_index:
            context_faf_expressions[pop] = table.faf[faf_index[published_key]].faf95
        else:
            context_faf_expressions[pop] = hl.if_else(
                hl.is_defined(frequency.AC)
                & hl.is_defined(frequency.AN)
                & (frequency.AN > 0),
                hl.experimental.filtering_allele_frequency(
                    frequency.AC, frequency.AN, FAF_CONFIDENCE
                ),
                hl.missing(hl.tfloat64),
            )
    selected = table.select(
        chrom=table.locus.contig,
        pos=table.locus.position,
        ref=table.alleles[0],
        alt=table.alleles[1],
        filter=hl.if_else(
            hl.len(table.filters) == 0,
            "PASS",
            hl.delimit(hl.sorted(hl.array(table.filters)), ";"),
        ),
        ac=overall.AC,
        an=overall.AN,
        af=overall.AF,
        nhomalt=overall.homozygote_count,
        **{f"ac_{pop}": populations[pop].AC for pop in ancestries},
        **{f"an_{pop}": populations[pop].AN for pop in ancestries},
        **{f"af_{pop}": populations[pop].AF for pop in ancestries},
        **{f"faf95_{pop}": faf_expressions[pop] for pop in ancestries},
        **{f"context_ac_{pop}": context_populations[pop].AC for pop in context_populations},
        **{f"context_an_{pop}": context_populations[pop].AN for pop in context_populations},
        **{f"context_af_{pop}": context_populations[pop].AF for pop in context_populations},
        **{f"context_faf95_{pop}": context_faf_expressions[pop] for pop in context_populations},
    )
    selected.key_by().export(_local_uri(destination))


def _export_coverage(hl, manifest: dict, dataset: dict, destination: Path) -> None:
    table = hl.read_table(dataset["coverage_hail_uri"])
    table = hl.filter_intervals(
        table, _hail_intervals(hl, manifest, dataset["assembly"])
    )
    row_fields = set(table.row_value.dtype.fields)
    median_field = "median" if "median" in row_fields else "median_approx"
    selected = table.select(
        chrom=table.locus.contig,
        pos=table.locus.position,
        mean_depth=table.mean,
        median_depth=table[median_field],
        over_20=table.over_20,
        over_25=table.over_25,
    )
    selected.key_by().export(_local_uri(destination))


def _optional_float(value: str | None) -> float | None:
    if value in (None, "", "NA", "null"):
        return None
    return float(value)


def _optional_int(value: str | None) -> int | None:
    number = _optional_float(value)
    return None if number is None else int(number)


def _gene_for_position(manifest: dict, assembly: str, chrom: str, pos: int) -> str:
    padding = int(manifest["region_padding_bp"])
    chrom = chrom.removeprefix("chr")
    matches = []
    for target in manifest["targets"]:
        interval = target["intervals"][assembly]
        if (
            str(interval["chrom"]).removeprefix("chr") == chrom
            and int(interval["start"]) - padding <= pos <= int(interval["end"]) + padding
        ):
            matches.append(target["gene"])
    if not matches:
        raise RuntimeError(
            f"position {assembly} {chrom}:{pos} does not map to a panel target"
        )
    # Gene intervals can legitimately overlap in larger panels.  The value is
    # informational only; frequency and coverage keys remain coordinate based.
    return ",".join(sorted(matches))


def _read_frequency_export(path: Path, manifest: dict, dataset: dict) -> dict:
    ancestries = tuple(_scored_ancestries(manifest))
    records = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row["chrom"].removeprefix("chr")
            pos = int(row["pos"])
            ancestry_faf = {
                pop: _optional_float(row[f"faf95_{pop}"]) for pop in ancestries
            }
            if any(value is None for value in ancestry_faf.values()):
                raise RuntimeError(
                    f"non-cancer FAF95 is incomplete for {dataset['dataset_key']} "
                    f"{chrom}:{pos}:{row['ref']}:{row['alt']}"
                )
            faf_pop, faf_max = max(ancestry_faf.items(), key=lambda item: item[1])
            ancestry_af = {
                pop: _optional_float(row[f"af_{pop}"]) for pop in ancestries
            }
            ancestry_ac = {
                pop: _optional_int(row[f"ac_{pop}"]) for pop in ancestries
            }
            ancestry_an = {
                pop: _optional_int(row[f"an_{pop}"]) for pop in ancestries
            }
            usable_af = {pop: value for pop, value in ancestry_af.items() if value is not None}
            popmax_pop, popmax_af = (
                max(usable_af.items(), key=lambda item: item[1])
                if usable_af
                else (None, None)
            )
            variant_id = f"{chrom}-{pos}-{row['ref']}-{row['alt']}"
            method = (
                "official_gnomad_hail_table_non_cancer_faf95"
                if dataset["faf95_mode"] == "published_hail_faf"
                else "hail.experimental.filtering_allele_frequency_from_official_non_cancer_ac_an"
            )
            excluded_population_context = {}
            for item in dataset.get("excluded_population_context", []):
                pop = item["code"]
                excluded_population_context[pop] = {
                    "label": item["label"],
                    "category": item["category"],
                    "ac": _optional_int(row[f"context_ac_{pop}"]),
                    "an": _optional_int(row[f"context_an_{pop}"]),
                    "af": _optional_float(row[f"context_af_{pop}"]),
                    "faf95": _optional_float(row[f"context_faf95_{pop}"]),
                    "used_for_ba1_bs1": False,
                    "used_for_pm2_presence": False,
                }
            non_founder_observed = any(
                (value or 0) > 0 for value in ancestry_ac.values()
            )
            record = {
                "variant_id": variant_id,
                "chrom": chrom,
                "pos": pos,
                "ref": row["ref"],
                "alt": row["alt"],
                "filter": row["filter"],
                "af": _optional_float(row["af"]),
                "ac": _optional_int(row["ac"]),
                "an": _optional_int(row["an"]),
                "nhomalt": _optional_int(row["nhomalt"]) or 0,
                "faf95_max": faf_max,
                "faf95_pop": faf_pop,
                "faf95_by_ancestry": ancestry_faf,
                "non_founder_ac_by_ancestry": ancestry_ac,
                "non_founder_an_by_ancestry": ancestry_an,
                "non_founder_observed": non_founder_observed,
                "faf95_scope": "non_cancer_non_founder_ancestries",
                "faf95_method": method,
                "excluded_population_context": excluded_population_context,
                "popmax_af": popmax_af,
                "popmax_pop": popmax_pop,
                "dataset": dataset["dataset_key"],
                "build": dataset["assembly"],
            }
            records.setdefault(variant_id, []).append(record)
    return records


def _read_coverage_export(path: Path, manifest: dict, dataset: dict) -> dict:
    records = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row["chrom"].removeprefix("chr")
            pos = int(row["pos"])
            mean = _optional_float(row["mean_depth"])
            key = f"{dataset['dataset_key']}|{dataset['assembly']}|{chrom}|{pos}"
            records[key] = {
                "chrom": chrom,
                "pos": pos,
                "mean_depth": mean,
                "median_depth": _optional_float(row["median_depth"]),
                "over_20": _optional_float(row["over_20"]),
                "over_25": _optional_float(row["over_25"]),
                "dataset": dataset["dataset_key"],
                "dataset_key": dataset["dataset_key"],
                "build": dataset["assembly"],
                "gene": _gene_for_position(
                    manifest, dataset["assembly"], chrom, pos
                ),
                "source": "official gnomAD Hail coverage table",
                "position_key": key,
            }
    return records


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, separators=(",", ":"))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def _regions_metadata(manifest: dict) -> dict:
    result = {}
    for target in manifest["targets"]:
        for assembly, interval in target["intervals"].items():
            result.setdefault(assembly, {})[target["gene"]] = interval
    return result


def _validate_built_payloads(variant_payload: dict, coverage_payload: dict, manifest: dict) -> None:
    variants = variant_payload.get("variants") or {}
    coverage = coverage_payload.get("coverage_by_position") or {}
    if not variants or not coverage:
        raise RuntimeError("refusing to publish an empty gnomAD snapshot")
    expected = {item["dataset_key"] for item in manifest["datasets"]}
    variant_datasets = {
        record["dataset"] for values in variants.values() for record in values
    }
    coverage_datasets = {item["dataset_key"] for item in coverage.values()}
    if variant_datasets != expected or coverage_datasets != expected:
        raise RuntimeError(
            "frequency or coverage exports do not contain exactly the manifest datasets"
        )
    for values in variants.values():
        for record in values:
            if record.get("faf95_max") is None or not record.get("faf95_method"):
                raise RuntimeError("frequency export contains unscorable FAF95")
            if set(record.get("non_founder_ac_by_ancestry") or {}) != set(
                _scored_ancestries(manifest)
            ):
                raise RuntimeError(
                    "frequency export lacks complete non-founder population counts"
                )
            dataset = next(
                item
                for item in manifest["datasets"]
                if item["dataset_key"] == record["dataset"]
            )
            expected_context = {
                item["code"] for item in dataset.get("excluded_population_context", [])
            }
            if set(record.get("excluded_population_context") or {}) != expected_context:
                raise RuntimeError(
                    "frequency export lacks excluded founder/population context"
                )
    if variant_payload["metadata"]["records_sha256"] != _canonical_sha256(variants):
        raise RuntimeError("frequency snapshot checksum validation failed")
    if coverage_payload["metadata"]["records_sha256"] != _canonical_sha256(coverage):
        raise RuntimeError("coverage snapshot checksum validation failed")


def refresh(
    manifest_path: Path,
    variants_path: Path,
    coverage_path: Path,
    staging_dir: Path | None = None,
    keep_staging: bool = False,
) -> None:
    manifest = _load_manifest(manifest_path)
    source_identities = _verify_sources(manifest)
    manifest_sha256 = _canonical_sha256(manifest)
    hl = _require_hail()

    own_staging = staging_dir is None
    staging = staging_dir or Path(tempfile.mkdtemp(prefix="ariane-gnomad-"))
    staging.mkdir(parents=True, exist_ok=True)
    hail_log = staging / "hail.log"
    hl.init(
        log=str(hail_log),
        quiet=True,
        # The pinned gnomAD bucket is public.  Explicit anonymous access keeps
        # local/container builds from probing the GCE metadata service or
        # depending on a developer's Google credentials.
        spark_conf={"spark.hadoop.fs.gs.auth.type": "UNAUTHENTICATED"},
    )
    all_variants: dict[str, list] = {}
    all_coverage: dict[str, dict] = {}
    dataset_metadata = []
    try:
        for dataset in manifest["datasets"]:
            frequency_tsv = staging / f"{dataset['dataset_key']}.frequency.tsv"
            coverage_tsv = staging / f"{dataset['dataset_key']}.coverage.tsv"
            _export_frequencies(hl, manifest, dataset, frequency_tsv)
            _export_coverage(hl, manifest, dataset, coverage_tsv)
            dataset_variants = _read_frequency_export(
                frequency_tsv, manifest, dataset
            )
            dataset_coverage = _read_coverage_export(
                coverage_tsv, manifest, dataset
            )
            for key, values in dataset_variants.items():
                all_variants.setdefault(key, []).extend(values)
            duplicate_coverage = set(all_coverage).intersection(dataset_coverage)
            if duplicate_coverage:
                raise RuntimeError("duplicate coverage keys across datasets")
            all_coverage.update(dataset_coverage)
            dataset_metadata.append(
                {
                    "dataset": dataset["dataset_key"],
                    "release": dataset["release"],
                    "assembly": dataset["assembly"],
                    "callset": dataset["callset"],
                    "subset": dataset["subset"],
                    "faf95_mode": dataset["faf95_mode"],
                    "sites_hail_uri": dataset["sites_hail_uri"],
                    "coverage_hail_uri": dataset["coverage_hail_uri"],
                    "sites_source_identity": source_identities[dataset["dataset_key"]]["sites"],
                    "coverage_source_identity": source_identities[dataset["dataset_key"]]["coverage"],
                    "variant_records": sum(len(v) for v in dataset_variants.values()),
                    "coverage_positions": len(dataset_coverage),
                }
            )
    finally:
        hl.stop()

    now = datetime.now(timezone.utc).isoformat()
    extraction_log = []
    for dataset in manifest["datasets"]:
        source_identity = source_identities[dataset["dataset_key"]]["sites"]
        for target in manifest["targets"]:
            interval = target["intervals"][dataset["assembly"]]
            extraction_log.append(
                {
                    "dataset": dataset["dataset_key"],
                    "status": "ok",
                    "chrom": interval["chrom"],
                    "gene": target["gene"],
                    "start": int(interval["start"]) - int(manifest["region_padding_bp"]),
                    "end": int(interval["end"]) + int(manifest["region_padding_bp"]),
                    "source": dataset["sites_hail_uri"],
                    "source_identity": source_identity,
                }
            )

    variant_payload = {
        "metadata": {
            "schema_version": 2,
            "source": "official gnomAD Hail Tables",
            "builder": "scripts/refresh_gnomad_panel_snapshot.py",
            "updated_utc": now,
            "manifest_path": str(manifest_path.relative_to(ROOT)).replace("\\", "/"),
            "manifest_sha256": manifest_sha256,
            "panel_id": manifest["panel_id"],
            "automatic_release_activation": False,
            "classification_policies": manifest["classification_policies"],
            "regions": _regions_metadata(manifest),
            "region_padding_bp": manifest["region_padding_bp"],
            "datasets": dataset_metadata,
            "extraction_log": extraction_log,
            "variants_count": sum(len(values) for values in all_variants.values()),
            "unique_ids": len(all_variants),
            "records_sha256": _canonical_sha256(all_variants),
            "v2_faf95": {
                "scope": (
                    "non-cancer subset; extracted populations: "
                    + ", ".join(_scored_ancestries(manifest))
                ),
                "method": "published gnomAD Hail Table FAF95",
                "raw_af_fallback_allowed": False,
            },
            "v3_faf95": {
                "scope": (
                    "non-cancer subset; extracted populations: "
                    + ", ".join(_scored_ancestries(manifest))
                ),
                "source_fields": "official Hail Table non-cancer AC and AN",
                "method": "hail.experimental.filtering_allele_frequency",
                "confidence": FAF_CONFIDENCE,
                "raw_af_fallback_allowed": False,
                "api_faf95_used": False,
            },
        },
        "variants": all_variants,
    }
    coverage_payload = {
        "metadata": {
            "schema_version": 2,
            "source": "official gnomAD Hail coverage Tables",
            "builder": "scripts/refresh_gnomad_panel_snapshot.py",
            "updated_utc": now,
            "manifest_sha256": manifest_sha256,
            "panel_id": manifest["panel_id"],
            "classification_policies": manifest["classification_policies"],
            "datasets": dataset_metadata,
            "records": len(all_coverage),
            "records_sha256": _canonical_sha256(all_coverage),
        },
        "coverage_by_position": all_coverage,
    }
    _validate_built_payloads(variant_payload, coverage_payload, manifest)
    _atomic_json(variants_path, variant_payload)
    _atomic_json(coverage_path, coverage_payload)

    print(f"frequency snapshot: {variants_path}")
    print(f"  records: {variant_payload['metadata']['variants_count']}")
    print(f"  sha256: {hashlib.sha256(variants_path.read_bytes()).hexdigest()}")
    print(f"coverage snapshot: {coverage_path}")
    print(f"  positions: {coverage_payload['metadata']['records']}")
    print(f"  sha256: {hashlib.sha256(coverage_path.read_bytes()).hexdigest()}")
    if own_staging and not keep_staging:
        shutil.rmtree(staging, ignore_errors=True)


def validate_existing(manifest_path: Path, variants_path: Path, coverage_path: Path) -> None:
    manifest = _load_manifest(manifest_path)
    variants = json.loads(variants_path.read_text(encoding="utf-8"))
    coverage = json.loads(coverage_path.read_text(encoding="utf-8"))
    _validate_built_payloads(variants, coverage, manifest)
    expected_manifest_hash = _canonical_sha256(manifest)
    for payload, label in ((variants, "frequency"), (coverage, "coverage")):
        if payload["metadata"].get("manifest_sha256") != expected_manifest_hash:
            raise RuntimeError(f"{label} snapshot was built from a different manifest")
    print("gnomAD panel snapshots are valid")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command", choices=("refresh", "check-updates", "validate")
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--variants-output", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--coverage-output", type=Path, default=DEFAULT_COVERAGE)
    parser.add_argument("--staging-dir", type=Path)
    parser.add_argument("--keep-staging", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    if args.command == "check-updates":
        report = check_updates(args.manifest)
        rendered = json.dumps(report, ensure_ascii=False, indent=2)
        print(rendered)
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(rendered + "\n", encoding="utf-8")
    elif args.command == "validate":
        validate_existing(args.manifest, args.variants_output, args.coverage_output)
    else:
        refresh(
            args.manifest,
            args.variants_output,
            args.coverage_output,
            staging_dir=args.staging_dir,
            keep_staging=args.keep_staging,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
