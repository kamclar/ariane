#!/usr/bin/env python3
"""Differentially audit sequence-derived HGVS against installed snapshots."""
from __future__ import annotations

import argparse
import json
from collections import Counter

from backend.lookups.indels import load_indel_snapshot
from backend.lookups.precomputed import load_classification_snapshot_index
from backend.modules.hgvs import normalize_protein_notation, protein_notations_compatible
from backend.modules.hgvs_engine import VariantNormalizationError, derive_protein_consequence


def records(include_indels: bool, only_indels: bool):
    if not only_indels:
        for key, record in load_classification_snapshot_index().items():
            gene, c_notation = key.split(":", 1)
            yield "coding_snv_snapshot", gene, c_notation, record.get("p_notation", "")
    if include_indels:
        indels, _aliases = load_indel_snapshot()
        for record in indels.values():
            yield (
                "normalized_indel_snapshot",
                record.get("gene", ""),
                record.get("canonical_c_notation") or record.get("c_notation", ""),
                record.get("p_notation", ""),
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-indels", action="store_true")
    parser.add_argument("--only-indels", action="store_true")
    parser.add_argument("--max-examples", type=int, default=25)
    args = parser.parse_args()
    counts = Counter()
    examples: dict[str, list[dict]] = {
        "conflict": [], "unsupported": [], "unknown_snapshot": []
    }
    for source, gene, c_notation, snapshot_p in records(
        args.include_indels or args.only_indels, args.only_indels
    ):
        counts["total"] += 1
        normalized_snapshot_p = normalize_protein_notation(str(snapshot_p or ""))
        if not normalized_snapshot_p or normalized_snapshot_p == "p.?":
            counts["unknown_snapshot"] += 1
            if len(examples["unknown_snapshot"]) < args.max_examples:
                examples["unknown_snapshot"].append(
                    {"source": source, "gene": gene, "c": c_notation, "snapshot_p": snapshot_p}
                )
            continue
        try:
            derived = derive_protein_consequence(gene, c_notation)
        except VariantNormalizationError as exc:
            counts[f"unsupported:{exc.code}"] += 1
            counts["unsupported"] += 1
            if len(examples["unsupported"]) < args.max_examples:
                examples["unsupported"].append(
                    {
                        "source": source, "gene": gene, "c": c_notation,
                        "snapshot_p": normalized_snapshot_p, "code": exc.code,
                        "reason": str(exc),
                    }
                )
            continue
        if derived.p_notation == "p.?":
            counts["derived_unknown"] += 1
            continue
        if protein_notations_compatible(normalized_snapshot_p, derived.p_notation):
            counts["compatible"] += 1
        else:
            counts["conflict"] += 1
            if len(examples["conflict"]) < args.max_examples:
                examples["conflict"].append(
                    {
                        "source": source, "gene": gene, "c": c_notation,
                        "snapshot_p": normalized_snapshot_p,
                        "derived_c": derived.canonical_c_notation,
                        "derived_p": derived.p_notation,
                    }
                )
    print(json.dumps({"counts": dict(sorted(counts.items())), "examples": examples}, indent=2))
    return 1 if counts["conflict"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
