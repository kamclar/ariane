"""Small constructors for auditable ENIGMA Figure 1A decision paths."""

from __future__ import annotations

from typing import Any


FIGURE_1A_SOURCE = {
    "source_id": "enigma-v1.2-specifications",
    "label": "Complete ENIGMA BRCA1/2 VCEP v1.2 specification",
    "url": "https://cspec.genome.network/cspec/File/id/11e62fec-23b0-4a3e-b2df-751855301746/data",
    "location": "Figure 1A",
}


def step(node_id: str, question: str, result: str, observed: str) -> dict[str, str]:
    return {
        "node_id": node_id,
        "question": question,
        "result": result,
        "observed": observed,
    }


def figure1a_path(
    *,
    branch_id: str,
    criterion: str,
    outcome_node: str,
    steps: list[dict[str, str]],
) -> dict[str, Any]:
    figure_url = (
        "/static/enigma/figure-1a-missense-inframe.png"
        if branch_id == "missense-inframe"
        else "/static/enigma/figure-1a-synonymous-intronic.png"
    )
    return {
        "tree_id": "figure-1a",
        "tree_version": "ENIGMA VCEP 1.2.0",
        "branch_id": branch_id,
        "criterion": criterion,
        "outcome": "applied",
        "outcome_node": outcome_node,
        "steps": steps,
        "sources": [{**FIGURE_1A_SOURCE, "figure_url": figure_url}],
    }
