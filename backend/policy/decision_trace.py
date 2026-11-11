"""Constructors for auditable ENIGMA Figure 1A decision paths."""

from __future__ import annotations

from typing import Any

from backend.policy.gene import policy_version, vcep_document_url


def step(node_id: str, question: str, result: str, observed: str) -> dict[str, str]:
    return {
        "node_id": node_id,
        "question": question,
        "result": result,
        "observed": observed,
    }


def figure1a_path(
    *,
    gene: str,
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
        "tree_version": f"ENIGMA VCEP {policy_version(gene)}",
        "branch_id": branch_id,
        "criterion": criterion,
        "outcome": "applied",
        "outcome_node": outcome_node,
        "steps": steps,
        "sources": [
            {
                "source_id": "enigma-v1.2-specifications",
                "label": "Complete ENIGMA BRCA1/2 VCEP v1.2 specification",
                "url": vcep_document_url("specifications", gene),
                "location": "Figure 1A",
                "figure_url": figure_url,
            }
        ],
    }
