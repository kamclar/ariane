"""Manual-review triage DAG nodes."""

from __future__ import annotations

from dataclasses import dataclass

from backend.classification_dag.types import NodeResult


@dataclass(frozen=True)
class ReviewTriageNode:
    id: str = "review.manual_triage"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {"classification_inputs", "retained_evidence", "classification_without_reviews"}
    )
    provides: frozenset[str] = frozenset({"unvalidated_classification_result"})

    def evaluate(self, context, inputs) -> NodeResult:
        from backend.modules.initiation_review import evaluate_initiation_review
        from backend.modules.protein_ps1_review import evaluate_protein_ps1_review
        from backend.modules.rna_review import evaluate_rna_review
        from backend.modules.splice_ps1_review import evaluate_splice_ps1_review

        ci = inputs["classification_inputs"]
        retained = inputs["retained_evidence"]
        result = dict(inputs["classification_without_reviews"])
        if retained.metadata.get("ba1_terminal"):
            return NodeResult.succeeded({"unvalidated_classification_result": result})
        result["rna_review"] = evaluate_rna_review(
            gene=ci.gene,
            variant_type=ci.variant_type,
            spliceai_score=ci.spliceai_score,
            pvs1_result=(
                retained.metadata["pvs1_rna"]
                if retained.metadata["pvs1_rna"].get("applies")
                else retained.metadata["pvs1"]
            ),
            criteria=result["criteria"],
        )
        result["splice_ps1_review"] = evaluate_splice_ps1_review(
            gene=ci.gene,
            variant_type=ci.variant_type,
            spliceai_score=ci.spliceai_score,
            ps1_result=retained.metadata["ps1"],
        )
        result["protein_ps1_review"] = evaluate_protein_ps1_review(
            retained.metadata["ps1"],
            gene=ci.gene,
        )
        result["initiation_review"] = evaluate_initiation_review(
            gene=ci.gene,
            variant_type=ci.variant_type
        )
        return NodeResult.succeeded({"unvalidated_classification_result": result})
