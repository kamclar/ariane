"""Small test adapter for the public synchronous classification DAG contract.

The adapter contains no clinical rules. It only turns explicit test evidence
into ``ClassificationInputs`` and returns the DAG result dictionary.
"""

from typing import Any, Mapping

from backend.classification_dag import ClassificationInputs, execute_classification
from backend.population_frequency.policy import classification_policy_for_gene


def classify_with_dag(
    gene: str,
    variant_type: str,
    p_notation: str,
    c_notation: str,
    spliceai_score: float | None = None,
    bayesdel_score: float | None = None,
    gnomad_data: Mapping[str, Any] | None = None,
    table9_result: Mapping[str, Any] | None = None,
    pp4_bp5_result: Mapping[str, Any] | None = None,
    ps1_result: Mapping[str, Any] | None = None,
    exon_cnv_result: Mapping[str, Any] | None = None,
    residue_info: Mapping[str, Any] | None = None,
    dup_type: str = "Unknown",
) -> dict[str, Any]:
    inputs = ClassificationInputs(
        gene=gene,
        variant_type=variant_type,
        p_notation=p_notation,
        c_notation=c_notation,
        spliceai_score=spliceai_score,
        bayesdel_score=bayesdel_score,
        gnomad_data=gnomad_data,
        frequency_policy=classification_policy_for_gene(gene),
        table9_result=table9_result,
        pp4_bp5_result=pp4_bp5_result,
        ps1_result=ps1_result,
        exon_cnv_result=exon_cnv_result,
        residue_info=residue_info,
        dup_type=dup_type,
    )
    return execute_classification(inputs, mode="dag").result
