from backend.review.initiation_review import evaluate_initiation_review
from backend.review.rna_review import evaluate_rna_review
from backend.review.splice_ps1_review import evaluate_splice_ps1_review


def test_empty_review_recommendations_have_one_shape():
    expected = {
        "recommended",
        "priority",
        "title",
        "summary",
        "reasons",
        "what_to_test",
        "potential_branches",
        "limitations",
        "reference_source",
        "source_url",
        "is_evidence_criterion",
        "manual_review_prefill",
    }
    results = (
        evaluate_initiation_review("BRCA1", "missense"),
        evaluate_splice_ps1_review("BRCA1", "missense", 0.0),
        evaluate_rna_review("BRCA1", "missense", 0.0),
    )
    assert all(set(result) == expected for result in results)
