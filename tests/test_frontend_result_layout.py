from pathlib import Path


FRONTEND_HTML = Path(__file__).resolve().parents[1] / "frontend" / "index.html"


def test_applied_criteria_immediately_follows_classification_header():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    classification = html.index("<!-- Classification badge -->")
    criteria = html.index("<!-- Criteria table -->")
    narrative = html.index("<!-- Narrative summary -->")

    assert classification < criteria < narrative
    assert html.count("<!-- Criteria table -->") == 1
