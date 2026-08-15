from pathlib import Path


FRONTEND_HTML = Path(__file__).resolve().parents[1] / "frontend" / "index.html"
FRONTEND_JS = Path(__file__).resolve().parents[1] / "frontend" / "static" / "js" / "app.js"


def test_applied_criteria_immediately_follows_classification_header():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    classification = html.index("<!-- Classification badge -->")
    criteria = html.index("<!-- Criteria table -->")
    narrative = html.index("<!-- Narrative summary -->")

    assert classification < criteria < narrative
    assert html.count("<!-- Criteria table -->") == 1


def test_point_thermometer_is_present_for_every_classified_result():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    thermometer = html.index('class="point-thermometer"')
    badge = html.index('class="class-badge"', thermometer)
    thermometer_markup = html[thermometer:badge]

    assert 'point-thermometer-mixed' in thermometer_markup
    assert 'point-thermometer-combination' in thermometer_markup
    assert 'x-show="result?.mixed_evidence"' not in thermometer_markup
    assert 'x-if="result?.mixed_evidence"' not in thermometer_markup
    assert ':style="{ top: pointMeterPosition(result?.total_points) + \'%\' }"' in thermometer_markup
    assert ">+10<" in thermometer_markup
    assert ">+6<" in thermometer_markup
    assert "&minus;2" in thermometer_markup
    assert "&minus;7" in thermometer_markup


def test_variant_input_can_synchronise_the_explicit_gene_selector():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    assert '@input="syncGeneFromVariantInput()"' in html
    assert '@change="syncGeneFromVariantInput()"' in html
    assert 'x-show="geneAutoSwitchNotice"' in html
    assert "explicitGeneFromVariantInput()" in javascript
    assert "syncGeneFromVariantInput()" in javascript
    assert 'accession === "NM_007294"' in javascript
    assert 'accession === "NM_000059"' in javascript
