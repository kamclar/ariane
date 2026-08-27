from pathlib import Path


FRONTEND_HTML = Path(__file__).resolve().parents[1] / "frontend" / "index.html"
FRONTEND_JS = Path(__file__).resolve().parents[1] / "frontend" / "static" / "js" / "app.js"


def test_application_version_is_loaded_from_backend_and_visible_in_header_and_footer():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    assert 'id="ariane-version"' in html
    assert 'id="ariane-footer-version"' in html
    assert "this.setAppVersion(resources.version)" in javascript
    assert 'getElementById("ariane-version")' in javascript
    assert 'getElementById("ariane-footer-version")' in javascript
    assert "ARIANE v1.8.0" not in html


def test_applied_criteria_immediately_follows_classification_header():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    classification = html.index("<!-- Classification badge -->")
    criteria = html.index("<!-- Criteria table -->")
    narrative = html.index("<!-- Narrative summary -->")

    assert classification < criteria < narrative
    assert html.count("<!-- Criteria table -->") == 1


def test_variant_specific_not_applicable_criteria_are_compact_and_backend_driven():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    assert 'class="not-applicable-criteria-details"' in html
    assert "Not applicable to this variant" in html
    assert "Show codes and reasons" in html
    assert "result?.not_applicable_criteria" in html
    assert "item.reason" in html
    assert "Criteria that were not met" in html
    assert "result.not_applicable_criteria = this.sortCriterionList" in javascript
    assert "currentNotUsedCriteria()" not in html
    assert "De novo evidence is not calibrated" not in html


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


def test_variant_specific_clinical_annotation_is_visible_but_not_scored():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    assert "result?.clinical_annotations" in html
    assert 'class="clinical-annotation-strip"' in html
    assert "annotation.affects_classification" not in html
    assert "annotation.source_url" in html
    assert "annotation.publications" in html


def test_external_comparison_remains_visible_for_not_found_or_failed_sources():
    html = FRONTEND_HTML.read_text(encoding="utf-8")

    assert 'class="external-section"' in html
    assert "result?.external?.clinvar_message" in html
    assert "result?.external?.clinvar_error" in html
    assert "result?.external?.clingen_message" in html
    assert "result?.external?.clingen_error" in html
    assert "result?.external?.erepo_evidence_codes" in html


def test_variant_input_can_synchronise_the_explicit_gene_selector():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    assert '@input="syncGeneFromVariantInput()"' in html
    assert '@change="syncGeneFromVariantInput()"' in html
    assert 'x-show="geneAutoSwitchNotice"' in html
    assert "explicitGeneFromVariantInput()" in javascript
    assert "syncGeneFromVariantInput()" in javascript
    assert "this.configuredGenes.find" in javascript
    assert "item.reference_transcript" in javascript


def test_unreviewed_splice_ps1_pilot_is_not_exposed_in_production_ui():
    project_root = FRONTEND_HTML.parents[1]
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")
    backend_main = (project_root / "backend" / "main.py").read_text(encoding="utf-8")

    assert not (project_root / "backend" / "data" / "splice_ps1_reference_set.json").exists()
    assert not (project_root / "backend" / "modules" / "splice_ps1_reference.py").exists()
    assert "splice_ps1_reference_candidates" not in backend_main
    assert "PS1_SPLICE_REFERENCE_PATH" not in backend_main
    assert "splice_ps1_reference_candidates" not in javascript
    assert "splicePs1ReferenceSet" not in javascript
    assert "applySplicePs1Candidate(item)" not in javascript
    assert "Candidate reference variant from pilot set" not in html
    assert "Candidate reference variant from ENIGMA ST2" in html
    assert "Why PS1 cannot be assigned from ST2 alone" in html
    assert "reference P/LP classification assigned using VCEP specifications" in html
    assert "precisely the same splice event" in html
    assert "similar to or stronger than the reference prediction" in html
    assert "positions of both variants within the donor or acceptor motif" in html
    assert "baseline PP3 or PVS1 result" in html
    assert "concurrent protein-level consequence" in html
    assert "applySplicePs1CandidateFacts(item)" in html
    safe_prefill_start = javascript.index("applySplicePs1CandidateFacts(item)")
    safe_prefill_end = javascript.index("async evaluateManualEvidence()", safe_prefill_start)
    safe_prefill = javascript[safe_prefill_start:safe_prefill_end]
    assert "curated_strength" not in safe_prefill
    assert "same_splice_event_confirmed" not in safe_prefill
    assert "prediction_strength_comparison" not in safe_prefill
    assert "PS1_SPLICE" in html
    assert "ARIANE does not infer this strength automatically" in html


def test_decision_path_stays_inside_applied_criteria_details():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    criteria = html.index("<!-- Criteria table -->")
    narrative = html.index("<!-- Narrative summary -->")
    fragment = html[criteria:narrative]

    assert "Show decision path" in fragment
    assert "Open full decision tree" in fragment
    assert "c.decision_path" in fragment
    assert 'class="criterion-composite-row"' in fragment
    assert '<td colspan="4">' in fragment
    assert 'class="criterion-summary-grid"' in fragment


def test_embedded_decision_path_uses_full_width_without_horizontal_scrollbar():
    css = (FRONTEND_HTML.parent / "static" / "css" / "style.css").read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    graph_start = css.index(".decision-path-graph {")
    svg_start = css.index(".decision-path-svg {", graph_start)
    graph_rule = css[graph_start:svg_start]
    svg_rule = css[svg_start:css.index(".decision-path-svg .path-edge-active", svg_start)]

    assert "overflow: hidden" in graph_rule
    assert "overflow-x: auto" not in graph_rule
    assert "width: 100%" in svg_rule
    assert "min-height: 250px" in svg_rule
    assert "const columns = Math.min(3, nodes.length)" in javascript
    assert "const verticalGap = 88" in javascript


def test_rules_navigation_labels_are_not_left_to_global_button_colours():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    css = (FRONTEND_HTML.parent / "static" / "css" / "style.css").read_text(encoding="utf-8")

    assert ">Decision trees</button>" in html
    assert ">Original figures</button>" in html
    assert ">Tables</button>" in html
    assert ">Sources and versions</button>" in html
    assert ".rules-subtabs button" in css
    assert "color: var(--color-text);" in css


def test_decision_tree_is_responsive_without_horizontal_scrollbar():
    css = (FRONTEND_HTML.parent / "static" / "css" / "style.css").read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    canvas_rule = css[css.index(".decision-tree-canvas"):css.index(".decision-tree-svg")]
    svg_rule = css[css.index(".decision-tree-svg {"):css.index(".decision-tree-svg .tree-connector")]

    assert "overflow: hidden" in canvas_rule
    assert "overflow-x: auto" not in canvas_rule
    assert "width: 100%" in svg_rule
    assert "width: max(" not in svg_rule
    assert 'preserveAspectRatio="xMidYMin meet"' in javascript
    assert "roundedOrthogonalPath(points, radius = 10)" in javascript
    assert "treeEdgeLabelPosition(points, edge)" in javascript
    assert 'markerUnits="userSpaceOnUse"' in javascript


def test_manual_review_ui_has_no_strength_override_control():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    assert "Reviewer-selected strength" not in html
    assert "reviewer override" not in html
    assert "override_strength" not in html
    assert "override_strength" not in javascript
    assert "reviewer_selected_strength" not in javascript
    assert "Criterion strength is calculated by the backend" in html


def test_frontend_does_not_calculate_manual_criterion_strength_or_eligibility():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    assert "suggestedManualStrength" not in javascript
    assert "numberOrNull" not in javascript
    assert "rule_derived_strength" not in javascript
    assert "Math.log10(2.08)" not in javascript
    assert "the submitted evidence does not meet an ENIGMA" not in javascript
    assert "Assessor and assessment date are required" not in javascript
    assert "evidence notes and at least one reference are required" not in javascript
    assert "Criterion strength is calculated by the backend" in html
    assert "criterion.suggested_strength" in html
    assert "criterion.selected_strength" in html


def test_protein_ps1_reference_facts_are_requested_from_backend():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    assert "Complete reference facts" in html
    assert "resolveProteinPs1Reference(item)" in html
    assert 'fetch("/api/manual-evidence/resolve-ps1-reference"' in javascript
    assert "resolved.reference.p_notation" in javascript
    assert "resolved.assessed.spliceai_score" in javascript
    assert "resolved.reference.spliceai_score" in javascript
    assert "resolved.classification_verification" in javascript
    assert "Using the ClinVar aggregate conclusion itself" in html
    assert "would be circular" in html
    resolver_start = javascript.index("async resolveProteinPs1Reference(item)")
    resolver_end = javascript.index("applySplicePs1CandidateFacts(item)", resolver_start)
    assert "item.references =" not in javascript[resolver_start:resolver_end]


def test_manual_review_ui_collects_required_enigma_stipulations():
    html = FRONTEND_HTML.read_text(encoding="utf-8")
    javascript = FRONTEND_JS.read_text(encoding="utf-8")

    for field in (
        "case_control_country_matched",
        "case_control_ethnicity_matched",
        "cooccurring_variant_classification_basis",
        "vua_benign_population_review",
        "very_strong_effect_basis",
    ):
        assert field in html
        assert field in javascript

    assert "Classified using VCEP specifications" in html
    assert "Does not meet a benign population evidence code" in html
    assert "Effect evidence relevant to PP1 Very Strong" in html
