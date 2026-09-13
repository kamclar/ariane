import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest


FRONTEND_JS_DIR = Path(__file__).resolve().parents[1] / "frontend" / "static" / "js"
APPLICATION_SCRIPTS = (
    "api.js",
    "composition.js",
    "core.js",
    "rules.js",
    "graphs.js",
    "classification.js",
    "formatters.js",
    "manual-review.js",
    "manual-review-api.js",
    "manual-review-ps1.js",
    "manual-review-persistence.js",
    "batch.js",
    "app.js",
)


def test_browser_shell_loads_every_application_script_in_dependency_order():
    shell = (FRONTEND_JS_DIR.parents[1] / "index.html").read_text(encoding="utf-8")
    positions = [shell.index(f"/static/js/{filename}") for filename in APPLICATION_SCRIPTS]

    assert positions == sorted(positions)
    loaded_names = set(re.findall(r'/static/js/([^"?]+)\?v=', shell))
    assert loaded_names == set(APPLICATION_SCRIPTS)


def run_javascript(body: str) -> dict:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is required to execute frontend behavior tests")
    prelude = f"""
const fs = require("fs");
const vm = require("vm");
global.window = {{
    ArianeFrontend: {{}},
    location: {{ href: "http://localhost/" }},
    requestAnimationFrame: callback => callback(),
}};
global.document = {{
    getElementById: () => null,
    createElement: () => ({{ click() {{}} }}),
}};
global.navigator = {{ userAgent: "frontend-test" }};
for (const filename of {json.dumps(APPLICATION_SCRIPTS)}) {{
    vm.runInThisContext(
        fs.readFileSync(`${{process.argv[1]}}/${{filename}}`, "utf8"),
        {{ filename }},
    );
}}
"""
    completed = subprocess.run(
        [node, "-e", prelude + body, str(FRONTEND_JS_DIR)],
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(completed.stdout)


def test_frontend_composition_rejects_duplicate_state_or_method_keys():
    result = run_javascript(
        """
let message = "";
try {
    window.ArianeFrontend.composeStrict([
        { name: "first", value: { shared: 1 } },
        { name: "second", value: { shared: 2 } },
    ]);
} catch (error) {
    message = error.message;
}
process.stdout.write(JSON.stringify({ message }));
"""
    )

    assert result["message"] == (
        "ARIANE frontend key collision: shared is defined by first and second"
    )


def test_complete_alpine_component_has_unique_expected_feature_methods():
    result = run_javascript(
        """
const component = ariane();
process.stdout.write(JSON.stringify({
    keyCount: Object.keys(component).length,
    methods: {
        classify: typeof component.classify,
        prefill: typeof component.prefillManualReviewFromResult,
        evaluate: typeof component.evaluateManualEvidence,
        resolvePs1: typeof component.resolveProteinPs1Reference,
        save: typeof component.saveManualReviewDraft,
    },
}));
"""
    )

    assert result["keyCount"] >= 160
    assert set(result["methods"].values()) == {"function"}


def test_classification_prefills_but_does_not_enable_recommended_manual_review():
    result = run_javascript(
        """
(async () => {
    const component = ariane();
    component.resetManualItems();
    component.configuredGenes = [{ symbol: "BRCA1" }];
    component.gene = "BRCA1";
    component.c_notation = "c.4185G>A";
    window.fetch = async (path) => {
        if (path === "/ui-api/classify") return {
            ok: true,
            json: async () => ({
                gene: "BRCA1",
                c_notation: "c.4185G>A",
                p_notation: "p.(Gln1395=)",
                criteria: [],
                rna_review: {
                    recommended: true,
                    title: "Review RNA evidence",
                    source_url: "https://example.test/evidence",
                    manual_review_prefill: {
                        source_citation: "ST2",
                        source_references: ["PMID:1"],
                        table4_context: "Curated RNA evidence requires review.",
                    },
                },
            }),
        };
        if (path === "/ui-api/manual-evidence/status") return {
            ok: true,
            json: async () => ({ criteria: [] }),
        };
        throw new Error(`Unexpected request: ${path}`);
    };
    await component.classify();
    await new Promise(resolve => setImmediate(resolve));
    const item = component.manualItems.find(value => value.code === "PVS1_RNA");
    process.stdout.write(JSON.stringify({
        classificationLoaded: component.result?.gene === "BRCA1",
        reviewOpen: component.manualReviewOpen,
        enabled: item.enabled,
        citation: item.evidence.source_citation,
        notes: item.notes,
        references: item.references,
    }));
})().catch(error => { console.error(error); process.exit(1); });
"""
    )

    assert result == {
        "classificationLoaded": True,
        "reviewOpen": True,
        "enabled": False,
        "citation": "ST2",
        "notes": "Curated RNA evidence requires review.",
        "references": "ST2\nPMID:1\nhttps://example.test/evidence",
    }


def test_manual_evidence_evaluation_uses_the_prefilled_form_payload():
    result = run_javascript(
        """
(async () => {
    const component = ariane();
    component.resetManualItems();
    component.result = {
        gene: "BRCA1",
        c_notation: "c.4185G>A",
        p_notation: "p.(Gln1395=)",
        criteria: [{ name: "PP4", applies: true }],
    };
    component.manualAssessor = "Reviewer 1";
    const pp1 = component.manualItems.find(value => value.code === "PP1");
    pp1.enabled = true;
    pp1.evidence.likelihood_ratio = 2.08;
    pp1.notes = "Segregation review";
    pp1.references = "PMID:1\\nPMID:2";
    let submitted = null;
    window.fetch = async (path, options) => {
        if (path !== "/ui-api/manual-evidence/evaluate") {
            throw new Error(`Unexpected request: ${path}`);
        }
        submitted = JSON.parse(options.body);
        return {
            ok: true,
            json: async () => ({ predicted_class: 4, total_points: 6 }),
        };
    };
    await component.evaluateManualEvidence();
    const submittedPp1 = submitted.manual_criteria.find(value => value.code === "PP1");
    process.stdout.write(JSON.stringify({
        resultClass: component.manualResult.predicted_class,
        assessor: submitted.assessor,
        variant: submitted.variant_context,
        criterion: submittedPp1,
    }));
})().catch(error => { console.error(error); process.exit(1); });
"""
    )

    assert result["resultClass"] == 4
    assert result["assessor"] == "Reviewer 1"
    assert result["variant"] == {
        "gene": "BRCA1",
        "c_notation": "c.4185G>A",
        "p_notation": "p.(Gln1395=)",
    }
    assert result["criterion"]["enabled"] is True
    assert result["criterion"]["evidence"]["likelihood_ratio"] == 2.08
    assert result["criterion"]["references"] == ["PMID:1", "PMID:2"]
