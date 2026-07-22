function ariane() {
    return {
        // shared
        mode: "single",

        // single mode
        gene: "BRCA1",
        c_notation: "",
        p_notation: "",
        dup_type: "Unknown",
        loading: false,
        progress: 0,
        progressText: "",
        error: "",
        result: null,

        // reference material and manually reviewed evidence
        manualDefinitions: {},
        resourceLinks: [],
        resourceError: "",
        splicePs1ReferenceSet: { candidates: [] },
        manualItems: [],
        manualAssessor: "",
        manualAssessedAt: new Date().toISOString().slice(0, 10),
        manualLoading: false,
        manualError: "",
        manualResult: null,

        // batch mode
        batchText: "",
        batchParsed: [],
        batchParseError: "",
        batchRunning: false,
        batchDone: 0,
        batchTotal: 0,
        batchResults: [],

        async init() {
            this.resetManualItems();
            try {
                const response = await fetch("/api/resources");
                if (response.ok) {
                    const resources = await response.json();
                    this.manualDefinitions = resources.manual_criteria || {};
                    this.resourceLinks = resources.links || [];
                    this.splicePs1ReferenceSet = resources.splice_ps1_reference_candidates || { candidates: [] };
                } else {
                    this.resourceError = `Reference materials could not be loaded (HTTP ${response.status}). Classification remains available, but manual-review guidance is incomplete.`;
                }
            } catch (e) {
                this.resourceError = `Reference materials could not be loaded: ${e?.message || e}. Classification remains available, but manual-review guidance is incomplete.`;
            }
        },

        resetManualItems() {
            this.manualItems = ["PS4", "PM3", "PP1", "BS2", "BS4", "PVS1_RNA", "BP7_RNA", "PVS1_INIT", "PS1_SPLICE"].map(code => ({
                code,
                enabled: false,
                evidence: {},
                override_strength: "",
                notes: "",
                references: "",
            }));
            this.manualResult = null;
            this.manualError = "";
        },

        manualDefinition(code) {
            return this.manualDefinitions[code] || {};
        },

        numberOrNull(value) {
            if (value === "" || value === null || value === undefined) return null;
            const parsed = Number(value);
            return Number.isFinite(parsed) ? parsed : null;
        },

        splicePs1CandidatesForCurrentGene() {
            const currentGene = this.result?.gene || this.gene;
            return (this.splicePs1ReferenceSet.candidates || [])
                .filter(candidate => candidate.gene === currentGene);
        },

        splicePs1CandidateLabel(candidate) {
            const protein = candidate.p_notation ? ` ${candidate.p_notation}` : "";
            return `${candidate.reference_variant}${protein} - ${candidate.classification}; ${candidate.splice_event_label}`;
        },

        batchReviewLabel(result) {
            if (!result) return "";
            const labels = [];
            if (result.rna_review && result.rna_review.recommended) {
                labels.push(`RNA ${result.rna_review.priority || ""}`.trim());
            }
            if (result.splice_ps1_review && result.splice_ps1_review.recommended) {
                labels.push(`PS1_SPLICE ${result.splice_ps1_review.priority || ""}`.trim());
            }
            if (result.initiation_review && result.initiation_review.recommended) {
                labels.push(`PVS1_INIT ${result.initiation_review.priority || ""}`.trim());
            }
            return labels.join("; ");
        },

        prefillManualReviewFromResult() {
            if (!this.result?.initiation_review?.recommended) return;
            const item = this.manualItems.find(value => value.code === "PVS1_INIT");
            if (!item || Object.keys(item.evidence || {}).length > 0) return;

            item.evidence.met1_loss_confirmed = true;
            item.evidence.initiation_flowchart_rationale =
                "Met1/start-loss variant flagged by ARIANE. Complete the ENIGMA initiation-codon flowchart review: alternative start assessment, upstream P/LP evidence, expected N-terminal impact, and curated PVS1_INIT strength.";
        },

        applySplicePs1Candidate(item) {
            const evidence = item.evidence || {};
            const candidate = (this.splicePs1ReferenceSet.candidates || [])
                .find(value => value.key === evidence.splice_ps1_candidate_key);
            if (!candidate) return;

            item.evidence.reference_variant = `${candidate.gene} ${candidate.reference_variant}`;
            item.evidence.reference_classification = candidate.classification;
            item.evidence.reference_classification_source =
                `${candidate.classification_source}; ${candidate.source_label}`;
            item.evidence.reference_splice_event = candidate.splice_event_label;
            item.evidence.reference_event_type = candidate.event_type;
            item.evidence.reference_assay_result_category = candidate.assay_result_category;
            item.evidence.reference_variant_context = candidate.variant_context;
            item.evidence.reference_curation_status = candidate.curation_status;
            item.evidence.reference_curation_note = candidate.curation_note;
            item.evidence.ps1_splice_prefill_suggestion = candidate.prefill_strength_suggestion;
            item.evidence.ps1_splice_prefill_basis = candidate.prefill_strength_basis;

            if (!item.evidence.curated_strength && candidate.prefill_strength_suggestion) {
                item.evidence.curated_strength = candidate.prefill_strength_suggestion;
            }

            if (!item.references.trim()) {
                item.references = [
                    candidate.source_label,
                    candidate.source_url,
                    "Pilot splice_ps1_reference_set.json; requires expert review before use",
                ].filter(Boolean).join("\n");
            }
        },

        suggestedManualStrength(item) {
            const evidence = item.evidence || {};
            if (item.code === "PVS1_RNA") {
                const hasRequired = evidence.assay_scope === "mrna_only"
                    && evidence.rna_conclusion === "damaging"
                    && ["absent_or_minimal", "reduced", "uncertain"].includes(evidence.functional_transcript_remaining)
                    && ["yes", "no", "not_applicable"].includes(evidence.nmd_assessed)
                    && Boolean((evidence.transcript_accession || "").trim())
                    && Boolean((evidence.tissue_or_cell_type || "").trim());
                return hasRequired ? (evidence.curated_strength || "") : "";
            }
            if (item.code === "BP7_RNA") {
                const hasRequired = evidence.assay_scope === "mrna_only"
                    && evidence.rna_conclusion === "no_damaging_effect"
                    && evidence.bp7_rna_eligible === true
                    && ["yes", "no", "not_applicable"].includes(evidence.nmd_assessed)
                    && Boolean((evidence.transcript_accession || "").trim())
                    && Boolean((evidence.tissue_or_cell_type || "").trim());
                return hasRequired ? "Strong" : "";
            }
            if (item.code === "PVS1_INIT") {
                const hasRequired = evidence.met1_loss_confirmed === true
                    && ["yes", "no"].includes(evidence.alternative_start_assessed)
                    && ["yes", "no", "not_applicable"].includes(evidence.upstream_pathogenic_evidence)
                    && ["yes", "no", "uncertain"].includes(evidence.functional_domain_impact)
                    && Boolean((evidence.nearest_alternative_start || "").trim())
                    && Boolean((evidence.initiation_flowchart_rationale || "").trim());
                return hasRequired ? (evidence.curated_strength || "") : "";
            }
            if (item.code === "PS1_SPLICE") {
                const hasRequired = Boolean((evidence.reference_variant || "").trim())
                    && ["Pathogenic", "Likely Pathogenic"].includes(evidence.reference_classification)
                    && Boolean((evidence.reference_classification_source || "").trim())
                    && evidence.same_splice_event_confirmed === true
                    && Boolean((evidence.vua_splice_event || "").trim())
                    && Boolean((evidence.reference_splice_event || "").trim())
                    && ["similar", "stronger"].includes(evidence.prediction_strength_comparison)
                    && Boolean((evidence.ps1_splice_rationale || "").trim());
                return hasRequired ? (evidence.curated_strength || "") : "";
            }
            if (item.code === "PS4") {
                const p = this.numberOrNull(evidence.p_value);
                const odds = this.numberOrNull(evidence.odds_ratio);
                const lower = this.numberOrNull(evidence.lower_ci);
                return p !== null && odds !== null && lower !== null
                    && p <= 0.05 && odds >= 4 && lower > 2 ? "Strong" : "";
            }
            if (item.code === "PM3" || item.code === "BS2") {
                const points = this.numberOrNull(evidence.evidence_points);
                if (points === null || points < 1) return "";
                if (points >= 4) return "Strong";
                if (points >= 2) return "Moderate";
                return "Supporting";
            }
            const lr = this.numberOrNull(evidence.likelihood_ratio);
            if (lr === null || lr < 0) return "";
            if (item.code === "PP1") {
                if (lr >= 350) return "Very Strong";
                if (lr >= 18.7) return "Strong";
                if (lr >= 4.3) return "Moderate";
                if (lr >= 2.08) return "Supporting";
            }
            if (item.code === "BS4") {
                if (lr <= 0.00285) return "Very Strong";
                if (lr <= 0.05) return "Strong";
                if (lr <= 0.23) return "Moderate";
                if (lr <= 0.48) return "Supporting";
            }
            return "";
        },

        async evaluateManualEvidence() {
            this.manualError = "";
            this.manualResult = null;
            if (!this.result) {
                this.manualError = "Classify a variant first.";
                return;
            }
            if (!this.manualAssessor.trim() || !this.manualAssessedAt) {
                this.manualError = "Assessor and assessment date are required.";
                return;
            }
            const enabled = this.manualItems.filter(item => item.enabled);
            if (enabled.length === 0) {
                this.manualError = "Select at least one manually reviewed criterion.";
                return;
            }
            for (const item of enabled) {
                if (!this.suggestedManualStrength(item) && !item.override_strength) {
                    this.manualError = `${item.code}: enter evidence meeting a threshold or select an allowed reviewer strength.`;
                    return;
                }
                if (!item.notes.trim() || !item.references.trim()) {
                    this.manualError = `${item.code}: evidence notes and at least one reference are required.`;
                    return;
                }
            }

            this.manualLoading = true;
            try {
                const response = await fetch("/api/manual-evidence/evaluate", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        base_criteria: this.result.criteria,
                        manual_criteria: this.manualItems.map(item => ({
                            code: item.code,
                            enabled: item.enabled,
                            evidence: item.evidence,
                            override_strength: item.override_strength || null,
                            notes: item.notes,
                            references: item.references
                                .split(/\r?\n/)
                                .map(value => value.trim())
                                .filter(Boolean),
                        })),
                        assessor: this.manualAssessor.trim(),
                        assessed_at: this.manualAssessedAt,
                    }),
                });
                if (!response.ok) {
                    const error = await response.json().catch(() => ({}));
                    this.manualError = this.formatApiError(error, response.status);
                    return;
                }
                this.manualResult = await response.json();
            } catch (e) {
                this.manualError = "Network error - amended result could not be calculated.";
            } finally {
                this.manualLoading = false;
            }
        },

        downloadManualAuditJson() {
            if (!this.manualResult || !this.result) return;
            const record = {
                schema_version: "1.0",
                exported_at: new Date().toISOString(),
                variant: {
                    gene: this.result.gene,
                    c_notation: this.result.c_notation,
                    p_notation: this.result.p_notation,
                },
                module1_result: {
                    predicted_class: this.result.predicted_class,
                    predicted_label: this.result.predicted_label,
                    total_points: this.result.total_points,
                    criteria: this.result.criteria,
                },
                amended_working_result: this.manualResult,
                submitted_manual_evidence: this.manualItems
                    .filter(item => item.enabled)
                    .map(item => ({
                        code: item.code,
                        evidence: item.evidence,
                        reviewer_selected_strength: item.override_strength || null,
                        notes: item.notes,
                        references: item.references
                            .split(/\r?\n/)
                            .map(value => value.trim())
                            .filter(Boolean),
                    })),
                disclaimer: "Audit support only; not a standalone clinical classification.",
            };
            const blob = new Blob(
                [JSON.stringify(record, null, 2)],
                { type: "application/json;charset=utf-8" }
            );
            const url = URL.createObjectURL(blob);
            const link = document.createElement("a");
            link.href = url;
            link.download = `${this.result.gene}_${this.result.c_notation.replaceAll(/[^\w.-]/g, "_")}_manual_evidence.json`;
            link.click();
            URL.revokeObjectURL(url);
        },

        logClientValidation(error, form = "single", input = null) {
            const submittedInput = input || {
                gene: this.gene,
                c_notation: this.c_notation.trim(),
                p_notation: this.p_notation.trim() || null,
                dup_type: this.dup_type,
            };
            fetch("/api/audit/client-validation", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                keepalive: true,
                body: JSON.stringify({
                    form,
                    input: submittedInput,
                    error,
                }),
            }).catch(() => {});
        },

        formatApiError(payload, status) {
            const detail = payload?.detail;
            if (Array.isArray(detail)) {
                const messages = detail
                    .map(item => String(item?.msg || "").replace(/^Value error,\s*/i, ""))
                    .filter(Boolean);
                if (messages.length > 0) return messages.join("; ");
            }
            if (typeof detail === "string" && detail.trim()) return detail;
            return `The input could not be processed. Use HGVS c. and p. notation, for example c.4185G>A and p.(Gln1395=). Error ${status}.`;
        },

        async classify() {
            this.error = "";
            this.result = null;
            this.resetManualItems();

            if (!this.c_notation.trim()) {
                this.error = "Please enter a c. notation.";
                this.logClientValidation(this.error);
                return;
            }
            if (!this.c_notation.trim().startsWith("c.")) {
                this.error = "Notation must start with 'c.' - e.g. c.4185G>A";
                this.logClientValidation(this.error);
                return;
            }

            this.loading = true;
            this.progress = 0;
            this.progressText = "Resolving coordinates...";

            // simulate progress steps while waiting for API
            const steps = [
                { pct: 15, text: "Resolving coordinates..." },
                { pct: 30, text: "Querying SpliceAI..." },
                { pct: 45, text: "Looking up gnomAD frequencies..." },
                { pct: 60, text: "Evaluating ACMG criteria..." },
                { pct: 75, text: "Checking ClinVar..." },
                { pct: 85, text: "Querying ClinGen ERepo..." },
                { pct: 95, text: "Finalising classification..." },
            ];

            let stepIdx = 0;
            const progressTimer = setInterval(() => {
                if (stepIdx < steps.length) {
                    this.progress = steps[stepIdx].pct;
                    this.progressText = steps[stepIdx].text;
                    stepIdx++;
                }
            }, 1500);

            try {
                const resp = await fetch("/api/classify", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        gene: this.gene,
                        c_notation: this.c_notation.trim(),
                        p_notation: this.p_notation.trim() || null,
                        dup_type: this.dup_type,
                    }),
                });

                clearInterval(progressTimer);

                if (!resp.ok) {
                    const err = await resp.json().catch(() => ({}));
                    this.error = this.formatApiError(err, resp.status);
                    this.loading = false;
                    return;
                }

                this.progress = 100;
                this.progressText = "Done.";
                this.result = await resp.json();
                this.prefillManualReviewFromResult();

            } catch (e) {
                clearInterval(progressTimer);
                this.error = "Network error - could not reach the server.";
            }

            this.loading = false;
        },

        // ── Batch: parse CSV text whenever batchText changes ──────────────
        parseBatch() {
            this.batchParseError = "";
            this.batchParsed = [];
            if (!this.batchText.trim()) return;

            const lines = this.batchText.trim().split("\n").filter(l => l.trim());
            const parsed = [];
            const errors = [];

            for (let i = 0; i < lines.length; i++) {
                // Skip header lines
                if (lines[i].toLowerCase().startsWith("gene") || lines[i].toLowerCase().startsWith("#")) continue;

                const parts = lines[i].split(",").map(s => s.trim());
                if (parts.length < 2) {
                    errors.push(`Line ${i + 1}: need at least gene and c. notation`);
                    continue;
                }
                const gene = parts[0].toUpperCase();
                if (!["BRCA1", "BRCA2"].includes(gene)) {
                    errors.push(`Line ${i + 1}: gene must be BRCA1 or BRCA2`);
                    continue;
                }

                let cRaw = parts[1];
                let pRaw = parts[2] || "";
                const dupRaw = parts[3] || "Unknown";

                // Handle appended protein notation, including "(p.Val2050del)".
                if (!pRaw) {
                    const combined = cRaw.match(/^(c\.\S+)\s+(\(?p\..+\)?)$/i);
                    if (combined) {
                        cRaw = combined[1];
                        pRaw = combined[2];
                    }
                }
                if (pRaw.startsWith("(p.") && pRaw.endsWith(")")) {
                    pRaw = pRaw.substring(1, pRaw.length - 1);
                }
                if (pRaw.startsWith("p.") && !pRaw.startsWith("p.(")) {
                    pRaw = `p.(${pRaw.substring(2)})`;
                }

                if (!cRaw.startsWith("c.")) {
                    errors.push(`Line ${i + 1}: c. notation must start with 'c.'`);
                    continue;
                }
                parsed.push({
                    gene,
                    c_notation: cRaw,
                    p_notation: pRaw,
                    dup_type: dupRaw,
                });
            }

            if (errors.length > 0) {
                this.batchParseError = errors.join("; ");
            }
            this.batchParsed = parsed;
        },

        // ── Batch: classify all parsed variants ───────────────────────────
        async classifyBatch() {
            this.parseBatch();
            if (this.batchParseError) {
                this.logClientValidation(
                    this.batchParseError,
                    "batch",
                    { batch_text: this.batchText.slice(0, 4000) },
                );
            }
            if (this.batchParsed.length === 0) return;

            this.batchRunning = true;
            this.batchDone = 0;
            this.batchTotal = this.batchParsed.length;
            this.batchResults = new Array(this.batchParsed.length).fill(null);

            // Run up to 3 classify calls concurrently for rate-limit safety
            const CONCURRENCY = 3;
            const queue = [...this.batchParsed.entries()];
            const active = new Set();

            const runNext = async () => {
                if (queue.length === 0) return;
                const [idx, item] = queue.shift();
                const task = (async () => {
                    try {
                        const resp = await fetch("/api/classify", {
                            method: "POST",
                            headers: { "Content-Type": "application/json" },
                            body: JSON.stringify({
                                gene: item.gene,
                                c_notation: item.c_notation,
                                p_notation: item.p_notation || null,
                                dup_type: item.dup_type || "Unknown",
                            }),
                        });
                        if (resp.ok) {
                            const data = await resp.json();
                            this.batchResults[idx] = {
                                status: "ok",
                                gene: item.gene,
                                c_notation: item.c_notation,
                                p_notation: item.p_notation,
                                result: data,
                            };
                        } else {
                            const err = await resp.json().catch(() => ({}));
                            this.batchResults[idx] = {
                                status: "error",
                                gene: item.gene,
                                c_notation: item.c_notation,
                                p_notation: item.p_notation,
                                error: this.formatApiError(err, resp.status),
                            };
                        }
                    } catch (e) {
                        this.batchResults[idx] = {
                            status: "error",
                            gene: item.gene,
                            c_notation: item.c_notation,
                            p_notation: item.p_notation,
                            error: "Network error",
                        };
                    }
                    this.batchDone++;
                    active.delete(task);
                    await runNext();
                })();
                active.add(task);
            };

            // Kick off initial workers
            const starters = [];
            for (let i = 0; i < Math.min(CONCURRENCY, this.batchParsed.length); i++) {
                starters.push(runNext());
            }
            await Promise.all(starters);

            // Wait for all active tasks to finish
            while (active.size > 0) {
                await Promise.race(active);
            }

            this.batchRunning = false;
        },

        // ── Batch: download results as CSV ────────────────────────────────
        downloadBatchCsv() {
            const header = [
                "Gene", "c_notation", "p_notation",
                "Class", "Label", "Points",
                "Criteria", "ClinVar", "ENIGMA_EP",
                "Classification_note", "VUS_category", "VUS_what_to_check",
                "RNA_review", "RNA_branches",
                "Splice_PS1_review", "Splice_PS1_branches",
                "PVS1_INIT_review", "PVS1_INIT_branches",
                "Warnings"
            ];
            const rows = this.batchResults.map(row => {
                if (!row || row.status === "error") {
                    return [
                        row ? row.gene : "",
                        row ? row.c_notation : "",
                        row ? (row.p_notation || "") : "",
                        "ERROR", "", "", "", "", "",
                        "", "", "",
                        "", "", "", "", "", "",
                        row ? (row.error || "") : ""
                    ];
                }
                const r = row.result;
                const criteria = r.criteria
                    .filter(c => c.applies)
                    .map(c => c.name + (c.strength ? "_" + c.strength.replace(" ", "") : ""))
                    .join("; ");
                const warnings = r.warnings.join(" | ");
                const rnaReview = r.rna_review && r.rna_review.recommended
                    ? `yes/${r.rna_review.priority || ""}`
                    : "";
                const rnaBranches = r.rna_review && r.rna_review.recommended
                    ? (r.rna_review.potential_branches || []).join("; ")
                    : "";
                const splicePs1Review = r.splice_ps1_review && r.splice_ps1_review.recommended
                    ? `yes/${r.splice_ps1_review.priority || ""}`
                    : "";
                const splicePs1Branches = r.splice_ps1_review && r.splice_ps1_review.recommended
                    ? (r.splice_ps1_review.potential_branches || []).join("; ")
                    : "";
                const initiationReview = r.initiation_review && r.initiation_review.recommended
                    ? `yes/${r.initiation_review.priority || ""}`
                    : "";
                const initiationBranches = r.initiation_review && r.initiation_review.recommended
                    ? (r.initiation_review.potential_branches || []).join("; ")
                    : "";
                return [
                    row.gene,
                    row.c_notation,
                    row.p_notation || "",
                    r.predicted_class,
                    r.predicted_label,
                    r.total_points,
                    criteria,
                    r.external ? r.external.clinvar_classification : "",
                    r.external ? r.external.enigma_ep_class : "",
                    r.classification_note,
                    r.vus_explanation ? r.vus_explanation.category : "",
                    r.vus_explanation ? r.vus_explanation.what_to_check : "",
                    rnaReview,
                    rnaBranches,
                    splicePs1Review,
                    splicePs1Branches,
                    initiationReview,
                    initiationBranches,
                    warnings,
                ];
            });

            const escape = v => '"' + String(v).replace(/"/g, '""') + '"';
            const csv = [header, ...rows]
                .map(row => row.map(escape).join(","))
                .join("\n");

            const blob = new Blob(["﻿" + csv], { type: "text/csv;charset=utf-8" });
            const url  = URL.createObjectURL(blob);
            const a    = document.createElement("a");
            a.href     = url;
            a.download = "ariane_batch_results.csv";
            a.click();
            URL.revokeObjectURL(url);
        },
    };
}
