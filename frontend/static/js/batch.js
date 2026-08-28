(function registerBatch(namespace) {
    "use strict";

    namespace.batchState = function batchState() {
        return {
        batchText: "",
        batchParsed: [],
        batchParseError: "",
        batchRunning: false,
        batchDone: 0,
        batchTotal: 0,
        batchResults: [],
        };
    };

    namespace.batchMethods = {
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
                const allowedGenes = this.configuredGenes.map(item => item.symbol);
                if (!allowedGenes.includes(gene)) {
                    errors.push(`Line ${i + 1}: gene must be one of ${allowedGenes.join(", ")}`);
                    continue;
                }

                let cRaw = parts[1];
                let pRaw = parts[2] || "";
                let assemblyRaw = "";
                if (/^(?:GRCh3[78]|hg(?:19|38))$/i.test(pRaw)) {
                    assemblyRaw = /^hg19$/i.test(pRaw) ? "GRCh37"
                        : /^hg38$/i.test(pRaw) ? "GRCh38"
                        : pRaw.replace(/^grch/i, "GRCh");
                    pRaw = "";
                }
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

                parsed.push({
                    gene,
                    c_notation: cRaw,
                    p_notation: pRaw,
                    assembly: assemblyRaw,
                    dup_type: dupRaw,
                });
            }

            if (errors.length > 0) {
                this.batchParseError = errors.join("; ");
            }
            this.batchParsed = parsed;
        },

        // â”€â”€ Batch: classify all parsed variants â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

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
                        const resp = await namespace.api.request("/api/classify", {
                            method: "POST",
                            headers: { "Content-Type": "application/json" },
                            body: JSON.stringify({
                                gene: item.gene,
                                c_notation: item.c_notation,
                                p_notation: item.p_notation || null,
                                assembly: item.assembly || null,
                                dup_type: item.dup_type || "Unknown",
                            }),
                        });
                        if (resp.ok) {
                            const data = this.normalizeCriterionOrder(await resp.json());
                            this.batchResults[idx] = {
                                status: "ok",
                                gene: item.gene,
                                c_notation: data.c_notation,
                                p_notation: data.p_notation,
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

        // â”€â”€ Batch: download results as CSV â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

        downloadBatchCsv() {
            const header = [
                "Gene", "c_notation", "p_notation",
                "Class", "Label", "Points",
                "Criteria", "Excluded_criteria", "ClinVar", "ENIGMA_EP",
                "Classification_note", "VUS_category", "VUS_what_to_check",
                "RNA_review", "RNA_branches",
                "Splice_PS1_review", "Splice_PS1_branches",
                "Protein_PS1_review", "Protein_PS1_branches",
                "PVS1_INIT_review", "PVS1_INIT_branches",
                "Warnings"
            ];
            const rows = this.batchResults.map(row => {
                if (!row || row.status === "error") {
                    return [
                        row ? row.gene : "",
                        row ? row.c_notation : "",
                        row ? (row.p_notation || "") : "",
                        "ERROR",
                        ...Array(header.length - 5).fill(""),
                        row ? (row.error || "") : "",
                    ];
                }
                const r = row.result;
                const criteria = r.criteria
                    .filter(c => c.applies)
                    .map(c => c.name + (c.strength ? "_" + c.strength.replace(" ", "") : ""))
                    .join("; ");
                const excludedCriteria = (r.excluded_criteria || [])
                    .map(c => `${c.name}${c.strength ? "_" + c.strength.replace(" ", "") : ""}: ${c.reason}`)
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
                const proteinPs1Review = r.protein_ps1_review && r.protein_ps1_review.recommended
                    ? `yes/${r.protein_ps1_review.priority || ""}`
                    : "";
                const proteinPs1Branches = r.protein_ps1_review && r.protein_ps1_review.recommended
                    ? (r.protein_ps1_review.potential_branches || []).join("; ")
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
                    excludedCriteria,
                    r.external ? r.external.clinvar_classification : "",
                    r.external ? r.external.enigma_ep_class : "",
                    r.classification_note,
                    r.vus_explanation ? r.vus_explanation.category : "",
                    r.vus_explanation ? r.vus_explanation.what_to_check : "",
                    rnaReview,
                    rnaBranches,
                    splicePs1Review,
                    splicePs1Branches,
                    proteinPs1Review,
                    proteinPs1Branches,
                    initiationReview,
                    initiationBranches,
                    warnings,
                ];
            });

            const escape = v => '"' + String(v).replace(/"/g, '""') + '"';
            const csv = [header, ...rows]
                .map(row => row.map(escape).join(","))
                .join("\n");

            const blob = new Blob(["ï»¿" + csv], { type: "text/csv;charset=utf-8" });
            const url  = URL.createObjectURL(blob);
            const a    = document.createElement("a");
            a.href     = url;
            a.download = "ariane_batch_results.csv";
            a.click();
            URL.revokeObjectURL(url);
        },
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});

