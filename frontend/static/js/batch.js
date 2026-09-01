(function registerBatch(namespace) {
    "use strict";

    const ASSEMBLY_PATTERN = /^(?:GRCh3[78]|hg(?:19|38))$/i;
    const PROTEIN_PATTERN = /^\(?p\..+\)?$/i;

    function normalizeAssembly(value) {
        if (/^hg19$/i.test(value)) return "GRCh37";
        if (/^hg38$/i.test(value)) return "GRCh38";
        return value.replace(/^grch/i, "GRCh");
    }

    function unwrapMarkdownLine(value) {
        let line = value.trim().replace(/\\$/, "").trim();
        if (line.startsWith("*") && line.endsWith("*") && line.length > 2) {
            line = line.substring(1, line.length - 1).trim();
        }
        return line.replace(/\\([_*])/g, "$1");
    }

    function splitBatchLine(value) {
        const line = unwrapMarkdownLine(value);
        if (line.includes(",")) {
            return line.split(",").map(item => item.trim());
        }
        return line.split(/\s+/).filter(Boolean);
    }

    namespace.parseBatchInput = function parseBatchInput(text, configuredGenes) {
        const lines = text.trim().split(/\r?\n/).filter(line => line.trim());
        const parsed = [];
        const errors = [];
        const allowedGenes = configuredGenes.map(item => item.symbol.toUpperCase());

        for (let i = 0; i < lines.length; i++) {
            const cleanedLine = unwrapMarkdownLine(lines[i]);
            if (/^(?:gene(?:\s|,|$)|#)/i.test(cleanedLine)) continue;

            const parts = splitBatchLine(lines[i]);
            if (parts.length < 2) {
                errors.push(`Line ${i + 1}: need at least gene and c. notation`);
                continue;
            }
            const gene = parts[0].toUpperCase();
            if (!allowedGenes.includes(gene)) {
                errors.push(`Line ${i + 1}: gene must be one of ${allowedGenes.join(", ")}`);
                continue;
            }

            let cRaw = parts[1];
            let pRaw = "";
            let assemblyRaw = "";
            let dupRaw = "Unknown";
            const remaining = parts.slice(2);

            while (remaining.length && !remaining[remaining.length - 1]) {
                remaining.pop();
            }
            if (remaining.length && !remaining[0]) {
                remaining.shift();
            }

            if (remaining.length && PROTEIN_PATTERN.test(remaining[0])) {
                pRaw = remaining.shift();
            }
            if (remaining.length && ASSEMBLY_PATTERN.test(remaining[0])) {
                assemblyRaw = normalizeAssembly(remaining.shift());
            }
            if (remaining.length) {
                dupRaw = remaining.join(" ");
            }

            // Preserve the established comma format where the protein notation
            // may be appended to the variant field instead of using its own column.
            if (!pRaw) {
                const combined = cRaw.match(/^(\S+)\s+(\(?p\..+\)?)$/i);
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

        return { parsed, errors };
    };

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

            const outcome = namespace.parseBatchInput(
                this.batchText,
                this.configuredGenes,
            );
            if (outcome.errors.length > 0) {
                this.batchParseError = outcome.errors.join("; ");
            }
            this.batchParsed = outcome.parsed;
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
