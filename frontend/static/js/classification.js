(function registerClassification(namespace) {
    "use strict";

    namespace.classificationState = function classificationState() {
        return {
        gene: "",
        c_notation: "",
        p_notation: "",
        assembly: "",
        dup_type: "Unknown",
        geneAutoSwitchNotice: "",
        geneInputConflict: "",
        loading: false,
        progress: 0,
        progressText: "",
        error: "",
        result: null,
        };
    };

    namespace.classificationMethods = {
        explicitGeneFromVariantInput() {
            let value = this.c_notation.trim();
            if (!value) return { gene: "", source: "", error: "" };

            let prefixGene = "";
            const prefix = value.match(/^([A-Za-z0-9-]+)\s*(?::\s*|\s+)(.+)$/i);
            if (
                prefix
                && this.configuredGenes.some(item => item.symbol === prefix[1].toUpperCase())
            ) {
                prefixGene = prefix[1].toUpperCase();
                value = prefix[2].trim();
            }

            let transcriptGene = "";
            let transcript = "";
            const transcriptMatch = value.match(/^(NM_\d+(?:\.\d+)?)\s*:/i);
            if (transcriptMatch) {
                transcript = transcriptMatch[1].toUpperCase();
                const accession = transcript.split(".", 1)[0];
                const match = this.configuredGenes.find(item =>
                    String(item.reference_transcript || "").split(".", 1)[0].toUpperCase() === accession
                );
                if (match) transcriptGene = match.symbol;
            }

            if (prefixGene && transcriptGene && prefixGene !== transcriptGene) {
                return {
                    gene: "",
                    source: "",
                    error: `Conflicting identifiers: ${prefixGene} does not match ${transcript} (${transcriptGene}).`,
                };
            }
            if (prefixGene) return { gene: prefixGene, source: "gene prefix", error: "" };
            if (transcriptGene) return { gene: transcriptGene, source: `transcript ${transcript}`, error: "" };
            return { gene: "", source: "", error: "" };
        },

        syncGeneFromVariantInput() {
            const explicit = this.explicitGeneFromVariantInput();
            this.geneInputConflict = explicit.error;
            if (explicit.error) {
                this.geneAutoSwitchNotice = "";
                return false;
            }
            if (!explicit.gene) {
                this.geneAutoSwitchNotice = "";
                return true;
            }
            if (this.gene !== explicit.gene) {
                this.gene = explicit.gene;
                this.geneAutoSwitchNotice = `Gene changed to ${explicit.gene} based on the ${explicit.source} in the variant input.`;
                this.loadManualDefinitions().catch(error => {
                    this.resourceError = `Manual-review guidance could not be loaded: ${error?.message || error}`;
                });
            }
            return true;
        },

        inputWithoutGenePrefix() {
            const value = this.c_notation.trim();
            const prefix = value.match(/^([A-Za-z0-9-]+)\s*(?::\s*|\s+)(.+)$/i);
            return (
                prefix
                && this.configuredGenes.some(item => item.symbol === prefix[1].toUpperCase())
            ) ? prefix[2].trim() : value;
        },

        isGenomicInput() {
            const value = this.inputWithoutGenePrefix();
            if (!value) return false;
            return /^(?:chr)?(?:13|17)[:\s-]+\d+/i.test(value);
        },

        logClientValidation(error, form = "single", input = null) {
            const submittedInput = input || {
                gene: this.gene,
                c_notation: this.c_notation.trim(),
                p_notation: this.p_notation.trim() || null,
                assembly: this.assembly || null,
                dup_type: this.dup_type,
            };
            namespace.api.request("/api/audit/client-validation", {
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
            if ([502, 503, 504].includes(status)) {
                return "The classification service did not finish because an upstream data source or the server timed out. The input was not rejected as invalid. Please retry; if the problem persists, report the time and variant to the administrator.";
            }
            return `The input could not be normalized. Enter c. HGVS or a genomic coordinate with GRCh37/GRCh38. Error ${status}.`;
        },

        async classify() {
            this.error = "";
            this.result = null;
            this.resetManualItems();

            if (!this.c_notation.trim()) {
                this.error = "Please enter a variant.";
                this.logClientValidation(this.error);
                return;
            }
            if (!this.syncGeneFromVariantInput()) {
                this.error = this.geneInputConflict;
                this.logClientValidation(this.error);
                return;
            }
            if (this.isGenomicInput() && !this.assembly) {
                this.error = "Select GRCh37 or GRCh38 for the genomic coordinate.";
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
                const resp = await namespace.api.request("/api/classify", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        gene: this.gene,
                        c_notation: this.c_notation.trim(),
                        p_notation: null,
                        assembly: this.isGenomicInput() ? this.assembly : null,
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
                this.result = this.normalizeCriterionOrder(await resp.json());
                this.gene = this.result.gene;
                this.prefillManualReviewFromResult();

            } catch (e) {
                clearInterval(progressTimer);
                this.error = "Network error - could not reach the server.";
            }

            this.loading = false;
        },

        // â”€â”€ Batch: parse CSV text whenever batchText changes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});

