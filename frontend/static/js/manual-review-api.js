(function registerManualReviewApi(namespace) {
    "use strict";

    namespace.manualReviewApiMethods = {
        manualCriteriaPayload() {
            return this.manualItems.map(item => ({
                code: item.code,
                enabled: item.enabled,
                evidence: item.evidence,
                notes: item.notes,
                references: item.references
                    .split(/\r?\n/)
                    .map(value => value.trim())
                    .filter(Boolean),
            }));
        },

        manualEvidenceRequestPayload() {
            return {
                base_criteria: this.result?.criteria || [],
                variant_context: this.result ? {
                    gene: this.result.gene,
                    c_notation: this.result.c_notation,
                    p_notation: this.result.p_notation,
                } : null,
                manual_criteria: this.manualCriteriaPayload(),
                assessor: this.manualAssessor.trim(),
                assessed_at: this.manualAssessedAt,
            };
        },

        async refreshManualFormStatuses() {
            if (!this.result) return;
            const requestId = ++this.manualStatusRequestId;
            this.manualStatusLoading = true;
            this.manualStatusError = "";
            try {
                const response = await namespace.api.request("/ui-api/manual-evidence/status", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        base_criteria: this.result.criteria,
                        variant_context: {
                            gene: this.result.gene,
                            c_notation: this.result.c_notation,
                            p_notation: this.result.p_notation,
                        },
                        manual_criteria: this.manualCriteriaPayload(),
                    }),
                });
                if (!response.ok) {
                    const error = await response.json().catch(() => ({}));
                    if (requestId === this.manualStatusRequestId) {
                        this.manualStatusError = this.formatApiError(error, response.status);
                    }
                    return;
                }
                const payload = await response.json();
                if (requestId === this.manualStatusRequestId) {
                    this.manualCriterionStatuses = Object.fromEntries(
                        (payload.criteria || []).map(item => [item.code, item])
                    );
                }
            } catch (e) {
                if (requestId === this.manualStatusRequestId) {
                    this.manualStatusError = "Form status could not be checked.";
                }
            } finally {
                if (requestId === this.manualStatusRequestId) {
                    this.manualStatusLoading = false;
                }
            }
        },

        async evaluateManualEvidence() {
            this.manualError = "";
            this.manualResult = null;
            this.manualSavedRecord = null;
            this.manualRecordError = "";
            if (!this.result) {
                this.manualError = "Classify a variant first.";
                return;
            }

            this.manualLoading = true;
            try {
                const response = await namespace.api.request("/ui-api/manual-evidence/evaluate", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify(this.manualEvidenceRequestPayload()),
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
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});
