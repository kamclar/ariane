(function registerManualReviewPersistence(namespace) {
    "use strict";

    namespace.manualReviewPersistenceMethods = {
        async saveManualReviewDraft() {
            if (!this.manualResult || !this.result) return;
            this.manualRecordLoading = true;
            this.manualRecordError = "";
            try {
                const response = await namespace.api.request("/api/review-records", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        module1_result: this.result,
                        manual_evidence: this.manualEvidenceRequestPayload(),
                        reviewer_role: this.manualReviewerRole.trim(),
                    }),
                });
                if (!response.ok) {
                    if (response.status === 401) {
                        this.manualRecordError = "Saving requires an authenticated reviewer account. Open the admin audit page, sign in, and try again.";
                        return;
                    }
                    const error = await response.json().catch(() => ({}));
                    this.manualRecordError = this.formatApiError(error, response.status);
                    return;
                }
                this.manualSavedRecord = await response.json();
                this.manualApprovalName = this.manualAssessor;
                this.manualApprovalRole = this.manualReviewerRole;
            } catch (e) {
                this.manualRecordError = "The review draft could not be saved.";
            } finally {
                this.manualRecordLoading = false;
            }
        },

        async approveManualReviewRecord() {
            if (!this.manualSavedRecord || this.manualSavedRecord.status !== "draft") return;
            this.manualRecordLoading = true;
            this.manualRecordError = "";
            try {
                const response = await namespace.api.request(
                    `/api/review-records/${encodeURIComponent(this.manualSavedRecord.record_id)}/approve`,
                    {
                        method: "POST",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify({
                            approver_name: this.manualApprovalName.trim(),
                            approver_role: this.manualApprovalRole.trim(),
                            approval_comment: this.manualApprovalComment.trim(),
                            attestation_confirmed: this.manualApprovalAttested,
                        }),
                    }
                );
                if (!response.ok) {
                    if (response.status === 401) {
                        this.manualRecordError = "Approval requires an authenticated reviewer account.";
                        return;
                    }
                    const error = await response.json().catch(() => ({}));
                    this.manualRecordError = this.formatApiError(error, response.status);
                    return;
                }
                this.manualSavedRecord = await response.json();
            } catch (e) {
                this.manualRecordError = "The review record could not be approved.";
            } finally {
                this.manualRecordLoading = false;
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
                persisted_review_record: this.manualSavedRecord,
                submitted_manual_evidence: this.manualItems
                    .filter(item => item.enabled)
                    .map(item => ({
                        code: item.code,
                        evidence: item.evidence,
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
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});
