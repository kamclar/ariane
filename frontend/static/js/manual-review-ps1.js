(function registerManualReviewPs1(namespace) {
    "use strict";

    namespace.manualReviewPs1Methods = {
        splicePs1CandidatesForCurrentGene() {
            const currentGene = this.result?.gene || this.gene;
            return (this.splicePs1Candidates.candidates || [])
                .filter(candidate => candidate.gene === currentGene);
        },

        splicePs1CandidateLabel(candidate) {
            const protein = candidate.p_notation ? ` ${candidate.p_notation}` : "";
            return `${candidate.reference_variant}${protein} - ${candidate.classification}; ${candidate.reference_splice_event}`;
        },

        async resolveProteinPs1Reference(item) {
            this.ps1ReferenceError = "";
            this.ps1ReferenceMessage = "";
            this.ps1ReferenceWarning = "";
            if (!this.result) {
                this.ps1ReferenceError = "Classify the assessed variant first.";
                return;
            }
            const referenceNotation = String(item?.evidence?.reference_variant || "").trim();
            if (!referenceNotation) {
                this.ps1ReferenceError = "Enter the reference c. notation first.";
                return;
            }
            this.ps1ReferenceLoading = true;
            try {
                const response = await namespace.api.request("/ui-api/manual-evidence/resolve-ps1-reference", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        gene: this.result.gene,
                        assessed_c_notation: this.result.c_notation,
                        reference_c_notation: referenceNotation,
                    }),
                });
                if (!response.ok) {
                    const error = await response.json().catch(() => ({}));
                    this.ps1ReferenceError = this.formatApiError(error, response.status);
                    return;
                }
                const resolved = await response.json();
                const hasAcceptedClassificationBasis = [
                    "external_vcep_assertion",
                    "enigma_st7_v1_2_reference_set",
                ].includes(resolved.classification_verification);
                item.evidence = {
                    ...item.evidence,
                    reference_variant: `${resolved.reference.gene} ${resolved.reference.c_notation}`,
                    reference_p_notation: resolved.reference.p_notation,
                    reference_classification: hasAcceptedClassificationBasis
                        ? resolved.classification
                        : item.evidence.reference_classification,
                    classification_verification: hasAcceptedClassificationBasis
                        ? resolved.classification_verification
                        : item.evidence.classification_verification,
                    classification_source: hasAcceptedClassificationBasis
                        ? resolved.classification_source
                        : item.evidence.classification_source,
                    same_missense_confirmed: resolved.same_missense_substitution === true,
                    different_nucleotide_change_confirmed: resolved.different_nucleotide_change === true,
                    vua_spliceai_score: resolved.assessed.spliceai_score ?? "",
                    reference_spliceai_score: resolved.reference.spliceai_score ?? "",
                    reference_confirmed_splice_status:
                        resolved.reference_confirmed_splice_status ||
                        item.evidence.reference_confirmed_splice_status,
                    reference_classification_used_ps1:
                        resolved.reference_classification_used_ps1 ||
                        item.evidence.reference_classification_used_ps1,
                    ps1_protein_rationale: hasAcceptedClassificationBasis
                        ? `ARIANE verified ${resolved.classification_source}. The variants have the same ` +
                          "normalized missense consequence and different nucleotide changes. Confirm any " +
                          "remaining assessed-variant RNA/splice facts before submission."
                        : item.evidence.ps1_protein_rationale,
                };
                this.ps1ReferenceMessage = resolved.review_message || "Reference facts resolved.";
                this.ps1ReferenceWarning = resolved.historical_expert_panel_warning || "";
            } catch (e) {
                this.ps1ReferenceError = "Network error - PS1 reference facts could not be resolved.";
            } finally {
                this.ps1ReferenceLoading = false;
            }
        },

        applySplicePs1CandidateFacts(item) {
            const evidence = item.evidence || {};
            const candidate = (this.splicePs1Candidates.candidates || [])
                .find(value => value.key === evidence.splice_ps1_candidate_key);
            if (!candidate) return;

            item.evidence.reference_variant = `${candidate.gene} ${candidate.reference_variant}`;
            item.evidence.reference_classification = candidate.classification;
            item.evidence.reference_classification_source =
                `${candidate.source_label}; ${candidate.classification_basis}`;
            item.evidence.reference_splice_event = candidate.reference_splice_event;
            item.evidence.reference_assay_result_category = candidate.assay_result_category;
            item.evidence.reference_variant_context = candidate.assay_context;
            item.evidence.candidate_source_status = candidate.eligibility_status;
            item.evidence.candidate_source_row = candidate.source_row;
            item.evidence.candidate_source_sha256 = candidate.source_file_sha256;

            if (!item.references.trim()) {
                item.references = [
                    candidate.source_label,
                    candidate.source_url,
                ].filter(Boolean).join("\n");
            }
        },
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});
