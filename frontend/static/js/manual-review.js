(function registerManualReview(namespace) {
    "use strict";

    namespace.manual_reviewState = function manual_reviewState() {
        return {
        manualItems: [],
        manualAssessor: "",
        manualAssessedAt: new Date().toISOString().slice(0, 10),
        manualLoading: false,
        manualError: "",
        manualResult: null,
        ps1ReferenceLoading: false,
        ps1ReferenceError: "",
        ps1ReferenceMessage: "",
        };
    };

    namespace.manual_reviewMethods = {
        resetManualItems() {
            this.manualItems = ["PVS1_INIT", "PVS1_RNA", "PS1_PROTEIN", "PS1_SPLICE", "PS4", "PM3", "PP1", "PP4", "BS2", "BS4", "BP7_RNA"].map(code => ({
                code,
                enabled: false,
                evidence: code === "PS4" ? {
                    case_control_country_matched: false,
                    case_control_ethnicity_matched: false,
                } : code === "PM3" ? {
                    cooccurring_variant_classification_basis: "not_assessed",
                    vua_benign_population_review: "not_assessed",
                } : code === "BS2" ? {
                    cooccurring_variant_classification_basis: "not_assessed",
                } : code === "PP1" ? {
                    very_strong_effect_basis: "",
                } : code === "BS4" ? {
                    likelihood_ratio_components: [],
                } : code === "PP4" ? {
                    clinical_lr_scale: "lr",
                    source_review_status: "appendix_b",
                } : code === "PS1_PROTEIN" ? {
                    reference_variant: "",
                    reference_p_notation: "",
                    reference_classification: "",
                    classification_verification: "",
                    classification_source: "",
                    same_missense_confirmed: false,
                    different_nucleotide_change_confirmed: false,
                    vua_spliceai_score: "",
                    reference_spliceai_score: "",
                    splice_source_check_completed: false,
                    splice_sources_checked: [
                        "ENIGMA Specifications Table 9 v1.2",
                        "ENIGMA Supplementary Table 2 v1.2",
                    ],
                    vua_confirmed_splice_status: "not_assessed",
                    reference_confirmed_splice_status: "not_assessed",
                    reference_classification_used_ps1: "unknown",
                    reference_ps1_dependency_reference: "",
                    direct_reciprocal_dependency_excluded: false,
                    ps1_protein_rationale: "",
                } : {},
                notes: "",
                references: "",
            }));
            this.manualResult = null;
            this.manualError = "";
            this.ps1ReferenceError = "";
            this.ps1ReferenceMessage = "";
        },

        manualDefinition(code) {
            return this.manualDefinitions[code] || {};
        },

        splicePs1CandidatesForCurrentGene() {
            const currentGene = this.result?.gene || this.gene;
            return (this.splicePs1Candidates.candidates || [])
                .filter(candidate => candidate.gene === currentGene);
        },

        splicePs1CandidateLabel(candidate) {
            const protein = candidate.p_notation ? ` ${candidate.p_notation}` : "";
            return `${candidate.reference_variant}${protein} - ${candidate.classification}; ${candidate.reference_splice_event}`;
        },

        prefillManualReviewFromResult() {
            if (this.result?.initiation_review?.recommended) {
                const item = this.manualItems.find(value => value.code === "PVS1_INIT");
                if (item && !item.evidence?.reference_variant) {
                    item.evidence.met1_loss_confirmed = true;
                    item.evidence.initiation_flowchart_rationale =
                        "Met1/start-loss variant flagged by ARIANE. Complete the ENIGMA initiation-codon flowchart review: alternative start assessment, upstream P/LP evidence, expected N-terminal impact, and curated PVS1_INIT strength.";
                }
            }

            const review = this.result?.protein_ps1_review;
            const candidate = review?.candidates?.[0];
            if (review?.recommended && candidate) {
                const item = this.manualItems.find(value => value.code === "PS1_PROTEIN");
                if (item && !item.evidence?.reference_variant) {
                    item.evidence = {
                        ...item.evidence,
                        ...(review.manual_review_prefill || {}),
                    };
                    item.references = [
                        candidate.source_dataset,
                        candidate.classification_source,
                        review.source_url,
                    ].filter(Boolean).join("\n");
                    // Check official ENIGMA/ClinGen assertion sources in the background.
                    // The result is already visible and an unavailable service does not block it.
                    void this.resolveProteinPs1Reference(item);
                }
            }
        },

        async resolveProteinPs1Reference(item) {
            this.ps1ReferenceError = "";
            this.ps1ReferenceMessage = "";
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
                const response = await namespace.api.request("/api/manual-evidence/resolve-ps1-reference", {
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
                const hasVerifiedAssertion =
                    resolved.classification_verification === "external_vcep_assertion";
                item.evidence = {
                    ...item.evidence,
                    reference_variant: `${resolved.reference.gene} ${resolved.reference.c_notation}`,
                    reference_p_notation: resolved.reference.p_notation,
                    reference_classification: hasVerifiedAssertion
                        ? resolved.classification
                        : item.evidence.reference_classification,
                    classification_verification: hasVerifiedAssertion
                        ? resolved.classification_verification
                        : item.evidence.classification_verification,
                    classification_source: hasVerifiedAssertion
                        ? resolved.classification_source
                        : item.evidence.classification_source,
                    same_missense_confirmed: resolved.same_missense_substitution === true,
                    different_nucleotide_change_confirmed: resolved.different_nucleotide_change === true,
                    vua_spliceai_score: resolved.assessed.spliceai_score ?? "",
                    reference_spliceai_score: resolved.reference.spliceai_score ?? "",
                    ps1_protein_rationale: hasVerifiedAssertion
                        ? `ARIANE verified ${resolved.classification_source}. The variants have the same ` +
                          "normalized missense consequence and different nucleotide changes. Complete or confirm " +
                          "the recorded RNA/splice source review and PS1 dependency review before submission."
                        : item.evidence.ps1_protein_rationale,
                };
                this.ps1ReferenceMessage = resolved.review_message || "Reference facts resolved.";
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

        async evaluateManualEvidence() {
            this.manualError = "";
            this.manualResult = null;
            if (!this.result) {
                this.manualError = "Classify a variant first.";
                return;
            }

            this.manualLoading = true;
            try {
                const response = await namespace.api.request("/api/manual-evidence/evaluate", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        base_criteria: this.result.criteria,
                        variant_context: {
                            gene: this.result.gene,
                            c_notation: this.result.c_notation,
                            p_notation: this.result.p_notation,
                        },
                        manual_criteria: this.manualItems.map(item => ({
                            code: item.code,
                            enabled: item.enabled,
                            evidence: item.evidence,
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

