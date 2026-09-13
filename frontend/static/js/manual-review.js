(function registerManualReview(namespace) {
    "use strict";

    namespace.manual_reviewState = function manual_reviewState() {
        return {
        manualItems: [],
        manualReviewOpen: false,
        manualAssessor: "",
        manualReviewerRole: "",
        manualAssessedAt: new Date().toISOString().slice(0, 10),
        manualLoading: false,
        manualError: "",
        manualResult: null,
        ps1ReferenceLoading: false,
        ps1ReferenceError: "",
        ps1ReferenceMessage: "",
        ps1ReferenceWarning: "",
        manualStatusLoading: false,
        manualStatusError: "",
        manualCriterionStatuses: {},
        manualStatusRequestId: 0,
        manualRecordLoading: false,
        manualRecordError: "",
        manualSavedRecord: null,
        manualApprovalName: "",
        manualApprovalRole: "",
        manualApprovalComment: "",
        manualApprovalAttested: false,
        };
    };

    namespace.manual_reviewMethods = {
        manualEvidenceDefaults(code) {
            if (code === "PS4") return {
                case_control_country_matched: false,
                case_control_ethnicity_matched: false,
            };
            if (code === "PM3") return {
                cooccurring_variant_classification_basis: "not_assessed",
                vua_benign_population_review: "not_assessed",
            };
            if (code === "BS2") return {
                cooccurring_variant_classification_basis: "not_assessed",
            };
            if (code === "PP1") return { very_strong_effect_basis: "" };
            if (code === "BS4") return { likelihood_ratio_components: [] };
            if (code === "PP4" || code === "BP5") return {
                clinical_lr_scale: "lr",
                source_review_status: "enigma_recognised",
                clinical_evidence_types: [],
                independence_review_confirmed: false,
            };
            if (code === "PS3" || code === "BS3") return {
                assay_scope: "",
                functional_conclusion: "",
                calibration_status: "not_reviewed",
                pathogenic_and_benign_controls_confirmed: false,
                curated_strength: "",
            };
            if (code === "PS1_PROTEIN") return {
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
                    "ENIGMA Supplementary Table 3 v1.2",
                ],
                vua_confirmed_splice_status: "not_assessed",
                reference_confirmed_splice_status: "not_assessed",
                reference_classification_used_ps1: "unknown",
                reference_ps1_dependency_reference: "",
                direct_reciprocal_dependency_excluded: false,
                ps1_protein_rationale: "",
            };
            return {};
        },

        resetManualItems() {
            this.manualItems = [
                "PVS1_INIT", "PVS1_RNA", "PS1_PROTEIN", "PS1_SPLICE", "PS3",
                "PS4", "PM3", "PP1", "PP4", "BS2", "BS3", "BS4", "BP5",
                "BP7_RNA",
            ].map(code => ({
                code,
                enabled: false,
                evidence: this.manualEvidenceDefaults(code),
                notes: "",
                references: "",
            }));
            this.manualResult = null;
            this.manualReviewOpen = false;
            this.manualError = "";
            this.ps1ReferenceError = "";
            this.ps1ReferenceMessage = "";
            this.ps1ReferenceWarning = "";
            this.manualStatusError = "";
            this.manualCriterionStatuses = {};
            this.manualStatusLoading = false;
            this.manualStatusRequestId += 1;
            this.manualRecordLoading = false;
            this.manualRecordError = "";
            this.manualSavedRecord = null;
            this.manualApprovalName = "";
            this.manualApprovalRole = "";
            this.manualApprovalComment = "";
            this.manualApprovalAttested = false;
        },

        manualDefinition(code) {
            return this.manualDefinitions[code] || {};
        },

        manualReviewRecommendations() {
            const recommendations = [];
            const add = (code, review, fallback) => {
                if (!review?.recommended || recommendations.some(item => item.code === code)) return;
                recommendations.push({
                    code,
                    priority: review.priority || "review",
                    reason: review.summary || review.title || fallback,
                });
            };
            add("PVS1_RNA", this.result?.rna_review, "RNA evidence is available for expert review.");
            add("PVS1_INIT", this.result?.initiation_review, "Complete the initiation-codon review.");
            add("PS1_SPLICE", this.result?.splice_ps1_review, "A splice PS1 candidate requires expert review.");
            add("PS1_PROTEIN", this.result?.protein_ps1_review, "A protein PS1 candidate requires expert review.");
            return recommendations;
        },

        isManualReviewRecommended(code) {
            return this.manualReviewRecommendations().some(item => item.code === code);
        },

        manualReviewGroups() {
            const groups = new Map();
            for (const item of this.manualItems) {
                const definition = this.manualDefinition(item.code);
                const id = definition.group_id || "other";
                if (!groups.has(id)) {
                    groups.set(id, {
                        id,
                        title: definition.group_title || "Other evidence",
                        order: Number(definition.group_order ?? 999),
                        items: [],
                    });
                }
                groups.get(id).items.push(item);
            }
            return Array.from(groups.values())
                .map(group => ({
                    ...group,
                    items: group.items.sort((left, right) =>
                        Number(this.manualDefinition(left.code).criterion_order ?? 999) -
                        Number(this.manualDefinition(right.code).criterion_order ?? 999)
                    ),
                }))
                .sort((left, right) => left.order - right.order);
        },

        manualGroupRecommendedCount(group) {
            return group.items.filter(item => this.isManualReviewRecommended(item.code)).length;
        },

        manualGroupEnabledCount(group) {
            return group.items.filter(item => item.enabled).length;
        },

        manualCriterionStatus(code) {
            return this.manualCriterionStatuses[code] || {
                status: "not_started",
                message: "Select this criterion to review its evidence.",
                suggested_strength: null,
            };
        },

        manualCriterionStatusLabel(code) {
            const status = this.manualCriterionStatus(code);
            if (status.status === "ready") {
                return status.suggested_strength
                    ? `Ready: ${status.suggested_strength}`
                    : "Ready";
            }
            if (status.status === "incomplete") return "Needs information";
            return "Not started";
        },

        enableManualCriterion(code) {
            const item = this.manualItems.find(value => value.code === code);
            if (!item) return;
            this.manualReviewOpen = true;
            item.enabled = true;
            void this.refreshManualFormStatuses();
            window.requestAnimationFrame(() => {
                document.getElementById(`manual-criterion-${code}`)?.scrollIntoView({
                    behavior: "smooth",
                    block: "start",
                });
            });
        },

        prefillManualReviewFromResult() {
            const rnaReview = this.result?.rna_review;
            if (rnaReview?.recommended && rnaReview?.manual_review_prefill) {
                const item = this.manualItems.find(value => value.code === "PVS1_RNA");
                const prefill = rnaReview.manual_review_prefill;
                if (item && Object.keys(prefill).length) {
                    item.evidence = {
                        ...item.evidence,
                        ...prefill,
                    };
                    item.notes = prefill.table4_context || item.notes;
                    item.references = [
                        prefill.source_citation,
                        ...(prefill.source_references || []),
                        rnaReview.source_url,
                    ].filter(Boolean).join("\n");
                }
            }

            for (const code of ["PS3", "BS3"]) {
                const item = this.manualItems.find(value => value.code === code);
                if (!item) continue;
                item.evidence = {
                    ...item.evidence,
                    assessed_gene: this.result?.gene || "",
                    assessed_c_notation: this.result?.c_notation || "",
                    assessed_p_notation: this.result?.p_notation || "",
                    assessed_variant_type: this.result?.variant_type || "",
                    reference_transcript: this.result?.reference_transcript || "",
                    spliceai_score: this.result?.spliceai_audit?.score ?? "",
                    spliceai_status: this.result?.spliceai_audit?.status || "unavailable",
                    table9_lookup_status: (this.result?.criteria || []).some(
                        criterion => criterion.name === code && criterion.applies
                    ) ? "already_applied" : "no_applied_record",
                };
            }

            const clinicalAudit = this.result?.clinical_lr_audit;
            if (clinicalAudit?.application_status === "review_required" &&
                clinicalAudit?.candidate_likelihood_ratio != null) {
                for (const code of ["PP4", "BP5"]) {
                    const item = this.manualItems.find(value => value.code === code);
                    if (!item) continue;
                    item.evidence = {
                        ...item.evidence,
                        clinical_lr_value: clinicalAudit.candidate_likelihood_ratio,
                        clinical_lr_scale: "lr",
                        source_review_status: "unreviewed",
                        source_citation: (clinicalAudit.source_bundle_ids || []).join(", "),
                        clinical_evidence_types: (clinicalAudit.clinical_evidence_types || [])
                            .map(value => String(value).trim().toLowerCase().replace(/[^a-z0-9]+/g, "_")),
                        clinical_data_summary: clinicalAudit.overlap_assessment_note ||
                            "Clinical LR candidate found by ARIANE. Review source overlap and independence before use.",
                    };
                    item.references = (clinicalAudit.source_bundle_ids || []).join("\n");
                }
            }

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
            if (this.manualReviewRecommendations().length) {
                this.manualReviewOpen = true;
            }
            void this.refreshManualFormStatuses();
        },

    };
})(window.ArianeFrontend = window.ArianeFrontend || {});
