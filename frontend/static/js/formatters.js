(function registerFormatters(namespace) {
    "use strict";

    namespace.formattersMethods = {
        signedPoints(value) {
            const points = Number(value);
            if (!Number.isFinite(points)) return "?";
            return points > 0 ? `+${points}` : String(points).replace("-", "\u2212");
        },

        pointMeterPosition(value) {
            const points = Number(value);
            if (!Number.isFinite(points)) return 50;
            const scaleMaximum = 15;
            const scaleMinimum = -12;
            const clamped = Math.min(scaleMaximum, Math.max(scaleMinimum, points));
            const rawPosition = ((scaleMaximum - clamped) / (scaleMaximum - scaleMinimum)) * 100;
            return Math.min(95, Math.max(5, rawPosition));
        },

        pointMeterSummary(value, mixedEvidence = false) {
            const points = Number(value);
            if (!Number.isFinite(points)) return "Point total is unavailable.";
            const formatted = this.signedPoints(points);
            let position;
            if (points >= 10) position = "Pathogenic point threshold reached.";
            else if (points >= 6) position = `${10 - points} point(s) from the Pathogenic threshold.`;
            else if (points >= -1) position = `${6 - points} point(s) from the Likely Pathogenic threshold.`;
            else if (points >= -6) position = `${points + 7} point(s) from the Benign threshold.`;
            else position = "Benign point threshold reached.";

            const method = mixedEvidence
                ? "Mixed evidence: the ENIGMA point-based classification method applies."
                : "One-direction evidence: ENIGMA criterion combinations determine the class; the point position is contextual only.";
            return `${formatted} points. ${position} ${method}`;
        },

        criterionSortKey(code) {
            const groupOrder = { PVS: 0, PS: 1, PM: 2, PP: 3, BA: 4, BS: 5, BP: 6 };
            const normalized = String(code || "").trim().toUpperCase();
            const match = normalized.match(/^(PVS|PS|PM|PP|BA|BS|BP)(\d+)/);
            if (!match) return [99, 99, 1, normalized];
            const base = `${match[1]}${match[2]}`;
            return [groupOrder[match[1]], Number(match[2]), normalized === base ? 0 : 1, normalized];
        },

        sortCriterionList(criteria) {
            return [...(criteria || [])].sort((left, right) => {
                const a = this.criterionSortKey(left.name || left.code);
                const b = this.criterionSortKey(right.name || right.code);
                for (let i = 0; i < 3; i++) {
                    if (a[i] !== b[i]) return a[i] - b[i];
                }
                return a[3].localeCompare(b[3]);
            });
        },

        normalizeCriterionOrder(result) {
            if (!result) return result;
            result.criteria = this.sortCriterionList(result.criteria);
            result.excluded_criteria = this.sortCriterionList(result.excluded_criteria);
            result.not_applicable_criteria = this.sortCriterionList(result.not_applicable_criteria);
            return result;
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
            if (result.protein_ps1_review && result.protein_ps1_review.recommended) {
                labels.push(`PS1_PROTEIN ${result.protein_ps1_review.priority || ""}`.trim());
            }
            if (result.initiation_review && result.initiation_review.recommended) {
                labels.push(`PVS1_INIT ${result.initiation_review.priority || ""}`.trim());
            }
            return labels.join("; ");
        },
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});

