(function registerCore(namespace) {
    "use strict";

    namespace.coreState = function coreState() {
        return {
        mode: "single",
        appVersion: "",
        buildRevision: "",
        issueTrackerUrl: "",
        configuredGenes: [],
        transcriptByGene: {},
        policyByGene: {},
        manualDefinitions: {},
        resourceLinks: [],
        resourceError: "",
        splicePs1Candidates: { status: "candidate_discovery_only", candidates: [] },
        };
    };

    namespace.coreMethods = {
        async init() {
            document.documentElement.dataset.arianeReady = "true";
            this.resetManualItems();
            try {
                const response = await namespace.api.request("/ui-api/resources");
                if (response.ok) {
                    const resources = await response.json();
                    this.setAppVersion(resources.version);
                    this.buildRevision = String(resources.build_revision || "").trim();
                    this.issueTrackerUrl = String(resources.issue_tracker_url || "").trim();
                    this.configuredGenes = resources.genes || [];
                    this.transcriptByGene = Object.fromEntries(
                        this.configuredGenes.map(item => [item.symbol, item.reference_transcript])
                    );
                    this.policyByGene = Object.fromEntries(
                        this.configuredGenes.map(item => [item.symbol, {
                            name: item.policy_name,
                            version: item.policy_version,
                            sourceUrl: item.policy_source_url,
                        }])
                    );
                    if (!this.gene && this.configuredGenes.length) {
                        this.gene = this.configuredGenes[0].symbol;
                    }
                    this.manualDefinitions = resources.manual_criteria || {};
                    this.resourceLinks = resources.links || [];
                    this.splicePs1Candidates = resources.splice_ps1_candidates || {
                        status: "candidate_discovery_only",
                        candidates: [],
                    };
                    if (this.gene && !Object.keys(this.manualDefinitions).length) {
                        await this.loadManualDefinitions();
                    }
                } else {
                    this.resourceError = `Reference materials could not be loaded (HTTP ${response.status}). Classification remains available, but manual-review guidance is incomplete.`;
                }
            } catch (e) {
                this.resourceError = `Reference materials could not be loaded: ${e?.message || e}. Classification remains available, but manual-review guidance is incomplete.`;
            }
        },

        setAppVersion(version) {
            const normalized = String(version || "").trim();
            if (!normalized) return;
            this.appVersion = normalized;
            const label = `v${normalized}`;
            const headerVersion = document.getElementById("ariane-version");
            const footerVersion = document.getElementById("ariane-footer-version");
            if (headerVersion) headerVersion.textContent = label;
            if (footerVersion) footerVersion.textContent = label;
        },

        currentPolicyInfo() {
            const symbol = String(this.result?.gene || this.gene || "").toUpperCase();
            return this.policyByGene[symbol] || null;
        },

        currentPolicyLabel() {
            const policy = this.currentPolicyInfo();
            if (!policy) return "the selected VCEP specification";
            return `${policy.name} v${policy.version}`;
        },

        issueReportUrl() {
            if (!this.issueTrackerUrl) return "#";
            const url = new URL(this.issueTrackerUrl, window.location.href);
            const lines = [
                "What happened and what did you expect?",
                "",
                "Technical context added by ARIANE",
                `ARIANE version: ${this.appVersion || "not reported"}`,
            ];
            if (this.buildRevision) lines.push(`Build revision: ${this.buildRevision}`);
            lines.push(`Page: ${this.mode === "batch" ? "Batch" : "Single variant"}`);
            const policy = this.currentPolicyInfo();
            if (policy) lines.push(`Policy: ${policy.name} v${policy.version}`);
            if (this.result) {
                lines.push(`Variant: ${this.result.variant || `${this.result.gene} ${this.result.c_notation}`}`);
                lines.push(`Module 1 result: Class ${this.result.predicted_class} ${this.result.predicted_label}; ${this.result.total_points} points`);
            } else if (this.gene || this.c_notation) {
                lines.push(`Current input: ${[this.gene, this.c_notation].filter(Boolean).join(" ")}`);
            }
            lines.push(`Browser: ${navigator.userAgent}`);
            lines.push(`Created: ${new Date().toISOString()}`);
            url.searchParams.set("description", lines.join("\n"));
            return url.toString();
        },

        async loadManualDefinitions() {
            if (!this.gene) return;
            const response = await namespace.api.request(`/ui-api/resources?gene=${encodeURIComponent(this.gene)}`);
            if (!response.ok) throw new Error(`Manual evidence resources HTTP ${response.status}`);
            const resources = await response.json();
            this.setAppVersion(resources.version);
            this.buildRevision = String(resources.build_revision || "").trim();
            this.issueTrackerUrl = String(resources.issue_tracker_url || "").trim();
            this.manualDefinitions = resources.manual_criteria || {};
            this.resetManualItems();
        },
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});
