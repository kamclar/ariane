(function registerCore(namespace) {
    "use strict";

    namespace.coreState = function coreState() {
        return {
        mode: "single",
        appVersion: "",
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
            this.resetManualItems();
            try {
                const response = await namespace.api.request("/api/resources");
                if (response.ok) {
                    const resources = await response.json();
                    this.setAppVersion(resources.version);
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

        async loadManualDefinitions() {
            if (!this.gene) return;
            const response = await namespace.api.request(`/api/resources?gene=${encodeURIComponent(this.gene)}`);
            if (!response.ok) throw new Error(`Manual evidence resources HTTP ${response.status}`);
            const resources = await response.json();
            this.setAppVersion(resources.version);
            this.manualDefinitions = resources.manual_criteria || {};
            this.resetManualItems();
        },
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});

