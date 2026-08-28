(function registerRules(namespace) {
    "use strict";

    namespace.rulesState = function rulesState() {
        return {
        rulesCatalog: null,
        rulesTree: null,
        rulesTreeCache: {},
        selectedRuleTreeId: "figure-1a",
        rulesLoading: false,
        rulesError: "",
        rulesView: "trees",
        highlightedRulePath: null,
        selectedRuleBranch: "missense-inframe",
        table9Gene: "",
        table9Query: "",
        table9Code: "",
        table9Page: 1,
        table9PageSize: 25,
        table9Total: 0,
        table9Items: [],
        selectedReferenceTableRole: "used_by_ariane",
        selectedReferenceTableId: "specification-table-3",
        selectedReferenceSectionId: "main",
        referenceTableQuery: "",
        referenceTablePage: 1,
        referenceTablePageSize: 25,
        referenceTableTotal: 0,
        referenceTableColumns: [],
        referenceTableItems: [],
        referenceTableLoading: false,
        referenceTableError: "",
        };
    };

    namespace.rulesMethods = {
        async openRules(decisionPath = null) {
            const resultPath = (this.result?.criteria || []).find(item => item.decision_path)?.decision_path || null;
            const selectedPath = decisionPath || resultPath || this.highlightedRulePath || null;
            this.mode = "rules";
            this.rulesView = selectedPath ? "trees" : this.rulesView;
            this.highlightedRulePath = selectedPath;
            if (selectedPath?.branch_id) this.selectedRuleBranch = selectedPath.branch_id;
            const requestedTreeId = selectedPath?.tree_id || this.selectedRuleTreeId || "figure-1a";
            if (this.rulesCatalog) {
                await this.selectRuleTree(requestedTreeId, selectedPath?.branch_id || null);
                return;
            }
            this.rulesLoading = true;
            this.rulesError = "";
            try {
                const catalogResponse = await namespace.api.request("/api/rules");
                if (!catalogResponse.ok) throw new Error(`HTTP ${catalogResponse.status}`);
                this.rulesCatalog = await catalogResponse.json();
                await this.selectRuleTree(requestedTreeId, selectedPath?.branch_id || null);
                await this.ensureReferenceTableSelected();
            } catch (error) {
                this.rulesError = `ENIGMA rules could not be loaded: ${error?.message || error}`;
            } finally {
                this.rulesLoading = false;
            }
        },

        async selectRuleTree(treeId, preferredBranch = null) {
            const available = (this.rulesCatalog?.decision_trees || []).some(tree => tree.id === treeId);
            const selectedId = available ? treeId : "figure-1a";
            this.selectedRuleTreeId = selectedId;
            if (!this.rulesTreeCache[selectedId]) {
                const response = await namespace.api.request(`/api/rules/trees/${encodeURIComponent(selectedId)}`);
                if (!response.ok) throw new Error(`Decision tree HTTP ${response.status}`);
                this.rulesTreeCache[selectedId] = await response.json();
            }
            this.rulesTree = this.rulesTreeCache[selectedId];
            const branchExists = (this.rulesTree?.branches || []).some(branch => branch.id === preferredBranch);
            this.selectedRuleBranch = branchExists
                ? preferredBranch
                : (this.rulesTree?.branches?.[0]?.id || "");
        },

        selectedRuleTreeSummary() {
            return (this.rulesCatalog?.decision_trees || []).find(tree => tree.id === this.selectedRuleTreeId) || null;
        },

        ruleTreeOptionLabel(tree) {
            const prefix = tree?.diagram_provenance === "official_redraw" ? "ENIGMA figure" : "ARIANE derived";
            return `${prefix}: ${tree?.title || ""}`;
        },

        selectedRuleOriginalFigures() {
            const ids = new Set(this.rulesTree?.original_figure_ids || []);
            return (this.rulesCatalog?.figures || []).filter(figure => ids.has(figure.id));
        },

        referenceTableRoles() {
            return this.rulesCatalog?.tables?.roles || [];
        },

        referenceTablesForRole(roleId = this.selectedReferenceTableRole) {
            return (this.rulesCatalog?.tables?.items || []).filter(table => table.role === roleId);
        },

        referenceTableOptionLabel(table) {
            const source = {
                specification: "Specification",
                appendix: "Appendix",
                supplementary: "Supplementary",
            }[table?.group] || "ENIGMA";
            return `${source} Table ${table?.number}`;
        },

        referenceTableUsageLabel(table = this.selectedReferenceTable()) {
            const labels = {
                rule_definition: "Rule definition used by ARIANE",
                runtime_rule: "Automated rule used by ARIANE",
                runtime_lookup: "Automated lookup used by ARIANE",
                candidate_registry: "Candidate registry used by ARIANE",
                expert_review: "Expert review",
                rule_support: "Supporting evidence",
                calibration_reference: "Calibration reference",
            };
            return labels[table?.usage] || String(table?.usage || "").replaceAll("_", " ");
        },

        selectedReferenceTable() {
            return (this.rulesCatalog?.tables?.items || []).find(
                table => table.id === this.selectedReferenceTableId
            ) || null;
        },

        selectedReferenceSection() {
            const table = this.selectedReferenceTable();
            return (table?.sections || []).find(
                section => section.id === this.selectedReferenceSectionId
            ) || table?.sections?.[0] || null;
        },

        referenceTableSourceUrl() {
            const sourceId = this.selectedReferenceTable()?.source_id;
            return (this.rulesCatalog?.sources || []).find(source => source.id === sourceId)?.official_url || "";
        },

        async ensureReferenceTableSelected() {
            if (!this.rulesCatalog?.tables?.items?.length) return;
            let table = this.selectedReferenceTable();
            if (!table) {
                table = this.referenceTablesForRole()[0] || this.rulesCatalog.tables.items[0];
                this.selectedReferenceTableId = table.id;
                this.selectedReferenceTableRole = table.role;
            }
            const sectionExists = (table.sections || []).some(
                section => section.id === this.selectedReferenceSectionId
            );
            this.selectedReferenceSectionId = sectionExists
                ? this.selectedReferenceSectionId
                : (table.sections?.[0]?.id || "main");
            if (!this.referenceTableItems.length) await this.searchReferenceTable(1);
        },

        async selectReferenceTableRole(roleId) {
            this.selectedReferenceTableRole = roleId;
            const table = this.referenceTablesForRole(roleId)[0];
            if (table) await this.selectReferenceTable(table.id);
        },

        async selectReferenceTable(tableId) {
            const table = (this.rulesCatalog?.tables?.items || []).find(item => item.id === tableId);
            if (!table) return;
            this.selectedReferenceTableId = table.id;
            this.selectedReferenceTableRole = table.role;
            this.selectedReferenceSectionId = table.sections?.[0]?.id || "main";
            this.referenceTableQuery = "";
            await this.searchReferenceTable(1);
        },

        async selectReferenceSection(sectionId) {
            this.selectedReferenceSectionId = sectionId;
            await this.searchReferenceTable(1);
        },

        async searchReferenceTable(page = 1) {
            if (!this.selectedReferenceTableId) return;
            this.referenceTableLoading = true;
            this.referenceTableError = "";
            this.referenceTablePage = Math.max(1, Number(page) || 1);
            const params = new URLSearchParams({
                page: String(this.referenceTablePage),
                page_size: String(this.referenceTablePageSize),
                section: this.selectedReferenceSectionId,
            });
            if (this.referenceTableQuery.trim()) params.set("query", this.referenceTableQuery.trim());
            try {
                const response = await namespace.api.request(
                    `/api/rules/tables/${encodeURIComponent(this.selectedReferenceTableId)}?${params.toString()}`
                );
                if (!response.ok) throw new Error(`HTTP ${response.status}`);
                const payload = await response.json();
                this.referenceTableItems = payload.items || [];
                this.referenceTableColumns = payload.columns || [];
                this.referenceTableTotal = payload.total || 0;
                this.referenceTablePage = payload.page || 1;
            } catch (error) {
                this.referenceTableItems = [];
                this.referenceTableColumns = [];
                this.referenceTableTotal = 0;
                this.referenceTableError = `ENIGMA table could not be loaded: ${error?.message || error}`;
            } finally {
                this.referenceTableLoading = false;
            }
        },

        referenceTableLastPage() {
            return Math.max(1, Math.ceil(this.referenceTableTotal / this.referenceTablePageSize));
        },

        referenceCellText(cell) {
            if (cell && typeof cell === "object" && !Array.isArray(cell)) {
                const value = cell.value;
                return String(value ?? cell.formula ?? "");
            }
            return String(cell ?? "");
        },

        referenceCellFormula(cell) {
            return cell && typeof cell === "object" && !Array.isArray(cell)
                ? String(cell.formula || "")
                : "";
        },

        async searchTable9(page = 1) {
            this.table9Page = Math.max(1, Number(page) || 1);
            const params = new URLSearchParams({
                page: String(this.table9Page),
                page_size: String(this.table9PageSize),
            });
            if (this.table9Gene) params.set("gene", this.table9Gene);
            if (this.table9Query.trim()) params.set("query", this.table9Query.trim());
            if (this.table9Code) params.set("code", this.table9Code);
            const response = await namespace.api.request(`/api/rules/tables/table9?${params.toString()}`);
            if (!response.ok) throw new Error(`Table 9 HTTP ${response.status}`);
            const payload = await response.json();
            this.table9Items = payload.items || [];
            this.table9Total = payload.total || 0;
            this.table9Page = payload.page || 1;
        },

        table9LastPage() {
            return Math.max(1, Math.ceil(this.table9Total / this.table9PageSize));
        },

        treeEdgesFrom(branch, nodeId) {
            return (branch?.edges || []).filter(edge => edge.from === nodeId);
        },

        treeNodeLabel(branch, nodeId) {
            return (branch?.nodes || []).find(node => node.id === nodeId)?.label || nodeId;
        },

        isHighlightedRuleNode(branchId, nodeId) {
            const path = this.highlightedRulePath;
            if (!path || path.branch_id !== branchId) return false;
            return path.outcome_node === nodeId || (path.steps || []).some(item => item.node_id === nodeId);
        },

        resultDecisionPaths() {
            return (this.result?.criteria || []).filter(item => item.decision_path);
        },

        selectResultDecisionPath(criterion) {
            if (!criterion?.decision_path) return;
            this.highlightedRulePath = criterion.decision_path;
            this.rulesView = "trees";
            this.selectRuleTree(
                criterion.decision_path.tree_id || "figure-1a",
                criterion.decision_path.branch_id,
            ).catch(error => {
                this.rulesError = `Decision tree could not be loaded: ${error?.message || error}`;
            });
        },

        selectedTreeBranch() {
            return (this.rulesTree?.branches || []).find(branch => branch.id === this.selectedRuleBranch)
                || this.rulesTree?.branches?.[0]
                || null;
        },

        ruleBranchTitle(branchId) {
            return (this.rulesTree?.branches || []).find(branch => branch.id === branchId)?.title || branchId;
        },
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});

