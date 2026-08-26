function ariane() {
    return {
        // shared
        mode: "single",
        appVersion: "",

        // single mode
        gene: "",
        configuredGenes: [],
        transcriptByGene: {},
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

        // reference material and manually reviewed evidence
        manualDefinitions: {},
        resourceLinks: [],
        resourceError: "",
        splicePs1Candidates: { status: "candidate_discovery_only", candidates: [] },
        manualItems: [],
        manualAssessor: "",
        manualAssessedAt: new Date().toISOString().slice(0, 10),
        manualLoading: false,
        manualError: "",
        manualResult: null,
        ps1ReferenceLoading: false,
        ps1ReferenceError: "",
        ps1ReferenceMessage: "",

        // batch mode
        batchText: "",
        batchParsed: [],
        batchParseError: "",
        batchRunning: false,
        batchDone: 0,
        batchTotal: 0,
        batchResults: [],

        // ENIGMA rules explorer
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

        async init() {
            this.resetManualItems();
            try {
                const response = await fetch("/api/resources");
                if (response.ok) {
                    const resources = await response.json();
                    this.setAppVersion(resources.version);
                    this.configuredGenes = resources.genes || [];
                    this.transcriptByGene = Object.fromEntries(
                        this.configuredGenes.map(item => [item.symbol, item.reference_transcript])
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

        async loadManualDefinitions() {
            if (!this.gene) return;
            const response = await fetch(`/api/resources?gene=${encodeURIComponent(this.gene)}`);
            if (!response.ok) throw new Error(`Manual evidence resources HTTP ${response.status}`);
            const resources = await response.json();
            this.setAppVersion(resources.version);
            this.manualDefinitions = resources.manual_criteria || {};
            this.resetManualItems();
        },

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
                const catalogResponse = await fetch("/api/rules");
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
                const response = await fetch(`/api/rules/trees/${encodeURIComponent(selectedId)}`);
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
                const response = await fetch(
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
            const response = await fetch(`/api/rules/tables/table9?${params.toString()}`);
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

        escapeXml(value) {
            return String(value ?? "")
                .replaceAll("&", "&amp;")
                .replaceAll("<", "&lt;")
                .replaceAll(">", "&gt;")
                .replaceAll('"', "&quot;")
                .replaceAll("'", "&#39;");
        },

        wrapGraphText(value, maximumCharacters = 32) {
            const words = String(value ?? "").split(/\s+/).filter(Boolean);
            const lines = [];
            let line = "";
            for (const word of words) {
                const candidate = line ? `${line} ${word}` : word;
                if (candidate.length > maximumCharacters && line) {
                    lines.push(line);
                    line = word;
                } else {
                    line = candidate;
                }
            }
            if (line) lines.push(line);
            return lines;
        },

        edgeDisplayLabel(value) {
            const labels = {
                start: "",
                yes: "Y",
                no: "N",
                impact: "IMPACT",
                no_impact: "NO IMPACT",
                not_informative: "NOT INFORMATIVE",
                no_impact_or_not_informative: "NO IMPACT / NOT INFORMATIVE",
            };
            return labels[value] ?? String(value ?? "").replaceAll("_", " ").toUpperCase();
        },

        edgeLabelLines(value) {
            const label = this.edgeDisplayLabel(value);
            if (label === "NO IMPACT / NOT INFORMATIVE") {
                return ["NO IMPACT OR", "NOT INFORMATIVE"];
            }
            return label ? [label] : [];
        },

        svgTextLines(value, x, y, maxCharacters = 27, lineHeight = 15, className = "", maximumLines = 4) {
            const lines = this.wrapGraphText(value, maxCharacters);
            const visibleLines = lines.slice(0, maximumLines);
            if (lines.length > maximumLines && visibleLines.length) {
                const last = visibleLines.length - 1;
                visibleLines[last] = `${visibleLines[last].replace(/[.,;:]$/, "")}...`;
            }
            return visibleLines.map((item, index) =>
                `<tspan x="${x}" y="${y + index * lineHeight}" class="${className}">${this.escapeXml(item)}</tspan>`
            ).join("");
        },

        decisionPathSvg(path) {
            if (!path) return "";
            const steps = path.steps || [];
            const entryTitles = {
                "missense-inframe": "Missense / in-frame",
                synonymous: "Synonymous (silent)",
                intronic: "Intronic",
            };
            const nodes = [
                {
                    id: `${path.branch_id}-entry`,
                    title: entryTitles[path.branch_id] || "Variant",
                    observed: "Variant type",
                    edgeResult: "start",
                    kind: "entry",
                },
                ...steps.map(item => ({
                    id: item.node_id,
                    title: item.question,
                    observed: item.observed,
                    edgeResult: item.result,
                    kind: "decision",
                })),
                {
                    id: path.outcome_node,
                    title: `Apply ${path.criterion}`,
                    observed: `${path.criterion} ${path.outcome}`,
                    result: "outcome",
                    kind: "outcome",
                },
            ];
            const nodeWidth = 230;
            const observedLinesByNode = Object.fromEntries(nodes.map(node => [
                node.id,
                this.wrapGraphText(node.observed, 35),
            ]));
            const titleLinesByNode = Object.fromEntries(nodes.map(node => [
                node.id,
                this.wrapGraphText(node.title, 29),
            ]));
            const nodeHeight = Math.max(118, ...nodes.map(node =>
                34 + titleLinesByNode[node.id].length * 16 + observedLinesByNode[node.id].length * 13
            ));
            const columns = Math.min(3, nodes.length);
            const rows = Math.ceil(nodes.length / columns);
            const horizontalGap = 92;
            const verticalGap = 88;
            const marginX = 28;
            const marginY = 34;
            const width = marginX * 2 + columns * nodeWidth + Math.max(0, columns - 1) * horizontalGap;
            const height = Math.max(
                250,
                marginY * 2 + rows * nodeHeight + Math.max(0, rows - 1) * verticalGap,
            );
            const nodePosition = index => {
                const row = Math.floor(index / columns);
                const offset = index % columns;
                const column = row % 2 === 0 ? offset : columns - 1 - offset;
                return {
                    x: marginX + column * (nodeWidth + horizontalGap),
                    y: marginY + row * (nodeHeight + verticalGap),
                    row,
                };
            };
            let edges = "";
            let cards = "";
            nodes.forEach((node, index) => {
                const { x, y, row } = nodePosition(index);
                if (index < nodes.length - 1) {
                    const next = nodePosition(index + 1);
                    const label = this.edgeDisplayLabel(node.edgeResult);
                    let route;
                    let labelX;
                    let labelY;
                    if (next.row === row) {
                        const movingRight = next.x > x;
                        const startX = movingRight ? x + nodeWidth : x;
                        const endX = movingRight ? next.x : next.x + nodeWidth;
                        const edgeY = y + nodeHeight / 2;
                        route = `M ${startX} ${edgeY} L ${endX} ${edgeY}`;
                        labelX = (startX + endX) / 2;
                        labelY = edgeY - 9;
                    } else {
                        const edgeX = x + nodeWidth / 2;
                        const startY = y + nodeHeight;
                        const endY = next.y;
                        route = `M ${edgeX} ${startY} L ${edgeX} ${endY}`;
                        labelX = edgeX + 13;
                        labelY = (startY + endY) / 2;
                    }
                    edges += `<path class="path-edge-active" d="${route}" marker-end="url(#path-arrow)"/>`;
                    if (label) {
                        edges += `<text class="path-edge-label" x="${labelX}" y="${labelY}" text-anchor="middle">${this.escapeXml(label)}</text>`;
                    }
                }
                const kindClass = node.kind === "entry"
                    ? "path-node-entry"
                    : (node.kind === "outcome" ? "path-node-outcome" : "path-node-decision");
                cards += `<g class="path-node ${kindClass}"><title>${this.escapeXml(node.observed)}</title>`;
                cards += `<rect x="${x}" y="${y}" width="${nodeWidth}" height="${nodeHeight}" rx="9"/>`;
                cards += `<text text-anchor="middle">${this.svgTextLines(node.title, x + nodeWidth / 2, y + 27, 29, 16, "path-node-title", titleLinesByNode[node.id].length)}</text>`;
                const observedStartY = y + nodeHeight - observedLinesByNode[node.id].length * 13 - 9;
                const observedMarkup = observedLinesByNode[node.id].map((line, lineIndex) =>
                    `<tspan x="${x + nodeWidth / 2}" y="${observedStartY + lineIndex * 13}">${this.escapeXml(line)}</tspan>`
                ).join("");
                cards += `<text x="${x + nodeWidth / 2}" y="${observedStartY}" text-anchor="middle" class="path-node-observed">${observedMarkup}</text>`;
                cards += `</g>`;
            });
            return `<svg class="decision-path-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMin meet" role="img" aria-label="Decision path for ${this.escapeXml(path.criterion)}">
                <defs><marker id="path-arrow" viewBox="0 0 10 10" refX="8.5" refY="5" markerUnits="userSpaceOnUse" markerWidth="8" markerHeight="8" orient="auto"><path d="M 0 1 L 8.5 5 L 0 9 z"/></marker></defs>
                ${edges}${cards}
            </svg>`;
        },

        treeEdgePoints(from, to) {
            const x1 = from.x + from.width;
            const y1 = from.y + from.height / 2;
            const x2 = to.x;
            const y2 = to.y + to.height / 2;
            const middle = x1 + Math.max(22, (x2 - x1) * 0.48);
            return [[x1, y1], [middle, y1], [middle, y2], [x2, y2]];
        },

        roundedOrthogonalPath(points, radius = 10) {
            const clean = [];
            for (const point of points || []) {
                const current = [Number(point[0]), Number(point[1])];
                const previous = clean[clean.length - 1];
                if (!previous || previous[0] !== current[0] || previous[1] !== current[1]) {
                    clean.push(current);
                }
            }
            if (clean.length < 2) return "";

            let route = `M ${clean[0][0]} ${clean[0][1]}`;
            for (let index = 1; index < clean.length - 1; index += 1) {
                const previous = clean[index - 1];
                const corner = clean[index];
                const next = clean[index + 1];
                const incomingHorizontal = previous[1] === corner[1];
                const outgoingHorizontal = corner[1] === next[1];
                const collinear = (incomingHorizontal && outgoingHorizontal)
                    || (!incomingHorizontal && !outgoingHorizontal);
                if (collinear) {
                    route += ` L ${corner[0]} ${corner[1]}`;
                    continue;
                }

                const incomingLength = Math.abs(corner[0] - previous[0]) + Math.abs(corner[1] - previous[1]);
                const outgoingLength = Math.abs(next[0] - corner[0]) + Math.abs(next[1] - corner[1]);
                const bend = Math.min(radius, incomingLength / 2, outgoingLength / 2);
                const before = [
                    corner[0] - Math.sign(corner[0] - previous[0]) * bend,
                    corner[1] - Math.sign(corner[1] - previous[1]) * bend,
                ];
                const after = [
                    corner[0] + Math.sign(next[0] - corner[0]) * bend,
                    corner[1] + Math.sign(next[1] - corner[1]) * bend,
                ];
                route += ` L ${before[0]} ${before[1]} Q ${corner[0]} ${corner[1]} ${after[0]} ${after[1]}`;
            }
            const last = clean[clean.length - 1];
            return `${route} L ${last[0]} ${last[1]}`;
        },

        treeEdgeLabelPosition(points, edge) {
            if (edge.label_at?.length >= 2) {
                return { x: edge.label_at[0], y: edge.label_at[1] };
            }
            const last = points[points.length - 1];
            const beforeLast = points[points.length - 2];
            if (!last || !beforeLast) return { x: 0, y: 0 };
            if (last[1] === beforeLast[1]) {
                return { x: (last[0] + beforeLast[0]) / 2, y: last[1] - 8 };
            }
            return { x: last[0] + 9, y: (last[1] + beforeLast[1]) / 2 };
        },

        decisionTreeSvg(branch, path = null) {
            if (!branch?.layout) return "";
            const layout = branch.layout;
            const width = layout.width;
            const height = layout.height;
            const defaultNodeWidth = layout.node_width;
            const defaultNodeHeight = layout.node_height;
            const positions = layout.positions || {};
            const nodeBox = (nodeId) => {
                const position = positions[nodeId] || [0, 0];
                return {
                    x: position[0],
                    y: position[1],
                    width: position[2] || defaultNodeWidth,
                    height: position[3] || defaultNodeHeight,
                };
            };
            const pathIds = path && path.branch_id === branch.id
                ? [branch.entry_node, ...(path.steps || []).map(item => item.node_id), path.outcome_node]
                : [];
            const observed = Object.fromEntries((path?.steps || []).map(item => [item.node_id, item.observed]));
            const observedAnnotations = Object.fromEntries(
                Object.entries(observed).map(([nodeId, value]) => {
                    const box = nodeBox(nodeId);
                    const maximumCharacters = Math.max(24, Math.floor((box.width + 58) / 5.4));
                    const lines = this.wrapGraphText(value, maximumCharacters);
                    const longestLine = Math.max(1, ...lines.map(line => line.length));
                    return [nodeId, {
                        lines,
                        width: Math.min(box.width + 58, Math.max(82, longestLine * 5.4 + 18)),
                        height: Math.max(20, lines.length * 12 + 9),
                    }];
                })
            );
            const topInset = Math.max(0, ...Object.entries(observedAnnotations).map(([nodeId, annotation]) => {
                const box = nodeBox(nodeId);
                return annotation.height + 7 - box.y;
            })) + (Object.keys(observedAnnotations).length ? 8 : 0);
            const activeEdge = (edge) => pathIds.some((id, index) => id === edge.from && pathIds[index + 1] === edge.to);
            let edgeMarkup = "";
            for (const edge of branch.edges || []) {
                if (!positions[edge.from] || !positions[edge.to]) continue;
                const from = nodeBox(edge.from);
                const to = nodeBox(edge.to);
                const active = activeEdge(edge);
                const labelLines = this.edgeLabelLines(edge.result);
                const routePoints = edge.points?.length
                    ? edge.points
                    : this.treeEdgePoints(from, to);
                const route = this.roundedOrthogonalPath(routePoints);
                edgeMarkup += `<path class="tree-connector ${active ? 'tree-connector-active' : ''}" d="${route}" marker-end="url(#tree-arrow${active ? '-active' : ''})"/>`;
                if (labelLines.length) {
                    const labelPosition = this.treeEdgeLabelPosition(routePoints, edge);
                    const labelX = labelPosition.x;
                    const labelY = labelPosition.y;
                    const labelWidth = Math.max(24, Math.max(...labelLines.map(line => line.length)) * 5.5 + 10);
                    const labelHeight = labelLines.length * 11 + 7;
                    const lineMarkup = labelLines.map((line, index) =>
                        `<tspan x="${labelX}" y="${labelY + index * 11}">${this.escapeXml(line)}</tspan>`
                    ).join("");
                    edgeMarkup += `<g class="tree-edge-label-group ${active ? 'tree-edge-label-active' : ''}">`;
                    edgeMarkup += `<rect class="tree-edge-label-bg" x="${labelX - labelWidth / 2}" y="${labelY - 10}" width="${labelWidth}" height="${labelHeight}" rx="6"/>`;
                    edgeMarkup += `<text class="tree-edge-label" x="${labelX}" y="${labelY}" text-anchor="middle">${lineMarkup}</text></g>`;
                }
            }
            let nodeMarkup = "";
            for (const node of branch.nodes || []) {
                const position = positions[node.id];
                if (!position) continue;
                const { x, y, width: nodeWidth, height: nodeHeight } = nodeBox(node.id);
                const active = pathIds.includes(node.id);
                const nodeClass = node.kind === "entry"
                    ? "graph-node-entry"
                    : (node.kind === "outcome" ? "graph-node-outcome" : "graph-node-decision");
                const visualClass = node.visual ? `graph-node-${this.escapeXml(node.visual)}` : "";
                nodeMarkup += `<g class="graph-node ${nodeClass} ${visualClass} ${active ? 'graph-node-active' : ''}"><title>${this.escapeXml(observed[node.id] || node.label)}</title>`;
                nodeMarkup += `<rect x="${x}" y="${y}" width="${nodeWidth}" height="${nodeHeight}" rx="10"/>`;
                const titleWidth = Math.max(18, Math.floor(nodeWidth / 7.2));
                const maximumTitleLines = Math.max(3, Math.floor((nodeHeight - 18) / 16));
                nodeMarkup += `<text text-anchor="middle">${this.svgTextLines(node.label, x + nodeWidth / 2, y + 27, titleWidth, 16, "graph-node-title", maximumTitleLines)}</text>`;
                if (observed[node.id]) {
                    const annotation = observedAnnotations[node.id];
                    const badgeX = x + nodeWidth / 2 - annotation.width / 2;
                    const badgeY = y - annotation.height - 7;
                    const textX = badgeX + annotation.width / 2;
                    const textY = badgeY + 13;
                    const textLines = annotation.lines.map((line, index) =>
                        `<tspan x="${textX}" y="${textY + index * 12}">${this.escapeXml(line)}</tspan>`
                    ).join("");
                    nodeMarkup += `<g class="graph-node-observed"><rect x="${badgeX}" y="${badgeY}" width="${annotation.width}" height="${annotation.height}" rx="7"/>`;
                    nodeMarkup += `<text x="${textX}" y="${textY}" text-anchor="middle">${textLines}</text></g>`;
                }
                nodeMarkup += `</g>`;
            }
            const renderedHeight = height + topInset;
            return `<svg class="decision-tree-svg ${pathIds.length ? 'tree-has-active-path' : ''}" viewBox="0 0 ${width} ${renderedHeight}" preserveAspectRatio="xMidYMin meet" style="aspect-ratio: ${width} / ${renderedHeight}" role="img" aria-label="${this.escapeXml(branch.title)} decision tree">
                <defs>
                    <marker id="tree-arrow" viewBox="0 0 8 8" refX="7" refY="4" markerUnits="userSpaceOnUse" markerWidth="5.5" markerHeight="5.5" orient="auto"><path d="M 0 1 L 7 4 L 0 7 z"/></marker>
                    <marker id="tree-arrow-active" viewBox="0 0 8 8" refX="7" refY="4" markerUnits="userSpaceOnUse" markerWidth="6" markerHeight="6" orient="auto"><path d="M 0 1 L 7 4 L 0 7 z"/></marker>
                </defs>
                <g transform="translate(0 ${topInset})">${edgeMarkup}${nodeMarkup}</g>
            </svg>`;
        },

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
            return result;
        },

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
                        reference_variant: `${candidate.gene} ${candidate.c_notation}`,
                        reference_p_notation: candidate.p_notation || "",
                        reference_classification: candidate.classification || "",
                        classification_source: candidate.classification_source || "",
                        classification_verification: "",
                        same_missense_confirmed: true,
                        different_nucleotide_change_confirmed: true,
                        vua_spliceai_score:
                            review.vua_spliceai_score ?? this.result?.spliceai_audit?.score ?? "",
                        reference_spliceai_score:
                            review.reference_spliceai_scores?.[candidate.c_notation] ?? "",
                        splice_source_check_completed: false,
                        splice_sources_checked: review.splice_sources_checked || [],
                        vua_confirmed_splice_status: review.vua_splice_evidence_status || "not_assessed",
                        reference_confirmed_splice_status: "not_assessed",
                        reference_classification_used_ps1: "unknown",
                        reference_ps1_dependency_reference: "",
                        direct_reciprocal_dependency_excluded: false,
                        ps1_protein_rationale: "",
                    };
                    item.references = [
                        candidate.source_dataset,
                        candidate.classification_source,
                        review.source_url,
                    ].filter(Boolean).join("\n");
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
                const response = await fetch("/api/manual-evidence/resolve-ps1-reference", {
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
                item.evidence = {
                    ...item.evidence,
                    reference_variant: `${resolved.reference.gene} ${resolved.reference.c_notation}`,
                    reference_p_notation: resolved.reference.p_notation,
                    reference_classification: resolved.classification || "",
                    classification_verification: resolved.classification_verification || "unresolved",
                    classification_source: resolved.classification_source || "",
                    same_missense_confirmed: resolved.same_missense_substitution === true,
                    different_nucleotide_change_confirmed: resolved.different_nucleotide_change === true,
                    vua_spliceai_score: resolved.assessed.spliceai_score ?? "",
                    reference_spliceai_score: resolved.reference.spliceai_score ?? "",
                };
                const references = [
                    ...(item.references || "").split(/\r?\n/),
                    ...(resolved.references || []),
                ].map(value => value.trim()).filter(Boolean);
                item.references = [...new Set(references)].join("\n");
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
                const response = await fetch("/api/manual-evidence/evaluate", {
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

        logClientValidation(error, form = "single", input = null) {
            const submittedInput = input || {
                gene: this.gene,
                c_notation: this.c_notation.trim(),
                p_notation: this.p_notation.trim() || null,
                assembly: this.assembly || null,
                dup_type: this.dup_type,
            };
            fetch("/api/audit/client-validation", {
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
                const resp = await fetch("/api/classify", {
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

        // ── Batch: parse CSV text whenever batchText changes ──────────────
        parseBatch() {
            this.batchParseError = "";
            this.batchParsed = [];
            if (!this.batchText.trim()) return;

            const lines = this.batchText.trim().split("\n").filter(l => l.trim());
            const parsed = [];
            const errors = [];

            for (let i = 0; i < lines.length; i++) {
                // Skip header lines
                if (lines[i].toLowerCase().startsWith("gene") || lines[i].toLowerCase().startsWith("#")) continue;

                const parts = lines[i].split(",").map(s => s.trim());
                if (parts.length < 2) {
                    errors.push(`Line ${i + 1}: need at least gene and c. notation`);
                    continue;
                }
                const gene = parts[0].toUpperCase();
                const allowedGenes = this.configuredGenes.map(item => item.symbol);
                if (!allowedGenes.includes(gene)) {
                    errors.push(`Line ${i + 1}: gene must be one of ${allowedGenes.join(", ")}`);
                    continue;
                }

                let cRaw = parts[1];
                let pRaw = parts[2] || "";
                let assemblyRaw = "";
                if (/^(?:GRCh3[78]|hg(?:19|38))$/i.test(pRaw)) {
                    assemblyRaw = /^hg19$/i.test(pRaw) ? "GRCh37"
                        : /^hg38$/i.test(pRaw) ? "GRCh38"
                        : pRaw.replace(/^grch/i, "GRCh");
                    pRaw = "";
                }
                const dupRaw = parts[3] || "Unknown";

                // Handle appended protein notation, including "(p.Val2050del)".
                if (!pRaw) {
                    const combined = cRaw.match(/^(c\.\S+)\s+(\(?p\..+\)?)$/i);
                    if (combined) {
                        cRaw = combined[1];
                        pRaw = combined[2];
                    }
                }
                if (pRaw.startsWith("(p.") && pRaw.endsWith(")")) {
                    pRaw = pRaw.substring(1, pRaw.length - 1);
                }
                if (pRaw.startsWith("p.") && !pRaw.startsWith("p.(")) {
                    pRaw = `p.(${pRaw.substring(2)})`;
                }

                parsed.push({
                    gene,
                    c_notation: cRaw,
                    p_notation: pRaw,
                    assembly: assemblyRaw,
                    dup_type: dupRaw,
                });
            }

            if (errors.length > 0) {
                this.batchParseError = errors.join("; ");
            }
            this.batchParsed = parsed;
        },

        // ── Batch: classify all parsed variants ───────────────────────────
        async classifyBatch() {
            this.parseBatch();
            if (this.batchParseError) {
                this.logClientValidation(
                    this.batchParseError,
                    "batch",
                    { batch_text: this.batchText.slice(0, 4000) },
                );
            }
            if (this.batchParsed.length === 0) return;

            this.batchRunning = true;
            this.batchDone = 0;
            this.batchTotal = this.batchParsed.length;
            this.batchResults = new Array(this.batchParsed.length).fill(null);

            // Run up to 3 classify calls concurrently for rate-limit safety
            const CONCURRENCY = 3;
            const queue = [...this.batchParsed.entries()];
            const active = new Set();

            const runNext = async () => {
                if (queue.length === 0) return;
                const [idx, item] = queue.shift();
                const task = (async () => {
                    try {
                        const resp = await fetch("/api/classify", {
                            method: "POST",
                            headers: { "Content-Type": "application/json" },
                            body: JSON.stringify({
                                gene: item.gene,
                                c_notation: item.c_notation,
                                p_notation: item.p_notation || null,
                                assembly: item.assembly || null,
                                dup_type: item.dup_type || "Unknown",
                            }),
                        });
                        if (resp.ok) {
                            const data = this.normalizeCriterionOrder(await resp.json());
                            this.batchResults[idx] = {
                                status: "ok",
                                gene: item.gene,
                                c_notation: data.c_notation,
                                p_notation: data.p_notation,
                                result: data,
                            };
                        } else {
                            const err = await resp.json().catch(() => ({}));
                            this.batchResults[idx] = {
                                status: "error",
                                gene: item.gene,
                                c_notation: item.c_notation,
                                p_notation: item.p_notation,
                                error: this.formatApiError(err, resp.status),
                            };
                        }
                    } catch (e) {
                        this.batchResults[idx] = {
                            status: "error",
                            gene: item.gene,
                            c_notation: item.c_notation,
                            p_notation: item.p_notation,
                            error: "Network error",
                        };
                    }
                    this.batchDone++;
                    active.delete(task);
                    await runNext();
                })();
                active.add(task);
            };

            // Kick off initial workers
            const starters = [];
            for (let i = 0; i < Math.min(CONCURRENCY, this.batchParsed.length); i++) {
                starters.push(runNext());
            }
            await Promise.all(starters);

            // Wait for all active tasks to finish
            while (active.size > 0) {
                await Promise.race(active);
            }

            this.batchRunning = false;
        },

        // ── Batch: download results as CSV ────────────────────────────────
        downloadBatchCsv() {
            const header = [
                "Gene", "c_notation", "p_notation",
                "Class", "Label", "Points",
                "Criteria", "Excluded_criteria", "ClinVar", "ENIGMA_EP",
                "Classification_note", "VUS_category", "VUS_what_to_check",
                "RNA_review", "RNA_branches",
                "Splice_PS1_review", "Splice_PS1_branches",
                "Protein_PS1_review", "Protein_PS1_branches",
                "PVS1_INIT_review", "PVS1_INIT_branches",
                "Warnings"
            ];
            const rows = this.batchResults.map(row => {
                if (!row || row.status === "error") {
                    return [
                        row ? row.gene : "",
                        row ? row.c_notation : "",
                        row ? (row.p_notation || "") : "",
                        "ERROR",
                        ...Array(header.length - 5).fill(""),
                        row ? (row.error || "") : "",
                    ];
                }
                const r = row.result;
                const criteria = r.criteria
                    .filter(c => c.applies)
                    .map(c => c.name + (c.strength ? "_" + c.strength.replace(" ", "") : ""))
                    .join("; ");
                const excludedCriteria = (r.excluded_criteria || [])
                    .map(c => `${c.name}${c.strength ? "_" + c.strength.replace(" ", "") : ""}: ${c.reason}`)
                    .join("; ");
                const warnings = r.warnings.join(" | ");
                const rnaReview = r.rna_review && r.rna_review.recommended
                    ? `yes/${r.rna_review.priority || ""}`
                    : "";
                const rnaBranches = r.rna_review && r.rna_review.recommended
                    ? (r.rna_review.potential_branches || []).join("; ")
                    : "";
                const splicePs1Review = r.splice_ps1_review && r.splice_ps1_review.recommended
                    ? `yes/${r.splice_ps1_review.priority || ""}`
                    : "";
                const splicePs1Branches = r.splice_ps1_review && r.splice_ps1_review.recommended
                    ? (r.splice_ps1_review.potential_branches || []).join("; ")
                    : "";
                const proteinPs1Review = r.protein_ps1_review && r.protein_ps1_review.recommended
                    ? `yes/${r.protein_ps1_review.priority || ""}`
                    : "";
                const proteinPs1Branches = r.protein_ps1_review && r.protein_ps1_review.recommended
                    ? (r.protein_ps1_review.potential_branches || []).join("; ")
                    : "";
                const initiationReview = r.initiation_review && r.initiation_review.recommended
                    ? `yes/${r.initiation_review.priority || ""}`
                    : "";
                const initiationBranches = r.initiation_review && r.initiation_review.recommended
                    ? (r.initiation_review.potential_branches || []).join("; ")
                    : "";
                return [
                    row.gene,
                    row.c_notation,
                    row.p_notation || "",
                    r.predicted_class,
                    r.predicted_label,
                    r.total_points,
                    criteria,
                    excludedCriteria,
                    r.external ? r.external.clinvar_classification : "",
                    r.external ? r.external.enigma_ep_class : "",
                    r.classification_note,
                    r.vus_explanation ? r.vus_explanation.category : "",
                    r.vus_explanation ? r.vus_explanation.what_to_check : "",
                    rnaReview,
                    rnaBranches,
                    splicePs1Review,
                    splicePs1Branches,
                    proteinPs1Review,
                    proteinPs1Branches,
                    initiationReview,
                    initiationBranches,
                    warnings,
                ];
            });

            const escape = v => '"' + String(v).replace(/"/g, '""') + '"';
            const csv = [header, ...rows]
                .map(row => row.map(escape).join(","))
                .join("\n");

            const blob = new Blob(["﻿" + csv], { type: "text/csv;charset=utf-8" });
            const url  = URL.createObjectURL(blob);
            const a    = document.createElement("a");
            a.href     = url;
            a.download = "ariane_batch_results.csv";
            a.click();
            URL.revokeObjectURL(url);
        },
    };
}
