(function registerGraphs(namespace) {
    "use strict";

    namespace.graphsMethods = {
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
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});

