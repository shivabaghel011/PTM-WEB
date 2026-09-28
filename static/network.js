/**
 * PTM-Net — Network Visualisation & Table Sync Script
 */

document.addEventListener('DOMContentLoaded', function() {
    
    // Globals injected by Flask template:
    // - NETWORK_DATA: { nodes: [...], edges: [...] }
    // - QUERY_GENES: ['GENE1', 'GENE2', ...]
    // - INIT_THRESHOLD: float
    
    var cy = null;

    // ──────────────────────────────────────────────────────────────
    // 1. HELPERS & FORMATTERS
    // ──────────────────────────────────────────────────────────────
    
    // Generates an SVG data-URI showing the PTM modification abbreviation inside the node
    function makeSvgBg(modLabel) {
        if (!modLabel) return 'none';
        var svgStr = '<svg xmlns="http://www.w3.org/2000/svg" width="50" height="50">'
                   + '<text x="50%" y="50%" font-family="Inter, system-ui, sans-serif" font-size="13" font-weight="900" '
                   + 'fill="#ffffff" text-anchor="middle" dominant-baseline="central">'
                   + modLabel + '</text></svg>';
        return 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgStr);
    }

    // Splits and groups modifications beside their respective Writer/Eraser annotations
    function splitRoleMods(weVal, modVal) {
        if (!weVal || weVal === 'NA') return null;
        var roles = weVal.split(';').map(function(r) { return r.trim().toUpperCase(); });
        var mods = modVal && modVal !== 'NA' ? modVal.split(';').map(function(m) { return m.trim(); }) : [];
        
        var writerMods = [];
        var eraserMods = [];
        mods.forEach(function(m) {
            // Heuristic: modifications starting with "De" are associated with Erasers
            if (/^de/i.test(m) || /removal/i.test(m)) {
                eraserMods.push(m);
            } else {
                writerMods.push(m);
            }
        });
        
        var isW = roles.indexOf('W') >= 0 || roles.indexOf('WRITER') >= 0 || roles.indexOf('W/E') >= 0;
        var isE = roles.indexOf('E') >= 0 || roles.indexOf('ERASER') >= 0 || roles.indexOf('W/E') >= 0;
        
        return {
            isW: isW,
            isE: isE,
            writerMods: writerMods,
            eraserMods: eraserMods
        };
    }

    // Formats PTM string cleanly by removing pipes and placing each PTM on a separate sub-row
    function formatPtmsHtml(rawPtmStr) {
        if (!rawPtmStr || rawPtmStr === 'NA' || rawPtmStr === '—' || rawPtmStr.trim() === '') {
            return '<span class="text-dim">—</span>';
        }

        var items = [];
        if (rawPtmStr.indexOf(' | ') >= 0) {
            items = rawPtmStr.split(' | ');
        } else if (rawPtmStr.indexOf(';') >= 0) {
            items = rawPtmStr.split(';');
        } else {
            items = [rawPtmStr];
        }

        var subRows = [];
        items.forEach(function(item) {
            item = item.trim();
            if (!item || item === 'NA' || item === '—') return;

            if (item.indexOf('|') >= 0) {
                var parts = item.split('|').map(function(p){ return p.trim(); }).filter(Boolean);
                if (parts.length >= 3) {
                    var ptmName = parts[0];
                    var res = parts[1];
                    var pos = parts[2];
                    subRows.push(`<div class="ptm-subrow"><span class="ptm-name">${ptmName}</span> <span class="ptm-pos">(${res} ${pos})</span></div>`);
                } else if (parts.length === 2) {
                    subRows.push(`<div class="ptm-subrow"><span class="ptm-name">${parts[0]}</span> <span class="ptm-pos">(${parts[1]})</span></div>`);
                } else {
                    subRows.push(`<div class="ptm-subrow"><span class="ptm-name">${parts[0]}</span></div>`);
                }
            } else {
                subRows.push(`<div class="ptm-subrow"><span class="ptm-name">${item}</span></div>`);
            }
        });

        if (subRows.length === 0) {
            return '<span class="text-dim">—</span>';
        }
        return `<div class="ptm-list-container">${subRows.join('')}</div>`;
    }

    // Highlight a row in the active table and scroll it into view
    function highlightTableRow(geneId) {
        document.querySelectorAll('.node-row').forEach(function(r) {
            r.classList.remove('row-highlight');
        });
        var rows = document.querySelectorAll(`.node-row[data-gene="${geneId}"]`);
        rows.forEach(function(row) {
            row.classList.add('row-highlight');
            var parentTab = row.closest('.table-tab-content');
            if (parentTab && parentTab.classList.contains('active')) {
                row.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
            }
        });
    }

    // Scroll chaining: forwards wheel scroll to parent table-wrapper when cell hits boundary
    function attachScrollChaining(cellContent, tableWrapper) {
        if (!cellContent || !tableWrapper) return;
        cellContent.addEventListener('wheel', function(e) {
            var scrollTop = cellContent.scrollTop;
            var scrollHeight = cellContent.scrollHeight;
            var clientHeight = cellContent.clientHeight;
            var delta = e.deltaY;

            var isAtTop = (scrollTop <= 0);
            var isAtBottom = (scrollTop + clientHeight >= scrollHeight - 1);
            var cannotScroll = (scrollHeight <= clientHeight);

            if (cannotScroll || (delta < 0 && isAtTop) || (delta > 0 && isAtBottom)) {
                tableWrapper.scrollTop += delta;
            }
        }, { passive: true });
    }

    // Handles expanding / collapsing rows on click
    function handleRowClick(row, geneId) {
        var isExpanded = row.classList.contains('expanded');
        var table = row.closest('table');
        if (table) {
            table.querySelectorAll('.node-row.expanded').forEach(function(r) {
                if (r !== row) r.classList.remove('expanded');
            });
        }

        if (isExpanded) {
            row.classList.remove('expanded');
        } else {
            row.classList.add('expanded');
            var tableWrapper = row.closest('.table-wrapper');
            if (tableWrapper) {
                var rowTopRelative = row.offsetTop - tableWrapper.offsetTop;
                tableWrapper.scrollTo({
                    top: Math.max(0, rowTopRelative - 4),
                    behavior: 'smooth'
                });
            }
        }

        if (cy) {
            cy.elements().unselect();
            var node = cy.$id(geneId);
            if (node.length > 0) {
                node.select();
                cy.animate({
                    center: { eles: node },
                    zoom: 1.3
                }, { duration: 300 });
            }
        }
        highlightTableRow(geneId);
    }

    // Column resizer initialization for draggable table column borders
    function initTableColumnResizers(table) {
        if (!table) return;
        var headers = table.querySelectorAll('thead th');
        headers.forEach(function(th) {
            var existingResizer = th.querySelector('.col-resizer');
            if (existingResizer) {
                existingResizer.remove();
            }

            var resizer = document.createElement('div');
            resizer.className = 'col-resizer';
            th.appendChild(resizer);

            var startX, startWidth;

            resizer.addEventListener('mousedown', function(e) {
                e.stopPropagation();
                e.preventDefault();
                startX = e.pageX;
                startWidth = th.offsetWidth;
                resizer.classList.add('resizing');
                document.body.style.cursor = 'col-resize';
                document.body.style.userSelect = 'none';

                function onMouseMove(e) {
                    var newWidth = Math.max(80, startWidth + (e.pageX - startX));
                    th.style.width = newWidth + 'px';
                    th.style.minWidth = newWidth + 'px';
                }

                function onMouseUp() {
                    resizer.classList.remove('resizing');
                    document.body.style.cursor = '';
                    document.body.style.userSelect = '';
                    document.removeEventListener('mousemove', onMouseMove);
                    document.removeEventListener('mouseup', onMouseUp);
                }

                document.addEventListener('mousemove', onMouseMove);
                document.addEventListener('mouseup', onMouseUp);
            });
        });
    }

    // Rebuild the HTML tables dynamically after filter updates
    function rebuildTable(nodes) {
        var nodeCountBadge = document.querySelector('.topbar-stats .stat-chip:nth-child(1)');
        if (nodeCountBadge) {
            nodeCountBadge.textContent = `${nodes.length} nodes`;
        }

        // 1. General table
        var tableGeneral = document.getElementById('protein-table');
        var tbodyGeneral = tableGeneral ? tableGeneral.querySelector('tbody') : null;
        var wrapperGeneral = tableGeneral ? tableGeneral.closest('.table-wrapper') : null;

        if (tbodyGeneral) {
            tbodyGeneral.innerHTML = '';
            nodes.forEach(function(n) {
                var row = document.createElement('tr');
                row.className = 'node-row';
                row.setAttribute('data-gene', n.id);
                row.setAttribute('data-role', n.role);

                var weColHtml = '<span class="we-badge we-interactor">INTERACTOR</span>';
                var modColHtml = '—';
                var parseInfo = splitRoleMods(n.weAnnotation, n.modification);
                if (parseInfo) {
                    if (parseInfo.isW && parseInfo.isE) {
                        weColHtml = `
                            <div class="cell-split"><span class="we-badge-sub we-writer">W</span></div>
                            <div class="cell-split" style="margin-top: 4px;"><span class="we-badge-sub we-eraser">E</span></div>
                        `;
                        var wModText = parseInfo.writerMods.length ? parseInfo.writerMods.join(', ') : '—';
                        var eModText = parseInfo.eraserMods.length ? parseInfo.eraserMods.join(', ') : '—';
                        modColHtml = `
                            <div class="cell-split"><span class="mod-sub">${wModText}</span></div>
                            <div class="cell-split" style="margin-top: 4px;"><span class="mod-sub">${eModText}</span></div>
                        `;
                    } else if (parseInfo.isW) {
                        weColHtml = `<div class="cell-split"><span class="we-badge-sub we-writer">W</span></div>`;
                        var wModText = parseInfo.writerMods.length ? parseInfo.writerMods.join(', ') : (parseInfo.eraserMods.join(', ') || '—');
                        modColHtml = `<div class="cell-split"><span class="mod-sub">${wModText}</span></div>`;
                    } else if (parseInfo.isE) {
                        weColHtml = `<div class="cell-split"><span class="we-badge-sub we-eraser">E</span></div>`;
                        var eModText = parseInfo.eraserMods.length ? parseInfo.eraserMods.join(', ') : (parseInfo.writerMods.join(', ') || '—');
                        modColHtml = `<div class="cell-split"><span class="mod-sub">${eModText}</span></div>`;
                    } else {
                        weColHtml = n.weAnnotation ? `<span>${n.weAnnotation}</span>` : `<span class="we-badge we-interactor">INTERACTOR</span>`;
                        modColHtml = `<span>${n.modification && n.modification !== 'NA' ? n.modification : '—'}</span>`;
                    }
                }

                var ptmsFormatted = formatPtmsHtml(n.ptms_raw);

                row.innerHTML = `
                    <td class="gene-name-cell">
                        <div class="cell-content gene-cell-content">
                            <span class="expand-indicator">▸</span>
                            <span class="role-dot role-${n.role}"></span>
                            <span class="gene-label-text">${n.label}</span>
                        </div>
                    </td>
                    <td><div class="cell-content mono-text">${n.accession || '—'}</div></td>
                    <td><div class="cell-content">${n.taxon || '—'}</div></td>
                    <td class="table-we-col"><div class="cell-content">${weColHtml}</div></td>
                    <td class="table-mod-col"><div class="cell-content">${modColHtml}</div></td>
                    <td><div class="cell-content">${(n.enzymeClass && n.enzymeClass !== 'NA' && n.enzymeClass !== '—') ? n.enzymeClass : 'NA'}</div></td>
                    <td><div class="cell-content">${(n.ecNumber && n.ecNumber !== 'NA' && n.ecNumber !== '—') ? n.ecNumber : 'NA'}</div></td>
                    <td class="activity-cell"><div class="cell-content">${n.activity || '—'}</div></td>
                    <td><div class="cell-content">${n.existence || '—'}</div></td>
                    <td><div class="cell-content">${n.regulation || '—'}</div></td>
                    <td class="ptm-cell"><div class="cell-content">${ptmsFormatted}</div></td>
                    <td><div class="cell-content">${n.subunit || '—'}</div></td>
                    <td><div class="cell-content">${n.tissue || '—'}</div></td>
                    <td><div class="cell-content">${n.string_ids || '—'}</div></td>
                `;

                row.querySelectorAll('.cell-content').forEach(function(cc) {
                    attachScrollChaining(cc, wrapperGeneral);
                });

                row.addEventListener('click', function() { handleRowClick(row, n.id); });

                // Clicking a protein's PTMs cell jumps straight to that protein's own
                // Domain Graph (replacing whichever protein's graph was showing), rather
                // than just highlighting the row like a normal row click does.
                var ptmCellEl = row.querySelector('.ptm-cell');
                var hasPtmData = n.ptms_raw && n.ptms_raw !== 'NA' && n.ptms_raw.trim() !== '';
                if (ptmCellEl && hasPtmData && n.accession) {
                    ptmCellEl.classList.add('ptm-cell-clickable');
                    ptmCellEl.title = 'View ' + n.label + '\u2019s domain graph';
                    ptmCellEl.addEventListener('click', function(evt) {
                        evt.stopPropagation();
                        if (window.PTMDomainGraph) window.PTMDomainGraph.showGene(n.id);
                    });
                }

                tbodyGeneral.appendChild(row);
            });
            initTableColumnResizers(tableGeneral);
        }

        // 2. Kinetic parameters table
        var tableKinetic = document.getElementById('kinetic-table');
        var tbodyKinetic = tableKinetic ? tableKinetic.querySelector('tbody') : null;
        var wrapperKinetic = tableKinetic ? tableKinetic.closest('.table-wrapper') : null;

        function formatKineticVal(v) {
            if (!v || v === 'NA' || v === '—' || String(v).trim() === '') return '—';
            return v;
        }

        if (tbodyKinetic) {
            tbodyKinetic.innerHTML = '';
            nodes.forEach(function(n) {
                var row = document.createElement('tr');
                row.className = 'node-row';
                row.setAttribute('data-gene', n.id);
                row.setAttribute('data-role', n.role);

                var sourceVal = n.kinetic_source || 'NA';
                var sourceHtml = '<span class="source-badge source-na">—</span>';
                if (sourceVal === 'BRENDA') {
                    sourceHtml = '<span class="source-badge source-brenda">BRENDA</span>';
                } else if (sourceVal === 'TRANSFERRED') {
                    sourceHtml = '<span class="source-badge source-transferred">TRANSFERRED</span>';
                } else if (sourceVal !== 'NA' && sourceVal !== '—') {
                    sourceHtml = `<span class="source-badge">${sourceVal}</span>`;
                }

                var ecText = (n.kinetic_ec && n.kinetic_ec !== 'NA' && n.kinetic_ec !== '—') 
                    ? n.kinetic_ec 
                    : ((n.ecNumber && n.ecNumber !== 'NA' && n.ecNumber !== '—') ? n.ecNumber : '—');

                row.innerHTML = `
                    <td class="gene-name-cell">
                        <div class="cell-content gene-cell-content">
                            <span class="expand-indicator">▸</span>
                            <span class="role-dot role-${n.role}"></span>
                            <span class="gene-label-text">${n.label}</span>
                        </div>
                    </td>
                    <td><div class="cell-content mono-text">${n.accession || '—'}</div></td>
                    <td><div class="cell-content">${sourceHtml}</div></td>
                    <td><div class="cell-content">${ecText}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_rec_name)}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_sys_name)}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_reaction_type)}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_turnover_no)}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_km)}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_kcat_km)}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_inhibitors)}</div></td>
                    <td><div class="cell-content">${formatKineticVal(n.kinetic_ki)}</div></td>
                `;
                row.querySelectorAll('.cell-content').forEach(function(cc) {
                    attachScrollChaining(cc, wrapperKinetic);
                });
                row.addEventListener('click', function() { handleRowClick(row, n.id); });
                tbodyKinetic.appendChild(row);
            });
            initTableColumnResizers(tableKinetic);
        }

        // 3. Pathway table
        var tablePathway = document.getElementById('pathway-table');
        var tbodyPathway = tablePathway ? tablePathway.querySelector('tbody') : null;
        var wrapperPathway = tablePathway ? tablePathway.closest('.table-wrapper') : null;

        if (tbodyPathway) {
            tbodyPathway.innerHTML = '';
            nodes.forEach(function(n) {
                var row = document.createElement('tr');
                row.className = 'node-row';
                row.setAttribute('data-gene', n.id);
                row.setAttribute('data-role', n.role);
                row.innerHTML = `
                    <td class="gene-name-cell">
                        <div class="cell-content gene-cell-content">
                            <span class="expand-indicator">▸</span>
                            <span class="role-dot role-${n.role}"></span>
                            <span class="gene-label-text">${n.label}</span>
                        </div>
                    </td>
                    <td><div class="cell-content mono-text">${n.accession || '—'}</div></td>
                    <td><div class="cell-content">—</div></td>
                    <td><div class="cell-content">—</div></td>
                    <td><div class="cell-content">—</div></td>
                    <td><div class="cell-content" style="color: var(--text-dim); font-size: 0.76rem;">No pathway annotations available for this network. Data will be populated from curated TSVs in a future update.</div></td>
                `;
                row.querySelectorAll('.cell-content').forEach(function(cc) {
                    attachScrollChaining(cc, wrapperPathway);
                });
                row.addEventListener('click', function() { handleRowClick(row, n.id); });
                tbodyPathway.appendChild(row);
            });
            initTableColumnResizers(tablePathway);
        }

        // 4. Domain table
        var tableDomain = document.getElementById('domain-table');
        var tbodyDomain = tableDomain ? tableDomain.querySelector('tbody') : null;
        var wrapperDomain = tableDomain ? tableDomain.closest('.table-wrapper') : null;

        if (tbodyDomain) {
            tbodyDomain.innerHTML = '';
            nodes.forEach(function(n) {
                var row = document.createElement('tr');
                row.className = 'node-row';
                row.setAttribute('data-gene', n.id);
                row.setAttribute('data-role', n.role);
                row.innerHTML = `
                    <td class="gene-name-cell">
                        <div class="cell-content gene-cell-content">
                            <span class="expand-indicator">▸</span>
                            <span class="role-dot role-${n.role}"></span>
                            <span class="gene-label-text">${n.label}</span>
                        </div>
                    </td>
                    <td><div class="cell-content mono-text">${n.accession || '—'}</div></td>
                    <td><div class="cell-content">—</div></td>
                    <td><div class="cell-content">—</div></td>
                    <td><div class="cell-content">—</div></td>
                    <td><div class="cell-content" style="color: var(--text-dim); font-size: 0.76rem;">No domain annotations available for this network. Data will be populated from curated TSVs in a future update.</div></td>
                `;
                row.querySelectorAll('.cell-content').forEach(function(cc) {
                    attachScrollChaining(cc, wrapperDomain);
                });
                row.addEventListener('click', function() { handleRowClick(row, n.id); });
                tbodyDomain.appendChild(row);
            });
            initTableColumnResizers(tableDomain);
        }
    }

    // Dialogue Tooltip hover bubble positioning
    function showTooltipBubble(node) {
        var tooltip = document.getElementById('network-tooltip');
        if (!tooltip) return;

        var role = node.data('role').toUpperCase();
        var annotation = node.data('weAnnotation') || 'N/A';
        var modification = node.data('modification') || 'N/A';
        
        var bodyHtml = '';
        var parseInfo = splitRoleMods(annotation, modification);
        if (parseInfo) {
            if (parseInfo.isW && parseInfo.isE) {
                var wModText = parseInfo.writerMods.length ? parseInfo.writerMods.join(', ') : '—';
                var eModText = parseInfo.eraserMods.length ? parseInfo.eraserMods.join(', ') : '—';
                bodyHtml = `
                    <div style="margin-bottom: 6px; font-size: 0.78rem;"><strong>Role:</strong> Writer &amp; Eraser</div>
                    <div style="display:flex; align-items:flex-start; margin-bottom: 4px;">
                        <span class="we-badge-sub we-writer" style="margin-top:2px;">W</span> <span class="mod-sub">${wModText}</span>
                    </div>
                    <div style="display:flex; align-items:flex-start;">
                        <span class="we-badge-sub we-eraser" style="margin-top:2px;">E</span> <span class="mod-sub">${eModText}</span>
                    </div>
                `;
            } else if (parseInfo.isW) {
                var wModText = parseInfo.writerMods.length ? parseInfo.writerMods.join(', ') : (parseInfo.eraserMods.join(', ') || '—');
                bodyHtml = `
                    <div style="margin-bottom: 6px; font-size: 0.78rem;"><strong>Role:</strong> Writer</div>
                    <div style="display:flex; align-items:flex-start;">
                        <span class="we-badge-sub we-writer" style="margin-top:2px;">W</span> <span class="mod-sub">${wModText}</span>
                    </div>
                `;
            } else if (parseInfo.isE) {
                var eModText = parseInfo.eraserMods.length ? parseInfo.eraserMods.join(', ') : (parseInfo.writerMods.join(', ') || '—');
                bodyHtml = `
                    <div style="margin-bottom: 6px; font-size: 0.78rem;"><strong>Role:</strong> Eraser</div>
                    <div style="display:flex; align-items:flex-start;">
                        <span class="we-badge-sub we-eraser" style="margin-top:2px;">E</span> <span class="mod-sub">${eModText}</span>
                    </div>
                `;
            }
        } else {
            bodyHtml = `
                <div style="margin-bottom: 4px; font-size: 0.78rem;"><strong>Role:</strong> ${role}</div>
                <div><strong>Modifications:</strong> ${modification !== 'NA' ? modification : '—'}</div>
            `;
        }
        
        tooltip.innerHTML = `
            <div style="font-weight: 700; font-size: 0.9rem; margin-bottom: 6px; color: var(--blue); border-bottom: 1px solid var(--border); padding-bottom: 4px;">${node.data('id')}</div>
            ${bodyHtml}
        `;
        
        tooltip.style.display = 'block';
        
        // Position relative to .network-wrapper
        var pos = node.renderedPosition();
        var nodeWidth = node.renderedWidth();
        
        var x = pos.x + (nodeWidth / 2) + 8;
        var y = pos.y - (tooltip.offsetHeight / 2);
        
        // Prevent going out of boundaries
        var wrapper = document.getElementById('network-wrapper');
        if (wrapper) {
            var wrapperWidth = wrapper.offsetWidth;
            if (x + tooltip.offsetWidth > wrapperWidth) {
                x = pos.x - (nodeWidth / 2) - tooltip.offsetWidth - 8;
                tooltip.classList.add('tooltip-left');
            } else {
                tooltip.classList.remove('tooltip-left');
            }
        }
        
        tooltip.style.left = x + 'px';
        tooltip.style.top = y + 'px';
    }

    function hideTooltipBubble() {
        var tooltip = document.getElementById('network-tooltip');
        if (tooltip) tooltip.style.display = 'none';
    }

    // ──────────────────────────────────────────────────────────────
    // 2. CYTOSCAPE INITIALIZATION
    // ──────────────────────────────────────────────────────────────
    
    function initCytoscape(networkData) {
        var elements = [];

        var COLORS = {
            writer:    '#e11d48', // Crimson Red
            eraser:    '#22c55e', // Bright Green
            both:      '#8b5cf6', // Purple
            interactor:'#cbd5e1'  // Light Grey
        };

        var SHAPES = {
            writer:    'round-triangle',
            eraser:    'ellipse',
            both:      'round-diamond',
            interactor:'round-rectangle'
        };

        // Format nodes
        networkData.nodes.forEach(function(n) {
            // For writers, erasers, and both, remove abbreviation inside node (make SVG background none)
            var bgImage = (n.role === 'writer' || n.role === 'eraser' || n.role === 'both') ? 'none' : makeSvgBg(n.modLabel);
            elements.push({
                group: 'nodes',
                data: {
                    id:           n.id,
                    label:        n.label,
                    role:         n.role,
                    shape:        SHAPES[n.role] || 'round-rectangle',
                    color:        COLORS[n.role] || COLORS.interactor,
                    bgImage:      bgImage,
                    modLabel:     n.modLabel,
                    accession:    n.accession,
                    activity:     n.activity,
                    weAnnotation: n.weAnnotation,
                    modification: n.modification,
                    enzymeClass:  n.enzymeClass,
                    ecNumber:     n.ecNumber,
                    isQuery:      n.isQuery
                }
            });
        });

        // Format edges
        networkData.edges.forEach(function(e) {
            var edgeId = e.source + '__' + e.target;
            elements.push({
                group: 'edges',
                data: {
                    id:     edgeId,
                    source: e.source,
                    target: e.target,
                    score:  e.score
                }
            });
        });

        if (cy) {
            cy.destroy();
        }

        cy = cytoscape({
            container: document.getElementById('cy'),
            elements:  elements,
            style: [
                {
                    selector: 'node',
                    style: {
                        'shape':               'data(shape)',
                        'background-color':    'data(color)',
                        'background-image':    'data(bgImage)',
                        'background-fit':      'none',
                        'background-position-x': '50%',
                        'background-position-y': '50%',
                        'label':               'data(label)',
                        'color':               '#1c1b18',
                        'font-size':           '11px',
                        'font-family':         'Inter, Arial, sans-serif',
                        'font-weight':         '600',
                        'text-valign':         'bottom',
                        'text-halign':         'center',
                        'text-margin-y':       6,
                        'text-outline-color':  '#faf9f6',
                        'text-outline-width':  2,
                        'width':               44,
                        'height':              44,
                        'border-width':        2,
                        'border-color':        'rgba(0,0,0,0.2)'
                    }
                },
                {
                    // Interactors: Squares size 42
                    selector: 'node[role = "interactor"]',
                    style: {
                        'width':               42,
                        'height':              42
                    }
                },
                {
                    selector: 'node[?isQuery]',
                    style: {
                        'border-width':        3.5,
                        'border-color':        '#3b82f6', // Glowing blue border for query genes
                        'border-opacity':      0.95
                    }
                },
                {
                    selector: 'node:selected',
                    style: {
                        'border-width':        4,
                        'border-color':        '#2563eb', // Solid blue border on select
                        'border-opacity':      1.0
                    }
                },
                {
                    selector: 'edge',
                    style: {
                        'width':              function(e){ return Math.max(1.0, e.data('score') * 4.5); },
                        'line-color':         '#cbd5e1',
                        'opacity':            0.6,
                        'curve-style':        'bezier'
                    }
                }
            ],
            layout: {
                name:           'cose',
                padding:        60,
                nodeRepulsion:  600000,
                idealEdgeLength: 110,
                edgeElasticity:  100,
                gravity:         140,
                numIter:         1000,
                animate:         true,
                animationDuration: 500
            }
        });

        // ──────────────────────────────────────────────────────────────
        // 3. EVENT BINDING (CYTOSCAPE)
        // ──────────────────────────────────────────────────────────────
        
        // Node select
        cy.on('tap', 'node', function(evt) {
            var node = evt.target;
            var geneId = node.id();
            highlightTableRow(geneId);
            
            var activeTab = document.querySelector('.table-tab-content.active');
            if (activeTab) {
                var targetRow = activeTab.querySelector(`.node-row[data-gene="${geneId}"]`);
                if (targetRow) {
                    handleRowClick(targetRow, geneId);
                }
            }
        });

        // Hover tooltip bubble handlers
        cy.on('mouseover', 'node', function(evt) {
            var node = evt.target;
            var role = node.data('role');
            if (role === 'writer' || role === 'eraser' || role === 'both') {
                showTooltipBubble(node);
            }
        });

        cy.on('mouseout', 'node', function(evt) {
            hideTooltipBubble();
        });

        // Canvas deselect
        cy.on('tap', function(evt) {
            if (evt.target === cy) {
                var detailSection = document.getElementById('node-detail-section');
                if (detailSection) detailSection.style.display = 'none';
                document.querySelectorAll('.node-row').forEach(function(r) {
                    r.classList.remove('row-highlight');
                });
            }
        });
    }

    // Initialize cytoscape and table on startup
    if (typeof NETWORK_DATA !== 'undefined' && NETWORK_DATA) {
        initCytoscape(NETWORK_DATA);
        rebuildTable(NETWORK_DATA.nodes);
    }

    // ──────────────────────────────────────────────────────────────
    // 4. VIEWPORT CONTROLS
    // ──────────────────────────────────────────────────────────────
    
    var fitBtn = document.getElementById('fit-btn');
    if (fitBtn) {
        fitBtn.addEventListener('click', function() {
            if (cy) cy.fit(50);
        });
    }

    var pngBtn = document.getElementById('png-btn');
    if (pngBtn) {
        pngBtn.addEventListener('click', function() {
            if (!cy) return;
            var pngContent = cy.png({ full: true, scale: 2 });
            var link = document.createElement('a');
            link.href = pngContent;
            link.download = 'ptm_network.png';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        });
    }

    // ──────────────────────────────────────────────────────────────
    // 5. LIVE AJAX INTERACTOR LIMIT FILTER (DEFAULT 10 INTERACTORS)
    // ──────────────────────────────────────────────────────────────
    
    var slider       = document.getElementById('live-limit');
    var badge        = document.getElementById('live-limit-badge');
    var applyBtn     = document.getElementById('apply-limit-btn');
    var edgeBadge    = document.querySelector('.topbar-stats .stat-chip:nth-child(2)');

    if (slider && badge) {
        var initialVal = parseInt(slider.value);
        if (initialVal === MAX_LIMIT && MAX_LIMIT > 10) {
            badge.textContent = 'All (' + initialVal + ')';
        } else {
            badge.textContent = initialVal;
        }

        slider.addEventListener('input', function() {
            var val = parseInt(slider.value);
            if (val === MAX_LIMIT && MAX_LIMIT > 10) {
                badge.textContent = 'All (' + val + ')';
            } else {
                badge.textContent = val;
            }
        });
    }

    if (applyBtn) {
        applyBtn.addEventListener('click', function() {
            var val = parseInt(slider.value);
            var limitVal = (val === MAX_LIMIT && MAX_LIMIT > 10) ? 'all' : val;
            applyBtn.textContent = 'Filtering...';
            applyBtn.disabled = true;

            fetch('/api/network', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    genes: QUERY_GENES,
                    limit: limitVal,
                    species: SPECIES_VAL
                })
            })
            .then(function(res) {
                if (!res.ok) throw new Error('API request failed');
                return res.json();
            })
            .then(function(newNetwork) {
                // Re-draw Cytoscape network
                initCytoscape(newNetwork);
                
                // Re-draw bottom table rows
                rebuildTable(newNetwork.nodes);

                // Update Edge count badge in top bar
                if (edgeBadge) {
                    edgeBadge.textContent = `${newNetwork.edges.length} edges`;
                }

                var detailSection = document.getElementById('node-detail-section');
                if (detailSection) detailSection.style.display = 'none';

                applyBtn.textContent = 'Apply Filter';
                applyBtn.disabled = false;
            })
            .catch(function(err) {
                console.error(err);
                alert('Failed to apply filter: ' + err.message);
                applyBtn.textContent = 'Apply Filter';
                applyBtn.disabled = false;
            });
        });
    }

    // ──────────────────────────────────────────────────────────────
    // 6. PANEL RESIZER & COLLAPSIBLE SIDEBAR
    // ──────────────────────────────────────────────────────────────
    
    var resizer = document.getElementById('h-resizer');
    var wrapper = document.getElementById('network-wrapper');
    var tableSec = document.getElementById('table-section');
    var layout = document.getElementById('results-layout');

    if (resizer && wrapper && tableSec && layout) {
        var isDragging = false;
        resizer.addEventListener('mousedown', function(e) {
            isDragging = true;
            resizer.classList.add('dragging');
            document.body.style.cursor = 'row-resize';
            document.body.style.userSelect = 'none';
        });

        document.addEventListener('mousemove', function(e) {
            if (!isDragging) return;
            var layoutRect = layout.getBoundingClientRect();
            var relativeY = e.clientY - layoutRect.top;
            
            if (relativeY > 150 && relativeY < layoutRect.height - 150) {
                wrapper.style.height = relativeY + 'px';
                
                var newTableHeight = layoutRect.height - relativeY - 6;
                document.documentElement.style.setProperty('--table-height', newTableHeight + 'px');
                var newExpandedHeight = Math.max(120, newTableHeight - 100);
                document.documentElement.style.setProperty('--expanded-row-height', newExpandedHeight + 'px');
                
                if (cy) {
                    cy.resize();
                }
            }
        });

        document.addEventListener('mouseup', function() {
            if (isDragging) {
                isDragging = false;
                resizer.classList.remove('dragging');
                document.body.style.cursor = '';
                document.body.style.userSelect = '';
            }
        });
    }

    // Sidebar toggle (Show/Hide menu)
    var toggleSidebarBtn = document.getElementById('toggle-sidebar-btn');
    var toggleSidebarText = document.getElementById('toggle-sidebar-text');
    var layoutContainer = document.getElementById('results-layout');

    if (toggleSidebarBtn && layoutContainer) {
        toggleSidebarBtn.addEventListener('click', function() {
            var collapsed = layoutContainer.classList.toggle('sidebar-collapsed');
            if (collapsed) {
                if (toggleSidebarText) toggleSidebarText.textContent = 'Show Menu';
                toggleSidebarBtn.innerHTML = `
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" style="margin-right: 4px; vertical-align: middle;"><path d="M15 3h4a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2h-4M9 17l-5-5 5-5M4 12h12"/></svg>
                    <span>Show Menu</span>
                `;
            } else {
                if (toggleSidebarText) toggleSidebarText.textContent = 'Hide Menu';
                toggleSidebarBtn.innerHTML = `
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" style="margin-right: 4px; vertical-align: middle;"><path d="M9 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h4M16 17l5-5-5-5M21 12H9"/></svg>
                    <span>Hide Menu</span>
                `;
            }
            setTimeout(function() {
                if (cy) cy.resize();
            }, 250);
        });
    }

    // ──────────────────────────────────────────────────────────────
    // 7. TABLE MENU TAB SWITCHER
    // ──────────────────────────────────────────────────────────────
    
    var menuButtons = document.querySelectorAll('.table-menu-btn');
    var tabContents = document.querySelectorAll('.table-tab-content');

    menuButtons.forEach(function(btn) {
        btn.addEventListener('click', function() {
            var targetTabId = btn.getAttribute('data-target');
            
            menuButtons.forEach(function(b) { b.classList.remove('active'); });
            tabContents.forEach(function(t) { t.classList.remove('active'); });
            
            btn.classList.add('active');
            var targetTab = document.getElementById(targetTabId);
            if (targetTab) {
                targetTab.classList.add('active');
                var activeTable = targetTab.querySelector('table');
                if (activeTable) {
                    initTableColumnResizers(activeTable);
                }
            }
        });
    });

    // ──────────────────────────────────────────────────────────────
    // 8. DOMAIN GRAPH VIEW — D3.js linear track
    //    PTM sites: from the site's OWN curated data (node.ptms_raw, already embedded
    //               in NETWORK_DATA by app.py from UniProt_25Jun2026.tsv). No network call,
    //               instant, and consistent with the rest of the site (tables, etc.).
    //    Domain boundaries: not present in the curated TSVs, so fetched from the backend
    //               /api/domains/<accession> (UniProt "Domain" feature + InterPro/Pfam),
    //               with a direct-UniProt-from-browser fallback if that route isn't reachable.
    // ──────────────────────────────────────────────────────────────

    (function () {
        var container  = document.getElementById('domain-cy');
        var statusEl    = document.getElementById('domain-status');
        var selectWrap  = document.getElementById('domain-protein-select-wrap');
        var selectEl    = document.getElementById('domain-protein-select');
        var tooltip     = document.getElementById('domain-tooltip');

        if (!container || typeof d3 === 'undefined' || !NETWORK_DATA || !NETWORK_DATA.nodes) return;

        // node.accession can be a semicolon-joined list of merged UniProt accessions
        // (e.g. "Q92830;Q8N1A2;Q9UCW1" — a primary plus historical/secondary IDs from
        // UniProt entry merges). Only the FIRST one is a fetchable accession; the rest
        // are bookkeeping and aren't needed for the graph.
        function primaryAccession(raw) {
            return (raw || '').split(';')[0].trim();
        }

        var allNodes = NETWORK_DATA.nodes;
        // The dropdown (for multi-gene queries) only lists the queried gene(s)...
        var queryNodes = allNodes.filter(function (n) {
            return QUERY_GENES.indexOf(n.id) >= 0 && n.accession;
        });
        // ...but ANY node with an accession (including interactors) can be shown here,
        // e.g. when clicked from the PTMs column of the General Information table.
        var hasAnyAccession = allNodes.some(function (n) { return !!n.accession; });

        var loadedFor = null; // accession currently drawn, so we don't refetch on every toggle
        var zoomBehavior = null;
        var svg;

        if (!hasAnyAccession) {
            if (statusEl) statusEl.textContent = 'No accession available for any protein in this network — cannot build a domain graph.';
            return;
        }

        if (queryNodes.length > 1) {
            selectWrap.style.display = '';
            queryNodes.forEach(function (n) {
                var opt = document.createElement('option');
                opt.value = n.id;
                opt.textContent = n.label + ' (' + primaryAccession(n.accession) + ')';
                selectEl.appendChild(opt);
            });
            selectEl.addEventListener('change', function () { loadForGene(selectEl.value); });
        }

        var purple = getComputedStyle(document.documentElement).getPropertyValue('--purple').trim() || '#7c3aed';
        function roleColor(role) {
            if (role === 'writer') return '#e11d48';
            if (role === 'eraser') return '#22c55e';
            if (role === 'both')   return purple;
            return '#94a3b8';
        }

        // PTM chemical sub-type → shape, independent of writer/eraser role (which is color).
        // So a marker's COLOR tells you the protein's role, its SHAPE tells you the modification type.
        var PTM_CATEGORIES = {
            phospho:    { label: 'Phosphorylation',    symbol: d3.symbolTriangle, color: '#3b82f6' }, // Blue
            acetyl:     { label: 'Acetylation',         symbol: d3.symbolCircle, color: '#ef4444' }, // Red
            acyl:       { label: 'Acylation (succinyl/glutaryl/malonyl…)', symbol: d3.symbolSquare, color: '#f97316' }, // Orange
            methyl:     { label: 'Methylation',         symbol: d3.symbolSquare, color: '#22c55e' }, // Green
            ubiquitin:  { label: 'Ubiquitin / SUMO',    symbol: d3.symbolDiamond, color: '#8b5cf6' }, // Purple
            glyco:      { label: 'Glycosylation',       symbol: d3.symbolWye, color: '#ec4899' }, // Pink
            lipid:      { label: 'Lipidation',          symbol: d3.symbolCross, color: '#14b8a6' }, // Teal
            hydroxyl:   { label: 'Hydroxylation',       symbol: d3.symbolStar, color: '#06b6d4' }, // Cyan
            nitration:  { label: 'Nitration',           symbol: d3.symbolWye, color: '#6366f1' }, // Indigo
            adpRibosyl: { label: 'ADP-ribosylation',    symbol: d3.symbolDiamond, color: '#84cc16' }, // Lime
            citrullin:  { label: 'Citrullination',      symbol: d3.symbolCross, color: '#10b981' }, // Emerald
            other:      { label: 'Other modification',  symbol: d3.symbolStar, color: '#94a3b8' } // Grey
        };
        // Classifies by the short modification name as it appears in the curated PTMs column
        // (e.g. "Acetylation", "Phosphoprotein", "Glutarylation", "Nitration"...).
        function classifyPtmCategory(name, featureType) {
            var d = (name || '').toLowerCase();
            if (featureType === 'Glycosylation') return 'glyco';
            if (featureType === 'Lipidation') return 'lipid';
            if (d.indexOf('phospho') >= 0) return 'phospho';
            if (d.indexOf('acetyl') >= 0) return 'acetyl';
            if (/(glutaryl|succinyl|malonyl|crotonyl|propionyl|butyryl)/.test(d)) return 'acyl';
            if (d.indexOf('methyl') >= 0) return 'methyl';
            if (d.indexOf('ubiquitin') >= 0 || d.indexOf('sumo') >= 0 || d.indexOf('isopeptide') >= 0) return 'ubiquitin';
            if (d.indexOf('glyco') >= 0) return 'glyco';
            if (/(palmitoyl|myristoyl|prenyl|farnesyl|geranylgeranyl)/.test(d)) return 'lipid';
            if (d.indexOf('hydroxy') >= 0) return 'hydroxyl';
            if (d.indexOf('nitrat') >= 0) return 'nitration';
            if (d.indexOf('adp-ribosyl') >= 0 || d.indexOf('adp ribosyl') >= 0) return 'adpRibosyl';
            if (d.indexOf('citrullin') >= 0) return 'citrullin';
            return 'other';
        }

        // Parses PTM features straight from a live UniProtKB JSON entry (used by the direct-
        // fetch fallback, when the backend route can't be reached). Mirrors what the backend's
        // fetch_uniprot_domain_data() does server-side.
        function parseLiveUniProtPtms(entry, role) {
            var ptms = [];
            (entry.features || []).forEach(function (f) {
                if (['Modified residue', 'Cross-link', 'Lipidation', 'Glycosylation'].indexOf(f.type) < 0) return;
                var start = f.location && f.location.start && f.location.start.value;
                if (!start) return;
                var desc = f.description || f.type;
                ptms.push({ position: start, type: desc, role: role, category: classifyPtmCategory(desc, f.type) });
            });
            return ptms;
        }

        // Last-resort fallback only: the site's own curated "PTMs" column (node.ptms_raw),
        // used only when live UniProt is unreachable through both the backend and direct fetch.
        // Format: "Acetylation|K|549 | Phosphoprotein|T|735 | ..."
        function parsePtmsRaw(raw, role) {
            if (!raw || raw === 'NA') return [];
            return raw.split(' | ').map(function (part) {
                var bits = part.split('|');
                var name = (bits[0] || '').trim();
                var pos = parseInt((bits[2] || '').trim(), 10);
                if (!name || isNaN(pos)) return null;
                return { position: pos, type: name, role: role, category: classifyPtmCategory(name) };
            }).filter(Boolean);
        }

        function showTip(evt, html) {
            var rect = container.getBoundingClientRect();
            tooltip.style.display = 'block';
            tooltip.style.left = (evt.clientX - rect.left + 14) + 'px';
            tooltip.style.top = (evt.clientY - rect.top - 10) + 'px';
            tooltip.innerHTML = html;
        }
        function hideTip() { tooltip.style.display = 'none'; }

        // Click-to-expand residue sequence panel (matches UniProt's own Family & Domains UX)
        var seqPanel  = document.getElementById('domain-sequence-panel');
        var seqTitle  = document.getElementById('domain-sequence-title');
        var seqBody   = document.getElementById('domain-sequence-body');
        var seqCloseBtn = document.getElementById('domain-sequence-close');
        var openDomainKey = null; // "start-end" of the currently expanded domain, so a second click collapses it

        function formatSequence(seq) {
            // Wrap into 60-char lines, UniProt-style, in a monospace block.
            var lines = [];
            for (var i = 0; i < seq.length; i += 60) lines.push(seq.slice(i, i + 60));
            return lines.join('\n');
        }

        function closeSeqPanel() {
            if (seqPanel) seqPanel.style.display = 'none';
            openDomainKey = null;
        }

        function toggleDomainSequence(d, sequence) {
            if (!seqPanel) return;
            var key = d.start + '-' + d.end + ':' + d.name;
            if (openDomainKey === key) { closeSeqPanel(); return; }
            openDomainKey = key;

            seqTitle.textContent = d.name + '  ·  residues ' + d.start + '–' + d.end + ' (' + (d.end - d.start + 1) + ' aa)';
            if (sequence) {
                var slice = sequence.substring(d.start - 1, d.end);
                seqBody.textContent = formatSequence(slice);
                seqBody.classList.remove('domain-sequence-empty');
            } else {
                seqBody.textContent = 'Sequence unavailable for this domain (backend/UniProt request for the full sequence failed).';
                seqBody.classList.add('domain-sequence-empty');
            }
            seqPanel.style.display = '';
        }

        if (seqCloseBtn) seqCloseBtn.addEventListener('click', closeSeqPanel);

        // Draws one UniProt-style track: a slim light-gray backbone bar spanning the full
        // protein, with colored domain segments layered on top. Segment labels are drawn
        // INSIDE the segment (white text) only when there's room — otherwise rely on the
        // hover tooltip / click-to-expand-sequence for the name. No separate label row,
        // which is what kept pushing the old layout's height past its box.
        function drawDomainTrack(gRoot, xScaleFn, domains, trackY, fillColor, sequence, length, sourceLabel) {
            var g = gRoot.append('g').attr('class', 'dg-track');

            g.append('rect').attr('class', 'dg-track-backbone')
                .attr('x', xScaleFn(0)).attr('y', trackY - 6)
                .attr('width', Math.max(xScaleFn(length) - xScaleFn(0), 1))
                .attr('height', 12).attr('rx', 6);

            var seg = g.selectAll('.dg-domain-rect').data(domains).enter().append('rect')
                .attr('class', 'dg-domain-rect')
                .attr('x', function (d) { return xScaleFn(d.start); })
                .attr('y', trackY - 6)
                .attr('width', function (d) { return Math.max(xScaleFn(d.end) - xScaleFn(d.start), 3); })
                .attr('height', 12).attr('rx', 3)
                .attr('fill', fillColor)
                .on('mousemove', function (evt, d) { showTip(evt, '<strong>' + d.name + '</strong> <span class="tip-dim">(' + sourceLabel + ')</span><br>Residues ' + d.start + '–' + d.end + '<br><em>Click to view sequence</em>'); })
                .on('mouseleave', hideTip)
                .on('click', function (evt, d) { toggleDomainSequence(d, sequence); });

            var labels = g.selectAll('.dg-domain-inline-label').data(domains).enter().append('text')
                .attr('class', 'dg-domain-inline-label')
                .attr('x', function (d) { return (xScaleFn(d.start) + xScaleFn(d.end)) / 2; })
                .attr('y', trackY + 4)
                .attr('text-anchor', 'middle')
                .style('pointer-events', 'none')
                .text(function (d) {
                    var w = xScaleFn(d.end) - xScaleFn(d.start);
                    return (d.name.length * 5.6) < (w - 8) ? d.name : '';
                });

            return { seg: seg, labels: labels };
        }

        function renderGraph(node, data) {
            // data = { length, sequence, uniprotDomains, pfamDomains, ptms, source }
            container.innerHTML = '';
            closeSeqPanel(); // switching protein/redrawing — any previously-open sequence is stale
            var hasPfam = data.pfamDomains && data.pfamDomains.length > 0;
            var margin = { top: 22, right: 20, bottom: 14, left: 24 };
            var width  = Math.max(container.clientWidth, 320);
            var ptmZoneH  = 36;
            var trackGap  = 22;
            var height = margin.top + ptmZoneH + trackGap * (hasPfam ? 1 : 0.4) + margin.bottom + 6;
            var length = data.length || (function () {
                var maxPos = 0;
                data.ptms.forEach(function (p) { if (p.position > maxPos) maxPos = p.position; });
                data.uniprotDomains.concat(data.pfamDomains).forEach(function (d) { if (d.end > maxPos) maxPos = d.end; });
                return maxPos ? maxPos + 50 : 500;
            })();

            var uniY  = margin.top + ptmZoneH;
            var pfamY = uniY + trackGap;

            svg = d3.select(container).append('svg')
                .attr('width', width).attr('height', height)
                .style('display', 'block');

            var clipId = 'dg-clip-' + node.id;
            svg.append('clipPath').attr('id', clipId).append('rect')
                .attr('x', margin.left).attr('y', 0)
                .attr('width', width - margin.left - margin.right).attr('height', height);

            var gRoot = svg.append('g').attr('clip-path', 'url(#' + clipId + ')');

            var baseX = d3.scaleLinear().domain([0, length]).range([margin.left, width - margin.right]);

            // Ruler on top, UniProt-style
            var gAxis = svg.append('g')
                .attr('class', 'dg-axis')
                .attr('transform', 'translate(0,' + margin.top + ')')
                .call(d3.axisTop(baseX).ticks(Math.max(Math.floor(width / 90), 4)).tickFormat(d3.format('d')));

            var gUni  = drawDomainTrack(gRoot, baseX, data.uniprotDomains, uniY, '#3b6fd4', data.sequence, length, 'UniProt');

            var gPfam = null;
            if (hasPfam) {
                gPfam = drawDomainTrack(gRoot, baseX, data.pfamDomains, pfamY, '#2f9e7a', data.sequence, length, 'Pfam');
            }

            if (data.uniprotDomains.length === 0 && !hasPfam) {
                gRoot.append('text').attr('class', 'dg-no-domains')
                    .attr('x', (baseX(0) + baseX(length)) / 2).attr('y', uniY + 4)
                    .attr('text-anchor', 'middle').text('No domain boundaries available for this accession.');
            }

            // PTM lollipops — thin stem from the UniProt backbone up to a shape+color marker.
            // SHAPE = modification type, COLOR = protein's writer/eraser/both role.
            var gPtms = gRoot.append('g').attr('class', 'dg-ptms');
            gPtms.selectAll('.dg-ptm-stem').data(data.ptms).enter().append('line')
                .attr('class', 'dg-ptm-stem')
                .attr('x1', function (d) { return baseX(d.position); })
                .attr('x2', function (d) { return baseX(d.position); })
                .attr('y1', uniY - 7).attr('y2', margin.top + 10);
            gPtms.selectAll('.dg-ptm-mark').data(data.ptms).enter().append('path')
                .attr('class', 'dg-ptm-mark')
                .each(function (d) {
                    var cat = PTM_CATEGORIES[d.category] || PTM_CATEGORIES.other;
                    d3.select(this).attr('d', d3.symbol().type(cat.symbol).size(60));
                })
                .attr('fill', function (d) { 
                    var cat = PTM_CATEGORIES[d.category] || PTM_CATEGORIES.other;
                    return cat.color; 
                })
                .attr('transform', function (d) { return 'translate(' + baseX(d.position) + ',' + (margin.top + 8) + ')'; })
                .on('mousemove', function (evt, d) {
                    var cat = PTM_CATEGORIES[d.category] || PTM_CATEGORIES.other;
                    showTip(evt, '<strong>' + d.type + '</strong><br>' + cat.label + ' · Position ' + d.position);
                })
                .on('mouseleave', hideTip);

            // Dynamic PTM-type legend — only shapes actually present on this protein
            var typeLegendEl = document.getElementById('domain-type-legend');
            if (typeLegendEl) {
                var seen = {};
                data.ptms.forEach(function (p) { seen[p.category] = true; });
                var cats = Object.keys(seen);
                typeLegendEl.innerHTML = '';
                if (cats.length) {
                    var headingEl = document.createElement('span');
                    headingEl.className = 'dg-type-legend-heading';
                    headingEl.textContent = 'Modification type:';
                    typeLegendEl.appendChild(headingEl);
                    cats.forEach(function (key) {
                        var cat = PTM_CATEGORIES[key] || PTM_CATEGORIES.other;
                        var item = document.createElement('span');
                        item.className = 'dg-type-legend-item';
                        var iconSvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
                        iconSvg.setAttribute('width', 14); iconSvg.setAttribute('height', 14);
                        var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                        path.setAttribute('d', d3.symbol().type(cat.symbol).size(60)());
                        path.setAttribute('transform', 'translate(7,7)');
                        path.setAttribute('fill', cat.color);
                        iconSvg.appendChild(path);
                        item.appendChild(iconSvg);
                        item.appendChild(document.createTextNode(' ' + cat.label));
                        typeLegendEl.appendChild(item);
                    });
                }
            }

            // Horizontal pan/zoom
            zoomBehavior = d3.zoom()
                .scaleExtent([1, 25])
                .translateExtent([[margin.left, 0], [width - margin.right, height]])
                .extent([[margin.left, 0], [width - margin.right, height]])
                .on('zoom', function (evt) {
                    var zx = evt.transform.rescaleX(baseX);
                    gAxis.call(d3.axisTop(zx).ticks(Math.max(Math.floor(width / 90), 4)).tickFormat(d3.format('d')));
                    [gUni, gPfam].forEach(function (trk) {
                        if (!trk) return;
                        trk.seg
                            .attr('x', function (d) { return zx(d.start); })
                            .attr('width', function (d) { return Math.max(zx(d.end) - zx(d.start), 3); });
                        trk.labels.attr('x', function (d) { return (zx(d.start) + zx(d.end)) / 2; });
                    });
                    gRoot.selectAll('.dg-track-backbone').attr('x', zx(0)).attr('width', Math.max(zx(length) - zx(0), 1));
                    gPtms.selectAll('.dg-ptm-stem')
                        .attr('x1', function (d) { return zx(d.position); })
                        .attr('x2', function (d) { return zx(d.position); });
                    gPtms.selectAll('.dg-ptm-mark')
                        .attr('transform', function (d) { return 'translate(' + zx(d.position) + ',' + (margin.top + 8) + ')'; });
                });
            svg.call(zoomBehavior);

            var degradedNote = data.source === 'ptm-only-curated'
                ? ' — showing offline curated data (UniProt unreachable right now)'
                : '';
            console.log('[Domain Graph] ' + node.id + ': source=' + data.source + ', accession=' + primaryAccession(node.accession));
            statusEl.textContent = node.label + ' (' + primaryAccession(node.accession) + ')' +
                (data.length ? ' — ' + data.length + ' aa' : '') + ' — ' +
                (data.uniprotDomains.length + (data.pfamDomains ? data.pfamDomains.length : 0)) + ' domain regions, ' +
                data.ptms.length + ' PTM sites' + degradedNote + '. Scroll to zoom, drag to pan.';
        }

        function fallbackDirectUniProt(node, primaryAccession) {
            statusEl.textContent = 'Fetching domain boundaries & PTM sites for ' + node.label + ' (' + primaryAccession + ') from UniProt…';
            fetch('https://rest.uniprot.org/uniprotkb/' + encodeURIComponent(primaryAccession) + '.json?fields=sequence,ft_domain,ft_mod_res,ft_carbohyd,ft_lipid,ft_crosslnk')
                .then(function (res) { if (!res.ok) throw new Error('UniProt request failed (' + res.status + ')'); return res.json(); })
                .then(function (entry) {
                    var length = entry.sequence ? entry.sequence.length : null;
                    var sequence = entry.sequence ? entry.sequence.value : null;
                    var uniprotDomains = [];
                    (entry.features || []).forEach(function (f) {
                        if (f.type !== 'Domain') return;
                        var start = f.location && f.location.start && f.location.start.value;
                        var end = f.location && f.location.end && f.location.end.value;
                        if (start) uniprotDomains.push({ start: start, end: end || start, name: f.description || 'Domain' });
                    });
                    var ptms = parseLiveUniProtPtms(entry, node.role);
                    renderGraph(node, { length: length, sequence: sequence, uniprotDomains: uniprotDomains, pfamDomains: [], ptms: ptms, source: 'uniprot-direct' });
                })
                .catch(function (err) {
                    console.error(err);
                    // Total network failure — last resort: the site's own curated PTM data, so the
                    // graph isn't empty. Clearly labelled as offline/curated in the status line.
                    var curatedPtms = parsePtmsRaw(node.ptms_raw, node.role);
                    renderGraph(node, { length: null, sequence: null, uniprotDomains: [], pfamDomains: [], ptms: curatedPtms, source: 'ptm-only-curated' });
                });
        }

        function loadForGene(geneId) {
            var node = allNodes.filter(function (n) { return n.id === geneId; })[0]
                     || queryNodes.filter(function (n) { return n.id === geneId; })[0]
                     || queryNodes[0];
            if (loadedFor === node.accession) return;
            loadedFor = node.accession;

            // node.accession can be a semicolon-list of merged/secondary accessions
            // (e.g. "Q92830;Q8N1A2;Q9UCW1") — shown in full to the user for reference,
            // but only the PRIMARY one is a valid id for UniProt/InterPro API calls.
            var primaryAccession = (node.accession || '').split(';')[0].trim();

            statusEl.textContent = 'Fetching domain boundaries & PTM sites for ' + node.label + ' (' + primaryAccession + ') from UniProt…';

            fetch('/api/domains/' + encodeURIComponent(primaryAccession), { cache: 'no-store' })
                .then(function (res) {
                    if (!res.ok) throw new Error('Backend request failed (' + res.status + ')');
                    return res.json();
                })
                .then(function (payload) {
                    renderGraph(node, {
                        length: payload.length,
                        sequence: payload.sequence || null,
                        uniprotDomains: (payload.domain_tracks && payload.domain_tracks.uniprot) || [],
                        pfamDomains: (payload.domain_tracks && payload.domain_tracks.pfam) || [],
                        ptms: payload.ptms || [], // backend already fetched live from UniProt + assigned our curated role
                        source: 'backend'
                    });
                })
                .catch(function (err) {
                    console.warn('Backend /api/domains unavailable, falling back:', err);
                    fallbackDirectUniProt(node, primaryAccession);
                });
        }

        window.PTMDomainGraph = {
            activate: function () {
                var geneId = selectEl.value || (queryNodes[0] && queryNodes[0].id) || (allNodes.filter(function (n) { return n.accession; })[0] || {}).id;
                if (geneId) loadForGene(geneId);
                if (svg) svg.attr('width', Math.max(container.clientWidth, 320));
            },
            // Called when the user clicks a protein's PTMs cell in the General Information
            // table: switches the whole page to the Domain Graph view and swaps in THAT
            // protein's graph, replacing whatever was showing before.
            showGene: function (geneId) {
                if (window.switchResultsView) window.switchResultsView('domain');
                loadForGene(geneId);
                if (svg) svg.attr('width', Math.max(container.clientWidth, 320));
            }
        };
    })();

});