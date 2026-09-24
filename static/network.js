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

});
