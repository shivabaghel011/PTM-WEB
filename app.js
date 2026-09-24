document.addEventListener('DOMContentLoaded', function () {

    // ──────────────────────────────────────────────────────────────
    // 1. STATE
    // ──────────────────────────────────────────────────────────────
    var masterData  = {};   // { GENE_NAME: { accession, activity, weAnnotation, modification, enzymeClass } }
    var stringEdges = [];   // [ { source, target, score } ]
    var queryGenes  = [];   // [ 'TP53', 'KAT2A', ... ]
    var cy          = null;

    // ──────────────────────────────────────────────────────────────
    // 2. DOM REFS  (all inside DOMContentLoaded — always available)
    // ──────────────────────────────────────────────────────────────
    var elMasterFile   = document.getElementById('master-file');
    var elStringFile   = document.getElementById('string-file');
    var elQueryFile    = document.getElementById('query-file');
    var elScoreSlider  = document.getElementById('score-threshold');
    var elScoreDisplay = document.getElementById('score-display');
    var elGenerateBtn  = document.getElementById('generate-btn');
    var elNodeDetails  = document.getElementById('node-details');
    var elDetailsBody  = document.getElementById('details-content');

    // ──────────────────────────────────────────────────────────────
    // 3. COLORS / LOOKUP
    // ──────────────────────────────────────────────────────────────
    var COLORS = {
        writer:    '#ec4899',
        eraser:    '#10b981',
        both:      '#8b5cf6',
        interactor:'#64748b'
    };

    // ──────────────────────────────────────────────────────────────
    // 4. BASIC EVENT WIRING
    // ──────────────────────────────────────────────────────────────
    elScoreSlider.addEventListener('input', function () {
        elScoreDisplay.textContent = parseFloat(elScoreSlider.value).toFixed(3);
    });

    [elMasterFile, elStringFile, elQueryFile].forEach(function (input) {
        input.addEventListener('change', checkReady);
    });

    function checkReady() {
        var ready = elMasterFile.files.length > 0 &&
                    elStringFile.files.length > 0 &&
                    elQueryFile.files.length  > 0;
        elGenerateBtn.disabled = !ready;
    }

    elGenerateBtn.addEventListener('click', function () {
        elGenerateBtn.disabled  = true;
        elGenerateBtn.textContent = 'Processing…';

        // Reset state
        masterData  = {};
        stringEdges = [];
        queryGenes  = [];

        // Chain: master → string → query → build
        parseMasterFile(elMasterFile.files[0], function (err) {
            if (err) return showError('Master TSV', err);
            parseStringFile(elStringFile.files[0], function (err) {
                if (err) return showError('STRING TSV', err);
                parseQueryFile(elQueryFile.files[0], function (err) {
                    if (err) return showError('Query file', err);
                    try {
                        buildNetwork();
                    } catch (e) {
                        showError('buildNetwork', e);
                    }
                    elGenerateBtn.disabled   = false;
                    elGenerateBtn.textContent = 'Generate Network';
                });
            });
        });
    });

    function showError(stage, err) {
        elGenerateBtn.disabled   = false;
        elGenerateBtn.textContent = 'Generate Network';
        alert('Error at [' + stage + ']:\n' + (err.message || err) + '\n\nLine: ' + (err.stack || ''));
        console.error(stage, err);
    }

    // ──────────────────────────────────────────────────────────────
    // 5. PARSERS
    // ──────────────────────────────────────────────────────────────

    // 5a. Master 6-column TSV  (tab-separated, header row)
    function parseMasterFile(file, cb) {
        Papa.parse(file, {
            header: true,
            delimiter: '\t',
            skipEmptyLines: true,
            complete: function (results) {
                if (!results.data || results.data.length === 0) {
                    return cb(new Error('Master TSV appears empty or header is missing.'));
                }

                var keys    = Object.keys(results.data[0]);
                // Strict positional fallback: col0=Gene, col1=Acc, col2=Activity, col3=W/E, col4=Mod, col5=Class
                var kGene   = keys.find(function(k){ return k.toLowerCase().includes('gene'); })        || keys[0];
                var kAcc    = keys.find(function(k){ return k.toLowerCase().includes('accession'); })   || keys[1];
                var kAct    = keys.find(function(k){ return k.toLowerCase().includes('catalytic'); })   || keys[2];
                var kWE     = keys.find(function(k){ return k.toLowerCase().includes('w/e') || k.toLowerCase().includes('annotation'); }) || keys[3];
                var kMod    = keys.find(function(k){ return k.toLowerCase().includes('modif'); })       || keys[4];
                var kClass  = keys.find(function(k){ return k.toLowerCase().includes('class') || k.toLowerCase().includes('enzyme'); })   || keys[5];

                results.data.forEach(function (row) {
                    var gene = (row[kGene] || '').trim().toUpperCase();
                    if (!gene) return;
                    masterData[gene] = {
                        accession:    (row[kAcc]   || '').trim(),
                        activity:     (row[kAct]   || '').trim(),
                        weAnnotation: (row[kWE]    || '').trim(),
                        modification: (row[kMod]   || '').trim(),
                        enzymeClass:  (row[kClass] || '').trim()
                    };
                });

                console.log('Master proteins loaded:', Object.keys(masterData).length);
                cb(null);
            },
            error: function (err) { cb(err); }
        });
    }

    // 5b. STRING interactors TSV
    //     Columns: Accession | Gene Name | Interactors (;-sep) | Scores (;-sep)
    function parseStringFile(file, cb) {
        Papa.parse(file, {
            header: true,
            delimiter: '\t',
            skipEmptyLines: true,
            complete: function (results) {
                if (!results.data || results.data.length === 0) {
                    return cb(new Error('STRING TSV appears empty.'));
                }

                var keys         = Object.keys(results.data[0]);
                var kGene        = keys[1]; // col 2: gene name of the source protein
                var kInteractors = keys[2]; // col 3: semicolon-list of interactor gene names
                var kScores      = keys[3]; // col 4: semicolon-list of scores

                results.data.forEach(function (row) {
                    var src  = (row[kGene]        || '').trim().toUpperCase();
                    var ints = (row[kInteractors]  || '').trim();
                    var scrs = (row[kScores]       || '').trim();
                    if (!src || !ints || !scrs) return;

                    var intArr = ints.split(';');
                    var scrArr = scrs.split(';');
                    var len    = Math.min(intArr.length, scrArr.length);

                    for (var i = 0; i < len; i++) {
                        var tgt   = intArr[i].trim().toUpperCase();
                        var score = parseFloat(scrArr[i].trim());
                        if (!tgt || isNaN(score)) continue;
                        if (score > 1) score = score / 1000; // normalise 0-1000 → 0-1
                        stringEdges.push({ source: src, target: tgt, score: score });
                    }
                });

                console.log('STRING edges loaded:', stringEdges.length);
                cb(null);
            },
            error: function (err) { cb(err); }
        });
    }

    // 5c. Query file: one gene name per line
    function parseQueryFile(file, cb) {
        var reader = new FileReader();
        reader.onload = function (ev) {
            var text  = ev.target.result;
            var lines = text.split(/\r?\n/);
            queryGenes = [];
            lines.forEach(function (line) {
                var g = line.trim().toUpperCase();
                if (g) queryGenes.push(g);
            });
            console.log('Query genes:', queryGenes.length);
            cb(null);
        };
        reader.onerror = function (ev) { cb(new Error('FileReader error: ' + ev.target.error)); };
        reader.readAsText(file);
    }

    // ──────────────────────────────────────────────────────────────
    // 6. BUILD NETWORK
    // ──────────────────────────────────────────────────────────────
    function buildNetwork() {
        var threshold = parseFloat(elScoreSlider.value);

        // Query gene lookup set
        var querySet = {};
        queryGenes.forEach(function (g) { querySet[g] = true; });

        // Filter STRING edges: keep edges where at least one end is a query gene and score passes
        var validEdges = stringEdges.filter(function (e) {
            return e.score >= threshold && (querySet[e.source] || querySet[e.target]);
        });

        // Collect all unique gene nodes
        var nodeSet = {};
        queryGenes.forEach(function (g) { nodeSet[g] = true; });
        validEdges.forEach(function (e) {
            nodeSet[e.source] = true;
            nodeSet[e.target] = true;
        });

        // Build Cytoscape element array
        var elements = [];

        Object.keys(nodeSet).forEach(function (gene) {
            var info = masterData[gene] || null;
            var we   = info ? (info.weAnnotation || '').toUpperCase() : '';

            var role  = 'interactor';
            var shape = 'rectangle';
            var color = COLORS.interactor;

            // Check for W/E annotation
            var isW = we === 'W' || we === 'W/E' || we.split(';').some(function(t){ return t.trim() === 'W'; });
            var isE = we === 'E' || we === 'W/E' || we.split(';').some(function(t){ return t.trim() === 'E'; });

            if (isW && isE) {
                role = 'both'; shape = 'diamond'; color = COLORS.both;
            } else if (isW) {
                role = 'writer'; shape = 'triangle'; color = COLORS.writer;
            } else if (isE) {
                role = 'eraser'; shape = 'ellipse'; color = COLORS.eraser;
            }

            // Modification abbreviation (first 2 chars, skip NA)
            var modLabel = '';
            if (info && info.modification && info.modification.toUpperCase() !== 'NA') {
                modLabel = info.modification.trim().substring(0, 2);
            }

            // SVG background image to render modLabel *inside* the node
            var bgImage = 'none';
            if (modLabel) {
                var svgStr = '<svg xmlns="http://www.w3.org/2000/svg" width="40" height="40">'
                           + '<text x="50%" y="55%" font-family="Arial" font-size="14" font-weight="bold" '
                           + 'fill="white" text-anchor="middle" dominant-baseline="middle">'
                           + modLabel + '</text></svg>';
                bgImage = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgStr);
            }

            elements.push({
                group: 'nodes',
                data: {
                    id:       gene,
                    label:    gene,
                    role:     role,
                    shape:    shape,
                    color:    color,
                    bgImage:  bgImage,
                    modLabel: modLabel,
                    info:     info
                }
            });
        });

        validEdges.forEach(function (e) {
            var edgeId = e.source + '__' + e.target;
            elements.push({
                group: 'edges',
                data: { id: edgeId, source: e.source, target: e.target, score: e.score }
            });
        });

        // Destroy previous instance
        if (cy) { cy.destroy(); cy = null; }

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
                        'color':               '#f8fafc',
                        'font-size':           '11px',
                        'font-family':         'Inter, Arial, sans-serif',
                        'text-valign':         'bottom',
                        'text-halign':         'center',
                        'text-margin-y':       6,
                        'text-outline-color':  '#0f172a',
                        'text-outline-width':  2,
                        'width':               44,
                        'height':              44
                    }
                },
                {
                    selector: 'edge',
                    style: {
                        'width':              function(e){ return Math.max(1, e.data('score') * 6); },
                        'line-color':         '#475569',
                        'opacity':            0.7,
                        'curve-style':        'bezier'
                    }
                },
                {
                    selector: 'node:selected',
                    style: {
                        'border-width': 3,
                        'border-color': '#f8fafc'
                    }
                }
            ],
            layout: {
                name:           'cose',
                padding:        60,
                nodeRepulsion:  500000,
                idealEdgeLength: 120,
                edgeElasticity:  100,
                gravity:         200,
                numIter:         1000,
                animate:         true,
                animationDuration: 600
            }
        });

        // Node click → side panel
        cy.on('tap', 'node', function (evt) {
            renderDetails(evt.target.data());
        });
        cy.on('tap', function (evt) {
            if (evt.target === cy) {
                elNodeDetails.classList.add('hidden');
                elDetailsBody.innerHTML = '<p>Select a node to view details.</p>';
            }
        });
    }

    // ──────────────────────────────────────────────────────────────
    // 7. SIDE PANEL
    // ──────────────────────────────────────────────────────────────
    function renderDetails(data) {
        elNodeDetails.classList.remove('hidden');
        var info = data.info || {};
        var html = '<div class="detail-row">'
                 + '<span class="detail-label">Gene Name:</span> '
                 + '<span class="detail-value">' + (data.label || '—') + '</span>'
                 + '</div>';

        if (info && info.accession) {
            html += row('Accession',         info.accession   || 'N/A');
            html += row('Type',              info.weAnnotation|| data.role);
            html += row('Modification',      info.modification|| 'N/A');
            html += row('Catalytic Activity',info.activity    || 'N/A');
            html += row('Enzyme Class',      info.enzymeClass || 'N/A');
        } else {
            html += '<p style="color:#94a3b8;font-size:0.8rem">Pure interactor — not in master TSV.</p>';
        }
        elDetailsBody.innerHTML = html;
    }

    function row(label, value) {
        return '<div class="detail-row">'
             + '<span class="detail-label">' + label + ':</span> '
             + '<span class="detail-value">' + value + '</span>'
             + '</div>';
    }

}); // end DOMContentLoaded
