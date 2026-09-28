import os
import csv
import json
import re
import time
import requests
from flask import Flask, render_template, request, jsonify

app = Flask(__name__)

# ─────────────────────────────────────────────────────────────────
# DATA PATHS  — adjust if your files are elsewhere
# ─────────────────────────────────────────────────────────────────
BASE_DIR       = os.path.dirname(os.path.abspath(__file__))

def find_data_file(filename, fallback_parent_path=None):
    # Check 1: in BASE_DIR/data/filename
    p1 = os.path.join(BASE_DIR, "data", filename)
    if os.path.exists(p1):
        return p1
    # Check 2: in BASE_DIR/filename
    p2 = os.path.join(BASE_DIR, filename)
    if os.path.exists(p2):
        return p2
    # Check 3: in BASE_DIR/../filename
    p3 = os.path.join(BASE_DIR, "..", filename)
    if os.path.exists(p3):
        return p3
    # Check 4: fallback_parent_path
    if fallback_parent_path:
        p4 = os.path.join(fallback_parent_path, filename)
        if os.path.exists(p4):
            return p4
    return p1 # default

MASTER_TSV     = find_data_file("Master_Writer_Eraser_25Jun2026.tsv", "E:/PTM_WEB")
if not os.path.exists(MASTER_TSV):
    MASTER_TSV = find_data_file("Master_Writers_Erasers.tsv", "E:/PTM_WEB")
if not os.path.exists(MASTER_TSV):
    MASTER_TSV = find_data_file("Master_Writers_Erasers_Human.tsv", "E:/PTM_WEB")

STRING_TSV     = find_data_file("Interactors_25Jun2026.tsv", "E:/PTM_WEB")
if not os.path.exists(STRING_TSV):
    STRING_TSV = find_data_file("Interactors.tsv", "E:/PTM_WEB")
if not os.path.exists(STRING_TSV):
    STRING_TSV = find_data_file("string_interactors_9606.tsv", "E:/PTM_WEB")

PROTEOMES_TSV  = find_data_file("UniProt_Proteomes.tsv", "E:/PTM_WEB")
VIRUSES_TSV    = find_data_file("UniProt_Proteomes_Viruses.tsv", "E:/PTM_WEB")

# ─────────────────────────────────────────────────────────────────
# IN-MEMORY DATA STORES
# ─────────────────────────────────────────────────────────────────
master_writers_erasers = {} # { "KAT2A": { 'accessions': [...], 'taxon': '9606' } }
gene_to_accessions     = {} # { ("KAT2A", "9606"): ["Q92830"] }
accession_to_gene      = {} # { "Q92830": ("KAT2A", "9606") }
uniprot_data           = {} # { "Q92830": { 'fullname', 'existence', 'ec', 'activity', 'regulation', 'ptms', 'tissue', 'similarity', 'subunit' } }
string_lookup          = {} # { "KAT2A": [ { "target": "TP53", "score": 0.95 }, ... ] }
kinetic_lookup         = {} # { "KAT2A": { 'kcat', 'km', 'kcat_km', 'substrate', 'refs' } }
kinetic_by_acc         = {} # { "Q92830": { ... } }
kinetic_by_gene_tax    = {} # { ("KAT2A", "homo sapiens"): { ... } }
kinetic_by_gene        = {} # { "KAT2A": { ... } }
species_list           = [] # [ { 'id': '9606', 'name': 'Homo sapiens (Human)' }, ... ]
gene_to_species        = {} # { "KAT2A": "9606" }
gene_to_uniprot_accs   = {} # { "TP53": ["P04637"] }
species_id_to_name     = {} # { "9606": "Homo sapiens (Human)" }

# Helper to extract taxon ID from Organism ID string
def extract_taxon_id(org_id_str):
    if not org_id_str:
        return ''
    match = re.search(r'(?:ID|TaxID|TaxonID|Taxon|ID:)?\s*(\d+)', org_id_str, re.IGNORECASE)
    if match:
        return match.group(1)
    match_num = re.search(r'\d+', org_id_str)
    if match_num:
        return match_num.group(0)
    return ''

# ─────────────────────────────────────────────────────────────────
# STARTUP: load TSV files
# ─────────────────────────────────────────────────────────────────
def load_master_writers_erasers():
    global master_writers_erasers, gene_to_accessions, accession_to_gene, gene_to_species
    master_writers_erasers = {}
    gene_to_accessions = {}
    accession_to_gene = {}
    gene_to_species = {}
    
    path = MASTER_TSV
    if not os.path.exists(path):
        print(f"[WARN] Master Writers Erasers TSV not found: {path}")
        return
        
    print(f"[INFO] Loading master W/E from: {path}")
    with open(path, newline='', encoding='utf-8', errors='ignore') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            keys = list(row.keys())
            if not keys:
                continue
            gene = (row.get(keys[0]) or '').strip().upper()
            accessions_raw = (row.get(keys[1]) or '').strip()
            org_id = (row.get(keys[2]) or '').strip()
            
            if not gene or not accessions_raw:
                continue
                
            taxon = extract_taxon_id(org_id)
            if not taxon:
                taxon = '9606' # default
                
            accessions = [acc.strip() for acc in re.split(r'[;,|]', accessions_raw) if acc.strip()]
            
            catalytic = (row.get('Catalytic Activity') or row.get('Catalytic activity') or '').strip()
            annotation = (row.get('Annotation') or '').strip()
            modification = (row.get('Modification') or '').strip()
            enzyme_class = (row.get('Enzyme Class') or row.get('Enzyme class') or '').strip()
            
            key = (gene, taxon)
            if key not in master_writers_erasers:
                master_writers_erasers[key] = {
                    'accessions': [],
                    'taxon': taxon,
                    'catalytic': catalytic,
                    'annotation': annotation,
                    'modification': modification,
                    'enzyme_class': enzyme_class
                }
            for acc in accessions:
                if acc not in master_writers_erasers[key]['accessions']:
                    master_writers_erasers[key]['accessions'].append(acc)
            
            if key not in gene_to_accessions:
                gene_to_accessions[key] = []
            for acc in accessions:
                if acc not in gene_to_accessions[key]:
                    gene_to_accessions[key].append(acc)
            
            for acc in accessions:
                accession_to_gene[acc] = (gene, taxon)
                
            if gene not in gene_to_species:
                gene_to_species[gene] = set()
            gene_to_species[gene].add(taxon)
                
    print(f"[INFO] Loaded {len(master_writers_erasers)} master writers/erasers.")


def load_uniprot_info():
    global uniprot_data, gene_to_uniprot_accs, gene_to_species, gene_to_accessions
    uniprot_data = {}
    gene_to_uniprot_accs = {}
    
    path = find_data_file("UniProt_25Jun2026_Updated.tsv", "E:/PTM_WEB")
    if not os.path.exists(path):
        path = find_data_file("UniProt_25Jun2026.tsv", "E:/PTM_WEB")
    if not os.path.exists(path):
        path = find_data_file("UniProt.tsv", "E:/PTM_WEB")
    if not os.path.exists(path):
        path = find_data_file("130Unique_ManuallyCurated_UniProtInformation.tsv", "E:/PTM_WEB")
        
    if not os.path.exists(path):
        print(f"[WARN] UniProt TSV not found: {path}")
        return
        
    print(f"[INFO] Loading UniProt info from: {path}")
    with open(path, newline='', encoding='utf-8', errors='ignore') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            gene_name = (row.get('Gene Name') or row.get('Gene name') or '').strip().upper()
            accs_raw = (row.get('Protein Accessions') or row.get('Protein accessions') or '').strip()
            if not accs_raw:
                keys = list(row.keys())
                if keys:
                    accs_raw = (row.get(keys[1]) or row.get(keys[0]) or '').strip()
            if not accs_raw:
                continue
                
            accs = [a.strip() for a in re.split(r'[;,|]', accs_raw) if a.strip()]
            
            org_raw = (row.get('Organism ID') or row.get('Organism id') or row.get('organism id') or '').strip()
            taxon = extract_taxon_id(org_raw)
            
            if gene_name:
                gene_to_uniprot_accs[gene_name] = accs
                if taxon:
                    key = (gene_name, taxon)
                    if key not in gene_to_accessions:
                        gene_to_accessions[key] = []
                    for a in accs:
                        if a not in gene_to_accessions[key]:
                            gene_to_accessions[key].append(a)
            
            existence = (row.get('Protein existence') or row.get('Protein Existence') or '').strip()
            ec = (row.get('Updated EC Number') or '').strip()
            if not ec or ec == 'NA':
                ec = (row.get('EC number') or row.get('EC Number') or '').strip()
            if not ec:
                ec = 'NA'
            activity = (row.get('Catalytic activity') or row.get('Catalytic Activity') or '').strip()
            regulation = (row.get('Activity regulation') or row.get('Activity Regulation') or '').strip()
            ptms_raw = (row.get('PTMs') or row.get('PTM') or '').strip()
            tissue = (row.get('Tissue specificity') or row.get('Tissue Specificity') or row.get('tissue specifity') or '').strip()
            subunit = (row.get('Subunit') or row.get('subunit information') or row.get('Subunit information') or '').strip()
            string_ids = (row.get('STRING IDs') or row.get('STRING ID') or row.get('string ids') or '').strip()
            
            if gene_name and taxon:
                if gene_name not in gene_to_species:
                    gene_to_species[gene_name] = set()
                gene_to_species[gene_name].add(taxon)
            
            ptms = []
            if ptms_raw and ptms_raw != 'NA':
                parts = [p.strip() for p in ptms_raw.split(' | ') if p.strip()]
                for part in parts:
                    subparts = part.split('|')
                    if len(subparts) >= 1:
                        mod_name = subparts[0].strip()
                        residue = subparts[1].strip() if len(subparts) >= 2 else ''
                        pos = subparts[2].strip() if len(subparts) >= 3 else ''
                        ptms.append({
                            'name': mod_name,
                            'residue': residue,
                            'pos': pos
                        })
                        
            for acc in accs:
                uniprot_data[(acc, taxon)] = {
                    'fullname': '—',
                    'existence': existence or '—',
                    'ec': ec or '—',
                    'activity': activity or '—',
                    'regulation': regulation or '—',
                    'ptms': ptms,
                    'tissue': tissue or '—',
                    'similarity': '—',
                    'subunit': subunit or '—',
                    'string_ids': string_ids or '—',
                    'taxon': taxon or '—'
                }
    print(f"[INFO] Loaded UniProt info for {len(uniprot_data)} accessions.")


def load_interactors():
    global string_lookup
    string_lookup = {}
    
    path = STRING_TSV
    if not os.path.exists(path):
        print(f"[WARN] Interactors TSV not found: {path}")
        return
        
    print(f"[INFO] Loading interactors from: {path}")
    with open(path, newline='', encoding='utf-8', errors='ignore') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            keys = list(row.keys())
            if not keys or len(keys) < 3:
                # Handle old layout if fallback
                gene = (row.get(keys[0]) or '').strip().upper() if keys else ''
                if not gene:
                    continue
                int_raw = (row.get('Gene names of interactors') or '').strip()
                scores_raw = (row.get('string scores') or '').strip()
                uniprot_raw = (row.get('Interactors from UniProt') or '').strip()
                
                edges_dict = {}
                if int_raw and int_raw != 'NA':
                    int_list = [x.strip().upper() for x in int_raw.split(';') if x.strip()]
                    score_list = [x.strip() for x in scores_raw.split(';') if x.strip()]
                    for i, target in enumerate(int_list):
                        if i >= len(score_list):
                            break
                        try:
                            score = float(score_list[i])
                            if score > 1.0:
                                score = score / 1000.0
                            if score >= 0.40:
                                edges_dict[target] = round(score, 4)
                        except ValueError:
                            continue
                if uniprot_raw and uniprot_raw != 'NA':
                    uniprot_list = [x.strip().upper() for x in re.split(r'[;,]', uniprot_raw) if x.strip()]
                    for target in uniprot_list:
                        edges_dict[target] = 1.0
                edges = [{'target': target, 'score': score} for target, score in edges_dict.items() if target != gene]
                if edges:
                    edges.sort(key=lambda x: x.get('score', 0.0), reverse=True)
                    string_lookup[gene] = edges
                continue

            gene = (row.get(keys[0]) or '').strip().upper()
            if not gene:
                continue
                
            string_raw = (row.get(keys[1]) or '').strip()
            uniprot_raw = (row.get(keys[2]) or '').strip()
            
            edges = []
            seen_targets = set()
            
            if string_raw and string_raw != 'NA':
                string_list = [x.strip().upper() for x in string_raw.split(';') if x.strip()]
                for target in string_list:
                    if target != gene and target not in seen_targets:
                        seen_targets.add(target)
                        # Sourced from pre-sorted STRING list, assign dummy score
                        edges.append({'target': target, 'score': 1.0})
                        
            if uniprot_raw and uniprot_raw != 'NA':
                uniprot_list = [x.strip().upper() for x in re.split(r'[;,]', uniprot_raw) if x.strip()]
                for target in uniprot_list:
                    if target != gene and target not in seen_targets:
                        seen_targets.add(target)
                        edges.append({'target': target, 'score': 1.0})
                        
            if edges:
                edges.sort(key=lambda x: x.get('score', 0.0), reverse=True)
                string_lookup[gene] = edges
                
    print(f"[INFO] Loaded interactors for {len(string_lookup)} source proteins.")


def load_kinetic():
    global kinetic_lookup, kinetic_by_acc, kinetic_by_gene_tax, kinetic_by_gene
    kinetic_lookup = {}
    kinetic_by_acc = {}
    kinetic_by_gene_tax = {}
    kinetic_by_gene = {}
    
    path = find_data_file("UniProt_Brenda_Kinetic_Mapped.tsv", "E:/PTM_WEB")
    if not os.path.exists(path):
        path = find_data_file("Brenda_Kinetic_Paramters_PTM_WEB.tsv", "E:/PTM_WEB")
    if not os.path.exists(path):
        print(f"[WARN] Kinetic parameters TSV not found: {path}")
        return
        
    print(f"[INFO] Loading kinetic parameters from: {path}")
    loaded_rows = 0
    with open(path, newline='', encoding='utf-8', errors='ignore') as fh:
        # Header: Gene name \t Protein accession(s) \t Species \t EC number \t Annotation \t Source \t Recommended Name \t Systematic Name \t Reaction Type \t Turnover Number \t Km Value \t Kcat/Km Value \t Inhibitors \t Ki Value
        hdr = fh.readline()
        for line in fh:
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) < 14:
                continue
            loaded_rows += 1
            gene_name = parts[0].strip()
            gene_uc   = gene_name.upper()
            accs_raw  = parts[1].strip()
            species   = parts[2].strip()
            ec_num    = parts[3].strip()
            role      = parts[4].strip()
            source    = parts[5].strip()
            rec_name  = parts[6].strip()
            sys_name  = parts[7].strip()
            rxn_type  = parts[8].strip()
            turnover  = parts[9].strip()
            km_val    = parts[10].strip()
            kcat_km   = parts[11].strip()
            inhib     = parts[12].strip()
            ki_val    = parts[13].strip()
            
            record = {
                'gene': gene_name,
                'accession': accs_raw,
                'species': species,
                'ec': ec_num,
                'role': role,
                'source': source,
                'rec_name': rec_name,
                'sys_name': sys_name,
                'reaction_type': rxn_type,
                'turnover_no': turnover,
                'km': km_val,
                'kcat_km': kcat_km,
                'inhibitors': inhib,
                'ki': ki_val
            }
            
            if accs_raw and accs_raw != 'NA':
                for a in re.split(r'[;,|]', accs_raw):
                    a_clean = a.strip()
                    if a_clean and a_clean != 'NA':
                        kinetic_by_acc[a_clean] = record
                        
            if gene_uc and gene_uc != 'NA':
                if gene_uc not in kinetic_by_gene:
                    kinetic_by_gene[gene_uc] = record
                if species and species != 'NA':
                    sp_clean = species.lower()
                    kinetic_by_gene_tax[(gene_uc, sp_clean)] = record
                    # Extract binomial if different
                    m_binom = re.match(r'^([A-Z][a-z0-9_.-]+(?:\s+[a-z0-9_.-]+)+)', species, re.I)
                    if m_binom:
                        kinetic_by_gene_tax[(gene_uc, m_binom.group(1).lower())] = record
                        
            # Keep legacy lookup populated for backwards compatibility
            kinetic_lookup[gene_uc] = {
                'kcat': turnover if turnover != 'NA' else '—',
                'km': km_val if km_val != 'NA' else '—',
                'kcat_km': kcat_km if kcat_km != 'NA' else '—',
                'substrate': rec_name if rec_name != 'NA' else '—',
                'refs': source if source != 'NA' else '—'
            }
            
    print(f"[INFO] Loaded {loaded_rows} kinetic records ({len(kinetic_by_acc)} accessions, {len(kinetic_by_gene)} genes).")


def load_species():
    """Load and filter species from UniProt_Proteomes.tsv using BUSCO >= 90% filters and excluding viruses."""
    global species_list, species_id_to_name
    import re
    species_id_to_name = {}
    
    viral_ids = set()
    if os.path.exists(VIRUSES_TSV):
        try:
            with open(VIRUSES_TSV, newline='', encoding='utf-8') as fh:
                header = fh.readline()
                for line in fh:
                    parts = line.strip('\r\n').split('\t')
                    if len(parts) >= 3:
                        v_id = parts[2].strip()
                        if v_id:
                            viral_ids.add(v_id)
            print(f"[INFO] Loaded {len(viral_ids)} viral taxonomy IDs.")
        except Exception as e:
            print(f"[ERROR] Failed to load viral proteomes file: {e}")

    if not os.path.exists(PROTEOMES_TSV):
        print(f"[WARN] Proteomes TSV not found: {PROTEOMES_TSV}")
        species_list = [
            {'id': '9606', 'name': 'Homo sapiens (Human)'},
            {'id': '10090', 'name': 'Mus musculus (Mouse)'},
            {'id': '10116', 'name': 'Rattus norvegicus (Rat)'},
            {'id': '559292', 'name': "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)"},
            {'id': '83333', 'name': 'Escherichia coli (strain K12)'},
            {'id': '3702', 'name': 'Arabidopsis thaliana (Mouse-ear cress)'}
        ]
        for sp in species_list:
            species_id_to_name[sp['id']] = sp['name']
        return

    temp_list = []
    try:
        with open(PROTEOMES_TSV, newline='', encoding='utf-8') as fh:
            header = fh.readline()
            for line in fh:
                parts = line.strip('\r\n').split('\t')
                if len(parts) < 6:
                    continue
                org_name = parts[1].strip()
                org_id = parts[2].strip()
                busco_str = parts[4].strip()

                if org_id in viral_ids:
                    continue

                comp_m = re.search(r'C:([\.\d]+)%', busco_str)
                frag_m = re.search(r',F:([\.\d]+)%', busco_str)
                miss_m = re.search(r',M:([\.\d]+)%', busco_str)

                if not (comp_m and frag_m and miss_m):
                    continue

                comp = float(comp_m.group(1))
                frag = float(frag_m.group(1))
                miss = float(miss_m.group(1))

                if comp >= 90.0 and frag <= 5.0 and miss <= 5.0:
                    temp_list.append({
                        'id': org_id,
                        'name': org_name
                    })
                    species_id_to_name[org_id] = org_name

        temp_list.sort(key=lambda x: x['name'].lower())
        species_list = temp_list
        print(f"[INFO] Loaded {len(species_list)} filtered species.")
    except Exception as e:
        print(f"[ERROR] Failed to parse species file: {e}")
        species_list = [{'id': '9606', 'name': 'Homo sapiens (Human)'}]


# Classify Writer, Eraser, or Both dynamically
def classify_enzyme(ec_number, catalytic_activity):
    is_w = False
    is_e = False
    
    act_lower = catalytic_activity.lower()
    
    ecs = [ec.strip() for ec in ec_number.split(';') if ec.strip()]
    for ec in ecs:
        if ec.startswith('2.') or ec.startswith('6.'):
            is_w = True
        elif ec.startswith('3.') or ec.startswith('1.14.11'):
            is_e = True
            
    eraser_kws = ["deacetyl", "demethyl", "dephospho", "deubiquitin", "depalmitoyl", "deribosyl", "removal of", "hydrolysis", "deimination", "deglycosyl"]
    for kw in eraser_kws:
        if kw in act_lower:
            is_e = True
            
    writer_kws = ["acetyl", "methyl", "phospho", "ubiquitin", "palmitoyl", "ribosyl", "sulfotran", "myristoyl", "glucosyl", "sumoyl"]
    for kw in writer_kws:
        if kw in act_lower:
            idx = act_lower.find(kw)
            if idx >= 2 and act_lower[idx-2:idx] == "de":
                continue
            is_w = True
            
    if is_w and is_e:
        return 'both'
    elif is_w:
        return 'writer'
    elif is_e:
        return 'eraser'
    elif is_w or is_e:
        return 'writer' if is_w else 'eraser'
    return 'writer'

# Detect modification type dynamically
def detect_modification(catalytic_activity, ec_number):
    act_lower = (catalytic_activity + " " + ec_number).lower()
    if "phospho" in act_lower:
        return "Phosphorylation"
    elif "acetyl" in act_lower:
        return "Acetylation"
    elif "methyl" in act_lower:
        return "Methylation"
    elif "ubiquitin" in act_lower:
        return "Ubiquitination"
    elif "sumo" in act_lower:
        return "Sumoylation"
    elif "sulf" in act_lower:
        return "Sulfonation"
    elif "myristoyl" in act_lower:
        return "Myristoylation"
    elif "palmitoyl" in act_lower:
        return "Palmitoylation"
    elif "glycosyl" in act_lower:
        return "Glycosylation"
    elif "adp-ribosyl" in act_lower or "ribosyl" in act_lower:
        return "ADP-riboylation"
    return "NA"


# Run loaders at startup
load_master_writers_erasers()
load_uniprot_info()
load_interactors()
load_kinetic()
load_species()


# ─────────────────────────────────────────────────────────────────
# NETWORK BUILDER
# ─────────────────────────────────────────────────────────────────
def build_network(query_genes, query_taxon, limit=None):
    """
    Given a list of query genes, a target taxon ID, and a maximum interactor limit,
    returns { nodes: [...], edges: [...] } filtered to only contain genes of that taxon.
    """
    valid_query_genes = []
    for g in query_genes:
        if (g, query_taxon) in gene_to_accessions:
            valid_query_genes.append(g)
            
    query_set = set(valid_query_genes)
    node_set  = set(valid_query_genes)
    
    # 1. Collect top N interactors for each query gene that exist in query_taxon
    for gene in valid_query_genes:
        edges = string_lookup.get(gene, [])
        valid_edges = [edge for edge in edges if (edge['target'], query_taxon) in gene_to_accessions]
        valid_edges.sort(key=lambda x: x.get('score', 0.0), reverse=True)
        
        if limit is not None and isinstance(limit, int):
            selected_edges = valid_edges[:limit]
        else:
            selected_edges = valid_edges
            
        for edge in selected_edges:
            node_set.add(edge['target'])
            
    # 2. Add connections between ANY nodes that are in node_set
    edges_out = []
    for source in node_set:
        for edge in string_lookup.get(source, []):
            target = edge['target']
            if target in node_set:
                edges_out.append({'source': source, 'target': target, 'score': edge.get('score', 1.0)})
                
    seen = set()
    unique_edges = []
    for e in edges_out:
        key = tuple(sorted([e['source'], e['target']]))
        if key not in seen:
            seen.add(key)
            unique_edges.append(e)
            
    # 3. Build node objects
    nodes_out = []
    for gene in node_set:
        in_master = (gene, query_taxon) in master_writers_erasers
        
        # Determine accessions
        if in_master:
            master_info = master_writers_erasers[(gene, query_taxon)]
            accessions = master_info['accessions']
        else:
            accessions = gene_to_accessions.get((gene, query_taxon), [])
            if not accessions:
                accessions = gene_to_uniprot_accs.get(gene, [])
            
        accession_str = ";".join(accessions) if accessions else '—'
        first_acc = accessions[0] if accessions else ''
        info = uniprot_data.get((first_acc, query_taxon)) if first_acc else None
        
        if in_master:
            role = master_info['annotation'] or '—'
            role_lower = role.lower()
            if role_lower == 'writer-eraser' or role_lower == 'both' or role_lower == 'w/e':
                role_cytoscape = 'both'
            elif role_lower == 'writer' or role_lower == 'w':
                role_cytoscape = 'writer'
            elif role_lower == 'eraser' or role_lower == 'e':
                role_cytoscape = 'eraser'
            else:
                role_cytoscape = role_lower
                
            mod_type = master_info['modification'] or '—'
            mod_label = mod_type[:2] if (mod_type and mod_type != 'NA') else ''
            enzyme_class = master_info['enzyme_class'] or 'NA'
            if not enzyme_class or enzyme_class in ('—', 'NA'):
                enzyme_class = 'NA'
            activity = master_info['catalytic'] or '—'
            
            # Additional info from uniprot_data if available
            existence = info.get('existence', '—') if info else '—'
            regulation = info.get('regulation', '—') if info else '—'
            ptms_list = info.get('ptms', []) if info else []
            ptm_str = " | ".join([f"{p['name']}|{p['residue']}|{p['pos']}" for p in ptms_list]) if ptms_list else 'NA'
            subunit = info.get('subunit', '—') if info else '—'
            tissue = info.get('tissue', '—') if info else '—'
            string_ids = info.get('string_ids', '—') if info else '—'
            
            # EC Number from UniProt, fallback to Master catalytic activity regex
            ec_number = info.get('ec', 'NA') if info else 'NA'
            if (not ec_number or ec_number in ('—', 'NA')) and master_info.get('catalytic'):
                m_ec = re.search(r'EC=([0-9.-]+)', master_info['catalytic'])
                if m_ec:
                    ec_number = m_ec.group(1)
            if not ec_number or ec_number == '—':
                ec_number = 'NA'
                
            org_id = query_taxon
        else:
            role_cytoscape = 'interactor'
            role = 'Interactor'
            mod_type = 'NA'
            mod_label = ''
            ptm_str = 'NA'
            enzyme_class = 'NA'
            ec_number = 'NA'
            activity = '—'
            existence = '—'
            regulation = '—'
            subunit = '—'
            tissue = '—'
            string_ids = '—'
            org_id = query_taxon
            
            if info:
                ec_number = info.get('ec', 'NA')
                if not ec_number or ec_number == '—':
                    ec_number = 'NA'
                activity = info.get('activity', '—')
                existence = info.get('existence', '—')
                regulation = info.get('regulation', '—')
                ptms_list = info.get('ptms', [])
                ptm_str = " | ".join([f"{p['name']}|{p['residue']}|{p['pos']}" for p in ptms_list]) if ptms_list else 'NA'
                subunit = info.get('subunit', '—')
                tissue = info.get('tissue', '—')
                string_ids = info.get('string_ids', '—')
                
        # Resolve full species name
        species_name = species_id_to_name.get(org_id, f"TaxID: {org_id}") if org_id and org_id != '—' else '—'
        
        # Resolve kinetic parameters
        kin = None
        for a in accessions:
            if a in kinetic_by_acc:
                kin = kinetic_by_acc[a]
                break
        if not kin and species_name and species_name != '—':
            sp_lower = species_name.lower()
            m_sp = re.match(r'^([A-Z][a-z0-9_.-]+(?:\s+[a-z0-9_.-]+)+)', species_name, re.I)
            if (gene, sp_lower) in kinetic_by_gene_tax:
                kin = kinetic_by_gene_tax[(gene, sp_lower)]
            elif m_sp and (gene, m_sp.group(1).lower()) in kinetic_by_gene_tax:
                kin = kinetic_by_gene_tax[(gene, m_sp.group(1).lower())]
        if not kin:
            kin = kinetic_by_gene.get(gene)
            
        if not kin:
            kin = {
                'source': 'NA',
                'ec': ec_number if ec_number and ec_number != '—' else 'NA',
                'rec_name': 'NA',
                'sys_name': 'NA',
                'reaction_type': 'NA',
                'turnover_no': 'NA',
                'km': 'NA',
                'kcat_km': 'NA',
                'inhibitors': 'NA',
                'ki': 'NA'
            }
        
        nodes_out.append({
            'id':           gene,
            'label':        gene,
            'role':         role_cytoscape,
            'modLabel':     mod_label,
            'accession':    accession_str,
            'activity':     activity,
            'weAnnotation': role.upper() if in_master else '',
            'modification': mod_type,
            'enzymeClass':  enzyme_class,
            'ecNumber':     ec_number,
            'inMaster':     in_master,
            'isQuery':      gene in query_set,
            'ptms_raw':     ptm_str,
            # Kinetic Parameters
            'kinetic_source':        kin.get('source', 'NA') or 'NA',
            'kinetic_ec':            kin.get('ec', 'NA') or ec_number or 'NA',
            'kinetic_rec_name':      kin.get('rec_name', 'NA') or 'NA',
            'kinetic_sys_name':      kin.get('sys_name', 'NA') or 'NA',
            'kinetic_reaction_type': kin.get('reaction_type', 'NA') or 'NA',
            'kinetic_turnover_no':   kin.get('turnover_no', 'NA') or 'NA',
            'kinetic_km':            kin.get('km', 'NA') or 'NA',
            'kinetic_kcat_km':       kin.get('kcat_km', 'NA') or 'NA',
            'kinetic_inhibitors':    kin.get('inhibitors', 'NA') or 'NA',
            'kinetic_ki':            kin.get('ki', 'NA') or 'NA',
            # Legacy fields
            'km':           kin.get('km', '—'),
            'kcat':         kin.get('turnover_no', '—'),
            'kcat_km':      kin.get('kcat_km', '—'),
            'substrate':    kin.get('rec_name', '—'),
            'kinetic_refs': kin.get('source', '—'),
            # Details
            'fullname':     info.get('fullname', '—') if info else '—',
            'existence':    existence,
            'regulation':   regulation,
            'subunit':      subunit,
            'tissue':       tissue,
            'string_ids':   string_ids,
            'taxon':        species_name,
        })
        
    return {'nodes': nodes_out, 'edges': unique_edges}


# ─────────────────────────────────────────────────────────────────
# ROUTES
# ─────────────────────────────────────────────────────────────────
@app.route('/', methods=['GET'])
def home():
    return render_template('home.html', species_list=species_list, active_page='home')


@app.route('/about', methods=['GET'])
def about():
    return render_template('about.html', species_list=species_list, active_page='about')


@app.route('/download', methods=['GET'])
def download():
    return render_template('download.html', species_list=species_list, active_page='download')


@app.route('/contact', methods=['GET'])
def contact():
    return render_template('contact.html', species_list=species_list, active_page='contact')


@app.route('/faq', methods=['GET'])
def faq():
    return render_template('faq.html', species_list=species_list, active_page='faq')


@app.route('/query', methods=['POST'])
def query():
    raw_input  = request.form.get('proteins', '')
    species_val = request.form.get('species', '').strip()

    query_genes = [
        line.strip().upper()
        for line in raw_input.splitlines()
        if line.strip()
    ]

    if not query_genes:
        return render_template('home.html', species_list=species_list, error="Please enter at least one gene name.", active_page='home')

    taxon = extract_taxon_id(species_val)
    if not taxon:
        taxon = '9606'

    # Compute max limit for the interactor slider
    max_limit = 0
    for gene in query_genes:
        edges = string_lookup.get(gene, [])
        valid_edges = [edge for edge in edges if (edge['target'], taxon) in gene_to_accessions]
        num_i = len(valid_edges)
        if num_i > max_limit:
            max_limit = num_i
    if max_limit < 1:
        max_limit = 1

    # Initially fetch network with default 10 interactors (or max available)
    default_limit = min(10, max_limit) if max_limit >= 1 else 1
    network = build_network(query_genes, taxon, limit=default_limit)
    
    # Sort nodes so PTM Writers/Erasers/Both are first, then interactors
    network['nodes'].sort(key=lambda x: 0 if x['role'] in ['writer', 'eraser', 'both'] else 1)

    return render_template(
        'results.html',
        network_json  = json.dumps(network),
        query_genes   = query_genes,
        max_limit     = max_limit,
        default_limit = default_limit,
        node_count    = len(network['nodes']),
        edge_count    = len(network['edges']),
        nodes         = network['nodes'],
        selected_species = species_val,
    )


@app.route('/api/network', methods=['POST'])
def api_network():
    """AJAX endpoint — re-compute network with updated limit."""
    payload     = request.get_json(force=True)
    genes       = [g.strip().upper() for g in payload.get('genes', []) if g.strip()]
    limit_val   = payload.get('limit')
    species_val = payload.get('species', '').strip()
    taxon = extract_taxon_id(species_val)
    if not taxon:
        taxon = '9606'
    
    if limit_val == 'all':
        limit = None
    else:
        try:
            limit = int(limit_val)
        except:
            limit = None
            
    network = build_network(genes, taxon, limit)
    
    # Sort nodes
    network['nodes'].sort(key=lambda x: 0 if x['role'] in ['writer', 'eraser', 'both'] else 1)
    
    return jsonify(network)


# ─────────────────────────────────────────────────────────────────
# DOMAIN GRAPH — UniProt + InterPro (Pfam/SMART) domain data, merged
# with our own curated Master Writer/Eraser role, cached in memory.
# ─────────────────────────────────────────────────────────────────
DOMAIN_CACHE = {}                # { accession: {'ts': float, 'data': dict} }
DOMAIN_CACHE_TTL_SECONDS = 60 * 60 * 24   # domain/PTM annotations barely change — cache a day


def classify_ptm_category(feature_type, description):
    """Buckets a UniProt PTM feature into a chemistry category so the frontend can
    give it a distinct marker shape (phospho=triangle, acetyl=circle, etc). Mirrors
    the frontend's classifyPtmCategory() so backend and client-fallback paths agree."""
    d = (description or '').lower()
    if feature_type == 'Glycosylation':
        return 'glyco'
    if feature_type == 'Lipidation':
        return 'lipid'
    if 'phospho' in d:
        return 'phospho'
    if 'acetyl' in d:
        return 'acetyl'
    if any(k in d for k in ('glutaryl', 'succinyl', 'malonyl', 'crotonyl', 'propionyl', 'butyryl')):
        return 'acyl'
    if 'methyl' in d:
        return 'methyl'
    if 'ubiquitin' in d or 'sumo' in d or 'isopeptide' in d:
        return 'ubiquitin'
    if any(k in d for k in ('palmitoyl', 'myristoyl', 'prenyl', 'farnesyl', 'geranylgeranyl')):
        return 'lipid'
    if 'hydroxy' in d:
        return 'hydroxyl'
    if 'nitrat' in d:
        return 'nitration'
    if 'adp-ribosyl' in d or 'adp ribosyl' in d:
        return 'adpRibosyl'
    if 'citrullin' in d:
        return 'citrullin'
    return 'other'


def fetch_uniprot_domain_data(accession):
    """Returns (length, sequence, uniprot_domains, ptms). Raises on network/HTTP failure —
    UniProt is the primary source, so we don't want to silently hide a failure here."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    params = {"fields": "sequence,ft_domain,ft_mod_res,ft_carbohyd,ft_lipid,ft_crosslnk"}
    resp = requests.get(url, params=params, timeout=8)
    resp.raise_for_status()
    entry = resp.json()

    seq_obj = entry.get('sequence') or {}
    length = seq_obj.get('length')
    sequence = seq_obj.get('value')
    domains = []
    ptms = []
    for f in entry.get('features', []):
        loc = f.get('location') or {}
        start = ((loc.get('start') or {}).get('value'))
        end = ((loc.get('end') or {}).get('value'))
        if not start:
            continue
        ftype = f.get('type')
        desc = f.get('description') or ftype
        if ftype == 'Domain':
            domains.append({'start': start, 'end': end or start, 'name': desc})
        elif ftype in ('Modified residue', 'Cross-link', 'Lipidation', 'Glycosylation'):
            ptms.append({
                'position': start,
                'type': desc,
                'category': classify_ptm_category(ftype, desc),
            })
    return length, sequence, domains, ptms


def fetch_interpro_domain_tracks(accession):
    """Returns e.g. {'pfam': [...], 'smart': [...]}. Never raises — InterPro is a
    'nice to have' second opinion on domain boundaries, so any failure just means
    the frontend falls back to UniProt's own (coarser) Domain annotations."""
    tracks = {}
    try:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{accession}/"
        resp = requests.get(url, params={'page_size': 100}, timeout=8)
        resp.raise_for_status()
        data = resp.json()
        for result in data.get('results', []):
            meta = result.get('metadata') or {}
            db = (meta.get('source_database') or '').lower()
            if db not in ('pfam', 'smart'):
                continue  # keep the graph readable — just the two most common signature DBs
            name = meta.get('name') or meta.get('accession')
            for prot in result.get('proteins', []) or []:
                for loc in prot.get('entry_protein_locations') or []:
                    for frag in loc.get('fragments') or []:
                        start = frag.get('start')
                        end = frag.get('end')
                        if not start:
                            continue
                        tracks.setdefault(db, []).append({'start': start, 'end': end or start, 'name': name})
    except Exception as e:
        print(f"[WARN] InterPro lookup failed for {accession}: {e}")
    return tracks


@app.route('/api/domains/<accession>')
def api_domains(accession):
    accession = accession.strip().upper()
    if not accession:
        return jsonify({'error': 'No accession provided.'}), 400

    cached = DOMAIN_CACHE.get(accession)
    if cached and (time.time() - cached['ts'] < DOMAIN_CACHE_TTL_SECONDS):
        return jsonify(cached['data'])

    try:
        length, sequence, uniprot_domains, ptms = fetch_uniprot_domain_data(accession)
    except Exception as e:
        return jsonify({'error': f'Could not fetch UniProt data for {accession}: {e}'}), 502

    warnings = []
    # === PTMWEB_FIX_27SEP: Pfam/InterPro intentionally disabled (see chat) ===
    interpro_tracks = {}

    # Cross-reference against our own curated Master Writer/Eraser table so the
    # PTM markers are colored by OUR curated role, not a guess made on the frontend.
    gene_name, role = '', ''
    gene_taxon = accession_to_gene.get(accession)
    if gene_taxon:
        gene_name, taxon = gene_taxon
        master_info = master_writers_erasers.get((gene_name, taxon))
        if master_info:
            raw_role = (master_info.get('annotation') or '').strip().lower()
            if raw_role in ('writer-eraser', 'both', 'w/e'):
                role = 'both'
            elif raw_role in ('writer', 'w'):
                role = 'writer'
            elif raw_role in ('eraser', 'e'):
                role = 'eraser'
    for p in ptms:
        p['role'] = role or 'other'

    domain_tracks = {'uniprot': uniprot_domains}
    domain_tracks.update(interpro_tracks)

    result = {
        'accession': accession,
        'gene': gene_name,
        'length': length,
        'sequence': sequence,
        'ptms': ptms,
        'domain_tracks': domain_tracks,
        'warnings': warnings,
    }
    DOMAIN_CACHE[accession] = {'ts': time.time(), 'data': result}
    resp = jsonify(result)
    resp.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    return resp


@app.route('/api/detect_species', methods=['POST'])
def api_detect_species():
    """AJAX endpoint — detect all species containing the queried gene(s)."""
    data = request.get_json(silent=True) or {}
    genes_raw = data.get('genes') or data.get('gene', '')
    if isinstance(genes_raw, str):
        genes = [g.strip().upper() for g in re.split(r'[\s,;\n]+', genes_raw) if g.strip()]
    elif isinstance(genes_raw, list):
        genes = [str(g).strip().upper() for g in genes_raw if str(g).strip()]
    else:
        genes = []
        
    matched_taxons = set()
    for g in genes:
        tax_set = gene_to_species.get(g, set())
        if isinstance(tax_set, (set, list)):
            matched_taxons.update(tax_set)
        elif isinstance(tax_set, str) and tax_set:
            matched_taxons.add(tax_set)
        
    matched_species = []
    # Prioritize Human (9606) first, then alphabetical by species name
    sorted_taxons = sorted(
        list(matched_taxons),
        key=lambda tid: (0 if tid == '9606' else 1, species_id_to_name.get(tid, tid).lower())
    )
    
    for tid in sorted_taxons:
        name = species_id_to_name.get(tid)
        if name:
            matched_species.append({'id': tid, 'name': name})
        else:
            matched_species.append({'id': tid, 'name': f"Organism (TaxID: {tid})"})
            
    default_species = matched_species[0] if matched_species else None
    
    return jsonify({
        'species': matched_species,
        'default': default_species,
        'count': len(matched_species),
        'species_id': default_species['id'] if default_species else ''
    })


# ─────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True, use_reloader=True)