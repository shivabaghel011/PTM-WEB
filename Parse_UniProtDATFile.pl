#!/usr/bin/perl
# =============================================================================
# parse_uniprot_human.pl
#
# Parses a UniProt/Swiss-Prot flat file (.dat) and extracts a defined set of
# annotation fields for human proteins (or any specified taxon) into a TSV.
#
# INPUT FILES:
#   ARG1 : UniProt/Swiss-Prot flat file  (e.g. uniprot_sprot.dat)
#   ARG2 : UniProt PTM controlled vocabulary (e.g. ptmlist.dat)
#           Download from: https://ftp.uniprot.org/pub/databases/uniprot/
#                          current_release/knowledgebase/complete/
#   ARG3 : Output file prefix (default: Swissprot)
#   ARG4 : UniProt Proteome TSV (6-column, downloaded from UniProt Proteomes)
#           Used to filter entries by CPD classification and BUSCO scores.
#           Column 3 = Organism ID (taxon), Column 5 = BUSCO scores,
#           Column 6 = CPD classification.
#
# PROTEOME FILTERS (applied in order):
#   1. CPD classification must be "Standard" or "Close to standard"
#   2. BUSCO Completeness (C)  >= 100%
#      BUSCO Fragmented  (F)   <=   5%
#      BUSCO Missing     (M)   <=   5%
#
# OUTPUT:
#   <prefix>_All.tsv            - All proteins from passing proteomes
#   <prefix>_Writers_Erasers.tsv- Subset with writer/eraser reactions
#   species_statistics.tsv      - Per-organism entry counts
#
# USAGE EXAMPLE:
#   perl Parse_UniProtDATFile.pl \
#       uniprot_sprot.dat \
#       ptmlist.dat \
#       Swissprot \
#       proteomes.tsv
#
# OUTPUT COLUMNS (23 total):
#   1.  UniProt_Accession
#   2.  Organism_Species_And_ID
#   3.  Protein_Existence
#   4.  Biophysicochemical_Properties       5.  Biophysicochemical_Properties_Refs
#   6.  Catalytic_Activity                  7.  Catalytic_Activity_Refs
#   8.  Disease                             9.  Disease_Refs
#   10. Activity_Regulation                 11. Activity_Regulation_Refs
#   12. Interaction_Partners                13. Interaction_Partners_Refs
#   14. Metabolic_Pathways                  15. Metabolic_Pathways_Refs
#   16. PTMs                                17. PTMs_Refs
#   18. Subunit                             19. Subunit_Refs
#   20. Tissue_Specificity                  21. Tissue_Specificity_Refs
#   22. Similarity                          23. Similarity_Refs
#
# NOTES ON MULTI-ENTRY COLUMNS:
#   - Within a single column, multiple annotation items are separated by " | "
#   - The parallel *_Refs column follows the same positional order;
#     positions with no PubMed reference carry "NA"
#   - PTM entries are formatted as:  ModificationType|AminoAcidLetter|Position
#   - Biophysicochemical sub-properties are labelled, e.g.:
#       "Absorption: Abs(max)=450 nm | Kinetic parameters: KM=0.5 mM for ATP"
#
# =============================================================================

use strict;
use warnings;

# ---------------------------------------------------------------------------
# 0.  Command-line argument handling
# ---------------------------------------------------------------------------
my $dat_file        = shift;
my $ptm_file        = shift;
my $proteome_file   = shift;   # Optional: 6-column UniProt proteome TSV

my $next_arg = shift;
my ($virus_file, $annot_rxn_file);
if (defined $next_arg) {
    if ($next_arg =~ /virus/i) {
        $virus_file = $next_arg;
        $annot_rxn_file = shift;
    } else {
        $annot_rxn_file = $next_arg;
    }
}

my $string_dir      = shift || 'String_Networks'; # STRING network directory
my $out_prefix      = shift || 'E:/PTM_WEB/';

my $out_master_file;
my $out_uniprot_file;
my $out_interactors_file;

if ($out_prefix =~ /[\\\/]$/ || -d $out_prefix) {
    my $dir = $out_prefix;
    $dir .= '/' unless $dir =~ /[\\\/]$/;
    $out_master_file      = $dir . "Master_Writers_Erasers_23Jun.tsv";
    $out_uniprot_file     = $dir . "UniProt_23Jun.tsv";
    $out_interactors_file = $dir . "Interactors_23Jun.tsv";
} else {
    $out_master_file      = "${out_prefix}_Master_Writers_Erasers_23Jun.tsv";
    $out_uniprot_file     = "${out_prefix}_UniProt_23Jun.tsv";
    $out_interactors_file = "${out_prefix}_Interactors_23Jun.tsv";
}


# ---------------------------------------------------------------------------
# 1.  Amino acid full-name -> single-letter code table
#     Used when the ptmlist.dat TG line gives a full name, e.g. "Lysine" -> K
# ---------------------------------------------------------------------------
my %AA_LETTER = (
    'alanine'        => 'A',  'arginine'       => 'R',
    'asparagine'     => 'N',  'aspartate'      => 'D',
    'aspartic acid'  => 'D',  'cysteine'       => 'C',
    'glutamine'      => 'Q',  'glutamate'      => 'E',
    'glutamic acid'  => 'E',  'glycine'        => 'G',
    'histidine'      => 'H',  'isoleucine'     => 'I',
    'leucine'        => 'L',  'lysine'         => 'K',
    'methionine'     => 'M',  'phenylalanine'  => 'F',
    'proline'        => 'P',  'serine'         => 'S',
    'threonine'      => 'T',  'tryptophan'     => 'W',
    'tyrosine'       => 'Y',  'valine'         => 'V',
    # Plural forms sometimes appear
    'serines'        => 'S',  'threonines'     => 'T',
    'tyrosines'      => 'Y',  'lysines'        => 'K',
    'arginines'      => 'R',  'histidines'     => 'H',
    'cysteines'      => 'C',  'tryptophans'    => 'W',
);

# ---------------------------------------------------------------------------
# 2.  FT feature types that represent post-translational modifications
# ---------------------------------------------------------------------------
my %PTM_FT_TYPES = map { $_ => 1 } qw(MOD_RES);

# ---------------------------------------------------------------------------
# 3.  Named sub-properties within BIOPHYSICOCHEMICAL PROPERTIES CC blocks
# ---------------------------------------------------------------------------
my @BIOPHYS_SUBPROPS = (
    'Absorption',
    'Kinetic parameters',
    'pH dependence',
    'Redox potential',
    'Temperature dependence',
);

# ---------------------------------------------------------------------------
# 4.  Amino Acid Residue Hash for Writer/Eraser Filtering
# ---------------------------------------------------------------------------
my %AMINO_ACID_RESIDUES = map { $_ => 1 } qw(
    ala alanine alanyl
    arg arginine arginyl
    asn asparagine asparaginyl
    asp aspartate aspartyl aspartic
    cys cysteine cysteinyl
    gln glutamine glutaminyl
    glu glutamate glutamyl glutamic
    gly glycine glycyl
    his histidine histidyl
    ile isoleucine isoleucyl
    leu leucine leucyl
    lys lysine lysyl
    met methionine methionyl
    phe phenylalanine phenylalanyl
    pro proline prolyl
    ser serine seryl
    thr threonine threonyl
    trp tryptophan tryptophyl
    tyr tyrosine tyrosyl
    val valine valyl
    cit citrulline citrullyl
    orn ornithine ornithyl
    cystine cystyl
    homocysteine homocysteinyl
    isoaspartate isoaspartyl isoaspartic
);

# ---------------------------------------------------------------------------
# 4.  TSV column header definitions
# ---------------------------------------------------------------------------
my @TSV_HEADER = (
    'UniProt_Accession',
    'Organism_Species_And_ID',
    'Protein_Existence',
    'EC_Numbers',
    'Biophysicochemical_Properties',        'Biophysicochemical_Properties_Refs',
    'Catalytic_Activity',                   'Catalytic_Activity_Refs',
    'Disease',                              'Disease_Refs',
    'Activity_Regulation',                  'Activity_Regulation_Refs',
    'Interaction_Partners',                 'Interaction_Partners_Refs',
    'Metabolic_Pathways',                   'Metabolic_Pathways_Refs',
    'PTMs',                                 'PTMs_Refs',
    'Subunit',                              'Subunit_Refs',
    'Tissue_Specificity',                   'Tissue_Specificity_Refs',
    'Similarity',                           'Similarity_Refs',
);

# ===========================================================================
# MAIN EXECUTION
# ===========================================================================

print STDERR "[INFO] Loading PTM controlled vocabulary: $ptm_file\n";
my %PTM_VOCAB = parse_ptm_vocabulary($ptm_file);
printf STDERR "[INFO] Loaded %d PTM vocabulary entries.\n", scalar keys %PTM_VOCAB;

# Define the 6 target species taxon IDs
my %TARGET_TAXA = map { $_ => 1 } qw(83333 9606 10116 10090 3702 559292);
my %VALID_TAXA  = %TARGET_TAXA;

# Fallback organism names
my %FALLBACK_ORGANISM_NAMES = (
    '83333'  => 'Escherichia coli',
    '9606'   => 'Homo sapiens (Human)',
    '10116'  => 'Rattus norvegicus (Rat)',
    '10090'  => 'Mus musculus (Mouse)',
    '3702'   => 'Arabidopsis thaliana (Mouse-ear cress)',
    '559292' => 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker\'s yeast)',
);

# Load exact organism names from UniProt proteomes file
my %ORGANISM_NAMES = %FALLBACK_ORGANISM_NAMES;
if (defined $proteome_file && -e $proteome_file) {
    print STDERR "[INFO] Loading organism names from: $proteome_file\n";
    open(my $pfh, '<:encoding(UTF-8)', $proteome_file) or warn "[WARN] Cannot open proteome file '$proteome_file': $!\n";
    if ($pfh) {
        my $p_header = <$pfh>;
        while (my $line = <$pfh>) {
            chomp $line;
            my @cols = split(/\t/, $line);
            next if scalar @cols < 3;
            my $org_name = $cols[1];
            my $org_id   = $cols[2];
            $org_name =~ s/^\s+|\s+$//g;
            $org_id   =~ s/^\s+|\s+$//g;
            if (exists $TARGET_TAXA{$org_id}) {
                $ORGANISM_NAMES{$org_id} = $org_name;
            }
        }
        close($pfh);
    }
}

if (!defined $virus_file || $virus_file eq '' || !-e $virus_file) {
    if (defined $proteome_file) {
        my ($dir) = $proteome_file =~ m{^(.*[\\/])};
        $dir //= '';
        $virus_file = $dir . "UniProt_Proteomes_Viruses.tsv";
        if (!-e $virus_file) {
            $virus_file = "UniProt_Proteomes_Viruses.tsv";
        }
    }
}

my %viral_taxa;
if (defined $proteome_file && $proteome_file ne '') {
    if (defined $virus_file && -e $virus_file) {
        print STDERR "[INFO] Loading viral proteomes from: $virus_file\n";
        open(my $vh, '<', $virus_file) or warn "[WARN] Cannot open viral proteomes file '$virus_file': $!\n";
        if ($vh) {
            my $v_header = <$vh>;
            while (my $line = <$vh>) {
                chomp $line;
                next unless $line =~ /\S/;
                my @cols = split(/\t/, $line, -1);
                next if scalar @cols < 3;
                my $v_org_id = $cols[2];
                $v_org_id =~ s/^\s+|\s+$//g;
                if ($v_org_id =~ /^\d+$/) {
                    $viral_taxa{$v_org_id} = 1;
                }
            }
            close($vh);
            printf STDERR "[INFO] Loaded %d viral organism IDs.\n", scalar keys %viral_taxa;
        }
    }
}

if (defined $proteome_file && $proteome_file ne '') {
    print STDERR "[INFO] Loading proteome quality filter: $proteome_file\n";
    my %passed_taxa = load_proteome_filter($proteome_file, \%viral_taxa);
    printf STDERR "[INFO] %d organism IDs passed proteome quality filters.\n", scalar keys %passed_taxa;
    # Replace valid taxa with quality-filtered taxa
    %VALID_TAXA = %passed_taxa;
    # Ensure target taxa are always included
    for my $tx (keys %TARGET_TAXA) {
        $VALID_TAXA{$tx} = 1;
    }
    
    # Populate ORGANISM_NAMES for all valid taxa
    open(my $pfh, '<:encoding(UTF-8)', $proteome_file) or warn "[WARN] Cannot open proteome file '$proteome_file': $!\n";
    if ($pfh) {
        my $p_header = <$pfh>;
        while (my $line = <$pfh>) {
            chomp $line;
            my @cols = split(/\t/, $line);
            next if scalar @cols < 3;
            my $org_name = $cols[1];
            my $org_id   = $cols[2];
            $org_name =~ s/^\s+|\s+$//g;
            $org_id   =~ s/^\s+|\s+$//g;
            if (exists $VALID_TAXA{$org_id}) {
                $ORGANISM_NAMES{$org_id} = $org_name;
            }
        }
        close($pfh);
    }
}

print STDERR "[INFO] Parsing Swiss-Prot flat file: $dat_file\n";
process_swissprot($dat_file, $out_master_file, $out_uniprot_file, $out_interactors_file, \%PTM_VOCAB, \%VALID_TAXA, \%ORGANISM_NAMES, $string_dir, $annot_rxn_file);
print STDERR "[INFO] Files successfully compiled and written.\n";

# ===========================================================================
# SUBROUTINES
# ===========================================================================

# ---------------------------------------------------------------------------
# parse_ptm_vocabulary($file)
#
# Parses UniProt ptmlist.dat.
# Each entry in that file has the form:
#   ID   <description>        <- PTM name as it appears in FT /note fields
#   AC   PTM-NNNN
#   FT   MOD_RES | LIPID | ...
#   ...
#   KW   <Keyword>.           <- first KW used as modification type label
#   TG   <AminoAcid>.         <- target residue (full name)
#   //
#
# Returns %vocab:  lc(description) => { mod_type => STR, residue => LETTER }
# ---------------------------------------------------------------------------
sub parse_ptm_vocabulary {
    my ($file) = @_;
    my %vocab;

    open(my $fh, '<', $file)
        or die "[ERROR] Cannot open PTM vocabulary file '$file': $!\n";

    my ($id, $mod_type, $target_aa, $kw_seen);

    while (my $line = <$fh>) {
        chomp $line;

        if ($line =~ /^ID\s+(.+?)\s*$/) {
            # ---- New entry ----
            $id        = $1;
            $mod_type  = undef;
            $target_aa = undef;
            $kw_seen   = 0;
        }
        elsif ($line =~ /^KW\s+(.+?)\s*$/ && !$kw_seen) {
            # First KW line = modification type (e.g. "Acetylation.", "Phosphorylation.")
            ($mod_type  = $1) =~ s/\.\s*$//;
            $kw_seen = 1;
        }
        elsif ($line =~ /^TG\s+(.+?)\s*$/) {
            # Target amino acid residue
            my $aa_raw = lc($1);
            $aa_raw =~ s/\.\s*$//;
            # Strip possible alternative forms joined by " or "
            $aa_raw =~ s/\s+or\s+.+$//;
            $target_aa = $AA_LETTER{$aa_raw} // '';
        }
        elsif ($line =~ m{^//}) {
            # ---- End of entry ----
            if (defined $id && defined $mod_type) {
                # If TG was absent, try inferring residue from the ID name
                if (!defined $target_aa || $target_aa eq '') {
                    $target_aa = infer_residue_from_name($id);
                }
                $vocab{ lc($id) } = {
                    mod_type => $mod_type,
                    residue  => $target_aa // '',
                };
            }
            $id = $mod_type = $target_aa = undef;
            $kw_seen = 0;
        }
    }
    close($fh);
    return %vocab;
}

# ---------------------------------------------------------------------------
# infer_residue_from_name($ptm_name)
#
# Fallback: tries to identify the target amino acid from the PTM description
# by pattern matching common naming conventions (e.g. "acetyllysine" -> K).
# Returns a single-letter code or empty string.
# ---------------------------------------------------------------------------
sub infer_residue_from_name {
    my ($name) = @_;
    my $lc = lc($name);

    # Ordered list of amino acid name fragments to search for
    my @patterns = (
        [ qr/lysine|lys\b/,          'K' ],
        [ qr/serine|ser\b/,           'S' ],
        [ qr/threonine|thr\b/,        'T' ],
        [ qr/tyrosine|tyr\b/,         'Y' ],
        [ qr/cysteine|cys\b/,         'C' ],
        [ qr/arginine|arg\b/,         'R' ],
        [ qr/histidine|his\b/,        'H' ],
        [ qr/aspartate|asp\b/,        'D' ],
        [ qr/glutamate|glu\b/,        'E' ],
        [ qr/asparagine|asn\b/,       'N' ],
        [ qr/glutamine|gln\b/,        'Q' ],
        [ qr/methionine|met\b/,       'M' ],
        [ qr/tryptophan|trp\b/,       'W' ],
        [ qr/phenylalanine|phe\b/,    'F' ],
        [ qr/proline|pro\b/,          'P' ],
        [ qr/valine|val\b/,           'V' ],
        [ qr/leucine|leu\b/,          'L' ],
        [ qr/isoleucine|ile\b/,       'I' ],
        [ qr/alanine|ala\b/,          'A' ],
        [ qr/glycine|gly\b/,          'G' ],
    );

    for my $pair (@patterns) {
        my ($pat, $letter) = @{$pair};
        return $letter if $lc =~ $pat;
    }
    return '';
}

# ---------------------------------------------------------------------------
# process_swissprot($dat_file, $out_file, $target_taxon, \%ptm_vocab)
#
# Streams through the Swiss-Prot flat file one entry at a time.
# Each entry is terminated by a "//" line.
# Matching entries are formatted and written as TSV rows.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# load_proteome_filter($file)
#
# Parses the 6-column UniProt proteome TSV.
# Column layout (1-indexed, tab-separated):
#   1: Proteome ID
#   2: Organism name
#   3: Organism ID  <-- used as hash key
#   4: Protein count
#   5: BUSCO score  (format: C:99.9%[S:99.9%,D:0.0%],F:0.1%,M:0.0%,n:688)
#   6: CPD classification
#
# FILTER 1: CPD must be "Standard" or "Close to standard" (case-insensitive)
# FILTER 2 (applied only if FILTER 1 passes):
#   Completeness (C) >= 100%  AND  Fragmented (F) <= 5%  AND  Missing (M) <= 5%
#
# Returns %valid_taxa: { organism_id => 1 }
# ---------------------------------------------------------------------------
sub load_proteome_filter {
    my ($file, $viral_ref) = @_;
    my %valid;

    open(my $fh, '<', $file)
        or die "[ERROR] Cannot open proteome filter file '$file': $!\n";

    my $header = <$fh>;   # Skip header line
    my ($scanned, $passed, $skipped_viral, $skipped_busco) = (0, 0, 0, 0);

    while (my $line = <$fh>) {
        chomp $line;
        next unless $line =~ /\S/;   # Skip blank lines

        my @cols = split(/\t/, $line, -1);

        # Need at least 6 columns
        next if scalar @cols < 6;

        $scanned++;

        my $organism_id = $cols[2];   # Column 3 (0-indexed: 2)
        my $busco_str   = $cols[4];   # Column 5 (0-indexed: 4)

        # Trim whitespace
        $organism_id =~ s/^\s+|\s+$//g;
        $busco_str   =~ s/^\s+|\s+$//g;

        # Skip if viral
        if ($viral_ref && $viral_ref->{$organism_id}) {
            $skipped_viral++;
            next;
        }

        # ---- FILTER: Parse BUSCO scores ----
        # Format: C:99.9%[S:99.9%,D:0.0%],F:0.1%,M:0.0%,n:688
        my ($completeness, $fragmented, $missing);

        if ($busco_str =~ /C:([\.\d]+)%/)   { $completeness = $1; }
        if ($busco_str =~ /,F:([\.\d]+)%/)  { $fragmented   = $1; }
        if ($busco_str =~ /,M:([\.\d]+)%/)  { $missing      = $1; }
        
        #print "$completeness\t$fragmented\t$missing";<>;

        # If BUSCO string is missing or unparseable, skip this proteome
        unless (defined $completeness && defined $fragmented && defined $missing) {
            $skipped_busco++;
            next;
        }

        # Apply thresholds:
        #   Completeness >= 90%   Fragmented <= 5%   Missing <= 5%
        if ($completeness >= 90 && $fragmented <= 5 && $missing <= 5) {
            $valid{$organism_id} = 1;
            $passed++;
        } else {
            $skipped_busco++;
        }
    }
    close($fh);

    printf STDERR "[INFO]   Proteomes scanned      : %d\n",  $scanned;
    printf STDERR "[INFO]   Skipped (Viral)         : %d\n",  $skipped_viral;
    printf STDERR "[INFO]   Skipped (BUSCO)         : %d\n",  $skipped_busco;
    printf STDERR "[INFO]   Passed quality filter   : %d\n",  $passed;

    return %valid;
}

sub process_swissprot {
    my ($dat_file, $out_master, $out_uniprot, $out_interactors, $ptm_vocab_ref, $valid_taxa_ref, $organism_names_ref, $string_dir, $annot_rxn_file) = @_;

    # Load reaction mappings
    my %reaction_map;
    if (defined $annot_rxn_file && -e $annot_rxn_file) {
        print STDERR "[INFO] Loading annotated reactions from: $annot_rxn_file\n";
        open(my $amfh, '<:encoding(UTF-8)', $annot_rxn_file) or die "[ERROR] Cannot open annotated reaction file '$annot_rxn_file': $!\n";
        my $hdr = <$amfh>; # skip header
        while (my $line = <$amfh>) {
            chomp $line;
            next unless $line =~ /\S/;
            my @cols = split(/\t/, $line);
            next if scalar @cols < 5;
            my $rxn        = $cols[0];
            my $annotation = $cols[2] // '';
            my $mod        = $cols[3] // '';
            my $ec         = $cols[4] // '';
            next unless defined $rxn && $rxn ne '';
            my $norm_rxn = normalize_reaction($rxn);
            
            my $annot_mapped = '';
            if ($annotation =~ /^W$/i || $annotation =~ /writer/i) {
                $annot_mapped = 'writer';
            } elsif ($annotation =~ /^E$/i || $annotation =~ /eraser/i) {
                $annot_mapped = 'eraser';
            } elsif ($annotation =~ /^W\/E$/i || $annotation =~ /writer-eraser/i || $annotation =~ /both/i) {
                $annot_mapped = 'writer-eraser';
            }
            
            if ($annot_mapped eq 'writer-eraser') {
                $reaction_map{$norm_rxn}{Annotation}{writer} = 1;
                $reaction_map{$norm_rxn}{Annotation}{eraser} = 1;
            } elsif ($annot_mapped ne '') {
                $reaction_map{$norm_rxn}{Annotation}{$annot_mapped} = 1;
            }
            $reaction_map{$norm_rxn}{Modification}{$mod} = 1 if defined $mod && $mod ne '' && $mod ne 'NA';
            $reaction_map{$norm_rxn}{EnzymeClass}{$ec} = 1 if defined $ec && $ec ne '' && $ec ne 'NA';
        }
        close($amfh);
        printf STDERR "[INFO] Loaded %d annotated reactions.\n", scalar keys %reaction_map;
    } else {
        warn "[WARN] Annotated reaction file not provided or not found.\n";
    }

    print STDERR "[INFO] Pass 1: Scanning Swiss-Prot file to map genes by taxon...\n";
    my %genes_by_taxon;
    
    open(my $in_pass1, '<', $dat_file)
        or die "[ERROR] Cannot open '$dat_file' for first pass: $!\n";
        
    my ($cur_entry_id, $cur_taxon, $cur_gene_name, @cur_accs);
    
    while (my $line = <$in_pass1>) {
        chomp $line;
        my $tag = substr($line, 0, 2);
        my $content = (length($line) >= 5) ? substr($line, 5) : '';
        
        if ($tag eq 'ID') {
            if ($content =~ /^(\S+)/) {
                $cur_entry_id = $1;
            }
        }
        elsif ($tag eq 'AC') {
            my @parts = split(/;/, $content);
            for my $p (@parts) {
                $p =~ s/^\s+|\s+$//g;
                push @cur_accs, $p if $p;
            }
        }
        elsif ($tag eq 'OX') {
            if ($content =~ /NCBI_TaxID=(\d+)/) {
                $cur_taxon = $1;
            }
        }
        elsif ($tag eq 'GN') {
            if ($content =~ /Name=([^;{\s]+)/) {
                $cur_gene_name = $1;
            }
        }
        elsif ($line =~ m{^//}) {
            if (defined $cur_taxon && exists $valid_taxa_ref->{$cur_taxon}) {
                my $g = $cur_gene_name;
                if (!defined $g || $g eq '') {
                    $g = $cur_entry_id || $cur_accs[0] || 'NA';
                }
                $g =~ s/^\s+|\s+$//g;
                $g = uc($g);
                $genes_by_taxon{$cur_taxon}{$g} = 1;
            }
            $cur_entry_id = $cur_taxon = $cur_gene_name = undef;
            @cur_accs = ();
        }
    }
    close($in_pass1);
    
    my $mapped_taxa_count = scalar keys %genes_by_taxon;
    my $total_mapped_genes = 0;
    for my $tx (keys %genes_by_taxon) {
        $total_mapped_genes += scalar keys %{$genes_by_taxon{$tx}};
    }
    print STDERR "[INFO] Pass 1 complete. Mapped $total_mapped_genes genes across $mapped_taxa_count taxa.\n";

    open(my $in,  '<', $dat_file)
        or die "[ERROR] Cannot open '$dat_file': $!\n";
    open(my $fh_master, '>', $out_master)
        or die "[ERROR] Cannot open '$out_master' for writing: $!\n";
    open(my $fh_uniprot, '>', $out_uniprot)
        or die "[ERROR] Cannot open '$out_uniprot' for writing: $!\n";

    # Write headers
    print $fh_master "Gene Name\tProtein Accessions\tOrganism ID\tCatalytic Activity\tAnnotation\tModification\tEnzyme Class\n";
    print $fh_uniprot join("\t", 'Gene Name', 'Protein Accessions', 'Organism ID', 'Protein existence', 'EC number', 'Catalytic activity', 'Activity regulation', 'PTMs', 'Subunit', 'Tissue specificity', 'STRING IDs'), "\n";

    my (@entry_buf);
    my ($total, $written_master, $written_uniprot) = (0, 0, 0);

    my %seen_genes;
    my %seen_taxa;
    my %uniprot_interactors_map;
    my %string_id_to_gene_override;
    my %taxa_with_we;

    while (my $raw = <$in>) {
        chomp $raw;
        if ($raw =~ m{^//}) {
            if (@entry_buf) {
                $total++;
                
                my %entry_data = process_entry(\@entry_buf, $ptm_vocab_ref, $valid_taxa_ref, $organism_names_ref, \%reaction_map);
                if (%entry_data) {
                    my $gene_name   = $entry_data{gene_name};
                    my $accessions  = $entry_data{accessions};
                    my $organism    = $entry_data{organism};
                    my $taxon_id    = $entry_data{taxon_id};
                    my $pe_line     = $entry_data{pe_line};
                    my $ec_column   = $entry_data{ec_column};
                    my $catalytic   = $entry_data{catalytic};
                    my $regulation  = $entry_data{regulation};
                    my $ptms        = $entry_data{ptms};
                    my $subunit     = $entry_data{subunit};
                    my $tissue      = $entry_data{tissue};
                    my $is_we       = $entry_data{is_we};
                    my $uniprot_int = $entry_data{uniprot_int};
                    my $string_ids  = $entry_data{string_ids};
                    
                    # Write to UniProt for all valid quality-filtered species
                    if (defined $taxon_id && exists $valid_taxa_ref->{$taxon_id}) {
                        # But ONLY build interactors maps and track seen genes for the 6 target species
                        if (exists $TARGET_TAXA{$taxon_id}) {
                            $seen_genes{$gene_name} = 1;
                            $seen_taxa{$taxon_id} = 1;
                            
                            # Filter CC interactors against same species's proteome
                            if (exists $genes_by_taxon{$taxon_id}) {
                                my $species_genes_ref = $genes_by_taxon{$taxon_id};
                                for my $int (@$uniprot_int) {
                                    if (exists $species_genes_ref->{$int} && $int ne $gene_name) {
                                        push @{$uniprot_interactors_map{$gene_name}}, $int;
                                    }
                                }
                            }

                            if (defined $string_ids && $string_ids ne 'NA') {
                                foreach my $sid (split(/;/, $string_ids)) {
                                    $string_id_to_gene_override{$sid} = $gene_name;
                                }
                            }
                        }

                        print $fh_uniprot join("\t", $gene_name, $accessions, $organism, $pe_line, $ec_column, $catalytic, $regulation, $ptms, $subunit, $tissue, $string_ids), "\n";
                        $written_uniprot++;
                    }

                    # Write PTM writers and erasers for ALL valid quality-filtered species
                    if ($is_we) {
                        print $fh_master join("\t", $gene_name, $accessions, $organism, $catalytic, $entry_data{annotation}, $entry_data{modification}, $entry_data{enzyme_class}), "\n";
                        $written_master++;
                        if (defined $taxon_id) {
                            $taxa_with_we{$taxon_id} = 1;
                        }
                    }
                    
                    print STDERR "[INFO] Scanned: $total | UniProt: $written_uniprot | Master: $written_master\r" if $total % 1000 == 0;
                }
                @entry_buf = ();
            }
        }
        else {
            push @entry_buf, $raw;
        }
    }

    close($in);
    close($fh_uniprot);
    close($fh_master);
    print STDERR "\n[INFO] Swiss-Prot parsing finished. Scanned $total entries, wrote $written_uniprot to UniProt.tsv, $written_master to Master_Writers_Erasers.tsv.\n";

    # Now parse STRING networks
    print STDERR "[INFO] Parsing STRING networks...\n";
    my %string_interactors = parse_string_networks($string_dir, \%seen_genes, \%seen_taxa, \%string_id_to_gene_override);

    # Write Interactors.tsv
    print STDERR "[INFO] Writing Interactors.tsv...\n";
    open(my $fh_int, '>', $out_interactors)
        or die "[ERROR] Cannot open '$out_interactors' for writing: $!\n";
    print $fh_int "Gene Name\tString Interactors\tUniProt Interactors\n";

    my $written_int = 0;
    foreach my $gene (sort keys %seen_genes) {
        my @str_ints = @{$string_interactors{$gene} // []};
        my @uni_ints = @{$uniprot_interactors_map{$gene} // []};

        # Clean/deduplicate lists
        my %seen_str;
        @str_ints = grep { $_ ne $gene && !$seen_str{$_}++ } @str_ints;

        my %seen_uni;
        @uni_ints = grep { $_ ne $gene && !$seen_uni{$_}++ } @uni_ints;

        my $str_joined = @str_ints ? join(';', @str_ints) : 'NA';
        my $uni_joined = @uni_ints ? join(';', @uni_ints) : 'NA';

        print $fh_int "$gene\t$str_joined\t$uni_joined\n";
        $written_int++;
    }
    close($fh_int);
    print STDERR "[INFO] Wrote $written_int rows to Interactors.tsv.\n";

    # Give all species for which no PTM writers/erasers were found and those for which they were found
    my ($dir) = $out_master =~ m{^(.*[\\/])};
    $dir //= '';
    my $with_file = $dir . "species_with_writers_erasers.tsv";
    my $without_file = $dir . "species_without_writers_erasers.tsv";
    
    print STDERR "[INFO] Writing species summaries to:\n  - $with_file\n  - $without_file\n";
    
    open(my $wfh, '>:encoding(UTF-8)', $with_file) or warn "[WARN] Cannot open '$with_file' for writing: $!\n";
    open(my $wofh, '>:encoding(UTF-8)', $without_file) or warn "[WARN] Cannot open '$without_file' for writing: $!\n";
    
    if ($wfh) {
        print $wfh "Taxon ID\tSpecies Name\n";
    }
    if ($wofh) {
        print $wofh "Taxon ID\tSpecies Name\n";
    }
    
    my ($with_count, $without_count) = (0, 0);
    for my $tx (sort { ($organism_names_ref->{$a} // '') cmp ($organism_names_ref->{$b} // '') } keys %{$valid_taxa_ref}) {
        my $name = $organism_names_ref->{$tx} // "Unknown Species (TaxID:$tx)";
        if (exists $taxa_with_we{$tx}) {
            if ($wfh) {
                print $wfh "$tx\t$name\n";
            }
            $with_count++;
        } else {
            if ($wofh) {
                print $wofh "$tx\t$name\n";
            }
            $without_count++;
        }
    }
    
    close($wfh) if $wfh;
    close($wofh) if $wofh;
    
    print STDERR "[INFO] Species with PTM writers/erasers: $with_count\n";
    print STDERR "[INFO] Species without PTM writers/erasers: $without_count\n";
}

# ---------------------------------------------------------------------------
# process_entry(\@lines, $target_taxon, \%ptm_vocab)
#
# Parses one Swiss-Prot entry (lines between start and "//").
# Returns a TSV row string, or undef if the entry does not match the taxon.
# ---------------------------------------------------------------------------
sub process_entry {
    my ($lines_ref, $ptm_vocab_ref, $valid_taxa_ref, $organism_names_ref, $reaction_map_ref) = @_;

    # ---- Early taxon gate check ----
    if (defined $valid_taxa_ref && %{$valid_taxa_ref}) {
        my $found_ox = 0;
        my $entry_taxon = undef;
        for my $line (@{$lines_ref}) {
            if (substr($line, 0, 2) eq 'OX') {
                if (substr($line, 5) =~ /NCBI_TaxID=(\d+)/) {
                    $entry_taxon = $1;
                    $found_ox = 1;
                    last;
                }
            }
        }
        if ($found_ox && !exists $valid_taxa_ref->{$entry_taxon}) {
            return ();
        }
    }

    # ---- Declare accumulators ----
    my @accessions;
    my $taxon_num = '';
    my $pe_line   = '';
    my @os_parts;          # OS lines (organism species name, possibly multi-line)
    my @de_lines;          # DE lines
    my @cc_lines;          # CC line contents (after the 5-char "CC   " prefix)
    my @ft_lines;          # Full FT lines (kept intact for position parsing)
    my @ref_blocks;        # [{rn,rp,rx,ra,rt,rl}] – one per reference block
    my %cur_ref;
    my $in_ref = 0;
    my $gene_name = '';
    my $entry_id  = '';
    my @string_ids;

    # ---- Scan every line in the entry ----
    for my $line (@{$lines_ref}) {

        # The first 2 chars are the line-type code; content starts at position 5
        my $tag     = substr($line, 0, 2);
        my $content = (length($line) >= 5) ? substr($line, 5) : '';

        if ($tag eq 'ID') {
            if ($content =~ /^(\S+)/) {
                $entry_id = $1;
            }
        }
        elsif ($tag eq 'AC') {
            my @parts = split(/;/, $content);
            foreach my $p (@parts) {
                $p =~ s/^\s+|\s+$//g;
                push @accessions, $p if $p;
            }
        }
        elsif ($tag eq 'OS') {
            push @os_parts, $content;
        }
        elsif ($tag eq 'GN') {
            if ($content =~ /Name=([^;{\s]+)/) {
                $gene_name = $1;
            }
        }
        elsif ($tag eq 'DE') {
            push @de_lines, $content;
        }
        elsif ($tag eq 'OX') {
            # e.g.  "NCBI_TaxID=10090 {ECO:0000304|Ref.6};"
            if ($content =~ /NCBI_TaxID=(\d+)/) {
                $taxon_num = $1;
            }
        }
        elsif ($tag eq 'DR') {
            # e.g. "DR   STRING; 10090.ENSMUSP00000149664; -."
            if ($content =~ /^STRING;\s*([^;]+);/) {
                my $sid = $1;
                $sid =~ s/^\s+|\s+$//g;
                push @string_ids, $sid if $sid;
            }
        }
        elsif ($tag eq 'PE') {
            # e.g.  "1: Evidence at protein level;"
            ($pe_line = $content) =~ s/;\s*$//;
            $pe_line =~ s/^\s+|\s+$//g;
        }
        elsif ($tag eq 'CC') {
            # Strip leading spaces that UniProt uses for CC continuation lines
            push @cc_lines, $content;
        }
        elsif ($tag eq 'FT') {
            push @ft_lines, $line;   # Keep the full line; parse positions from it
        }
        elsif ($tag eq 'RN') {
            # Start of a new reference block
            if ($in_ref && %cur_ref) {
                push @ref_blocks, {%cur_ref};
            }
            %cur_ref = (rn => $content, rp => '', rx => '',
                        ra => '', rt => '', rl => '');
            $in_ref = 1;
        }
        elsif ($in_ref) {
            if    ($tag eq 'RP') { $cur_ref{rp} .= $content . ' '; }
            elsif ($tag eq 'RX') { $cur_ref{rx} .= $content . ' '; }
            elsif ($tag eq 'RA') { $cur_ref{ra} .= $content . ' '; }
            elsif ($tag eq 'RT') { $cur_ref{rt} .= $content . ' '; }
            elsif ($tag eq 'RL') { $cur_ref{rl} .= $content . ' '; }
        }
    }

    # Save the last reference block (no RN follows it)
    if ($in_ref && %cur_ref) {
        push @ref_blocks, {%cur_ref};
    }

    # ---- Taxon filter ----
    if (defined $valid_taxa_ref && %{$valid_taxa_ref}) {
        return () unless (defined $taxon_num && exists $valid_taxa_ref->{$taxon_num});
    }

    # ---- Gene name fallback ----
    if (!defined $gene_name || $gene_name eq '') {
        $gene_name = $entry_id || $accessions[0] || 'NA';
    }
    $gene_name =~ s/^\s+|\s+$//g;
    $gene_name = uc($gene_name);

    # ---- Organism string ----
    my $os_raw = join(' ', @os_parts);
    $os_raw =~ s/\s+/ /g;
    my $org_name = (defined $organism_names_ref && exists $organism_names_ref->{$taxon_num}) 
                   ? $organism_names_ref->{$taxon_num} 
                   : $os_raw;
    $org_name =~ s/^\s+|\s+$//g;
    $org_name =~ s/\.\s*$//;
    my $organism = "$org_name (TaxID:$taxon_num)";

    # ---- Parse EC Numbers from DE Lines ----
    my @ecs;
    foreach my $de (@de_lines) {
        if ($de =~ /EC=(\d+\.\d+\.\d+\.\d+)/) {
            push @ecs, $1;
        }
    }
    my $ec_column = @ecs ? join(';', @ecs) : 'NA';

    # ---- Join CC lines into one searchable string ----
    my $cc_text = join("\n", @cc_lines);

    # ---- Parse each annotation type from CC and FT lines ----
    my @catalytic   = parse_cc_topic($cc_text, 'CATALYTIC ACTIVITY');
    my @regulation  = parse_cc_topic($cc_text, 'ACTIVITY REGULATION');
    my @subunit     = parse_cc_topic($cc_text, 'SUBUNIT');
    my @tissue      = parse_cc_topic($cc_text, 'TISSUE SPECIFICITY');
    my @ptms        = parse_ft_ptms(\@ft_lines, $ptm_vocab_ref);

    # ---- Parse UniProt Interactors from CC Interaction lines ----
    my @uniprot_interactors;
    if ($cc_text =~ /-!-\s+INTERACTION:\s*(.*?)(?=\s*-!-|\z)/si) {
        my $block = $1;
        while ($block =~ /^\s*([A-Z0-9]{6,10}(?:-\d+)?);\s+([A-Z0-9]{6,10}(?:-\d+)?)(?::\s*([a-zA-Z0-9_-]+(?::[a-zA-Z0-9_-]+)?))?;/gm) {
            my $int_acc  = $2;
            my $int_gene = $3 // $2;
            push @uniprot_interactors, $int_gene if defined $int_gene && $int_gene ne '';
        }
    }
    my %seen_int;
    my @unique_uniprot_interactors;
    foreach my $int (@uniprot_interactors) {
        my $u = $int;
        $u =~ s/^\s+|\s+$//g;
        $u = uc($u);
        next if $u eq '' || $u eq $gene_name;
        push @unique_uniprot_interactors, $u unless $seen_int{$u}++;
    }

    # ---- Writer/Eraser Filter ----
    my $cat_text_tmp = join(' | ', map { $_->{text} // '' } @catalytic);
    my $is_writer_eraser = 0;
    
    # Extract all Reaction strings
    while ($cat_text_tmp =~ /Reaction=([^;\|]+)/g) {
        my $rxn = $1;

        # Split reaction by ' + ' or ' = ' to get individual participants
        my @parts = split(/\s+\+\s+|\s+=\s+/, $rxn);
        
        # Check if there are normal brackets inside the square brackets in any participant
        my $has_normal_brackets_inside = 0;
        foreach my $p (@parts) {
            my @bracket_contents = ($p =~ /\[(.*?)\]/g);
            foreach my $content (@bracket_contents) {
                if ($content =~ /[()]/) {
                    $has_normal_brackets_inside = 1;
                    last;
                }
            }
            last if $has_normal_brackets_inside;
        }
        
        # If any participant has normal brackets inside its square brackets, skip this reaction
        next if $has_normal_brackets_inside;

        my $reaction_is_we = 0;
        foreach my $p (@parts) {
            $p =~ s/^\s+|\s+$//g; # Trim whitespace
            
            # Check if it contains a square bracket block and is connected to outside text with a hyphen
            if ($p =~ /\[.*?\]/ && ($p =~ /-\s*\[/ || $p =~ /\]\s*-/)) {
                # Extract everything outside the brackets
                my $outside = $p;
                $outside =~ s/\[.*?\]//g;
                
                # Split the outside part by non-alphabet characters (e.g. hyphens, parentheses, spaces)
                my @words = split(/[^a-zA-Z]+/, $outside);
                foreach my $w (@words) {
                    next unless $w;
                    my $part = lc($w);
                    
                    # Lookup in amino acid hash
                    if (exists $AMINO_ACID_RESIDUES{$part}) {
                        $reaction_is_we = 1;
                        last;
                    }
                }
            }
            last if $reaction_is_we;
        }
        if ($reaction_is_we) {
            $is_writer_eraser = 1;
            last;
        }
    }

    my $catalytic_joined  = join_clean_text(\@catalytic);
    my $regulation_joined = join_clean_text(\@regulation);
    my $ptms_joined       = join_clean_text(\@ptms);
    my $subunit_joined    = join_clean_text(\@subunit);
    my $tissue_joined     = join_clean_text(\@tissue);

    my $accession_str = join(';', @accessions) || 'NA';

    my %entry_annots;
    my %entry_mods;
    my %entry_ecs;
    
    if ($is_writer_eraser && defined $reaction_map_ref) {
        if ($catalytic_joined ne 'NA' && $catalytic_joined ne '') {
            my @blocks = split(/\s*\|\s*/, $catalytic_joined);
            for my $block (@blocks) {
                if ($block =~ /Reaction=(.*?)(?:;\s*Xref=|\z)/i) {
                    my $rxn_eqn = $1;
                    my $norm = normalize_reaction($rxn_eqn);
                    if (exists $reaction_map_ref->{$norm}) {
                        for my $a (keys %{$reaction_map_ref->{$norm}{Annotation}}) {
                            $entry_annots{$a} = 1;
                        }
                        for my $m (keys %{$reaction_map_ref->{$norm}{Modification}}) {
                            $entry_mods{$m} = 1;
                        }
                        for my $e (keys %{$reaction_map_ref->{$norm}{EnzymeClass}}) {
                            $entry_ecs{$e} = 1;
                        }
                    }
                }
            }
        }
    }
    
    my $annot_val = %entry_annots ? join('|', sort keys %entry_annots) : 'NA';
    my $mod_val   = %entry_mods   ? join('|', sort keys %entry_mods)   : 'NA';
    my $ec_val    = %entry_ecs    ? join('|', sort keys %entry_ecs)    : 'NA';

    return (
        gene_name     => $gene_name,
        accessions    => $accession_str,
        organism      => $organism,
        taxon_id      => $taxon_num,
        pe_line       => $pe_line,
        ec_column     => $ec_column,
        catalytic     => $catalytic_joined,
        regulation    => $regulation_joined,
        ptms          => $ptms_joined,
        subunit       => $subunit_joined,
        tissue        => $tissue_joined,
        is_we         => $is_writer_eraser,
        uniprot_int   => \@unique_uniprot_interactors,
        string_ids    => join(';', @string_ids) || 'NA',
        annotation    => $annot_val,
        modification  => $mod_val,
        enzyme_class  => $ec_val,
    );
}

# ===========================================================================
# HELPER SUBROUTINES
# ===========================================================================

# ---------------------------------------------------------------------------
# build_pubmed_map(\@ref_blocks)
#
# Builds a hash mapping each PubMed ID found in the RX lines to a
# human-readable reference string: "PubMed:NNNNN | Authors | Title | Journal"
# ---------------------------------------------------------------------------
sub build_pubmed_map {
    my ($ref_list_ref) = @_;
    my %map;

    for my $ref (@{$ref_list_ref}) {
        my $rx = $ref->{rx} // '';

        while ($rx =~ /PubMed=(\d+)/g) {
            my $pmid = $1;
            next if exists $map{$pmid};   # first occurrence wins

            # Clean up each sub-field
            my $ra = _clean_ref_field($ref->{ra} // '');
            my $rt = _clean_ref_field($ref->{rt} // '');
            $rt =~ s/^"|"$//g;           # Remove surrounding quotes from title
            my $rl = _clean_ref_field($ref->{rl} // '');

            my $ref_str = "PubMed:$pmid";
            $ref_str .= " | Authors: $ra"   if $ra;
            $ref_str .= " | Title: $rt"     if $rt;
            $ref_str .= " | Journal: $rl"   if $rl;

            $map{$pmid} = $ref_str;
        }
    }
    return %map;
}

sub _clean_ref_field {
    my ($s) = @_;
    $s =~ s/\s+/ /g;
    $s =~ s/^\s+|\s+$//g;
    $s =~ s/;\s*$//;
    return $s;
}

# ---------------------------------------------------------------------------
# extract_pubmed_ids($text)
#
# Extracts all unique PubMed IDs from inline ECO evidence tags of the form:
#   {ECO:0000269|PubMed:12345678}
# Returns a list of numeric PubMed ID strings, in order of appearance.
# ---------------------------------------------------------------------------
sub extract_pubmed_ids {
    my ($text) = @_;
    my (%seen, @ids);
    while ($text =~ /PubMed:(\d+)/g) {
        push @ids, $1 unless $seen{$1}++;
    }
    return @ids;
}

# ---------------------------------------------------------------------------
# clean_text($text)
#
# Removes ECO evidence tags, normalises whitespace and trims trailing
# semicolons / punctuation that would clutter the output.
# ---------------------------------------------------------------------------
sub clean_text {
    my ($text) = @_;
    $text =~ s/\{ECO:[^}]+\}//g;   # Remove {ECO:...|PubMed:...} evidence tags
    $text =~ s/\s+/ /g;             # Collapse all whitespace to single space
    $text =~ s/^\s+|\s+$//g;        # Trim leading / trailing whitespace
    $text =~ s/\s*;\s*$//;          # Remove trailing semicolon
    $text =~ s/\s*\.\s*$//;         # Remove trailing period
    return $text;
}

# ---------------------------------------------------------------------------
# parse_cc_topic($cc_text, $topic)
#
# Finds every occurrence of "-!- TOPIC:" in the joined CC text and returns
# an array of annotation item hashrefs:
#   { text => <cleaned annotation text>, refs => [<PubMed IDs>] }
#
# Most topics (CATALYTIC ACTIVITY, DISEASE, etc.) can appear multiple times
# per entry; each block is a separate item.
# ---------------------------------------------------------------------------
sub parse_cc_topic {
    my ($cc_text, $topic) = @_;
    my @entries;

    # The regex:
    #   -!- TOPIC:    <- topic header
    #   \s*           <- optional whitespace/newlines following the colon
    #   (.*?)         <- lazily captured block content  (. matches \n with /s)
    #   (?=\s*-!-|\z) <- stops at the next topic header or end of CC text
    while ($cc_text =~ /-!-\s+\Q$topic\E:\s*(.*?)(?=\s*-!-|\z)/gsi) {
        my $block = $1;

        # Normalise to single line for ECO extraction before cleaning
        $block =~ s/\s+/ /g;
        $block =~ s/^\s+|\s+$//g;
        next unless $block;

        my @pmids = extract_pubmed_ids($block);
        my $clean = clean_text($block);
        next unless $clean;

        push @entries, { text => $clean, refs => \@pmids };
    }
    return @entries;
}

# ---------------------------------------------------------------------------
# parse_biophys_props($cc_text)
#
# Special-cases the BIOPHYSICOCHEMICAL PROPERTIES block, which contains
# named sub-sections (Absorption, Kinetic parameters, etc.).
# Each sub-section becomes a separate item, labelled with its name.
# Returns [{text, refs}] items.
# ---------------------------------------------------------------------------
sub parse_biophys_props {
    my ($cc_text) = @_;
    my @entries;

    return @entries
        unless $cc_text =~ /-!-\s+BIOPHYSICOCHEMICAL PROPERTIES:\s*(.*?)(?=\s*-!-|\z)/si;

    my $block = $1;

    # Build an alternation pattern for all known sub-property names
    my $alt = join('|', map { quotemeta($_) } @BIOPHYS_SUBPROPS);

    # Match each sub-property section: runs from the sub-property header until
    # the next sub-property header or the end of the biophys block
    my $found_any = 0;

    while ($block =~ /(?:^|\n)\s*($alt):\s*(.*?)(?=\n\s*(?:$alt):|\z)/gsi) {
        my ($sp_name, $sp_content) = ($1, $2);

        $sp_content =~ s/\s+/ /g;
        $sp_content =~ s/^\s+|\s+$//g;
        next unless $sp_content;

        my @pmids = extract_pubmed_ids($sp_content);
        my $clean = clean_text($sp_content);
        next unless $clean;

        push @entries, {
            text => "$sp_name: $clean",
            refs => \@pmids,
        };
        $found_any = 1;
    }

    # Fallback: if none of the named sub-sections were detected, treat the
    # entire block as a single item
    unless ($found_any) {
        $block =~ s/\s+/ /g;
        $block =~ s/^\s+|\s+$//g;
        if ($block) {
            my @pmids = extract_pubmed_ids($block);
            my $clean = clean_text($block);
            push @entries, { text => $clean, refs => \@pmids } if $clean;
        }
    }

    return @entries;
}

# ---------------------------------------------------------------------------
# parse_ft_ptms(\@ft_lines, \%ptm_vocab)
#
# Parses FT (Feature Table) lines to extract PTM features.
# Handles both the new (post-2021_02) and old UniProt FT line formats.
#
# New format:
#   FT   MOD_RES         100
#   FT                   /note="N6-acetyllysine"
#   FT                   /evidence="ECO:0000269|PubMed:12345678"
#
# Old format:
#   FT   MOD_RES        100    100       N6-acetyllysine.
#
# Returns [{text => "ModType|Residue|Position", refs => [PubMed IDs]}]
# ---------------------------------------------------------------------------
sub parse_ft_ptms {
    my ($ft_lines_ref, $ptm_vocab_ref) = @_;
    my @entries;

    my ($feat_type, $feat_pos, $feat_note, $feat_evidence);

    for my $line (@{$ft_lines_ref}) {

        # ---------- New format: feature-key line ----------
        # Pattern: FT   FEATUREKEY     <position_spec>
        # Position spec examples: 100   100..200   <1..200   100..>200
        if ($line =~ /^FT   ([A-Z_]+)\s+([\d<>?.]+(?:\.\.[\d<>?.]+)?)\s*$/) {

            # Before starting a new feature, save the previous one if it is a PTM
            if (defined $feat_type && $PTM_FT_TYPES{$feat_type}) {
                my $entry = make_ptm_entry(
                    $feat_type, $feat_pos, $feat_note, $feat_evidence, $ptm_vocab_ref
                );
                push @entries, $entry if defined $entry;
            }

            $feat_type     = $1;
            $feat_pos      = $2;
            $feat_note     = '';
            $feat_evidence = '';
        }

        # ---------- Old format: key + positions + description on one line ----------
        # Pattern: FT   FEATUREKEY   start   end   description text
        elsif ($line =~ /^FT   ([A-Z_]+)\s+(\d+)\s+\d+\s+(.+)$/) {

            if (defined $feat_type && $PTM_FT_TYPES{$feat_type}) {
                my $entry = make_ptm_entry(
                    $feat_type, $feat_pos, $feat_note, $feat_evidence, $ptm_vocab_ref
                );
                push @entries, $entry if defined $entry;
            }

            $feat_type     = $1;
            $feat_pos      = $2;
            ($feat_note    = $3) =~ s/\.\s*$//;
            $feat_evidence = '';
        }

        # ---------- New format: /note= qualifier ----------
        elsif ($line =~ /^FT\s+\/note="([^"]+)"/ && defined $feat_type) {
            $feat_note = $1;
            $feat_note =~ s/\.\s*$//;
        }

        # ---------- New format: /evidence= qualifier ----------
        elsif ($line =~ /^FT\s+\/evidence="([^"]+)"/ && defined $feat_type) {
            $feat_evidence = $1;
        }

        # ---------- Old-format continuation lines (description overflow) ----------
        # These start with FT followed by only whitespace then more description text
        elsif ($line =~ /^FT\s{20,}([^\/].+)$/ && defined $feat_type && $feat_note) {
            # Append to current note (rare in MOD_RES but possible)
            my $extra = $1;
            $extra =~ s/^\s+|\s+$//g;
            $feat_note .= ' ' . $extra if $extra;
        }
    }

    # ---- Save the last feature ----
    if (defined $feat_type && $PTM_FT_TYPES{$feat_type}) {
        my $entry = make_ptm_entry(
            $feat_type, $feat_pos, $feat_note, $feat_evidence, $ptm_vocab_ref
        );
        push @entries, $entry if defined $entry;
    }

    return @entries;
}

# ---------------------------------------------------------------------------
# make_ptm_entry($feat_type, $position, $note, $evidence, \%ptm_vocab)
#
# Converts a raw FT feature into a standardised PTM annotation hashref.
# Output text format: "ModificationType|AminoAcidLetter|Position"
#
# For DISULFID:  "Disulfide bond|C|position"
# For CROSSLNK:  "Cross-link (note)|?|position"
# For others:    vocabulary lookup -> "Phosphorylation|S|100"
#                fallback          -> "original note text|?|100"
# ---------------------------------------------------------------------------
sub make_ptm_entry {
    my ($feat_type, $position, $note, $evidence, $ptm_vocab_ref) = @_;
    return undef unless defined $note && $note ne '';

    my @pmids = extract_pubmed_ids($evidence // '');

    my $ptm_text;

    if ($feat_type eq 'DISULFID') {
        # Disulfide bonds always involve Cys; position may be "start..end"
        $ptm_text = "Disulfide bond|C|$position";
    }
    elsif ($feat_type eq 'CROSSLNK') {
        my $note_clean = clean_text($note);
        $ptm_text = "Cross-link ($note_clean)|?|$position";
    }
    else {
        # MOD_RES, LIPID, CARBOHYD
        # ---- 1. Exact lookup in PTM vocabulary (case-insensitive) ----
        my $note_lc = lc($note);
        $note_lc =~ s/[.;]\s*$//;

        if (exists $ptm_vocab_ref->{$note_lc}) {
            my $mod = $ptm_vocab_ref->{$note_lc}{mod_type};
            my $res = $ptm_vocab_ref->{$note_lc}{residue} || '?';
            $ptm_text = "$mod|$res|$position";
        }
        else {
            # ---- 2. Substring / partial match ----
            my $matched = 0;
            for my $key (sort keys %{$ptm_vocab_ref}) {
                if (   $note_lc eq $key
                    || index($note_lc, $key) >= 0
                    || index($key, $note_lc) >= 0 )
                {
                    my $mod = $ptm_vocab_ref->{$key}{mod_type};
                    my $res = $ptm_vocab_ref->{$key}{residue} || '?';
                    $ptm_text = "$mod|$res|$position";
                    $matched  = 1;
                    last;
                }
            }

            unless ($matched) {
                # ---- 3. Ultimate fallback: keep original note ----
                my $note_clean = clean_text($note);
                $ptm_text = "$note_clean|?|$position";
            }
        }
    }

    return { text => $ptm_text, refs => \@pmids };
}

# ---------------------------------------------------------------------------
# format_col(\@entries, \%pmid_ref_map)
#
# Converts an array of annotation items into two parallel column strings.
# Items are separated by " | " within each column.
#
# Returns a two-element list:  ($data_string, $refs_string)
#
# The ref string mirrors the data string position-for-position:
#   - If an item has PubMed refs  -> comma-separated ref strings
#   - If an item has NO refs      -> "NA"
#
# Example:
#   data : "Absorption: Abs(max)=450 nm | Kinetic parameters: KM=0.5 mM"
#   refs : "PubMed:12345 | NA"                 (second item lacked a ref)
# ---------------------------------------------------------------------------
sub format_col {
    my ($entries_ref, $pmid_ref_map_ref) = @_;
    return ('NA', 'NA') unless @{$entries_ref};

    my (@data_parts, @ref_parts);

    for my $item (@{$entries_ref}) {
        next unless defined $item && defined $item->{text} && $item->{text} ne '';

        push @data_parts, $item->{text};

        if ( $item->{refs} && @{ $item->{refs} } ) {
            my @ref_strs;
            for my $pmid ( @{ $item->{refs} } ) {
                if ( exists $pmid_ref_map_ref->{$pmid} ) {
                    push @ref_strs, $pmid_ref_map_ref->{$pmid};
                }
                else {
                    # PubMed ID found in ECO tag but not in the entry's R-lines
                    push @ref_strs, "PubMed:$pmid";
                }
            }
            push @ref_parts, join(', ', @ref_strs);
        }
        else {
            push @ref_parts, 'NA';
        }
    }

    return ('NA', 'NA') unless @data_parts;

    return (join(' | ', @data_parts), join(' | ', @ref_parts));
}

sub join_clean_text {
    my ($entries_ref) = @_;
    return 'NA' unless defined $entries_ref && @{$entries_ref};
    my @parts;
    for my $item (@{$entries_ref}) {
        if (defined $item && defined $item->{text} && $item->{text} ne '') {
            push @parts, $item->{text};
        }
    }
    return @parts ? join(' | ', @parts) : 'NA';
}

sub parse_string_networks {
    my ($string_dir, $seen_genes_ref, $seen_taxa_ref, $string_id_to_gene_override_ref) = @_;
    $string_id_to_gene_override_ref //= {};
    my %string_interactors;

    my %STRING_TO_UNIPROT_TAXON = (
        '511145' => '83333',
        '4932'   => '559292',
    );

    # Dynamically scan the String_Networks directory to identify available taxons
    my @string_taxa;
    if (opendir(my $dh, $string_dir)) {
        while (my $entry = readdir($dh)) {
            if ($entry =~ /^(\d+)\.protein\.links\.detailed\.v12\.0(?:\.txt)?$/) {
                push @string_taxa, $1;
            }
        }
        closedir($dh);
    } else {
        warn "[WARN] Cannot open STRING directory '$string_dir': $!\n";
    }

    # Iterate over each taxon
    foreach my $string_taxon (@string_taxa) {
        my $uniprot_taxon = $STRING_TO_UNIPROT_TAXON{$string_taxon} // $string_taxon;
        if (defined $seen_taxa_ref && %{$seen_taxa_ref}) {
            next unless exists $seen_taxa_ref->{$uniprot_taxon};
        }

        # 1. Find info mapping file
        my $info_file = "$string_dir/$string_taxon.protein.info.v12.0.txt";
        if (!-e $info_file) {
            $info_file = "$string_dir/$string_taxon.protein.info.v12.0";
        }
        
        my %id_to_gene;
        if (-e $info_file) {
            print STDERR "[INFO] Loading STRING info for taxon $string_taxon from: $info_file\n";
            open(my $infh, '<:encoding(UTF-8)', $info_file) or warn "[WARN] Cannot open info file '$info_file': $!\n";
            if ($infh) {
                my $dummy = <$infh>; # skip header
                while (my $line = <$infh>) {
                    chomp $line;
                    my @cols = split(/\t/, $line);
                    next if scalar @cols < 2;
                    my $str_id = trim_whitespace($cols[0]);
                    my $gene   = uc(trim_whitespace($cols[1]));
                    $id_to_gene{$str_id} = $gene;
                }
                close($infh);
            }
        } else {
            warn "[WARN] STRING info file not found for taxon $string_taxon: $info_file\n";
            next;
        }

        # 2. Find links detailed file
        my $links_file = "$string_dir/$string_taxon.protein.links.detailed.v12.0.txt";
        if (!-e $links_file) {
            $links_file = "$string_dir/$string_taxon.protein.links.detailed.v12.0";
        }

        if (-e $links_file) {
            print STDERR "[INFO] Parsing STRING links for taxon $string_taxon from: $links_file\n";
            open(my $lfh, '<', $links_file) or warn "[WARN] Cannot open links file '$links_file': $!\n";
            if ($lfh) {
                my $dummy = <$lfh>; # skip header
                
                # Find STRING protein IDs matching seen genes for fast filtering
                my %seen_string_ids;
                foreach my $str_id (keys %{$string_id_to_gene_override_ref}) {
                    my $gene = $string_id_to_gene_override_ref->{$str_id};
                    if (exists $seen_genes_ref->{$gene}) {
                        $seen_string_ids{$str_id} = 1;
                    }
                }
                foreach my $str_id (keys %id_to_gene) {
                    my $gene = $id_to_gene{$str_id};
                    if (exists $seen_genes_ref->{$gene}) {
                        $seen_string_ids{$str_id} = 1;
                    }
                }

                # Temp storage: gene_name -> { target_gene => score }
                my %temp_links;

                while (my $line = <$lfh>) {
                    # Fast extract first column (p1) to see if we care about it
                    my ($p1) = split(/\s+/, $line, 2);
                    next unless exists $seen_string_ids{$p1};

                    chomp $line;
                    my @cols = split(/\s+/, $line);
                    next if scalar @cols < 10;

                    my $score = $cols[-1];
                    next if $score < 400; # filter combined_score >= 400

                    my $p2 = $cols[1];

                    my $g1 = $string_id_to_gene_override_ref->{$p1} // $id_to_gene{$p1};
                    my $g2 = $string_id_to_gene_override_ref->{$p2} // $id_to_gene{$p2};
                    next unless defined $g1 && defined $g2;

                    $temp_links{$g1}{$g2} = $score;
                }
                close($lfh);

                # Sort by score descending and save
                foreach my $g1 (keys %temp_links) {
                    my $targets_ref = $temp_links{$g1};
                    my @sorted = sort { $targets_ref->{$b} <=> $targets_ref->{$a} } keys %{$targets_ref};
                    
                    my @cleaned = grep { $_ ne $g1 } @sorted;
                    $string_interactors{$g1} = \@cleaned;
                }
            }
        } else {
            warn "[WARN] STRING links file not found for taxon $string_taxon: $links_file\n";
        }
    }

    return %string_interactors;
}

sub trim_whitespace {
    my ($s) = @_;
    return '' unless defined $s;
    $s =~ s/^\s+|\s+$//g;
    return $s;
}

sub normalize_reaction {
    my ($r) = @_;
    return '' unless defined $r;
    $r =~ s/^\s+|\s+$//g;
    $r = lc($r);
    $r =~ s/\s+//g; # remove all spaces
    return $r;
}

# ---------------------------------------------------------------------------
# END OF SCRIPT
# ---------------------------------------------------------------------------
__END__

=head1 NAME

parse_uniprot_human.pl - Extract selected Swiss-Prot annotations for human
proteins into a structured TSV file.

=head1 SYNOPSIS

  perl parse_uniprot_human.pl \
      --dat   uniprot_sprot.dat \
      --ptm   ptmlist.dat \
      --out   human_swissprot.tsv \
      [--taxon 10090]

=head1 DESCRIPTION

Reads the UniProt/Swiss-Prot flat-file format (.dat) and writes a TSV where
each row represents one human protein (filtered by NCBI Taxon ID, default
10090).  Columns captured include:

  UniProt accession, organism, protein existence level, biophysicochemical
  properties, catalytic activity, associated diseases, activity regulation,
  interaction partners, metabolic pathways, PTMs (converted to
  ModType|AA|Position via ptmlist.dat), subunit structure, tissue
  specificity, and similarity annotations.

Every annotation column has a companion *_Refs column that mirrors it
position-for-position, listing PubMed references or "NA" where none exist.

=head1 INPUT FILES

=over 4

=item B<--dat>   uniprot_sprot.dat

Swiss-Prot flat-file.  Download from:
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/
    knowledgebase/complete/uniprot_sprot.dat.gz

=item B<--ptm>   ptmlist.dat

UniProt PTM controlled vocabulary.  Download from:
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/
    knowledgebase/complete/docs/ptmlist.txt

=back

=head1 OPTIONS

=over 4

=item B<--taxon> ID   (default: 10090)

NCBI Taxonomy ID to filter on.  Use 10090 for mouse, 7955 for zebrafish, etc.

=item B<--out>   FILE  (default: human_swissprot.tsv)

Output file path.

=back

=head1 AUTHOR

Generated for Swiss-Prot human proteome annotation extraction.

=cut