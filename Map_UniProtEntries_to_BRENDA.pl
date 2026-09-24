#!/usr/bin/env perl
use strict;
use warnings;
use utf8;

# ==============================================================================
# Script: Map_UniProtEntries_to_BRENDA.pl
# Description: Maps proteins from UniProt_25Jun2026.tsv to BRENDA kinetic parameters
#              using two prioritized hash lookups:
#                1. EC Number + Species hash
#                2. Protein Accession + Species hash (fallback)
#              Also annotates Writer/Eraser/Both/Interactor using Master file.
#              Entries without a match are included with 'NA' in kinetic columns.
#
# Output file: E:\PTM_WEB\InteractionNetworks_Antigravity\UniProt_Brenda_Kinetic_Mapped.tsv
# ==============================================================================

# Input & output file paths
my $uniprot_file = $ARGV[0] || (
    -e 'E:/PTM_WEB/InteractionNetworks_Antigravity/UniProt_25Jun2026_Updated.tsv'
    ? 'E:/PTM_WEB/InteractionNetworks_Antigravity/UniProt_25Jun2026_Updated.tsv'
    : 'E:/PTM_WEB/InteractionNetworks_Antigravity/UniProt_25Jun2026.tsv'
);

my $brenda_file = $ARGV[1] || (
    -e 'E:/PTM_WEB/InteractionNetworks_Antigravity/Brenda_Kinetic_Parameters.tsv'
    ? 'E:/PTM_WEB/InteractionNetworks_Antigravity/Brenda_Kinetic_Parameters.tsv'
    : 'E:/PTM_WEB/Brenda_Kinetic_Parameters.tsv'
);

my $master_file = $ARGV[2] || 'E:/PTM_WEB/InteractionNetworks_Antigravity/Master_Writer_Eraser_25Jun2026.tsv';
my $output_file = $ARGV[3] || 'E:/PTM_WEB/InteractionNetworks_Antigravity/UniProt_Brenda_Kinetic_Mapped.tsv';

print "===============================================================================\n";
print "Map UniProt Entries to BRENDA Kinetic Parameters\n";
print "===============================================================================\n";
print "UniProt Input:   $uniprot_file\n";
print "BRENDA Input:    $brenda_file\n";
print "Master Input:    $master_file\n";
print "Mapped Output:   $output_file\n\n";

# Subroutine to clean strings
sub clean_str {
    my ($s) = @_;
    return 'NA' unless defined $s;
    $s =~ s/[\r\n\t]+/ /g;
    $s =~ s/\s+/ /g;
    $s =~ s/^\s+|\s+$//g;
    return ($s eq '' || $s eq '—') ? 'NA' : $s;
}

# Subroutine to extract candidate normalized species names from an organism string
sub get_candidate_species {
    my ($org_raw) = @_;
    return () unless defined $org_raw && $org_raw ne '' && $org_raw ne 'NA';

    my @candidates;
    my $s = $org_raw;
    $s =~ s/\s*\(TaxID:\d+\)//gi;
    $s =~ s/^\s+|\s+$//g;

    if ($s ne '') {
        push @candidates, lc($s);
    }

    # Extract text before first opening parenthesis (e.g. "Homo sapiens (Human)" -> "Homo sapiens")
    if ($s =~ /^([^(]+)/) {
        my $pre = $1;
        $pre =~ s/^\s+|\s+$//g;
        if ($pre ne '' && lc($pre) ne lc($s)) {
            push @candidates, lc($pre);
        }
    }

    # Extract standard binomial (Genus species)
    if ($s =~ /^([A-Z][a-z0-9_.-]+\s+[a-z0-9_.-]+)/) {
        my $binom = $1;
        $binom =~ s/^\s+|\s+$//g;
        push @candidates, lc($binom);
    }

    my %seen;
    return grep { !$seen{$_}++ } @candidates;
}

# Subroutine to extract the binomial scientific species name
sub extract_binomial_species {
    my ($org_raw) = @_;
    return 'NA' unless defined $org_raw && $org_raw ne '' && $org_raw ne 'NA';

    my $s = $org_raw;
    $s =~ s/\s*\(TaxID:\d+\)//gi;
    $s =~ s/^\s+|\s+$//g;

    # Strip parenthetical annotations (e.g. "Homo sapiens (Human)" -> "Homo sapiens")
    if ($s =~ /^([^(]+)/) {
        $s = $1;
        $s =~ s/^\s+|\s+$//g;
    }

    # Match binomial/trinomial scientific name (e.g. "Homo sapiens", "Escherichia coli")
    if ($s =~ /^([A-Z][a-z0-9_.-]+(?:\s+[a-z0-9_.-]+)+)/) {
        return $1;
    }
    return ($s ne '') ? $s : 'NA';
}

# ------------------------------------------------------------------------------
# STEP 1: Load Master Writer/Eraser catalog for role annotations
# ------------------------------------------------------------------------------
print "[1/4] Loading Master Writer/Eraser annotations...\n";
my %master_by_gene_tax;
my %master_by_gene_sp;
my %master_by_acc;
my $total_master_entries = 0;

if (-e $master_file) {
    open my $mfh, '<:encoding(UTF-8)', $master_file or die "Cannot open $master_file: $!\n";
    my $m_hdr = <$mfh>;
    while (my $line = <$mfh>) {
        chomp($line);
        $line =~ s/\r$//;
        $total_master_entries++;
        my @p = split(/\t/, $line, -1);
        next if @p < 5;

        my $gene     = uc(clean_str($p[0]));
        my $accs_raw = $p[1];
        my $org_raw  = $p[2];
        my $annot    = lc(clean_str($p[4]));

        my $role = 'Writer';
        if ($annot =~ /writer.*eraser|eraser.*writer|both|w\/e/i) {
            $role = 'Both';
        } elsif ($annot =~ /eraser/i) {
            $role = 'Eraser';
        } elsif ($annot =~ /writer/i) {
            $role = 'Writer';
        }

        my ($taxid) = ($org_raw =~ /\(TaxID:(\d+)\)/i);
        $taxid //= '';

        if ($gene ne '' && $gene ne 'NA') {
            if ($taxid ne '') {
                $master_by_gene_tax{$gene}{$taxid} = $role;
            }
            for my $sp (get_candidate_species($org_raw)) {
                $master_by_gene_sp{$gene}{$sp} = $role;
            }
        }

        for my $a (split(/[;,|]/, $accs_raw)) {
            $a =~ s/^\s+|\s+$//g;
            if ($a ne '' && $a ne 'NA') {
                $master_by_acc{$a} = $role;
            }
        }
    }
    close $mfh;
    print "      Loaded master annotations.\n";
} else {
    print "      [WARN] Master file not found. Defaulting non-matched to Interactor.\n";
}

# ------------------------------------------------------------------------------
# STEP 2: Parse BRENDA Kinetics file and build the two hashes
# ------------------------------------------------------------------------------
print "[2/4] Parsing BRENDA kinetics file and creating hashes...\n";
open my $bfh, '<:encoding(UTF-8)', $brenda_file or die "Cannot open $brenda_file: $!\n";
my $b_hdr = <$bfh>; # EC Number \t Species (Accession) \t RN \t SN \t RT \t TN \t KM \t KKM \t IN \t KI

# Hash 1: $hash_ec{$ec_number}{$species} = [ @kinetic_parameters, $ec_number ]
my %hash_ec;

# Hash 2: $hash_acc{$accession}{$species} = [ @kinetic_parameters, $ec_number ]
my %hash_acc;

# Hash 3: $hash_ec_general{$ec_number} = [ @kinetic_parameters, $ec_number ] (for transferred parameters)
my %hash_ec_general;

my $brenda_ec_entries = 0;
my $brenda_acc_keys   = 0;
my $total_brenda_entries = 0;

while (my $line = <$bfh>) {
    chomp($line);
    $line =~ s/\r$//;
    $total_brenda_entries++;
    my @parts = split(/\t/, $line, -1);
    next if @parts < 2;

    my $ec_num  = clean_str($parts[0]);
    my $sp_accs = $parts[1];

    # Kinetic parameter columns (index 2 to 9)
    my $rec_name     = clean_str($parts[2]);
    my $sys_name     = clean_str($parts[3]);
    my $react_type   = clean_str($parts[4]);
    my $turnover_no  = clean_str($parts[5]);
    my $km_val       = clean_str($parts[6]);
    my $kcat_km_val  = clean_str($parts[7]);
    my $inhibitors   = clean_str($parts[8]);
    my $ki_val       = clean_str($parts[9]);

    my @param_vals = (
        $rec_name,
        $sys_name,
        $react_type,
        $turnover_no,
        $km_val,
        $kcat_km_val,
        $inhibitors,
        $ki_val,
        $ec_num # store EC number as well for accession matches lacking EC
    );

    $brenda_ec_entries++;

    # Store generic parameters for this EC number (for transferred parameters across species)
    if ($ec_num ne 'NA' && $ec_num ne '') {
        $hash_ec_general{$ec_num} //= \@param_vals;
    }

    # Parse column 2: "Species (Accession); Species (Accession); ..."
    if (defined $sp_accs && $sp_accs ne '' && $sp_accs ne 'NA') {
        my @items = split(/;\s*(?=[^;()]+\s*\([A-Za-z0-9_.-]+\))/, $sp_accs);
        for my $item (@items) {
            $item =~ s/^\s+|\s+$//g;
            if ($item =~ /^([^()]+?)\s*\(([A-Za-z0-9_.-]+)\)$/) {
                my ($species, $acc) = ($1, $2);
                $species =~ s/^\s+|\s+$//g;
                $acc     =~ s/^\s+|\s+$//g;

                for my $cand_sp (get_candidate_species($species)) {
                    # Populate Hash 1: EC number + Species
                    if ($ec_num ne 'NA' && $ec_num ne '') {
                        $hash_ec{$ec_num}{$cand_sp} //= \@param_vals;
                    }

                    # Populate Hash 2: Protein Accession + Species
                    if ($acc ne '' && $acc ne 'NA') {
                        $hash_acc{$acc}{$cand_sp} //= \@param_vals;
                        $brenda_acc_keys++;
                    }
                }
            }
        }
    }
}
close $bfh;
print "      Processed $brenda_ec_entries BRENDA EC entries.\n";
print "      Populated EC and Accession hashes.\n\n";

# ------------------------------------------------------------------------------
# STEP 3: Parse UniProt TSV, match against hashes, and write mapped TSV
# ------------------------------------------------------------------------------
print "[3/4] Parsing UniProt file and performing prioritized matching...\n";
open my $ufh, '<:encoding(UTF-8)', $uniprot_file or die "Cannot open $uniprot_file: $!\n";
open my $ofh, '>:encoding(UTF-8)', $output_file  or die "Cannot open $output_file: $!\n";

# Output TSV Header
print $ofh join("\t",
    'Gene name',
    'Protein accession(s)',
    'Species',
    'EC number',
    'Annotation',
    'Source',
    'Recommended Name',
    'Systematic Name',
    'Reaction Type',
    'Turnover Number',
    'Km Value',
    'Kcat/Km Value',
    'Inhibitors',
    'Ki Value'
) . "\n";

my $u_hdr = <$ufh>;
chomp($u_hdr);
$u_hdr =~ s/\r$//;
my @u_cols = split(/\t/, $u_hdr, -1);

# Find column indices dynamically
my ($idx_gene, $idx_acc, $idx_org, $idx_ec) = (0, 1, 2, 4);
for my $i (0 .. $#u_cols) {
    my $c = lc($u_cols[$i]);
    if ($c eq 'gene name' || $c eq 'gene') { $idx_gene = $i; }
    elsif ($c eq 'protein accessions' || $c eq 'protein accession' || $c eq 'accession') { $idx_acc = $i; }
    elsif ($c eq 'organism id' || $c eq 'organism' || $c eq 'species') { $idx_org = $i; }
    elsif ($c eq 'updated ec number') { $idx_ec = $i; } # Prioritize updated EC column if present
    elsif ($c eq 'ec number' && $idx_ec == 4) { $idx_ec = $i; }
}

my $total_uniprot_rows          = 0;
my $matched_by_ec               = 0;
my $matched_by_acc              = 0;
my $matched_transferred         = 0;
my $unmatched_count             = 0;
my $we_brenda_count             = 0;
my $we_transferred_count        = 0;
my $we_unmapped_count           = 0;
my $interactor_brenda_count     = 0;
my $interactor_transferred_count= 0;
my $interactor_unmapped_count   = 0;

while (my $line = <$ufh>) {
    chomp($line);
    $line =~ s/\r$//;
    my @parts = split(/\t/, $line, -1);
    next unless @parts;

    $total_uniprot_rows++;
    if ($total_uniprot_rows % 50000 == 0) {
        print "      Processed $total_uniprot_rows UniProt rows...\n";
    }

    my $gene     = defined $parts[$idx_gene] ? clean_str($parts[$idx_gene]) : 'NA';
    my $accs_raw = defined $parts[$idx_acc]  ? clean_str($parts[$idx_acc])  : 'NA';
    my $org_raw  = defined $parts[$idx_org]  ? clean_str($parts[$idx_org])  : 'NA';
    my $ec_raw   = defined $parts[$idx_ec]   ? clean_str($parts[$idx_ec])   : 'NA';

    my $species_name = extract_binomial_species($org_raw);

    # Fallback to column 4 if updated EC column was NA
    if (($ec_raw eq 'NA' || $ec_raw eq '') && @parts > 4 && $idx_ec != 4) {
        my $c4 = clean_str($parts[4]);
        if ($c4 ne 'NA' && $c4 ne '') {
            $ec_raw = $c4;
        }
    }

    # Extract candidate species names and TaxID
    my @cand_species = get_candidate_species($org_raw);
    my ($taxid) = ($org_raw =~ /\(TaxID:(\d+)\)/i);
    $taxid //= '';

    # Determine Writer/Eraser/Both/Interactor annotation
    my $role = 'Interactor';
    my $gene_uc = uc($gene);
    if ($taxid ne '' && exists $master_by_gene_tax{$gene_uc}{$taxid}) {
        $role = $master_by_gene_tax{$gene_uc}{$taxid};
    } else {
        for my $sp (@cand_species) {
            if (exists $master_by_gene_sp{$gene_uc}{$sp}) {
                $role = $master_by_gene_sp{$gene_uc}{$sp};
                last;
            }
        }
    }
    if ($role eq 'Interactor') {
        for my $a (split(/[;,|]/, $accs_raw)) {
            $a =~ s/^\s+|\s+$//g;
            if (exists $master_by_acc{$a}) {
                $role = $master_by_acc{$a};
                last;
            }
        }
    }

    # --------------------------------------------------------------------------
    # Prioritized Matching:
    #   Step A: Check in EC Number hash (species-specific -> BRENDA)
    #   Step B: Fallback check in Protein Accession hash (species-specific -> BRENDA)
    #   Step C: Fallback check in EC General hash (cross-species -> TRANSFERRED)
    # --------------------------------------------------------------------------
    my $matched_params = undef;
    my $final_ec       = $ec_raw;
    my $source         = 'NA';

    # Step A: EC check (species-specific)
    if ($ec_raw ne 'NA' && $ec_raw ne '') {
        my @ecs_to_check = split(/[;, ]+/, $ec_raw);
        EC_SEARCH:
        for my $ec_cand (@ecs_to_check) {
            $ec_cand =~ s/^\s+|\s+$//g;
            next if $ec_cand eq '' || $ec_cand eq 'NA';
            for my $sp_cand (@cand_species) {
                if (exists $hash_ec{$ec_cand}{$sp_cand}) {
                    $matched_params = $hash_ec{$ec_cand}{$sp_cand};
                    $source         = 'BRENDA';
                    $matched_by_ec++;
                    last EC_SEARCH;
                }
            }
        }
    }

    # Step B: Fallback to Protein Accession check (species-specific)
    if (!defined $matched_params && $accs_raw ne 'NA' && $accs_raw ne '') {
        my @accs_to_check = split(/[;,| ]+/, $accs_raw);
        ACC_SEARCH:
        for my $acc_cand (@accs_to_check) {
            $acc_cand =~ s/^\s+|\s+$//g;
            next if $acc_cand eq '' || $acc_cand eq 'NA';
            for my $sp_cand (@cand_species) {
                if (exists $hash_acc{$acc_cand}{$sp_cand}) {
                    $matched_params = $hash_acc{$acc_cand}{$sp_cand};
                    $source         = 'BRENDA';
                    $matched_by_acc++;
                    # Fill EC number from Brenda if UniProt EC was NA
                    if ($final_ec eq 'NA' || $final_ec eq '') {
                        $final_ec = $matched_params->[8]; # Brenda EC stored at index 8
                    }
                    last ACC_SEARCH;
                }
            }
        }
    }

    # Step C: Fallback to Transferred EC check (generic parameters across species)
    if (!defined $matched_params && $final_ec ne 'NA' && $final_ec ne '') {
        my @ecs_to_check = split(/[;, ]+/, $final_ec);
        EC_TRANS_SEARCH:
        for my $ec_cand (@ecs_to_check) {
            $ec_cand =~ s/^\s+|\s+$//g;
            next if $ec_cand eq '' || $ec_cand eq 'NA';
            if (exists $hash_ec_general{$ec_cand}) {
                $matched_params = $hash_ec_general{$ec_cand};
                $source         = 'TRANSFERRED';
                $matched_transferred++;
                last EC_TRANS_SEARCH;
            }
        }
    }

    # Assign output parameter values
    my ($rec_name, $sys_name, $react_type, $turnover_no, $km_val, $kcat_km_val, $inhibitors, $ki_val);
    my $is_we = ($role eq 'Writer' || $role eq 'Eraser' || $role eq 'Both') ? 1 : 0;

    if (defined $matched_params) {
        if ($is_we) {
            if ($source eq 'BRENDA') {
                $we_brenda_count++;
            } elsif ($source eq 'TRANSFERRED') {
                $we_transferred_count++;
            }
        } else {
            if ($source eq 'BRENDA') {
                $interactor_brenda_count++;
            } elsif ($source eq 'TRANSFERRED') {
                $interactor_transferred_count++;
            }
        }
        $rec_name    = $matched_params->[0];
        $sys_name    = $matched_params->[1];
        $react_type  = $matched_params->[2];
        $turnover_no = $matched_params->[3];
        $km_val      = $matched_params->[4];
        $kcat_km_val = $matched_params->[5];
        $inhibitors  = $matched_params->[6];
        $ki_val      = $matched_params->[7];
    } else {
        $source = 'NA';
        $unmatched_count++;
        if ($is_we) {
            $we_unmapped_count++;
        } else {
            $interactor_unmapped_count++;
        }
        ($rec_name, $sys_name, $react_type, $turnover_no, $km_val, $kcat_km_val, $inhibitors, $ki_val) =
            ('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA');
    }

    # Print record to output TSV
    print $ofh join("\t",
        $gene,
        $accs_raw,
        $species_name,
        $final_ec,
        $role,
        $source,
        $rec_name,
        $sys_name,
        $react_type,
        $turnover_no,
        $km_val,
        $kcat_km_val,
        $inhibitors,
        $ki_val
    ) . "\n";
}

close $ufh;
close $ofh;

# ------------------------------------------------------------------------------
# STEP 4: Summary Report & Statistics
# ------------------------------------------------------------------------------
my $total_direct_brenda = $matched_by_ec + $matched_by_acc;
my $total_mapped        = $total_direct_brenda + $matched_transferred;
my $we_mapped_count     = $we_brenda_count + $we_transferred_count;
my $total_we            = $we_mapped_count + $we_unmapped_count;
my $interactor_mapped   = $interactor_brenda_count + $interactor_transferred_count;
my $total_interactors   = $interactor_mapped + $interactor_unmapped_count;

print "\n[4/4] MAPPING STATISTICS:\n";
print "===============================================================================\n";
print "a) Total number of entries in each of the three files:\n";
print "   - Master Writer/Eraser file:                 $total_master_entries\n";
print "   - BRENDA Kinetics file:                      $total_brenda_entries\n";
print "   - UniProt reference file:                    $total_uniprot_rows\n\n";

print "b) Total number of entries in UniProt mapped to BRENDA:\n";
print "   - Direct BRENDA match (Species-specific):    $total_direct_brenda (out of $total_uniprot_rows)\n";
print "     * Mapped via EC Number:                    $matched_by_ec\n";
print "     * Mapped via Protein Accession:            $matched_by_acc\n";
print "   - Transferred match (EC-level parameters):   $matched_transferred\n";
print "   - Total entries with kinetic parameters:     $total_mapped\n";
print "   - Unmapped entries (NA):                     $unmatched_count\n\n";

print "c) Total number of Writers/Erasers mapped from BRENDA:\n";
print "   - Writers/Erasers direct BRENDA match:       $we_brenda_count\n";
print "   - Writers/Erasers transferred parameters:    $we_transferred_count\n";
print "   - Total Writers/Erasers with parameters:     $we_mapped_count\n\n";

print "d) Total number of Writers/Erasers with NO mapping (parameters written as NA):\n";
print "   - Writers/Erasers unmapped (NA):             $we_unmapped_count\n";
print "   - Total Writers/Erasers in UniProt dataset:  $total_we\n\n";

print "Additional Info (Interactors):\n";
print "   - Interactors direct BRENDA match:           $interactor_brenda_count\n";
print "   - Interactors transferred parameters:        $interactor_transferred_count\n";
print "   - Interactors unmapped (NA):                 $interactor_unmapped_count\n";
print "   - Total Interactors in UniProt dataset:      $total_interactors\n";
print "-------------------------------------------------------------------------------\n";
print "Output saved to: $output_file\n";
print "===============================================================================\n";
print "Completed successfully.\n";
