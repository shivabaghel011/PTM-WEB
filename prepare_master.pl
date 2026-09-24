#!/usr/bin/perl
use strict;
use warnings;

# ==============================================================================
# Configuration
# ==============================================================================
my $uniprot_dat   = "E:/PTM_WEB/uniprot_sprot.dat"; 
my $tsv1_proteins = "E:/PTM_WEB/Human_Writers_Erasers_MC_Parsing_Combined.tsv";
my $tsv2_reactions= "E:/PTM_WEB/Human_Writers_Erasers_Unique_Reactions_Annotated.tsv";
my $output_file   = "Master_Writers_Erasers_Human.tsv";
# ==============================================================================

print "1. Parsing UniProt DAT file to extract Gene Names...\n";
my %acc2gene;
open(my $fh_dat, '<', $uniprot_dat) or die "Cannot open $uniprot_dat: $!";
my @current_accs;
while (my $line = <$fh_dat>) {
    chomp $line;
    if ($line =~ /^AC\s+(.+)$/) {
        my $accs_str = $1;
        foreach my $acc (split(/;\s*/, $accs_str)) {
            push @current_accs, $acc if $acc;
        }
    } 
    elsif ($line =~ /^GN\s+Name=([^; {]+)/) {
        if (@current_accs) {
            my $gene = $1;
            foreach my $acc (@current_accs) {
                $acc2gene{$acc} = $gene;
            }
        }
    } 
    elsif ($line =~ /^\/\//) {
        @current_accs = (); 
    }
}
close($fh_dat);
my $gene_count = keys %acc2gene;
print "   -> Found gene names for $gene_count accessions.\n\n";

print "2. Parsing TSV 2 (Reactions Annotations) for Accessions...\n";
my %accession_annotations;
open(my $fh_tsv2, '<', $tsv2_reactions) or die "Cannot open $tsv2_reactions: $!";
my $header2 = <$fh_tsv2>; # Skip header
while (my $line = <$fh_tsv2>) {
    chomp $line;
    my @cols = split(/\t/, $line);
    
    # Columns based on: Square_bracket, Normal_bracket_inside, Amino_acid_outside, Reaction, Proteins, W/E, Modification, EnzymeClass
    my $proteins_col = $cols[4] || "";
    my $we           = $cols[5] || "";
    my $mod          = $cols[6] || "";
    my $eclass       = $cols[7] || "";
    
    # Check if W/E is actually valid
    if ($we && $we !~ /none|false|NA/i) {
        my @accessions = split(/;/, $proteins_col);
        foreach my $acc (@accessions) {
            $acc =~ s/^\s+|\s+$//g; # Trim
            if ($acc) {
                $accession_annotations{$acc}{we}{$we} = 1 if $we ne 'NA';
                $accession_annotations{$acc}{mod}{$mod} = 1 if $mod ne 'NA';
                $accession_annotations{$acc}{eclass}{$eclass} = 1 if $eclass ne 'NA';
            }
        }
    }
}
close($fh_tsv2);
my $annotated_acc_count = keys %accession_annotations;
print "   -> Loaded W/E annotations for $annotated_acc_count accessions.\n\n";

print "3. Parsing TSV 1 and generating output...\n";
open(my $fh_tsv1, '<', $tsv1_proteins) or die "Cannot open $tsv1_proteins: $!";
my $header1 = <$fh_tsv1>; # Skip header

open(my $out, '>', $output_file) or die "Cannot open $output_file: $!";
print $out "Gene name\tAccession\tCatalytic activity\tW/E annotation\tModification\tEnzyme class\n";

my $kept_count = 0;
my $skipped_count = 0;

while (my $line = <$fh_tsv1>) {
    chomp $line;
    my ($acc, $catalytic_activity) = split(/\t/, $line);
    
    if (!$acc) {
        next;
    }
    
    # Check if this accession has any W/E annotations mapped from TSV 2
    if (!exists $accession_annotations{$acc}) {
        $skipped_count++;
        next;
    }
    
    $kept_count++;
    
    my $gene = $acc2gene{$acc} || $acc;
    
    # Aggregate annotations (multiple separated by semicolon)
    my $all_we    = join(";", keys %{$accession_annotations{$acc}{we}});
    my $all_mod   = join(";", keys %{$accession_annotations{$acc}{mod}});
    my $all_class = join(";", keys %{$accession_annotations{$acc}{eclass}});
    
    print $out "$gene\t$acc\t$catalytic_activity\t$all_we\t$all_mod\t$all_class\n";
}
close($fh_tsv1);
close($out);

print "====================================================\n";
print "Done! Master file generated: $output_file\n";
print "Proteins kept (W/E annotated): $kept_count\n";
print "Proteins skipped (No W/E reaction): $skipped_count\n";
