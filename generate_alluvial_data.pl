#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;

# Paths
my $master_file = "Master_Writer_Eraser_25Jun2026.tsv";
my $interactors_file = "Interactors_25Jun2026.tsv";
my $output_file = "alluvial_data.tsv";

# # Target modifications mapping
my %target_mods = (
    'phosphorylation' => 'Phosphorylation',
    'methylation'     => 'Methylation',
    'acetylation'     => 'Acetylation',
    'ubiquitination'  => 'Ubiquitination',
    'glycosylation'   => 'Glycosylation',
    'glycation'       => 'Glycation'
);

# Target species and their TaxID mappings in the master file
my %target_taxons = (
    '9606'   => 'Human',
    '10090'  => 'Mouse',
    '10116'  => 'Rat',
    '3702'   => 'Arabidopsis',
    '559292' => 'Yeast',
    '285006' => 'Yeast',
    '83333'  => 'E.coli',
    '83334'  => 'E.coli',
    '331111' => 'E.coli',
    '199310' => 'E.coli',
    '585035' => 'E.coli',
    '574521' => 'E.coli',
    '585055' => 'E.coli'
);

print "Step 1: Reading master file $master_file...\n";
open(my $fh, "<:encoding(UTF-8)", $master_file) or die "Cannot open $master_file: $!";
my $header = <$fh>;
chomp $header;
my @headers = split(/\t/, $header);

# Find column indices
my ($gene_col, $annot_col, $mod_col, $org_col);
for (my $i = 0; $i < @headers; $i++) {
    my $h = $headers[$i];
    if ($h =~ /^Gene Name$/i) { $gene_col = $i; }
    elsif ($h =~ /^Annotation$/i) { $annot_col = $i; }
    elsif ($h =~ /^Modification$/i) { $mod_col = $i; }
    elsif ($h =~ /^Organism ID$/i || $h =~ /^Organism$/i) { $org_col = $i; }
}

if (!defined $gene_col || !defined $annot_col || !defined $mod_col || !defined $org_col) {
    die "Could not find all required columns (Gene Name, Annotation, Modification, Organism ID) in master file.";
}

my %is_writer;
my %is_eraser;
my %enzyme_mods; # Keep track of mods for each enzyme

while (my $line = <$fh>) {
    chomp $line;
    my @cols = split(/\t/, $line);
    next if @cols < 3;
    
    my $gene = uc(trim($cols[$gene_col]));
    my $annot = trim($cols[$annot_col]);
    my $mod_field = trim($cols[$mod_col]);
    my $org_field = trim($cols[$org_col]);
    
    next unless $gene && $annot && $mod_field && $org_field;
    
    # Extract Taxon ID
    my $taxon = "";
    if ($org_field =~ /TaxID:(\d+)/i) {
        $taxon = $1;
    } elsif ($org_field =~ /(\d+)/) {
        $taxon = $1;
    }
    
    next unless exists $target_taxons{$taxon};
    
    # Check which target modifications match
    foreach my $mod_key (keys %target_mods) {
        if ($mod_field =~ /$mod_key/i) {
            my $standard_mod = $target_mods{$mod_key};
            
            if ($annot =~ /writer/i) {
                $is_writer{$gene}{$standard_mod} = 1;
                $enzyme_mods{$gene}{$standard_mod} = 1;
            }
            if ($annot =~ /eraser/i) {
                $is_eraser{$gene}{$standard_mod} = 1;
                $enzyme_mods{$gene}{$standard_mod} = 1;
            }
        }
    }
}
close($fh);

print "Loaded " . scalar(keys %enzyme_mods) . " enzymes from master file.\n";

print "Step 2: Parsing interactors file $interactors_file...\n";
open(my $ifh, "<:encoding(UTF-8)", $interactors_file) or die "Cannot open $interactors_file: $!";
my $iheader = <$ifh>;
chomp $iheader;
my @iheaders = split(/\t/, $iheader);

my ($igene_col, $istring_col, $iuni_col);
for (my $i = 0; $i < @iheaders; $i++) {
    my $h = $iheaders[$i];
    if ($h =~ /^Gene Name$/i) { $igene_col = $i; }
    elsif ($h =~ /^String Interactors$/i) { $istring_col = $i; }
    elsif ($h =~ /^UniProt Interactors$/i) { $iuni_col = $i; }
}

if (!defined $igene_col || !defined $istring_col || !defined $iuni_col) {
    die "Could not find all required columns (Gene Name, String Interactors, UniProt Interactors) in interactors file.";
}

# Mapping: $interactor_data{$mod}{$interactor}{writers|erasers} = { $enzyme => 1 }
my %interactor_data;

while (my $line = <$ifh>) {
    chomp $line;
    my @cols = split(/\t/, $line);
    next if @cols < 1;
    
    my $enzyme = uc(trim($cols[$igene_col]));
    next unless $enzyme;
    next unless exists $enzyme_mods{$enzyme}; # Only care about target enzymes
    
    # Parse interactors
    my %seen_interactors;
    
    # String interactors
    if (@cols > $istring_col) {
        my $str_field = trim($cols[$istring_col]);
        if ($str_field && $str_field ne 'NA') {
            foreach my $int (split(/;/, $str_field)) {
                $int = uc(trim($int));
                if ($int && $int ne $enzyme) {
                    $seen_interactors{$int} = 1;
                }
            }
        }
    }
    
    # UniProt interactors
    if (@cols > $iuni_col) {
        my $uni_field = trim($cols[$iuni_col]);
        if ($uni_field && $uni_field ne 'NA') {
            foreach my $int (split(/[;,]/, $uni_field)) {
                $int = uc(trim($int));
                if ($int && $int ne $enzyme) {
                    $seen_interactors{$int} = 1;
                }
            }
        }
    }
    
    # Map to modifications
    foreach my $mod (keys %{$enzyme_mods{$enzyme}}) {
        my $role_writer = exists $is_writer{$enzyme}{$mod} ? 1 : 0;
        my $role_eraser = exists $is_eraser{$enzyme}{$mod} ? 1 : 0;
        
        foreach my $int (keys %seen_interactors) {
            if ($role_writer) {
                $interactor_data{$mod}{$int}{writers}{$enzyme} = 1;
            }
            if ($role_eraser) {
                $interactor_data{$mod}{$int}{erasers}{$enzyme} = 1;
            }
        }
    }
}
close($ifh);

print "Step 3: Generating alluvial relationships...\n";
open(my $ofh, ">:encoding(UTF-8)", $output_file) or die "Cannot open $output_file for writing: $!";
print $ofh "Writer\tInteractor\tEraser\tModification\tType\n";

my $total_flows = 0;
my %stats; # { $mod => { total_substrates => X, shared => Y, writer_only => Z, eraser_only => W } }

foreach my $mod (sort values %target_mods) {
    next unless exists $interactor_data{$mod};
    
    my $substrates_ref = $interactor_data{$mod};
    my $total_sub = scalar(keys %$substrates_ref);
    $stats{$mod}{total_substrates} = $total_sub;
    $stats{$mod}{shared} = 0;
    $stats{$mod}{writer_only} = 0;
    $stats{$mod}{eraser_only} = 0;
    
    foreach my $int (sort keys %$substrates_ref) {
        my @writers = sort keys %{$substrates_ref->{$int}{writers} || {}};
        my @erasers = sort keys %{$substrates_ref->{$int}{erasers} || {}};
        
        my $has_writers = @writers > 0 ? 1 : 0;
        my $has_erasers = @erasers > 0 ? 1 : 0;
        
        if ($has_writers && $has_erasers) {
            $stats{$mod}{shared}++;
            foreach my $w (@writers) {
                foreach my $e (@erasers) {
                    print $ofh "$w\t$int\t$e\t$mod\tShared\n";
                    $total_flows++;
                }
            }
        }
        elsif ($has_writers) {
            $stats{$mod}{writer_only}++;
            foreach my $w (@writers) {
                print $ofh "$w\t$int\tNone\t$mod\tWriter_Only\n";
                $total_flows++;
            }
        }
        elsif ($has_erasers) {
            $stats{$mod}{eraser_only}++;
            foreach my $e (@erasers) {
                print $ofh "None\t$int\t$e\t$mod\tEraser_Only\n";
                $total_flows++;
            }
        }
    }
}
close($ofh);

print "Alluvial data written to $output_file. Total flows: $total_flows\n";
print "\n--- STATISTICS BY MODIFICATION ---\n";
printf("%-18s | %-16s | %-12s | %-12s | %-12s | %-15s\n", 
       "Modification", "Total Substrates", "Shared", "Writer-Only", "Eraser-Only", "% Shared");
print "-" x 95 . "\n";
foreach my $mod (sort keys %stats) {
    my $s = $stats{$mod};
    my $pct_shared = $s->{total_substrates} > 0 ? ($s->{shared} / $s->{total_substrates}) * 100 : 0;
    printf("%-18s | %-16d | %-12d | %-12d | %-12d | %-14.2f%%\n",
           $mod, $s->{total_substrates}, $s->{shared}, $s->{writer_only}, $s->{eraser_only}, $pct_shared);
}

sub trim {
    my $str = shift;
    return "" unless defined $str;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    # Strip quotes if present
    $str =~ s/^"//;
    $str =~ s/"$//;
    return $str;
}
