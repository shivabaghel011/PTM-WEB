#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl Extract_Domain_Data.pl uniprot_sprot.dat
my $dat_file = shift @ARGV or die "Usage: perl Extract_Domain_Data.pl uniprot_sprot.dat\n";
my $out_file = "Domain_Graph_Data.tsv";

open(my $in, '<', $dat_file) or die "Cannot open $dat_file: $!\n";
open(my $out, '>', $out_file) or die "Cannot open $out_file: $!\n";

print $out "Gene_Name\tAccessions\tTaxID\tSpecies\tTotal Amino Acid\tDomains\tVariants\tSequence\tDomain Sequence\n";

my ($gene, $taxid, $species, $length, $sequence) = ("", "", "", "", "");
my @accs;
my (@domains, @variants);
my ($ft_type, $ft_loc, $ft_note) = ("", "", "");
my $in_sq = 0; 

sub save_feature {
    if ($ft_type eq 'DOMAIN' && $ft_note) {
        push @domains, "$ft_loc:$ft_note";
    } elsif (($ft_type eq 'VARIANT' || $ft_type eq 'VAR_SEQ') && $ft_note) {
        push @variants, "$ft_loc:$ft_note";
    }
    $ft_type = ""; $ft_loc = ""; $ft_note = "";
}

print "Parsing $dat_file. Please wait...\n";

while (my $line = <$in>) {
    chomp $line;
    if ($line =~ m{^//}) {
        save_feature();
        if (@accs && $length) {
            my $acc_str = join(";", @accs);
            $gene = $accs[0] unless $gene; 
            
            $species =~ s/\s+$//;
            $species =~ s/\.$//;
            $species = "NA" unless $species;

            my $dom_str = @domains ? join(" | ", @domains) : "NA";
            my $var_str = @variants ? join(" | ", @variants) : "NA";
            my $seq_str = $sequence ? $sequence : "NA";
            
            my @dom_seqs;
            if ($sequence) {
                foreach my $d (@domains) {
                    if ($d =~ /^(\d+)\.\.(\d+):/) {
                        my $start = $1;
                        my $end = $2;
                        if ($start > 0 && $end >= $start && $end <= length($sequence)) {
                            my $d_seq = substr($sequence, $start - 1, $end - $start + 1);
                            push @dom_seqs, "$start..$end:$d_seq";
                        }
                    }
                }
            }
            my $dom_seq_str = @dom_seqs ? join(" | ", @dom_seqs) : "NA";

            print $out "$gene\t$acc_str\t$taxid\t$species\t$length\t$dom_str\t$var_str\t$seq_str\t$dom_seq_str\n";
        }
        $gene = $taxid = $species = $length = $sequence = "";
        @accs = @domains = @variants = ();
        $in_sq = 0;
    } 
    # NAYA: ID line se "_" se pehle ka hissa Gene Name ke roop mein nikalna
    elsif ($line =~ /^ID\s+([^_]+)_/) {
        $gene = $1;
    } 
    elsif ($line =~ /^AC\s+(.*)/) {
        my $ac_line = $1;
        push @accs, grep { $_ } map { s/\s+//g; s/;//g; $_ } split(/;/, $ac_line);
    } 
    elsif ($line =~ /^OS\s+(.+)/) {
        $species .= $1 . " ";
    } 
    elsif ($line =~ /^OX\s+NCBI_TaxID=(\d+)/) {
        $taxid = $1;
    } elsif ($line =~ /^FT\s+([A-Z_]+)\s+([\d<>?.]+(?:\.\.[\d<>?.]+)?)/) {
        save_feature();
        $ft_type = $1;
        my $loc = $2;
        $loc =~ s/[<>?]//g; 
        $ft_loc = $loc;
    } elsif ($line =~ /^FT\s+\/note="([^"]+)"/ && $ft_type) {
        $ft_note = $1;
    } elsif ($line =~ /^FT\s{20,}([^\/].+)$/ && $ft_type && $ft_note) {
        my $extra = $1;
        $extra =~ s/^\s+|\s+$//g;
        $extra =~ s/"$//;
        $ft_note .= " " . $extra;
    } elsif ($line =~ /^SQ\s+SEQUENCE\s+(\d+)\s+AA/) {
        $length = $1;
        $in_sq = 1; 
    } elsif ($in_sq && $line !~ /^[A-Z]{2}\s/) {
        my $seq_part = $line;
        $seq_part =~ s/\s+//g;
        $sequence .= $seq_part;
    }
}
close($in); close($out);
print "Success: Domain_Graph_Data.tsv updated successfully!\n";