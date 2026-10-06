#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl Extract_Domain_Data.pl uniprot_sprot.dat
my $dat_file = shift @ARGV or die "Usage: perl Extract_Domain_Data.pl uniprot_sprot.dat\n";
my $out_file = "Domain_Graph_Data.tsv";

open(my $in, '<', $dat_file) or die "Cannot open $dat_file: $!\n";
open(my $out, '>', $out_file) or die "Cannot open $out_file: $!\n";

# TSV Header mein ek naya column 'Sequence' add kiya gaya hai
print $out "Gene_Name\tAccessions\tTaxID\tLength\tDomains\tVariants\tSequence\n";

my ($gene, $taxid, $length, $sequence) = ("", "", "", "");
my @accs;
my (@domains, @variants);
my ($ft_type, $ft_loc, $ft_note) = ("", "", "");
my $in_sq = 0; # Sequence block track karne ke liye

sub save_feature {
    if ($ft_type eq 'DOMAIN' && $ft_note) {
        push @domains, "$ft_loc:$ft_note";
    } elsif (($ft_type eq 'VARIANT' || $ft_type eq 'VAR_SEQ') && $ft_note) {
        push @variants, "$ft_loc:$ft_note";
    }
    $ft_type = ""; $ft_loc = ""; $ft_note = "";
}

print "Parsing $dat_file (4GB) for Domains and Sequences. Please wait...\n";

while (my $line = <$in>) {
    chomp $line;
    if ($line =~ m{^//}) {
        save_feature();
        if (@accs && $length) {
            my $acc_str = join(";", @accs);
            $gene = $accs[0] unless $gene;
            my $dom_str = @domains ? join(" | ", @domains) : "NA";
            my $var_str = @variants ? join(" | ", @variants) : "NA";
            my $seq_str = $sequence ? $sequence : "NA";
            
            print $out "$gene\t$acc_str\t$taxid\t$length\t$dom_str\t$var_str\t$seq_str\n";
        }
        # Reset variables
        $gene = $taxid = $length = $sequence = "";
        @accs = @domains = @variants = ();
        $in_sq = 0;
    } elsif ($line =~ /^AC\s+(.*)/) {
        my $ac_line = $1;
        push @accs, grep { $_ } map { s/\s+//g; s/;//g; $_ } split(/;/, $ac_line);
    } elsif ($line =~ /^GN\s+Name=([^;{\s]+)/) {
        $gene = $1 unless $gene;
    } elsif ($line =~ /^OX\s+NCBI_TaxID=(\d+)/) {
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
        $in_sq = 1; # Sequence block start
    } elsif ($in_sq && $line !~ /^[A-Z]{2}\s/) {
        # Sequence block ke andar ki spaces hata kar add karna
        my $seq_part = $line;
        $seq_part =~ s/\s+//g;
        $sequence .= $seq_part;
    }
}
close($in); close($out);
print "Success: Domain_Graph_Data.tsv with Sequences has been created!\n";