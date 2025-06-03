#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $file = $ARGV[0];
die "Usage: perl $0 input.fasta\n" unless $file;

my $prefix = basename($file);
$prefix =~ s/\.\w+$//;

open(my $IN,  "<", $file) or die "Cannot open $file: $!";
open(my $OUT, ">", "$prefix") or die "Cannot write fasta: $!";
open(my $LOG, ">", "$prefix.skip.log") or die "Cannot write log: $!";
open(my $MAP, ">", "$prefix.idmap.tsv") or die "Cannot write idmap: $!";

my $print_seq = 0;
my $new_id = "";

while (<$IN>) {
    chomp;
    if (/^>(.+)/) {
        my $line = $_;
        my $header = $1;

        if ($header =~ /chloroplast|mitochondrion|scaffold|plasmid|uncluster|unplaced/i) {
            print $LOG "Skipped non-chromosome header: $line\n";
            $print_seq = 0;
            next;
        }

        my $chr_id = "";

        if ($header =~ /\bChr([0-9]+[A-Za-z]*)\b/i) {
            $chr_id = "Chr$1";
        }
        elsif ($header =~ /\b(?:chr(?:omosome)?)0*([0-9]+[A-Za-z]*)\b/i) {
            my $val = $1;
            $chr_id = ($val =~ /^\d+$/) ? sprintf("Chr%02d", $val) : "Chr$val";
        }
        elsif ($header =~ /\b(?:chr(?:omosome)?)([XYM])\b/i) {
            $chr_id = "Chr$1";
        }
        elsif ($header =~ /\b([0-9]+[A-Za-z]*|[XYM])\b/i) {
            my $val = $1;
            $chr_id = ($val =~ /^\d+$/) ? sprintf("Chr%02d", $val) : "Chr$val";
        }

        if ($chr_id) {
            $new_id = "$prefix.$chr_id";
            print $OUT ">$new_id\n";
            print $MAP "$line\t$chr_id\n";
            $print_seq = 1;
        } else {
            print $LOG "Skipped unrecognized header: $line\n";
            $print_seq = 0;
        }
    } else {
        print $OUT "$_\n" if $print_seq;
    }
}

close $IN;
close $OUT;
close $LOG;
close $MAP;

print "✅ Done. Output:\n";
print " - Rename fasta: $prefix\n";
print " - Log of skipped headers: $prefix.skip.log\n";

