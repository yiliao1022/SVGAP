#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/pub";
use lib "$Bin/../lib";
use lib "$Bin/../pub";
use lib "$Bin/../../lib";
use lib "$Bin/../../pub";
use Getopt::Long qw(GetOptions);

# ============================================================================
# SVGAP standardized script: SVGAP_vcf_tools.pl
# Original author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Standardized date: 2026-06-12
#
# Utility functions for standardized SVGAP VCF outputs.
# Current mode:
#   merge: concatenate multiple VCFs and force sample columns to follow g.list
# ============================================================================

my %opt = (
    mode => 'merge',
    vcfs => '',
    out  => '',
    query_list => '',
    glist => '',
    refname => '',
    refsize => '',
    missing_gt => './.',
    help => 0,
);

GetOptions(
    'mode=s' => \$opt{mode},
    'vcfs=s' => \$opt{vcfs},
    'out=s' => \$opt{out},
    'query_list|query-list=s' => \$opt{query_list},
    'glist=s' => \$opt{glist},
    'refname|target=s' => \$opt{refname},
    'refsize=s' => \$opt{refsize},
    'missing_gt=s' => \$opt{missing_gt},
    'help|h' => \$opt{help},
) or die usage();

if ($opt{help} || $opt{mode} ne 'merge' || !$opt{vcfs} || !$opt{out}) {
    print usage();
    exit($opt{help} ? 0 : 1);
}

my @samples = read_sample_order($opt{query_list} || $opt{glist}, $opt{refname});
die "ERROR: --query_list or --glist is required for deterministic sample order\n" unless @samples;
my %want_sample = map { $_ => 1 } @samples;
my @vcfs = split(/,/, $opt{vcfs});
@vcfs = grep { defined $_ && $_ ne '' && -s $_ } @vcfs;
die "ERROR: no readable VCFs in --vcfs\n" unless @vcfs;
my %chr_rank = read_chr_rank($opt{refsize});

my (@meta, @records);
my %seen_meta;
for my $vcf (@vcfs) {
    open(my $IN, '<', $vcf) or die "ERROR: cannot read $vcf: $!\n";
    my @file_samples;
    my %sample_col;
    while (my $line = <$IN>) {
        chomp $line;
        next if $line eq '';
        if ($line =~ /^##/) {
            if (!$seen_meta{$line}++ && $line !~ /^##fileDate=/ && $line !~ /^##source=/) {
                push @meta, $line;
            }
            next;
        }
        if ($line =~ /^#CHROM\t/) {
            my @h = split(/\t/, $line);
            @file_samples = @h[9..$#h] if @h > 9;
            for (my $i=0; $i<@file_samples; $i++) {
                $sample_col{$file_samples[$i]} = 9 + $i;
            }
            next;
        }
        my @f = split(/\t/, $line, -1);
        next unless @f >= 8;
        my @prefix = @f[0..8];
        my @gt;
        for my $s (@samples) {
            if (exists $sample_col{$s} && defined $f[$sample_col{$s}]) {
                push @gt, $f[$sample_col{$s}];
            } else {
                push @gt, normalize_missing_gt($opt{missing_gt}, $prefix[8]);
            }
        }
        push @records, { chr=>$f[0], pos=>($f[1] =~ /^\d+$/ ? $f[1] : 0), line=>join("\t", @prefix, @gt) };
    }
    close $IN;
}

@records = sort {
    (exists $chr_rank{$a->{chr}} ? $chr_rank{$a->{chr}} : 1e12) <=> (exists $chr_rank{$b->{chr}} ? $chr_rank{$b->{chr}} : 1e12)
    || $a->{chr} cmp $b->{chr}
    || $a->{pos} <=> $b->{pos}
    || $a->{line} cmp $b->{line}
} @records;

open(my $OUT, '>', $opt{out}) or die "ERROR: cannot write $opt{out}: $!\n";
print $OUT "##fileformat=VCFv4.2\n";
print $OUT "##source=SVGAP_vcf_tools.pl merge\n";
for my $m (@meta) {
    next if $m eq '##fileformat=VCFv4.2';
    print $OUT "$m\n";
}
print $OUT join("\t", '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',@samples), "\n";
for my $r (@records) {
    print $OUT $r->{line}, "\n";
}
close $OUT;
print STDERR "Merged ".scalar(@records)." records into $opt{out}\n";

sub normalize_missing_gt {
    my ($gt, $format) = @_;
    return $gt if !defined $format || $format eq 'GT';
    my @keys = split(/:/, $format);
    return $gt if @keys <= 1;
    my @vals = ($gt, ('.') x (@keys-1));
    return join(':', @vals);
}

sub read_sample_order {
    my ($file, $refname) = @_;
    return () unless defined $file && $file ne '';
    open(my $IN, '<', $file) or die "ERROR: cannot read sample list $file: $!\n";
    my (@samples, %seen);
    while (<$IN>) {
        chomp;
        s/^\s+|\s+$//g;
        next if $_ eq '' || /^#/;
        next if defined($refname) && $refname ne '' && $_ eq $refname;
        die "ERROR: duplicated sample in $file: $_\n" if $seen{$_}++;
        push @samples, $_;
    }
    close $IN;
    return @samples;
}

sub read_chr_rank {
    my ($file) = @_;
    my %rank;
    return %rank unless defined $file && $file ne '' && -s $file;
    open(my $IN, '<', $file) or die "ERROR: cannot read --refsize $file: $!\n";
    my $i = 0;
    while (<$IN>) {
        chomp;
        next if /^\s*$/;
        my ($chr) = split(/\s+/);
        $rank{$chr} = ++$i if defined $chr && $chr ne '';
    }
    close $IN;
    return %rank;
}

sub usage {
    return <<'USAGE';
Usage:
  perl SVGAP_vcf_tools.pl --mode merge \
    --vcfs DEL.vcf,INS.vcf,INV.vcf \
    --query_list query.list \
    --refname Ref \
    --refsize Ref.sizes \
    --out All.SV.raw.vcf

Purpose:
  Merge several SVGAP VCF files while forcing sample columns to follow
  --query_list/--glist. Missing samples in a component VCF are filled as ./.
USAGE
}
