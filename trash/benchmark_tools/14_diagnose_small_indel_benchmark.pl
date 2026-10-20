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
use File::Basename qw(basename);

################################################################################
# diagnose_small_indel_vcf_benchmark.pl
#
# Purpose:
#   Diagnose why small InDel benchmarking between a simulation truth VCF and an
#   SVGAP prediction VCF gives low precision/recall.  This script reports how
#   many predicted small InDels can be matched to truth records under different
#   breakpoint tolerances, and summarizes the offset distribution.  It is useful
#   for distinguishing true genotyping/calling errors from VCF coordinate or
#   representation differences.
#
# Matching logic:
#   - Only INS/DEL with 1 <= abs(SVLEN or allele length difference) <= --max_len
#     are used.
#   - By default, records must have the same CHROM, same SVTYPE, and the same
#     absolute length.  Use --len_tol to allow a small length difference.
#   - Event-level benchmarking is repeated for multiple position tolerances
#     given by --tolerances, e.g. 0,1,2,5,10,20,50.
#   - Genotype-level benchmarking is also reported for each tolerance when sample
#     columns are available.
#
# Author: Yi Liao; diagnostic implementation generated with ChatGPT
################################################################################

my %opt = (
    max_len     => 50,
    min_len     => 1,
    len_tol     => 0,
    tolerances  => '0,1,2,5,10,20,50,100',
    out         => 'small_indel_diagnosis',
    samples     => '',
    require_type=> 1,
    help        => 0,
);

GetOptions(
    'truth=s'       => \$opt{truth},
    'pred=s'        => \$opt{pred},
    'samples=s'     => \$opt{samples},
    'out=s'         => \$opt{out},
    'min_len=i'     => \$opt{min_len},
    'max_len=i'     => \$opt{max_len},
    'len_tol=i'     => \$opt{len_tol},
    'tolerances=s'  => \$opt{tolerances},
    'require_type!' => \$opt{require_type},
    'help|h'        => \$opt{help},
) or die usage();

if ($opt{help} || !$opt{truth} || !$opt{pred}) {
    print usage();
    exit($opt{help} ? 0 : 1);
}

my @tol = sort { $a <=> $b } grep { /^\d+$/ } split(/,/, $opt{tolerances});
die "ERROR: --tolerances must contain at least one non-negative integer\n" unless @tol;

my @sample_order = read_sample_list($opt{samples});
my $truth = read_vcf($opt{truth}, 'TRUTH', \@sample_order, \%opt);
my $pred  = read_vcf($opt{pred},  'PRED',  \@sample_order, \%opt);

my @truth_records = @{ $truth->{records} };
my @pred_records  = @{ $pred->{records} };

open(my $SUM, '>', "$opt{out}.summary.tsv") or die "ERROR: cannot write $opt{out}.summary.tsv: $!\n";
print $SUM join("\t", qw(TOLERANCE LEVEL TRUTH_N PRED_N TP FP FN TN PRECISION RECALL F1)), "\n";

open(my $OFF, '>', "$opt{out}.nearest_offset.tsv") or die "ERROR: cannot write $opt{out}.nearest_offset.tsv: $!\n";
print $OFF join("\t", qw(PRED_ID PRED_CHR PRED_POS PRED_TYPE PRED_LEN NEAREST_TRUTH_ID NEAREST_TRUTH_POS OFFSET SAME_TYPE SAME_LEN ALLELE_EXACT TRUTH_GT_CARRIERS PRED_GT_CARRIERS)), "\n";
write_nearest_offsets($OFF, \@truth_records, \@pred_records, \%opt);
close $OFF;

open(my $MATCH, '>', "$opt{out}.matches_by_tolerance.tsv") or die "ERROR: cannot write $opt{out}.matches_by_tolerance.tsv: $!\n";
print $MATCH join("\t", qw(TOLERANCE TRUTH_ID PRED_ID CHROM TRUTH_POS PRED_POS OFFSET TYPE LEN ALLELE_EXACT TRUTH_CARRIERS PRED_CARRIERS)), "\n";

for my $tolerance (@tol) {
    my ($pairs, $event_stats) = greedy_match(\@truth_records, \@pred_records, $tolerance, \%opt);
    print $SUM join("\t", $tolerance, 'EVENT', scalar(@truth_records), scalar(@pred_records),
        $event_stats->{TP}, $event_stats->{FP}, $event_stats->{FN}, 'NA',
        fmt($event_stats->{precision}), fmt($event_stats->{recall}), fmt($event_stats->{f1})), "\n";

    my $gt_stats = genotype_stats($pairs, \@truth_records, \@pred_records, \@sample_order);
    print $SUM join("\t", $tolerance, 'GENOTYPE', scalar(@truth_records), scalar(@pred_records),
        $gt_stats->{TP}, $gt_stats->{FP}, $gt_stats->{FN}, $gt_stats->{TN},
        fmt($gt_stats->{precision}), fmt($gt_stats->{recall}), fmt($gt_stats->{f1})), "\n";

    for my $pair (@$pairs) {
        my ($ti, $pi) = @$pair;
        my $tr = $truth_records[$ti];
        my $pr = $pred_records[$pi];
        print $MATCH join("\t", $tolerance, $tr->{id}, $pr->{id}, $tr->{chr}, $tr->{pos}, $pr->{pos},
            $pr->{pos} - $tr->{pos}, $tr->{type}, $tr->{len}, allele_exact($tr, $pr),
            carrier_string($tr), carrier_string($pr)), "\n";
    }
}
close $MATCH;
close $SUM;

write_offset_histogram("$opt{out}.offset_histogram.tsv", \@truth_records, \@pred_records, \%opt);

print STDERR "Done. Output files:\n";
print STDERR "  $opt{out}.summary.tsv\n";
print STDERR "  $opt{out}.nearest_offset.tsv\n";
print STDERR "  $opt{out}.offset_histogram.tsv\n";
print STDERR "  $opt{out}.matches_by_tolerance.tsv\n";

sub usage {
    return <<'USAGE';
Usage:
  perl diagnose_small_indel_vcf_benchmark.pl \
    --truth SIM.population.truth.vcf \
    --pred All.smallIndel.vcf \
    --samples query.list \
    --out smallIndel.diagnosis

Options:
  --truth FILE          simulation truth VCF containing SNP, small InDel and SVs
  --pred FILE           SVGAP small InDel VCF
  --samples FILE        optional sample/query list; one sample per line
  --out PREFIX          output prefix [small_indel_diagnosis]
  --min_len INT         minimum InDel length [1]
  --max_len INT         maximum InDel length [50]
  --len_tol INT         allowed length difference when matching [0]
  --tolerances LIST     comma-separated breakpoint tolerances [0,1,2,5,10,20,50,100]
  --no-require_type     allow INS/DEL type mismatch when diagnosing

Interpretation:
  If F1 is low at tolerance=0 but high at tolerance=5/10/50, the main problem is
  VCF coordinate/normalization rather than failure to detect the InDels.
USAGE
}

sub read_sample_list {
    my ($file) = @_;
    return () unless defined $file && $file ne '';
    open(my $IN, '<', $file) or die "ERROR: cannot open sample list $file: $!\n";
    my (@s, %seen);
    while (my $line = <$IN>) {
        chomp $line;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '' || $line =~ /^#/;
        my ($x) = split(/\s+/, $line);
        next if $seen{$x}++;
        push @s, $x;
    }
    close $IN;
    return @s;
}

sub read_vcf {
    my ($file, $label, $sample_order_ref, $opt) = @_;
    open(my $IN, '<', $file) or die "ERROR: cannot open VCF $file: $!\n";
    my (@header_samples, %sample_col, @records);
    my %wanted = map { $_ => 1 } @$sample_order_ref;
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^##/;
        if ($line =~ /^#CHROM/) {
            my @h = split(/\t/, $line);
            @header_samples = @h[9..$#h] if @h > 9;
            for my $i (0 .. $#header_samples) {
                $sample_col{$header_samples[$i]} = 9 + $i;
            }
            if (!@$sample_order_ref) {
                @$sample_order_ref = @header_samples;
            }
            next;
        }
        next if $line =~ /^#/ || $line =~ /^\s*$/;
        my @f = split(/\t/, $line);
        next unless @f >= 8;
        my ($chr, $pos, $id, $ref, $alt, undef, undef, $info_str) = @f[0..7];
        next unless defined $pos && $pos =~ /^\d+$/;
        my %info = parse_info($info_str);
        my @alts = split(/,/, $alt);
        for my $alt_one (@alts) {
            my ($type, $len, $end) = infer_indel($pos, $ref, $alt_one, \%info);
            next unless $type eq 'INS' || $type eq 'DEL';
            next unless $len >= $opt->{min_len} && $len <= $opt->{max_len};
            my ($npos, $nref, $nalt) = trim_normalize($pos, $ref, $alt_one);
            my (%gt, %carrier);
            if (@$sample_order_ref && @f > 9) {
                my @fmt = split(/:/, $f[8] || 'GT');
                my %fmt_idx = map { $fmt[$_] => $_ } 0..$#fmt;
                my $gt_idx = exists $fmt_idx{GT} ? $fmt_idx{GT} : 0;
                for my $s (@$sample_order_ref) {
                    my $col = exists $sample_col{$s} ? $sample_col{$s} : undef;
                    my $field = defined($col) && defined($f[$col]) ? $f[$col] : './.';
                    my @v = split(/:/, $field);
                    my $g = defined $v[$gt_idx] ? $v[$gt_idx] : './.';
                    $gt{$s} = $g;
                    $carrier{$s} = gt_is_carrier($g);
                }
            }
            push @records, {
                label => $label,
                chr => $chr,
                pos => int($pos),
                id => ($id && $id ne '.') ? $id : "$label:$chr:$pos:$type:" . scalar(@records),
                ref => uc($ref),
                alt => uc($alt_one),
                npos => int($npos),
                nref => uc($nref),
                nalt => uc($nalt),
                type => $type,
                len => int($len),
                end => int($end),
                gt => \%gt,
                carrier => \%carrier,
                raw => $line,
            };
        }
    }
    close $IN;
    return { samples => \@header_samples, records => \@records };
}

sub parse_info {
    my ($s) = @_;
    my %h;
    for my $p (split(/;/, $s || '')) {
        next if $p eq '';
        my ($k, $v) = split(/=/, $p, 2);
        $h{$k} = defined $v ? $v : 1;
    }
    return %h;
}

sub infer_indel {
    my ($pos, $ref, $alt, $info) = @_;
    my $svtype = uc($info->{SVTYPE} || '');
    my $svlen  = defined $info->{SVLEN} ? $info->{SVLEN} : undef;
    my $end    = defined $info->{END} ? $info->{END} : $pos;
    if ($svtype eq 'INS' || $svtype eq 'DEL') {
        my $len = defined $svlen && $svlen =~ /^-?\d+/ ? abs(int($svlen)) : allele_indel_len($ref, $alt);
        return ($svtype, $len, $end);
    }
    my $d = length($alt || '') - length($ref || '');
    return ('INS', abs($d), $pos) if $d > 0;
    return ('DEL', abs($d), $pos + length($ref || '') - 1) if $d < 0;
    return ('OTHER', 0, $pos);
}

sub allele_indel_len {
    my ($ref, $alt) = @_;
    return 0 if !defined $ref || !defined $alt || $alt =~ /^<[^>]+>$/;
    my ($p, $r, $a) = trim_normalize(1, $ref, $alt);
    return abs(length($a) - length($r));
}

sub trim_normalize {
    my ($pos, $ref, $alt) = @_;
    $ref = uc($ref // '');
    $alt = uc($alt // '');
    return ($pos, $ref, $alt) if $alt =~ /^<[^>]+>$/;
    while (length($ref) > 1 && length($alt) > 1 && substr($ref, -1) eq substr($alt, -1)) {
        chop $ref; chop $alt;
    }
    while (length($ref) > 1 && length($alt) > 1 && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
        $ref = substr($ref, 1);
        $alt = substr($alt, 1);
        $pos++;
    }
    return ($pos, $ref, $alt);
}

sub gt_is_carrier {
    my ($gt) = @_;
    return undef if !defined $gt || $gt eq '' || $gt eq '.' || $gt =~ /^\.([\/|]\.)?$/;
    return 0 if $gt =~ /^0([\/|]0)*$/;
    return 1 if $gt =~ /[1-9]/;
    return undef;
}

sub greedy_match {
    my ($truth, $pred, $tol, $opt) = @_;
    my @candidates;
    for my $ti (0 .. $#$truth) {
        for my $pi (0 .. $#$pred) {
            my $t = $truth->[$ti];
            my $p = $pred->[$pi];
            next unless $t->{chr} eq $p->{chr};
            next if $opt->{require_type} && $t->{type} ne $p->{type};
            next if abs($t->{len} - $p->{len}) > $opt->{len_tol};
            my $dist = abs($t->{pos} - $p->{pos});
            next if $dist > $tol;
            my $allele_bonus = allele_exact($t, $p) ? 0 : 0.1;
            push @candidates, [$dist + $allele_bonus, $ti, $pi];
        }
    }
    @candidates = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @candidates;
    my (%used_t, %used_p, @pairs);
    for my $c (@candidates) {
        my ($score, $ti, $pi) = @$c;
        next if $used_t{$ti} || $used_p{$pi};
        $used_t{$ti} = 1;
        $used_p{$pi} = 1;
        push @pairs, [$ti, $pi];
    }
    my $tp = scalar @pairs;
    my $fp = scalar(@$pred) - $tp;
    my $fn = scalar(@$truth) - $tp;
    return (\@pairs, calc_stats($tp, $fp, $fn));
}

sub genotype_stats {
    my ($pairs, $truth, $pred, $samples) = @_;
    my ($tp, $fp, $fn, $tn) = (0, 0, 0, 0);
    return { TP=>0, FP=>0, FN=>0, TN=>'NA', precision=>'NA', recall=>'NA', f1=>'NA' } unless @$samples;
    for my $pair (@$pairs) {
        my ($ti, $pi) = @$pair;
        my $t = $truth->[$ti];
        my $p = $pred->[$pi];
        for my $s (@$samples) {
            my $tc = $t->{carrier}{$s};
            my $pc = $p->{carrier}{$s};
            next unless defined $tc && defined $pc;
            if ($tc && $pc) { $tp++; }
            elsif (!$tc && $pc) { $fp++; }
            elsif ($tc && !$pc) { $fn++; }
            else { $tn++; }
        }
    }
    my $st = calc_stats($tp, $fp, $fn);
    $st->{TN} = $tn;
    return $st;
}

sub calc_stats {
    my ($tp, $fp, $fn) = @_;
    my $precision = ($tp + $fp) ? $tp / ($tp + $fp) : 'NA';
    my $recall    = ($tp + $fn) ? $tp / ($tp + $fn) : 'NA';
    my $f1 = 'NA';
    if ($precision ne 'NA' && $recall ne 'NA' && ($precision + $recall) > 0) {
        $f1 = 2 * $precision * $recall / ($precision + $recall);
    }
    return { TP=>$tp, FP=>$fp, FN=>$fn, precision=>$precision, recall=>$recall, f1=>$f1 };
}

sub write_nearest_offsets {
    my ($OUT, $truth, $pred, $opt) = @_;
    for my $p (@$pred) {
        my ($best, $best_score);
        for my $t (@$truth) {
            next unless $t->{chr} eq $p->{chr};
            my $same_type = $t->{type} eq $p->{type} ? 1 : 0;
            my $same_len  = abs($t->{len} - $p->{len}) <= $opt->{len_tol} ? 1 : 0;
            my $score = abs($t->{pos} - $p->{pos}) + ($same_type ? 0 : 1000000) + ($same_len ? 0 : 100000);
            if (!defined $best_score || $score < $best_score) {
                $best_score = $score;
                $best = $t;
            }
        }
        if ($best) {
            print $OUT join("\t", $p->{id}, $p->{chr}, $p->{pos}, $p->{type}, $p->{len},
                $best->{id}, $best->{pos}, $p->{pos} - $best->{pos},
                ($p->{type} eq $best->{type} ? 1 : 0),
                (abs($p->{len} - $best->{len}) <= $opt->{len_tol} ? 1 : 0),
                allele_exact($best, $p), carrier_string($best), carrier_string($p)), "\n";
        } else {
            print $OUT join("\t", $p->{id}, $p->{chr}, $p->{pos}, $p->{type}, $p->{len}, '.', '.', '.', '.', '.', '.', '.', carrier_string($p)), "\n";
        }
    }
}

sub write_offset_histogram {
    my ($file, $truth, $pred, $opt) = @_;
    my %hist;
    for my $p (@$pred) {
        my $best_abs;
        for my $t (@$truth) {
            next unless $t->{chr} eq $p->{chr};
            next if $opt->{require_type} && $t->{type} ne $p->{type};
            next if abs($t->{len} - $p->{len}) > $opt->{len_tol};
            my $d = abs($p->{pos} - $t->{pos});
            $best_abs = $d if !defined($best_abs) || $d < $best_abs;
        }
        my $bin = !defined($best_abs) ? 'no_same_chr_type_len'
                : $best_abs == 0 ? '0'
                : $best_abs <= 1 ? '1'
                : $best_abs <= 2 ? '2'
                : $best_abs <= 5 ? '3-5'
                : $best_abs <= 10 ? '6-10'
                : $best_abs <= 20 ? '11-20'
                : $best_abs <= 50 ? '21-50'
                : $best_abs <= 100 ? '51-100'
                : '>100';
        $hist{$bin}++;
    }
    open(my $OUT, '>', $file) or die "ERROR: cannot write $file: $!\n";
    print $OUT "BIN\tCOUNT\n";
    for my $b (qw(0 1 2 3-5 6-10 11-20 21-50 51-100 >100 no_same_chr_type_len)) {
        print $OUT "$b\t", ($hist{$b} || 0), "\n";
    }
    close $OUT;
}

sub allele_exact {
    my ($a, $b) = @_;
    return 1 if $a->{npos} == $b->{npos} && $a->{nref} eq $b->{nref} && $a->{nalt} eq $b->{nalt};
    return 0;
}

sub carrier_string {
    my ($r) = @_;
    my @c = sort grep { defined $r->{carrier}{$_} && $r->{carrier}{$_} } keys %{ $r->{carrier} || {} };
    return @c ? join(',', @c) : '.';
}

sub fmt {
    my ($x) = @_;
    return 'NA' if !defined $x || $x eq 'NA';
    return sprintf('%.6f', $x);
}
