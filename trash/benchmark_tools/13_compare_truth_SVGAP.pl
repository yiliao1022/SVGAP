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
use File::Path qw(make_path);
use POSIX qw(strftime);

################################################################################
# compare_svgap_vcf_to_sim_truth.pl
#
# Purpose:
#   Compare SVGAP-called VCF files against one integrated simulated truth VCF.
#   The truth VCF can contain SNPs, 1-50 bp small InDels, and large SVs in one
#   file, whereas SVGAP outputs them as separate VCF files.  This script evaluates
#   each SVGAP VCF separately against the relevant records in the truth VCF.
#
# Author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Date: 2026-06-13
#
# Main outputs:
#   OUT.summary.tsv           event-level and genotype/carrier-level metrics
#   OUT.per_sample.tsv        genotype/carrier-level metrics for each sample
#   OUT.matches.tsv           matched predicted-truth record pairs
#   OUT.unmatched_pred.tsv    false-positive predicted records
#   OUT.unmatched_truth.tsv   false-negative truth records
#
# Definitions:
#   Event-level TP/FP/FN:
#     TP = predicted variant matched to a truth variant.
#     FP = predicted variant not matched to any truth variant.
#     FN = truth variant not matched by any predicted variant.
#
#   Genotype/carrier-level TP/FP/FN:
#     For matched records, compare carrier status per sample.
#     For unmatched predicted records, predicted carriers are FP.
#     For unmatched truth records, truth carriers are FN.
#     TN is reported only for matched records because the full variant universe is
#     not defined for VCF-vs-VCF benchmarking.
################################################################################

my %opt = (
    out               => 'SVGAP_vs_truth',
    samples           => '',
    pred              => [],
    large_vcf         => [],
    small_vcf         => [],
    snp_vcf           => [],

    large_min_len     => 51,     # used only when truth lacks VARTYPE
    small_max_len     => 50,
    sv_ov             => 0.50,
    sv_len_ratio      => 0.50,   # max relative difference; 0.50 allows 1.5-fold
    sv_bp_tol         => 1000,
    ins_tol           => 100,
    ins_len_ratio     => 0.50,
    small_pos_tol     => 20,   # v2 default: small InDels are matched within 20 bp to allow VCF normalization shifts
    snp_pos_tol       => 0,

    require_type      => 1,
    report_multi_alt  => 0,
    keep_missing_gt   => 0,
    quiet             => 0,
    help              => 0,
);

GetOptions(
    'truth=s'            => \$opt{truth},
    'pred=s@'            => \$opt{pred},
    'large_vcf=s@'       => \$opt{large_vcf},
    'sv_vcf=s@'          => \$opt{large_vcf},
    'small_vcf=s@'       => \$opt{small_vcf},
    'small_indel_vcf=s@' => \$opt{small_vcf},
    'snp_vcf=s@'         => \$opt{snp_vcf},
    'out=s'              => \$opt{out},
    'samples=s'          => \$opt{samples},
    'query_list=s'       => \$opt{samples},
    'glist=s'            => \$opt{samples},

    'large_min_len=i'    => \$opt{large_min_len},
    'small_max_len=i'    => \$opt{small_max_len},
    'sv_ov=f'            => \$opt{sv_ov},
    'sv_len_ratio=f'     => \$opt{sv_len_ratio},
    'sv_bp_tol=i'        => \$opt{sv_bp_tol},
    'ins_tol=i'          => \$opt{ins_tol},
    'ins_len_ratio=f'    => \$opt{ins_len_ratio},
    'small_pos_tol=i'    => \$opt{small_pos_tol},
    'snp_pos_tol=i'      => \$opt{snp_pos_tol},
    'require_type!'      => \$opt{require_type},
    'report_multi_alt!'  => \$opt{report_multi_alt},
    'keep_missing_gt!'   => \$opt{keep_missing_gt},
    'quiet!'             => \$opt{quiet},
    'help|h'             => \$opt{help},
) or die usage();

if ($opt{help} || !$opt{truth}) {
    print usage();
    exit($opt{help} ? 0 : 1);
}

my @evals = parse_eval_inputs(\%opt);
die "ERROR: no SVGAP prediction VCF was provided. Use --large_vcf, --small_vcf, --snp_vcf or --pred LABEL:CLASS:FILE.\n" unless @evals;

die "ERROR: truth VCF does not exist: $opt{truth}\n" unless -s $opt{truth};
for my $e (@evals) {
    die "ERROR: prediction VCF does not exist for $e->{label}: $e->{file}\n" unless -s $e->{file};
}

my $outdir = dirname_of_prefix($opt{out});
make_path($outdir) if $outdir ne '' && !-d $outdir;

logmsg("Reading truth VCF: $opt{truth}");
my $truth_vcf = read_vcf($opt{truth}, 'TRUTH', \%opt);
logmsg("Truth records loaded: " . scalar(@{ $truth_vcf->{records} }));

my @sample_order = determine_sample_order($truth_vcf, \@evals, \%opt);
logmsg("Samples used for genotype-level comparison: " . scalar(@sample_order));

open(my $SUM, '>', "$opt{out}.summary.tsv") or die "ERROR: cannot write $opt{out}.summary.tsv: $!\n";
open(my $PS,  '>', "$opt{out}.per_sample.tsv") or die "ERROR: cannot write $opt{out}.per_sample.tsv: $!\n";
open(my $MAT, '>', "$opt{out}.matches.tsv") or die "ERROR: cannot write $opt{out}.matches.tsv: $!\n";
open(my $UP,  '>', "$opt{out}.unmatched_pred.tsv") or die "ERROR: cannot write $opt{out}.unmatched_pred.tsv: $!\n";
open(my $UT,  '>', "$opt{out}.unmatched_truth.tsv") or die "ERROR: cannot write $opt{out}.unmatched_truth.tsv: $!\n";

print $SUM join("\t", qw/LABEL CLASS LEVEL TRUTH_N PRED_N TP FP FN TN PRECISION RECALL F1/), "\n";
print $PS  join("\t", qw/LABEL CLASS SAMPLE TRUTH_CARRIERS PRED_CARRIERS TP FP FN TN PRECISION RECALL F1/), "\n";
print $MAT join("\t", qw/LABEL CLASS MATCH_SCORE PRED_KEY TRUTH_KEY PRED_ID TRUTH_ID PRED_CHR PRED_POS PRED_END PRED_TYPE PRED_LEN TRUTH_CHR TRUTH_POS TRUTH_END TRUTH_TYPE TRUTH_LEN/), "\n";
print $UP  join("\t", qw/LABEL CLASS PRED_KEY PRED_ID CHR POS END TYPE LEN CARRIERS/), "\n";
print $UT  join("\t", qw/LABEL CLASS TRUTH_KEY TRUTH_ID CHR POS END TYPE LEN CARRIERS/), "\n";

for my $eval (@evals) {
    logmsg("Evaluating $eval->{label} [$eval->{class}] : $eval->{file}");
    my $pred_vcf = read_vcf($eval->{file}, 'PRED', \%opt);
    my @truth_records = grep { record_in_class($_, $eval->{class}, \%opt) } @{ $truth_vcf->{records} };
    my @pred_records  = grep { record_in_class($_, $eval->{class}, \%opt, 1) } @{ $pred_vcf->{records} };

    logmsg("  truth records for class $eval->{class}: " . scalar(@truth_records));
    logmsg("  predicted records for class $eval->{class}: " . scalar(@pred_records));

    my $res = compare_records(\@truth_records, \@pred_records, \@sample_order, $eval, \%opt);
    write_result($SUM, $PS, $MAT, $UP, $UT, $eval, $res, \@sample_order);
}

close $SUM; close $PS; close $MAT; close $UP; close $UT;

logmsg("Done.");
logmsg("Summary:          $opt{out}.summary.tsv");
logmsg("Per-sample:       $opt{out}.per_sample.tsv");
logmsg("Matches:          $opt{out}.matches.tsv");
logmsg("Unmatched pred:   $opt{out}.unmatched_pred.tsv");
logmsg("Unmatched truth:  $opt{out}.unmatched_truth.tsv");
exit 0;

################################################################################
# Input and option handling
################################################################################
sub usage {
    return <<'USAGE';
Usage:
  perl compare_svgap_vcf_to_sim_truth.pl \
    --truth SIM.population.truth.vcf \
    --large_vcf SVGAP.largeSV.vcf \
    --small_vcf SVGAP.smallInDel.vcf \
    --snp_vcf SVGAP.SNP.vcf \
    --samples query.list \
    --out benchmark/SVGAP_vs_truth

Alternative generic input form:
  --pred LABEL:CLASS:FILE

where CLASS is one of:
  LARGE      large SVs: DEL/INS/INV/CNV/DUP/TRL larger than small InDel range
  SMALL      small insertions/deletions, usually 1-50 bp
  SNP        SNPs
  ALL        all records in the VCF
  DEL/INS/INV/CNV/DUP/TRL  one specific SVTYPE

Examples:
  perl compare_svgap_vcf_to_sim_truth.pl \
    --truth SVsim_out/SIM.population.truth.vcf \
    --pred LargeSV:LARGE:SVGAP_run/08_vcf/All.largeSV.raw.vcf \
    --pred SmallInDel:SMALL:SVGAP_run/08_vcf/All.smallIndel.vcf \
    --pred SNP:SNP:SVGAP_run/08_vcf/All.sorted.SNP.vcf \
    --samples SVsim_out/query.list \
    --out SVGAP_benchmark

Recommended for your current SVGAP benchmark:
  1. Use the integrated simulated truth VCF as --truth.
  2. Compare SVGAP large SV, small InDel and SNP VCFs separately.
  3. Use --samples or --glist to force the sample order. If the list contains
     the reference genome name, it will be ignored unless that sample exists in VCF.

Matching defaults:
  SNP:          same CHROM, POS, REF and ALT after simple allele trimming.
  Small InDel: same type, near POS with --small_pos_tol; exact normalized allele gets priority, otherwise same length within tolerance can match
               when allele sequences are available.
  Large DEL/INV/CNV/DUP/TRL:
               same type, same chromosome, reciprocal overlap >= --sv_ov,
               length ratio within --sv_len_ratio, or both breakpoints within
               --sv_bp_tol.
  Large INS:   same type, same chromosome, POS distance <= --ins_tol and length
               ratio within --ins_len_ratio.

Options:
  --truth FILE             integrated simulated truth VCF [required]
  --large_vcf FILE         SVGAP large SV VCF; can be repeated
  --small_vcf FILE         SVGAP small InDel VCF; can be repeated
  --snp_vcf FILE           SVGAP SNP VCF; can be repeated
  --pred LABEL:CLASS:FILE  generic prediction VCF input; can be repeated
  --samples FILE           sample/query list, one name per line; fixes order
  --glist FILE             alias of --samples
  --out PREFIX             output prefix [SVGAP_vs_truth]

Thresholds:
  --small_max_len 50       truth records with abs(SVLEN) <= this are small InDels
  --large_min_len 51       truth records with abs(SVLEN) >= this are large if no VARTYPE
  --sv_ov 0.50             reciprocal overlap threshold for interval SVs
  --sv_len_ratio 0.50      max relative length difference for interval SVs
  --sv_bp_tol 1000         breakpoint tolerance for interval SVs
  --ins_tol 100            breakpoint tolerance for insertions
  --ins_len_ratio 0.50     max relative length difference for insertions
  --small_pos_tol 20       position tolerance for small InDels [v2 default; use 0 for strict exact-coordinate benchmark]
  --snp_pos_tol 0          position tolerance for SNPs, normally 0
  --no-require_type        allow matching without requiring the same SVTYPE
  --report_multi_alt       keep all ALT alleles; default uses first ALT allele only

Output metrics:
  Event-level:
    TP = matched predicted records
    FP = unmatched predicted records
    FN = unmatched truth records
  Genotype/carrier-level:
    TP/FP/FN are counted across samples using non-reference GT as carrier.
    TN is reported only across matched variant records.

USAGE
}

sub parse_eval_inputs {
    my ($opt) = @_;
    my @evals;

    for my $f (@{ $opt->{large_vcf} || [] }) {
        push @evals, { label => clean_label(basename($f)), class => 'LARGE', file => $f };
    }
    for my $f (@{ $opt->{small_vcf} || [] }) {
        push @evals, { label => clean_label(basename($f)), class => 'SMALL', file => $f };
    }
    for my $f (@{ $opt->{snp_vcf} || [] }) {
        push @evals, { label => clean_label(basename($f)), class => 'SNP', file => $f };
    }

    for my $p (@{ $opt->{pred} || [] }) {
        my ($label, $class, $file);
        if ($p =~ /^([^:]+):([^:]+):(.+)$/) {
            ($label, $class, $file) = ($1, $2, $3);
        } elsif ($p =~ /^([^=]+)=([^=]+)=(.+)$/) {
            ($label, $class, $file) = ($1, $2, $3);
        } elsif ($p =~ /^([^=]+)=(.+)$/) {
            ($label, $file) = ($1, $2);
            $class = infer_class_from_label($label);
        } else {
            $file = $p;
            $label = clean_label(basename($file));
            $class = infer_class_from_label($label);
        }
        $class = uc($class);
        die "ERROR: unknown class '$class' in --pred $p\n" unless $class =~ /^(LARGE|SMALL|SNP|ALL|DEL|INS|INV|CNV|DUP|TRL|TRA|BND)$/;
        push @evals, { label => $label, class => $class, file => $file };
    }
    return @evals;
}

sub clean_label {
    my ($x) = @_;
    $x =~ s/\.vcf(?:\.gz)?$//;
    $x =~ s/[^A-Za-z0-9_.+-]+/_/g;
    return $x;
}

sub infer_class_from_label {
    my ($label) = @_;
    return 'SNP'   if $label =~ /snp/i;
    return 'SMALL' if $label =~ /small|indel/i;
    return 'LARGE';
}

sub dirname_of_prefix {
    my ($p) = @_;
    return '' unless $p =~ m{/};
    $p =~ s{/[^/]+$}{};
    return $p;
}

sub logmsg {
    my ($msg) = @_;
    return if $opt{quiet};
    my $t = strftime('%Y-%m-%d %H:%M:%S', localtime);
    print STDERR "[$t] $msg\n";
}

################################################################################
# VCF parsing
################################################################################
sub read_vcf {
    my ($file, $source, $opt) = @_;
    my $fh = open_vcf($file);
    my @samples;
    my @records;
    my $line_no = 0;

    while (my $line = <$fh>) {
        $line_no++;
        chomp $line;
        next if $line =~ /^##/;
        if ($line =~ /^#CHROM\t/) {
            my @h = split(/\t/, $line);
            @samples = @h[9 .. $#h] if @h > 9;
            next;
        }
        next if $line =~ /^#/ || $line =~ /^\s*$/;
        my @f = split(/\t/, $line);
        next unless @f >= 8;
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info_str) = @f[0..7];
        next unless defined $pos && $pos =~ /^\d+$/;

        my @alts = split(/,/, $alt);
        @alts = ($alts[0]) unless $opt->{report_multi_alt};
        for my $alt_i (0 .. $#alts) {
            my $alt_one = $alts[$alt_i];
            my %info = parse_info($info_str);
            my ($type, $vartype, $svlen, $end) = infer_variant_properties($chr, $pos, $ref, $alt_one, \%info);
            my %carrier;
            my %gt_raw;
            if (@samples && @f > 9) {
                my @fmt = split(/:/, $f[8] || '');
                my %fmt_idx;
                for my $i (0 .. $#fmt) { $fmt_idx{$fmt[$i]} = $i; }
                my $gt_idx = exists $fmt_idx{GT} ? $fmt_idx{GT} : 0;
                for my $i (0 .. $#samples) {
                    my $sample = $samples[$i];
                    my $field = $f[9 + $i];
                    my @val = split(/:/, defined($field) ? $field : '.');
                    my $gt = defined $val[$gt_idx] ? $val[$gt_idx] : '.';
                    $gt_raw{$sample} = $gt;
                    $carrier{$sample} = gt_is_carrier($gt);
                }
            }
            my ($n_pos, $n_ref, $n_alt) = normalize_allele($pos, $ref, $alt_one);
            my $key = join(':', $chr, $n_pos, $n_ref, $n_alt);
            my $rec = {
                source   => $source,
                file     => $file,
                line_no  => $line_no,
                chr      => $chr,
                pos      => int($pos),
                end      => int($end),
                id       => $id && $id ne '.' ? $id : "$source:$chr:$pos:$type:$line_no" . ($alt_i ? ":ALT$alt_i" : ''),
                ref      => uc($ref),
                alt      => uc($alt_one),
                npos     => int($n_pos),
                nref     => uc($n_ref),
                nalt     => uc($n_alt),
                allele_key => $key,
                svtype   => $type,
                vartype  => $vartype,
                len      => abs($svlen || 0),
                signed_len => $svlen || 0,
                info     => \%info,
                carriers => \%carrier,
                gt_raw   => \%gt_raw,
                raw      => $line,
            };
            push @records, $rec;
        }
    }
    close $fh;
    return { file => $file, samples => \@samples, records => \@records };
}

sub open_vcf {
    my ($file) = @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gzip -dc '$file' |") or die "ERROR: cannot open gzipped VCF $file: $!\n";
    } else {
        open($fh, '<', $file) or die "ERROR: cannot open VCF $file: $!\n";
    }
    return $fh;
}

sub parse_info {
    my ($s) = @_;
    my %h;
    for my $item (split(/;/, $s || '')) {
        next if $item eq '';
        if ($item =~ /^([^=]+)=(.*)$/) {
            $h{$1} = $2;
        } else {
            $h{$item} = 1;
        }
    }
    return %h;
}

sub infer_variant_properties {
    my ($chr, $pos, $ref, $alt, $info) = @_;
    my $svtype  = $info->{SVTYPE}  // '';
    my $vartype = $info->{VARTYPE} // '';
    $svtype  = uc($svtype);
    $vartype = uc($vartype);

    if (!$svtype || $svtype eq '.') {
        if ($alt =~ /^<([^>]+)>$/) {
            $svtype = uc($1);
        } elsif (length($ref) == 1 && length($alt) == 1) {
            $svtype = 'SNP';
        } elsif (length($ref) > length($alt)) {
            $svtype = 'DEL';
        } elsif (length($alt) > length($ref)) {
            $svtype = 'INS';
        } else {
            $svtype = 'SUB';
        }
    }
    $svtype = 'CNV' if $svtype =~ /^(DUP|TD|TANDEM_DUP|TANDEM_DUPLICATION)$/;
    $svtype = 'TRL' if $svtype =~ /^(TRA|BND|TRANSLOCATION)$/;

    my $svlen = 0;
    if (defined $info->{SVLEN} && $info->{SVLEN} ne '.') {
        my ($v) = split(/,/, $info->{SVLEN});
        $svlen = $v + 0 if defined $v && $v =~ /^-?\d+(?:\.\d+)?$/;
    }
    if (!$svlen) {
        if ($svtype eq 'SNP') {
            $svlen = 1;
        } elsif ($svtype eq 'DEL') {
            $svlen = -(length($ref) - length($alt));
        } elsif ($svtype eq 'INS') {
            $svlen = length($alt) - length($ref);
        } elsif (defined $info->{END} && $info->{END} =~ /^\d+$/) {
            $svlen = $info->{END} - $pos + 1;
        } else {
            $svlen = abs(length($alt) - length($ref));
        }
    }

    my $end;
    if (defined $info->{END} && $info->{END} =~ /^\d+$/) {
        $end = int($info->{END});
    } elsif ($svtype eq 'INS' || $svtype eq 'SNP') {
        $end = int($pos);
    } else {
        $end = int($pos + abs($svlen) - 1);
    }

    if (!$vartype || $vartype eq '.') {
        $vartype = $svtype;
    }
    return ($svtype, $vartype, $svlen, $end);
}

sub gt_is_carrier {
    my ($gt) = @_;
    return undef if !defined($gt) || $gt eq '.' || $gt =~ /^\.([\/|]\.)?$/;
    return 0 if $gt =~ /^0([\/|]0)*$/;
    return 1 if $gt =~ /[1-9]/;
    return undef;
}

sub normalize_allele {
    my ($pos, $ref, $alt) = @_;
    $ref = uc($ref // '');
    $alt = uc($alt // '');
    return ($pos, $ref, $alt) if $alt =~ /^<[^>]+>$/;

    # Remove identical suffix while keeping both alleles non-empty.
    while (length($ref) > 1 && length($alt) > 1 && substr($ref, -1) eq substr($alt, -1)) {
        chop $ref;
        chop $alt;
    }
    # Remove identical prefix while keeping both alleles non-empty.
    while (length($ref) > 1 && length($alt) > 1 && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
        $ref = substr($ref, 1);
        $alt = substr($alt, 1);
        $pos++;
    }
    return ($pos, $ref, $alt);
}

################################################################################
# Sample order and class filters
################################################################################
sub determine_sample_order {
    my ($truth_vcf, $evals, $opt) = @_;
    my @order;
    my %seen;

    if ($opt->{samples}) {
        open(my $IN, '<', $opt->{samples}) or die "ERROR: cannot open sample list $opt->{samples}: $!\n";
        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/^\s+|\s+$//g;
            next if $line eq '' || $line =~ /^#/;
            my @f = split(/\s+/, $line);
            my $s = $f[0];
            next if $seen{$s}++;
            push @order, $s;
        }
        close $IN;
    } else {
        @order = @{ $truth_vcf->{samples} };
        %seen = map { $_ => 1 } @order;
    }

    # If no samples are in truth header, use sample list only; otherwise keep only
    # samples that appear in truth.  This automatically removes a reference genome
    # name from g.list when the truth VCF contains only query samples.
    my %truth_sample = map { $_ => 1 } @{ $truth_vcf->{samples} || [] };
    if (%truth_sample) {
        @order = grep { $truth_sample{$_} } @order;
    }

    return @order;
}

sub record_in_class {
    my ($r, $class, $opt, $is_pred) = @_;
    $class = uc($class);
    return 1 if $class eq 'ALL';

    my $type = uc($r->{svtype} || '');
    my $vt   = uc($r->{vartype} || '');
    my $len  = abs($r->{signed_len} || $r->{len} || 0);

    # Predictions are usually already separated by SVGAP into different files.
    # Still apply light filtering to protect against mixed VCFs.
    if ($class eq 'SNP') {
        return 1 if $type eq 'SNP' || $vt eq 'SNP';
        return 1 if length($r->{ref}) == 1 && length($r->{alt}) == 1 && $r->{alt} !~ /^</;
        return 0;
    }
    if ($class eq 'SMALL') {
        return 0 if $type eq 'SNP' || $vt eq 'SNP';
        return 1 if $vt =~ /^S(?:DEL|INS|INDEL)$/;
        return 1 if ($type eq 'DEL' || $type eq 'INS') && $len > 0 && $len <= $opt->{small_max_len};
        return 0;
    }
    if ($class eq 'LARGE') {
        return 0 if $type eq 'SNP' || $vt eq 'SNP';
        return 0 if $vt =~ /^S(?:DEL|INS|INDEL)$/;
        return 1 if $type =~ /^(DEL|INS|INV|CNV|DUP|TRL|TRA|BND)$/ && ($is_pred || $len >= $opt->{large_min_len});
        return 0;
    }

    $class = 'CNV' if $class eq 'DUP';
    $class = 'TRL' if $class eq 'TRA' || $class eq 'BND';
    return $type eq $class;
}

################################################################################
# Comparison
################################################################################
sub compare_records {
    my ($truth, $pred, $samples, $eval, $opt) = @_;

    my %truth_by_chr;
    for my $i (0 .. $#$truth) {
        push @{ $truth_by_chr{ $truth->[$i]{chr} } }, $i;
    }

    my %truth_matched;
    my %pred_matched;
    my @matches;
    my @unmatched_pred;

    for my $pi (0 .. $#$pred) {
        my $p = $pred->[$pi];
        my ($best_ti, $best_score) = (-1, -1);
        for my $ti (@{ $truth_by_chr{ $p->{chr} } || [] }) {
            next if $truth_matched{$ti};
            my $t = $truth->[$ti];
            my $score = match_score($t, $p, $eval->{class}, $opt);
            if ($score > $best_score) {
                $best_score = $score;
                $best_ti = $ti;
            }
        }
        if ($best_ti >= 0 && $best_score > 0) {
            $truth_matched{$best_ti} = 1;
            $pred_matched{$pi} = 1;
            push @matches, { pred_i => $pi, truth_i => $best_ti, pred => $p, truth => $truth->[$best_ti], score => $best_score };
        } else {
            push @unmatched_pred, $p;
        }
    }

    my @unmatched_truth;
    for my $ti (0 .. $#$truth) {
        push @unmatched_truth, $truth->[$ti] unless $truth_matched{$ti};
    }

    my %event = (
        truth_n => scalar(@$truth),
        pred_n  => scalar(@$pred),
        tp      => scalar(@matches),
        fp      => scalar(@unmatched_pred),
        fn      => scalar(@unmatched_truth),
        tn      => 'NA',
    );
    add_prf(\%event);

    my %gt_total = (truth_n => scalar(@$truth), pred_n => scalar(@$pred), tp => 0, fp => 0, fn => 0, tn => 0);
    my %per_sample;
    for my $s (@$samples) {
        $per_sample{$s} = { truth_carriers => 0, pred_carriers => 0, tp => 0, fp => 0, fn => 0, tn => 0 };
    }

    # Count carriers in all truth and prediction records.
    for my $t (@$truth) {
        for my $s (@$samples) {
            $per_sample{$s}{truth_carriers}++ if carrier_value($t, $s, $opt) == 1;
        }
    }
    for my $p (@$pred) {
        for my $s (@$samples) {
            $per_sample{$s}{pred_carriers}++ if carrier_value($p, $s, $opt) == 1;
        }
    }

    # Matched records: compare per-sample carrier status.
    for my $m (@matches) {
        my ($t, $p) = ($m->{truth}, $m->{pred});
        for my $s (@$samples) {
            my $tc = carrier_value($t, $s, $opt);
            my $pc = carrier_value($p, $s, $opt);
            if ($tc == 1 && $pc == 1) {
                $gt_total{tp}++; $per_sample{$s}{tp}++;
            } elsif ($tc == 0 && $pc == 1) {
                $gt_total{fp}++; $per_sample{$s}{fp}++;
            } elsif ($tc == 1 && $pc == 0) {
                $gt_total{fn}++; $per_sample{$s}{fn}++;
            } elsif ($tc == 0 && $pc == 0) {
                $gt_total{tn}++; $per_sample{$s}{tn}++;
            }
        }
    }

    # Unmatched predicted records: all predicted carriers are genotype-level FPs.
    for my $p (@unmatched_pred) {
        for my $s (@$samples) {
            if (carrier_value($p, $s, $opt) == 1) {
                $gt_total{fp}++; $per_sample{$s}{fp}++;
            }
        }
    }

    # Unmatched truth records: all truth carriers are genotype-level FNs.
    for my $t (@unmatched_truth) {
        for my $s (@$samples) {
            if (carrier_value($t, $s, $opt) == 1) {
                $gt_total{fn}++; $per_sample{$s}{fn}++;
            }
        }
    }

    add_prf(\%gt_total);
    for my $s (@$samples) {
        add_prf($per_sample{$s});
    }

    return {
        event => \%event,
        gt => \%gt_total,
        per_sample => \%per_sample,
        matches => \@matches,
        unmatched_pred => \@unmatched_pred,
        unmatched_truth => \@unmatched_truth,
    };
}

sub match_score {
    my ($t, $p, $class, $opt) = @_;
    return 0 unless $t->{chr} eq $p->{chr};

    my $tt = canonical_type($t->{svtype});
    my $pt = canonical_type($p->{svtype});
    if ($opt->{require_type}) {
        return 0 unless $tt eq $pt;
    }

    $class = uc($class);
    if ($class eq 'SNP' || $tt eq 'SNP' || $pt eq 'SNP') {
        return match_snp($t, $p, $opt);
    }
    if ($class eq 'SMALL') {
        return match_small_indel($t, $p, $opt);
    }

    if ($tt eq 'INS' || $pt eq 'INS') {
        return match_insertion($t, $p, $opt);
    }
    return match_interval_sv($t, $p, $opt);
}

sub canonical_type {
    my ($x) = @_;
    $x = uc($x || '');
    return 'CNV' if $x =~ /^(DUP|TD|TANDEM_DUP|TANDEM_DUPLICATION)$/;
    return 'TRL' if $x =~ /^(TRA|BND|TRANSLOCATION)$/;
    return $x;
}

sub match_snp {
    my ($t, $p, $opt) = @_;
    return 0 if abs($t->{pos} - $p->{pos}) > $opt->{snp_pos_tol};
    return 1.0 if $t->{nref} eq $p->{nref} && $t->{nalt} eq $p->{nalt};
    # Allow same position and ALT if REF is unavailable/ambiguous in one VCF.
    return 0.8 if $t->{nalt} eq $p->{nalt} && $t->{alt} !~ /^</ && $p->{alt} !~ /^</;
    return 0;
}

sub match_small_indel {
    my ($t, $p, $opt) = @_;
    return 0 if abs($t->{pos} - $p->{pos}) > $opt->{small_pos_tol};
    if ($t->{nref} eq $p->{nref} && $t->{nalt} eq $p->{nalt}) {
        return 1.0;
    }
    my $len1 = max2(1, $t->{len});
    my $len2 = max2(1, $p->{len});
    my $ratio = ratio2($len1, $len2);
    return 0 if $ratio > 1 + $opt->{ins_len_ratio};
    return 0.5;
}

sub match_insertion {
    my ($t, $p, $opt) = @_;
    return 0 if abs($t->{pos} - $p->{pos}) > $opt->{ins_tol};
    my $len1 = max2(1, $t->{len});
    my $len2 = max2(1, $p->{len});
    my $ratio = ratio2($len1, $len2);
    return 0 if $ratio > 1 + $opt->{ins_len_ratio};

    # If both exact inserted alleles are available, reward exact sequence match.
    my $ti = inserted_sequence($t);
    my $pi = inserted_sequence($p);
    if (is_real_seq($ti) && is_real_seq($pi)) {
        return 1.0 if $ti eq $pi;
        return 0.8 if kmer_similarity($ti, $pi, 10) >= 0.50;
    }
    my $dist_score = 1 - (abs($t->{pos} - $p->{pos}) / max2(1, $opt->{ins_tol}));
    my $len_score  = 1 / $ratio;
    return 0.50 + 0.49 * (($dist_score + $len_score) / 2);
}

sub match_interval_sv {
    my ($t, $p, $opt) = @_;
    my $len1 = max2(1, $t->{end} - $t->{pos} + 1);
    my $len2 = max2(1, $p->{end} - $p->{pos} + 1);
    my $ratio = ratio2($len1, $len2);
    return 0 if $ratio > 1 + $opt->{sv_len_ratio};

    my $ol = overlap($t->{pos}, $t->{end}, $p->{pos}, $p->{end});
    if ($ol > 0) {
        my $ro1 = $ol / $len1;
        my $ro2 = $ol / $len2;
        if ($ro1 >= $opt->{sv_ov} && $ro2 >= $opt->{sv_ov}) {
            return 0.60 + 0.40 * (($ro1 + $ro2) / 2) * (1 / $ratio);
        }
    }
    if (abs($t->{pos} - $p->{pos}) <= $opt->{sv_bp_tol} && abs($t->{end} - $p->{end}) <= $opt->{sv_bp_tol}) {
        return 0.55 * (1 / $ratio);
    }
    return 0;
}

sub carrier_value {
    my ($r, $sample, $opt) = @_;
    return 0 unless exists $r->{carriers}{$sample};
    my $v = $r->{carriers}{$sample};
    return 0 if !defined $v && !$opt->{keep_missing_gt};
    return $v ? 1 : 0;
}

sub add_prf {
    my ($h) = @_;
    my $tp = $h->{tp} || 0;
    my $fp = $h->{fp} || 0;
    my $fn = $h->{fn} || 0;
    my $p = ($tp + $fp) ? $tp / ($tp + $fp) : 0;
    my $r = ($tp + $fn) ? $tp / ($tp + $fn) : 0;
    my $f = ($p + $r) ? 2 * $p * $r / ($p + $r) : 0;
    $h->{precision} = $p;
    $h->{recall} = $r;
    $h->{f1} = $f;
}

################################################################################
# Output
################################################################################
sub write_result {
    my ($SUM, $PS, $MAT, $UP, $UT, $eval, $res, $samples) = @_;
    my ($label, $class) = ($eval->{label}, $eval->{class});

    print_metric_row($SUM, $label, $class, 'EVENT', $res->{event});
    print_metric_row($SUM, $label, $class, 'GENOTYPE', $res->{gt});

    for my $s (@$samples) {
        my $h = $res->{per_sample}{$s};
        print $PS join("\t",
            $label, $class, $s,
            $h->{truth_carriers} || 0, $h->{pred_carriers} || 0,
            $h->{tp} || 0, $h->{fp} || 0, $h->{fn} || 0, $h->{tn} || 0,
            fmt($h->{precision}), fmt($h->{recall}), fmt($h->{f1})
        ), "\n";
    }

    for my $m (@{ $res->{matches} }) {
        my ($p, $t) = ($m->{pred}, $m->{truth});
        print $MAT join("\t",
            $label, $class, fmt($m->{score}),
            record_key($p), record_key($t),
            $p->{id}, $t->{id},
            $p->{chr}, $p->{pos}, $p->{end}, $p->{svtype}, $p->{len},
            $t->{chr}, $t->{pos}, $t->{end}, $t->{svtype}, $t->{len}
        ), "\n";
    }
    for my $p (@{ $res->{unmatched_pred} }) {
        print $UP join("\t", $label, $class, record_key($p), $p->{id}, $p->{chr}, $p->{pos}, $p->{end}, $p->{svtype}, $p->{len}, carrier_string($p)), "\n";
    }
    for my $t (@{ $res->{unmatched_truth} }) {
        print $UT join("\t", $label, $class, record_key($t), $t->{id}, $t->{chr}, $t->{pos}, $t->{end}, $t->{svtype}, $t->{len}, carrier_string($t)), "\n";
    }
}

sub print_metric_row {
    my ($fh, $label, $class, $level, $h) = @_;
    print $fh join("\t",
        $label, $class, $level,
        $h->{truth_n} // 'NA', $h->{pred_n} // 'NA',
        $h->{tp} || 0, $h->{fp} || 0, $h->{fn} || 0, defined($h->{tn}) ? $h->{tn} : 'NA',
        fmt($h->{precision}), fmt($h->{recall}), fmt($h->{f1})
    ), "\n";
}

sub record_key {
    my ($r) = @_;
    return join(':', $r->{chr}, $r->{pos}, $r->{end}, $r->{svtype}, $r->{len}, $r->{nref}, $r->{nalt});
}

sub carrier_string {
    my ($r) = @_;
    my @c = sort grep { $r->{carriers}{$_} } keys %{ $r->{carriers} || {} };
    return @c ? join(',', @c) : '.';
}

sub fmt {
    my ($x) = @_;
    return 'NA' unless defined $x && $x ne 'NA';
    return sprintf('%.6f', $x);
}

################################################################################
# Utility
################################################################################
sub dirname_of_path {
    my ($p) = @_;
    $p =~ s{/[^/]+$}{};
    return $p;
}

sub inserted_sequence {
    my ($r) = @_;
    return '' if $r->{alt} =~ /^</;
    my ($pos, $ref, $alt) = normalize_allele($r->{pos}, $r->{ref}, $r->{alt});
    return length($alt) > length($ref) ? substr($alt, length($ref)) : '';
}

sub is_real_seq {
    my ($s) = @_;
    return 0 unless defined $s;
    return 0 if $s eq '' || $s eq '.';
    return 0 if $s =~ /^N+$/i;
    return 1;
}

sub kmer_similarity {
    my ($s1, $s2, $k) = @_;
    return 0 unless is_real_seq($s1) && is_real_seq($s2);
    $k = 4 if length($s1) < 10 || length($s2) < 10;
    return 0 if length($s1) < $k || length($s2) < $k;
    my ($short, $long) = length($s1) <= length($s2) ? ($s1, $s2) : ($s2, $s1);
    my %long;
    for (my $i = 0; $i <= length($long) - $k; $i++) {
        $long{substr($long, $i, $k)} = 1;
    }
    my ($hit, $total) = (0, 0);
    for (my $i = 0; $i <= length($short) - $k; $i++) {
        $total++;
        $hit++ if $long{substr($short, $i, $k)};
    }
    return $total ? $hit / $total : 0;
}

sub overlap {
    my ($s1, $e1, $s2, $e2) = @_;
    my $s = max2($s1, $s2);
    my $e = min2($e1, $e2);
    return $e >= $s ? $e - $s + 1 : 0;
}

sub ratio2 {
    my ($a, $b) = @_;
    $a = max2(1, $a);
    $b = max2(1, $b);
    return $a >= $b ? $a / $b : $b / $a;
}

sub min2 { $_[0] < $_[1] ? $_[0] : $_[1] }
sub max2 { $_[0] > $_[1] ? $_[0] : $_[1] }
