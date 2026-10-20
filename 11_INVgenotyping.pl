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
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename);
use POSIX qw(strftime WIFEXITED WEXITSTATUS);

################################################################################
# 11_INVgenotyping.v2.robust.pl
#
# Original author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Original date: 2023-12-31
# Standardized/optimized date: 2026-06-13
#
# Purpose:
#   Genotype merged inversion candidates across population-scale query assemblies
#   using target-oriented single-coverage MAF files.
#
# Main improvements over the older 11_INVgenotyping.pl:
#   1. Does not require mafSplit or per-chromosome MAF splitting.
#   2. Preserves query/sample order from --query_list / --glist.
#   3. Uses robust orientation evidence from left flank, inversion body and
#      right flank.  A local inversion is called when body orientation is
#      opposite to flanking orientation.
#   4. Handles chromosome names containing dots, e.g. ISO1.seq111.
#   5. Avoids the old foreach(%hash) bug and reduces excessive missing ./. calls.
#   6. Uses built-in fork for parallelism; Parallel::ForkManager is not required.
#
# Coordinate assumption:
#   Input INV BED is 0-based half-open by default, as produced by SVGAP merged
#   BED files.  VCF POS is therefore start+1 and INFO/END is the BED END.
#   Use --bed_coord 1 if input POS/END are already 1-based VCF-like coordinates.
################################################################################

my %opt = (
    inv            => '',
    maf_dir        => '',
    sing           => '',       # alias of maf_dir
    refsize        => '',
    refname        => '',
    speciesname    => 'Unknown',
    query_list     => '',
    glist          => '',
    out            => '',
    outdir         => '',
    maf_suffix     => '.sing.maf',
    threads        => 1,
    bed_coord      => 0,        # 0: 0-based half-open; 1: already VCF-like 1-based

    # Genotyping parameters.
    flank_ratio    => 0.50,     # flank length = inv_len * flank_ratio
    min_flank      => 10000,
    max_flank      => 200000,
    min_body_cov   => 0.30,
    min_flank_cov  => 0.15,
    dominant_frac  => 0.70,
    min_body_bp    => 1000,
    trust_carriers => 1,        # candidates' MERGE_SAMPLES are directly set as 1|1
    body_only      => 0,        # if flanks unavailable, allow body-only calls
    phased         => 1,
    keep_tmp       => 0,
    debug          => 0,
    help           => 0,
);

GetOptions(
    'inv=s'              => \$opt{inv},
    'maf_dir|maf-dir=s'  => \$opt{maf_dir},
    'sing|singmaf=s'     => \$opt{sing},
    'refsize=s'          => \$opt{refsize},
    'refname|target=s'   => \$opt{refname},
    'speciesname=s'      => \$opt{speciesname},
    'query_list|query-list=s' => \$opt{query_list},
    'glist=s'            => \$opt{glist},
    'out=s'              => \$opt{out},
    'outdir=s'           => \$opt{outdir},
    'maf_suffix|maf-suffix=s' => \$opt{maf_suffix},
    'threads|t=i'        => \$opt{threads},
    'bed_coord=i'        => \$opt{bed_coord},

    'flank_ratio|flank-ratio=f' => \$opt{flank_ratio},
    'min_flank|min-flank=i'     => \$opt{min_flank},
    'max_flank|max-flank=i'     => \$opt{max_flank},
    'min_body_cov|min-body-cov=f' => \$opt{min_body_cov},
    'min_flank_cov|min-flank-cov=f' => \$opt{min_flank_cov},
    'dominant_frac|dominant-frac=f' => \$opt{dominant_frac},
    'min_body_bp|min-body-bp=i' => \$opt{min_body_bp},
    'trust_carriers!'    => \$opt{trust_carriers},
    'body_only!'         => \$opt{body_only},
    'phased!'            => \$opt{phased},
    'keep_tmp!'          => \$opt{keep_tmp},
    'debug!'             => \$opt{debug},
    'help|h'             => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

$opt{maf_dir} ||= $opt{sing};
$opt{outdir}  ||= 'INV_genotype_out';
$opt{out}     ||= "$opt{outdir}/INV.genotyped.vcf";

validate_options(\%opt);
make_path($opt{outdir}) unless -d $opt{outdir};
my $tmpdir = "$opt{outdir}/tmp.inv_genotype.$$";
make_path($tmpdir) unless -d $tmpdir;

my ($variants_ref, $header_ref) = read_inv_file($opt{inv}, \%opt);
my @variants = @{$variants_ref};

my @queries = read_queries(\%opt);
die "ERROR: no query samples found. Provide --query_list/--glist or check --maf_dir.\n" unless @queries;

if (!@variants) {
    warn "[WARN] no inversion records were loaded from $opt{inv}; writing an empty INV VCF and exiting normally.\n";
    write_empty_inv_vcf($opt{out}, \@queries, \%opt);
    open(my $ST, '>', "$opt{out}.sample_stats.tsv") or die "ERROR: cannot write $opt{out}.sample_stats.tsv: $!\n";
    print $ST join("\t", qw/SAMPLE N_INV GT_11 GT_00 GT_MISSING/), "\n";
    for my $q (@queries) { print $ST join("\t", $q, 0, 0, 0, 0), "\n"; }
    close $ST;
    remove_tree($tmpdir) unless $opt{keep_tmp};
    print STDERR "[", timestamp(), "] Done. Empty INV VCF: $opt{out}\n";
    exit 0;
}

my %queries_wanted = map { $_ => 1 } @queries;
my %maf_for = discover_mafs(\%opt, \@queries);
for my $q (@queries) {
    warn "[WARN] no MAF file found for query $q; all genotypes will be missing.\n" unless $maf_for{$q};
}

print STDERR "[", timestamp(), "] Loaded ", scalar(@variants), " INV records\n";
print STDERR "[", timestamp(), "] Query samples: ", scalar(@queries), "\n";
print STDERR "[", timestamp(), "] Output: $opt{out}\n";

run_query_jobs(\@queries, \%maf_for, \@variants, \%opt, $tmpdir);
write_final_vcf($opt{out}, \@variants, \@queries, \%opt, $tmpdir);
write_stats("$opt{out}.sample_stats.tsv", \@queries, $tmpdir);

unless ($opt{keep_tmp}) {
    remove_tree($tmpdir);
}

print STDERR "[", timestamp(), "] Done. VCF: $opt{out}\n";
exit 0;

################################################################################
# Input parsing
################################################################################
sub usage {
    return <<'USAGE';
Usage:
  perl 11_INVgenotyping.v2.robust.pl \
    --inv CombinedSV/All.INVs.bed.combined.sorted.txt \
    --maf_dir SynNetFilter_out/filtered_sing \
    --refsize ISO1.sizes \
    --refname ISO1 \
    --query_list query.list \
    --out INV.genotyped.vcf \
    --threads 8

Required:
  --inv          Merged inversion BED/TSV from SVGAP CombineSV
  --maf_dir      Directory containing target-oriented single-coverage MAF files
                 named <refname>.<query>.sing.maf by default
  --refsize      Reference genome .sizes file
  --refname      Reference/target genome name

Strongly recommended:
  --query_list   Query sample list; one query per line. Do not include reference.
                 This fixes VCF sample order.

Optional:
  --glist        Genome list. If used, --refname is removed and the rest are queries.
  --out          Output VCF [outdir/INV.genotyped.vcf]
  --outdir       Output directory [INV_genotype_out]
  --threads      Parallel query jobs [1]
  --maf_suffix   MAF suffix after <refname>.<query> [.sing.maf]
  --bed_coord 0  Input INV coordinates are 0-based half-open [default]
  --bed_coord 1  Input INV POS/END are already 1-based VCF-like

Genotyping parameters:
  --flank_ratio 0.50     Flank length = INV length * this ratio
  --min_flank 10000      Minimum flank length
  --max_flank 200000     Maximum flank length
  --min_body_cov 0.30    Minimum body coverage to make a call
  --min_flank_cov 0.15   Minimum combined flank coverage to use flank orientation
  --dominant_frac 0.70   Dominant orientation fraction threshold
  --trust_carriers       Set MERGE_SAMPLES carriers directly to 1|1 [default]
  --no-trust_carriers    Genotype every sample from MAF evidence only
  --body_only            Allow body-only calls when flank coverage is low [off]
  --keep_tmp             Keep temporary per-query genotype files

Output:
  A VCF with REF=N, ALT=<INV>, FORMAT=GT:GQ.
  GT is 1|1 for inverted, 0|0 for reference-like, and .|. for uncertain.
USAGE
}

sub validate_options {
    my ($opt) = @_;
    for my $k (qw/inv maf_dir refsize refname/) {
        die "ERROR: --$k is required.\n" . usage() unless defined $opt->{$k} && $opt->{$k} ne '';
    }
    die "ERROR: cannot find --inv $opt->{inv}\n" unless -s $opt->{inv};
    die "ERROR: cannot find --maf_dir $opt->{maf_dir}\n" unless -d $opt->{maf_dir};
    die "ERROR: cannot find --refsize $opt->{refsize}\n" unless -s $opt->{refsize};
    die "ERROR: --threads must be >= 1\n" unless $opt->{threads} >= 1;
    die "ERROR: --bed_coord must be 0 or 1\n" unless $opt->{bed_coord} == 0 || $opt->{bed_coord} == 1;
    die "ERROR: --dominant_frac must be between 0.5 and 1\n" unless $opt->{dominant_frac} >= 0.5 && $opt->{dominant_frac} <= 1;
    die "ERROR: --min_body_cov must be 0-1\n" unless $opt->{min_body_cov} >= 0 && $opt->{min_body_cov} <= 1;
    die "ERROR: --min_flank_cov must be 0-1\n" unless $opt->{min_flank_cov} >= 0 && $opt->{min_flank_cov} <= 1;
}

sub read_inv_file {
    my ($file, $opt) = @_;
    open(my $IN, '<', $file) or die "ERROR: cannot read $file: $!\n";
    my @vars;
    my @header;
    my $idx = 0;
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        if ($line =~ /^#/) {
            @header = split(/\t/, $line);
            next;
        }
        my @f = split(/\t/, $line);
        next unless @f >= 6;
        my $svtype = uc($f[4] || 'INV');
        next unless $svtype eq 'INV' || $svtype =~ /INV/;
        my ($chrom, $start, $end, $id, $len) = @f[0,1,2,3,5];
        $id = "INV_$idx" if !defined($id) || $id eq '' || $id eq '.';
        $len = abs(num($len));
        $len = abs(num($end) - num($start)) if !$len;
        my $bed_start = num($start);
        my $bed_end   = num($end);
        if ($opt->{bed_coord} == 1) {
            # Convert VCF-like POS/END to internal 0-based half-open.
            $bed_start = num($start) - 1;
            $bed_end   = num($end);
        }
        next if $bed_end <= $bed_start;
        my %carriers;
        if (defined $f[7]) {
            for my $m (split(/,/, $f[7])) {
                my $s = sample_from_member($m);
                $carriers{$s} = 1 if $s ne '' && $s ne '.';
            }
        }
        push @vars, {
            index     => $idx++,
            chrom     => $chrom,
            start     => $bed_start,
            end       => $bed_end,
            orig_pos  => num($start),
            orig_end  => num($end),
            id        => $id,
            svtype    => 'INV',
            len       => $len,
            raw       => $line,
            carriers  => \%carriers,
        };
    }
    close $IN;
    return (\@vars, \@header);
}

sub sample_from_member {
    my ($m) = @_;
    return '' unless defined $m;
    $m =~ s/^\s+|\s+$//g;
    return '' if $m eq '' || $m eq '.';
    # Common SVGAP member ID: SIM001.Chr01_100_200_100 or SIM001.seq111_...
    if ($m =~ /^([^\.\s]+)\./) { return $1; }
    if ($m =~ /^([^_\s]+)_/)   { return $1; }
    return $m;
}

sub read_queries {
    my ($opt) = @_;
    my @q;
    my $file = $opt->{query_list} || $opt->{glist} || '';
    if ($file) {
        open(my $IN, '<', $file) or die "ERROR: cannot read query/g list $file: $!\n";
        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/^\s+|\s+$//g;
            next if $line eq '' || $line =~ /^#/;
            my ($id) = split(/\s+/, $line);
            next if defined $opt->{refname} && $id eq $opt->{refname};
            push @q, $id;
        }
        close $IN;
        return unique_ordered(@q);
    }

    opendir(my $dh, $opt->{maf_dir}) or die "ERROR: cannot open $opt->{maf_dir}: $!\n";
    while (my $f = readdir($dh)) {
        next if $f =~ /^\./;
        if ($f =~ /^\Q$opt->{refname}\E\.(.+)\Q$opt->{maf_suffix}\E(?:\.gz)?$/) {
            push @q, $1;
        }
    }
    closedir($dh);
    return sort unique_ordered(@q);
}

sub unique_ordered {
    my %seen;
    my @out;
    for my $x (@_) {
        next if $seen{$x}++;
        push @out, $x;
    }
    return @out;
}

sub discover_mafs {
    my ($opt, $queries) = @_;
    my %maf;
    for my $q (@{$queries}) {
        my @candidates = (
            "$opt->{maf_dir}/$opt->{refname}.$q$opt->{maf_suffix}",
            "$opt->{maf_dir}/$opt->{refname}.$q$opt->{maf_suffix}.gz",
            "$opt->{maf_dir}/$q$opt->{maf_suffix}",
            "$opt->{maf_dir}/$q$opt->{maf_suffix}.gz",
        );
        for my $f (@candidates) {
            if (-s $f) { $maf{$q} = $f; last; }
        }
    }
    return %maf;
}

################################################################################
# Genotyping
################################################################################
sub run_query_jobs {
    my ($queries, $maf_for, $vars, $opt, $tmpdir) = @_;
    my $threads = $opt->{threads} || 1;
    if ($threads <= 1 || @{$queries} <= 1) {
        for my $q (@{$queries}) {
            genotype_one_query($q, $maf_for->{$q}, $vars, $opt, $tmpdir);
        }
        return;
    }

    my %children;
    my $had_error = 0;
    for my $q (@{$queries}) {
        while (scalar(keys %children) >= $threads) {
            my $pid = wait();
            last if $pid == -1;
            my $status = $?;
            my $name = delete $children{$pid};
            if ($status != 0) {
                warn "ERROR: child failed for query $name, status=$status\n";
                $had_error = 1;
            }
        }
        my $pid = fork();
        die "ERROR: fork failed: $!\n" unless defined $pid;
        if ($pid == 0) {
            eval { genotype_one_query($q, $maf_for->{$q}, $vars, $opt, $tmpdir); };
            if ($@) { warn "ERROR in query $q: $@\n"; exit 1; }
            exit 0;
        }
        $children{$pid} = $q;
    }
    while (scalar keys %children) {
        my $pid = wait();
        last if $pid == -1;
        my $status = $?;
        my $name = delete $children{$pid};
        if ($status != 0) {
            warn "ERROR: child failed for query $name, status=$status\n";
            $had_error = 1;
        }
    }
    die "ERROR: at least one query job failed\n" if $had_error;
}

sub genotype_one_query {
    my ($query, $maf, $vars, $opt, $tmpdir) = @_;
    my $outfile = "$tmpdir/$query.gt.tsv";
    open(my $OUT, '>', $outfile) or die "ERROR: cannot write $outfile: $!\n";
    print $OUT join("\t", qw/INDEX GT GQ BODY_COV FLANK_COV BODY_STRAND FLANK_STRAND BODY_PLUS BODY_MINUS FLANK_PLUS FLANK_MINUS REASON/), "\n";

    my @gt = map { init_evidence($_, $opt) } @{$vars};

    if (!$maf || !-s $maf) {
        for my $i (0 .. $#{$vars}) {
            print_gt_line($OUT, $vars->[$i]{index}, missing_call('NO_MAF'));
        }
        close $OUT;
        return;
    }

    # Trust the samples that originally supported this merged INV candidate.
    if ($opt->{trust_carriers}) {
        for my $i (0 .. $#{$vars}) {
            if ($vars->[$i]{carriers}{$query}) {
                $gt[$i]{trusted} = 1;
            }
        }
    }

    my %by_chr;
    for my $i (0 .. $#{$vars}) {
        push @{$by_chr{$vars->[$i]{chrom}}}, $i;
    }
    for my $chr (keys %by_chr) {
        @{$by_chr{$chr}} = sort { $vars->[$a]{start} <=> $vars->[$b]{start} } @{$by_chr{$chr}};
    }

    my $fh = open_maybe_gz($maf);
    local $/ = "\n\n";
    while (my $block = <$fh>) {
        # A MAF header is commonly followed immediately by the first alignment
        # block.  Do not discard that block merely because the paragraph starts
        # with a comment/header line; selecting s-lines already ignores headers.
        my @s = grep { /^s\s+/ } split(/\n/, $block);
        next unless @s >= 2;
        my $t = parse_maf_sline($s[0]);
        my $q = parse_maf_sline($s[1]);
        next unless $t && $q;
        next unless exists $by_chr{$t->{src}};

        my ($ts, $te) = normalized_interval($t);
        next if $te <= $ts;
        my $rel_strand = relative_strand($t->{strand}, $q->{strand});

        my $idxs = $by_chr{$t->{src}};
        for my $idx (@{$idxs}) {
            my $v = $vars->[$idx];
            my $flank = flank_len($v, $opt);
            my $win_s = $v->{start} - $flank;
            my $win_e = $v->{end} + $flank;
            last if $v->{start} > $te + $flank;   # sorted by start
            next if $win_e < $ts;
            next if $win_s > $te;

            add_overlap_evidence($gt[$idx]{left},  $ts, $te, max2(0, $v->{start} - $flank), $v->{start}, $rel_strand);
            add_overlap_evidence($gt[$idx]{body},  $ts, $te, $v->{start}, $v->{end}, $rel_strand);
            add_overlap_evidence($gt[$idx]{right}, $ts, $te, $v->{end}, $v->{end} + $flank, $rel_strand);
        }
    }
    close $fh;
    local $/ = "\n";

    for my $i (0 .. $#{$vars}) {
        my $call = call_from_evidence($gt[$i], $vars->[$i], $opt);
        print_gt_line($OUT, $vars->[$i]{index}, $call);
    }
    close $OUT;
}

sub init_evidence {
    my ($v, $opt) = @_;
    return {
        left    => { plus => 0, minus => 0 },
        body    => { plus => 0, minus => 0 },
        right   => { plus => 0, minus => 0 },
        trusted => 0,
    };
}

sub add_overlap_evidence {
    my ($e, $a1, $a2, $b1, $b2, $strand) = @_;
    return if $b2 <= $b1;
    my $ov = overlap($a1, $a2, $b1, $b2);
    return unless $ov > 0;
    if ($strand eq '-') { $e->{minus} += $ov; }
    else                { $e->{plus}  += $ov; }
}

sub call_from_evidence {
    my ($e, $v, $opt) = @_;
    if ($e->{trusted}) {
        return {
            gt => gt_string(1, $opt), gq => 22,
            body_cov => 1, flank_cov => 1,
            body_strand => 'TRUST', flank_strand => 'TRUST',
            body_plus => 0, body_minus => 0, flank_plus => 0, flank_minus => 0,
            reason => 'MERGE_SAMPLE'
        };
    }

    my $body_len = max2(1, $v->{end} - $v->{start});
    my $flank    = flank_len($v, $opt);
    my $left_len = $v->{start} < $flank ? $v->{start} : $flank;
    my $flank_len = max2(1, $left_len + $flank);

    my $bp = $e->{body}{plus};
    my $bm = $e->{body}{minus};
    my $fp = $e->{left}{plus} + $e->{right}{plus};
    my $fm = $e->{left}{minus} + $e->{right}{minus};

    my $bt = $bp + $bm;
    my $ft = $fp + $fm;
    my $body_cov = $bt / $body_len;
    my $flank_cov = $ft / $flank_len;

    return missing_call('LOW_BODY_COV', $body_cov, $flank_cov, $bp, $bm, $fp, $fm)
        if $body_cov < $opt->{min_body_cov} || $bt < $opt->{min_body_bp};

    my ($bstrand, $bfrac) = dominant_strand($bp, $bm);
    my ($fstrand, $ffrac) = dominant_strand($fp, $fm);

    if ($flank_cov >= $opt->{min_flank_cov} && $bfrac >= $opt->{dominant_frac} && $ffrac >= $opt->{dominant_frac}) {
        if ($bstrand ne $fstrand) {
            return made_call(1, $opt, $body_cov, $flank_cov, $bstrand, $fstrand, $bp, $bm, $fp, $fm, 'BODY_OPPOSITE_FLANK');
        } else {
            return made_call(0, $opt, $body_cov, $flank_cov, $bstrand, $fstrand, $bp, $bm, $fp, $fm, 'BODY_SAME_AS_FLANK');
        }
    }

    if ($opt->{body_only} && $bfrac >= $opt->{dominant_frac}) {
        # This is useful when flank alignment is absent, but it can confuse whole-
        # chromosome orientation changes.  It is therefore off by default.
        if ($bstrand eq '-') {
            return made_call(1, $opt, $body_cov, $flank_cov, $bstrand, 'LOWCOV', $bp, $bm, $fp, $fm, 'BODY_ONLY_MINUS');
        } else {
            return made_call(0, $opt, $body_cov, $flank_cov, $bstrand, 'LOWCOV', $bp, $bm, $fp, $fm, 'BODY_ONLY_PLUS');
        }
    }

    return missing_call('AMBIGUOUS', $body_cov, $flank_cov, $bp, $bm, $fp, $fm, $bstrand, $fstrand);
}

sub missing_call {
    my ($reason, $bc, $fc, $bp, $bm, $fp, $fm, $bs, $fs) = @_;
    return {
        gt => '.|.', gq => 0,
        body_cov => defined($bc) ? sprintf('%.4f', $bc) : 0,
        flank_cov => defined($fc) ? sprintf('%.4f', $fc) : 0,
        body_strand => $bs || '.', flank_strand => $fs || '.',
        body_plus => $bp || 0, body_minus => $bm || 0,
        flank_plus => $fp || 0, flank_minus => $fm || 0,
        reason => $reason || 'MISSING',
    };
}

sub made_call {
    my ($allele, $opt, $bc, $fc, $bs, $fs, $bp, $bm, $fp, $fm, $reason) = @_;
    return {
        gt => gt_string($allele, $opt), gq => 22,
        body_cov => sprintf('%.4f', $bc), flank_cov => sprintf('%.4f', $fc),
        body_strand => $bs, flank_strand => $fs,
        body_plus => $bp || 0, body_minus => $bm || 0,
        flank_plus => $fp || 0, flank_minus => $fm || 0,
        reason => $reason,
    };
}

sub gt_string {
    my ($allele, $opt) = @_;
    my $sep = $opt->{phased} ? '|' : '/';
    return $allele ? "1${sep}1" : "0${sep}0";
}

sub print_gt_line {
    my ($OUT, $idx, $c) = @_;
    print $OUT join("\t",
        $idx, $c->{gt}, $c->{gq}, $c->{body_cov}, $c->{flank_cov},
        $c->{body_strand}, $c->{flank_strand}, $c->{body_plus}, $c->{body_minus},
        $c->{flank_plus}, $c->{flank_minus}, $c->{reason}
    ), "\n";
}

sub dominant_strand {
    my ($plus, $minus) = @_;
    my $tot = $plus + $minus;
    return ('.', 0) if $tot <= 0;
    return $plus >= $minus ? ('+', $plus / $tot) : ('-', $minus / $tot);
}

sub flank_len {
    my ($v, $opt) = @_;
    my $len = max2(1, $v->{end} - $v->{start});
    my $flank = int($len * $opt->{flank_ratio} + 0.5);
    $flank = $opt->{min_flank} if $flank < $opt->{min_flank};
    $flank = $opt->{max_flank} if $flank > $opt->{max_flank};
    return $flank;
}

################################################################################
# MAF utilities
################################################################################
sub open_maybe_gz {
    my ($file) = @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gzip -dc '$file' |") or die "ERROR: cannot open gzipped MAF $file: $!\n";
    } else {
        open($fh, '<', $file) or die "ERROR: cannot open MAF $file: $!\n";
    }
    return $fh;
}

sub parse_maf_sline {
    my ($line) = @_;
    my @f = split(/\s+/, $line);
    return undef unless @f >= 7 && $f[0] eq 's';
    return {
        src    => $f[1],
        start  => num($f[2]),
        size   => num($f[3]),
        strand => $f[4],
        srcSize=> num($f[5]),
        text   => $f[6],
    };
}

sub normalized_interval {
    my ($s) = @_;
    if ($s->{strand} && $s->{strand} eq '-') {
        my $end = $s->{srcSize} - $s->{start};
        my $start = $end - $s->{size};
        return ($start, $end);
    }
    return ($s->{start}, $s->{start} + $s->{size});
}

sub relative_strand {
    my ($tstrand, $qstrand) = @_;
    $tstrand ||= '+';
    $qstrand ||= '+';
    return $tstrand eq $qstrand ? '+' : '-';
}

################################################################################
# Output
################################################################################
sub write_empty_inv_vcf {
    my ($out, $queries, $opt) = @_;
    open(my $VCF, '>', $out) or die "ERROR: cannot write $out: $!\n";
    my $date = strftime('%Y%m%d', localtime);
    print $VCF "##fileformat=VCFv4.2\n";
    print $VCF "##fileDate=$date\n";
    print $VCF "##source=SVGAP_11_INVgenotyping.v2.robust_empty\n";
    print $VCF "##reference=$opt->{refname}\n";
    print_contigs($VCF, $opt);
    print $VCF "##ALT=<ID=INV,Description=\"Inversion\">\n";
    print $VCF "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    print $VCF "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n";
    print $VCF "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">\n";
    print $VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $VCF "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
    print $VCF join("\t", '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', @{$queries}), "\n";
    close $VCF;
}

sub write_final_vcf {
    my ($out, $vars, $queries, $opt, $tmpdir) = @_;
    my %gt;
    for my $q (@{$queries}) {
        my $file = "$tmpdir/$q.gt.tsv";
        next unless -s $file;
        open(my $IN, '<', $file) or die "ERROR: cannot read $file: $!\n";
        my $h = <$IN>;
        while (my $line = <$IN>) {
            chomp $line;
            my @f = split(/\t/, $line);
            next unless @f >= 3;
            $gt{$q}{$f[0]} = "$f[1]:$f[2]";
        }
        close $IN;
    }

    open(my $VCF, '>', $out) or die "ERROR: cannot write $out: $!\n";
    my $date = strftime('%Y%m%d', localtime);
    print $VCF "##fileformat=VCFv4.2\n";
    print $VCF "##fileDate=$date\n";
    print $VCF "##source=SVGAP_11_INVgenotyping.v2.robust\n";
    print $VCF "##reference=$opt->{refname}\n";
    print_contigs($VCF, $opt);
    print $VCF "##ALT=<ID=INV,Description=\"Inversion\">\n";
    print $VCF "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    print $VCF "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n";
    print $VCF "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">\n";
    print $VCF "##INFO=<ID=BED_START,Number=1,Type=Integer,Description=\"Original 0-based BED start used internally\">\n";
    print $VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $VCF "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
    print $VCF join("\t", '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', @{$queries}), "\n";

    for my $v (sort { $a->{chrom} cmp $b->{chrom} || $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end} } @{$vars}) {
        my $vcf_pos = $opt->{bed_coord} == 1 ? $v->{orig_pos} : $v->{start} + 1;
        my $vcf_end = $opt->{bed_coord} == 1 ? $v->{orig_end} : $v->{end};
        my $id = $v->{id};
        $id = "INV_$v->{chrom}_$vcf_pos" if !defined($id) || $id eq '' || $id eq '.';
        my $info = join(';', "SVTYPE=INV", "END=$vcf_end", "SVLEN=$v->{len}", "BED_START=$v->{start}");
        my @sample_gt = map { $gt{$_}{$v->{index}} // '.|.:0' } @{$queries};
        print $VCF join("\t", $v->{chrom}, $vcf_pos, $id, 'N', '<INV>', '.', 'PASS', $info, 'GT:GQ', @sample_gt), "\n";
    }
    close $VCF;
}

sub print_contigs {
    my ($VCF, $opt) = @_;
    open(my $IN, '<', $opt->{refsize}) or die "ERROR: cannot read $opt->{refsize}: $!\n";
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my @f = split(/\s+/, $line);
        next unless @f >= 2;
        print $VCF "##contig=<ID=$f[0],length=$f[1],assembly=$opt->{refname},species=\"$opt->{speciesname}\">\n";
    }
    close $IN;
}

sub write_stats {
    my ($out, $queries, $tmpdir) = @_;
    open(my $OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";
    print $OUT join("\t", qw/SAMPLE N_INV GT_11 GT_00 GT_MISSING/), "\n";
    for my $q (@{$queries}) {
        my ($n, $g1, $g0, $gm) = (0,0,0,0);
        my $file = "$tmpdir/$q.gt.tsv";
        if (-s $file) {
            open(my $IN, '<', $file) or die "ERROR: cannot read $file: $!\n";
            my $h = <$IN>;
            while (my $line = <$IN>) {
                chomp $line;
                my @f = split(/\t/, $line);
                $n++;
                if ($f[1] =~ /^1[\/|]1/) { $g1++; }
                elsif ($f[1] =~ /^0[\/|]0/) { $g0++; }
                else { $gm++; }
            }
            close $IN;
        }
        print $OUT join("\t", $q, $n, $g1, $g0, $gm), "\n";
    }
    close $OUT;
}

################################################################################
# Small utilities
################################################################################
sub overlap {
    my ($s1, $e1, $s2, $e2) = @_;
    my $s = max2($s1, $s2);
    my $e = min2($e1, $e2);
    return $e > $s ? $e - $s : 0;
}
sub num {
    my ($x) = @_;
    return 0 unless defined $x && $x ne '';
    $x =~ s/,//g;
    return $x + 0;
}
sub min2 { $_[0] < $_[1] ? $_[0] : $_[1] }
sub max2 { $_[0] > $_[1] ? $_[0] : $_[1] }
sub timestamp { strftime('%F %T', localtime) }
