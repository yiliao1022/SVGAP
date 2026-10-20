#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
# ============================================================================
# SVGAP standardized script: 09_SNPgenotyping.pl
# Original author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Standardized date: 2026-06-12
#
# Parameter convention used in this standardized package:
#   --refname      reference genome ID/prefix
#   --glist        genome list, one genome ID per line; --refname is excluded
#                 from query/sample order when appropriate
#   --query_list   query/sample list; if supplied, VCF samples follow this order
#   --threads / -t number of parallel workers
#   --outdir       output directory; --out is used only for output prefixes
#
# Coordinates follow the original SVGAP convention unless stated otherwise in
# each script: raw BED/SV tables are 0-based half-open; VCF POS is 1-based.
# ============================================================================
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
use POSIX qw(strftime);
use File::Basename qw(basename dirname);
use File::Path qw(make_path remove_tree);
use File::Spec;
use IO::File;

################################################################################
# 9_SNPgenotyping_fast.pl
#
# Fast SNP genotyping from multiple single-coverage pairwise MAF files.
#
# Major changes relative to the older 9_SNPgenotyping.pl:
#   1. Parse each query MAF only once.
#   2. Do not require mafSplit/faSplit and avoid thousands of intermediate VCFs.
#   3. Store compact per-query SNP and callable-interval files.
#   4. Use built-in fork for parallel jobs; Parallel::ForkManager is not required.
#   5. Distinguish reference genotype from missing genotype using callable intervals.
#   6. Use proper VCF INFO fields: NS, AN, AC, AF.
#   7. Skip ambiguous bases by default to reduce false SNPs.
#
# Input MAF names are expected to be:
#   <refname>.<query>.sing.maf
# or
#   <refname>.<query>.sing.maf.gz
# in --singmaf.
################################################################################

my ($singmaf, $refsize, $refseq, $refname, $speciesname, $outdir, $threads, $help);
my ($tmpdir, $keep_tmp, $allow_ambiguous, $min_ns, $phased, $legacy_format, $force, $query_list);

GetOptions(
    'singmaf=s'          => \$singmaf,
    'refsize=s'          => \$refsize,
    'refseq=s'           => \$refseq,      # kept for compatibility; not required by this version
    'refname|target=s'   => \$refname,
    'query_list|query-list|glist=s' => \$query_list,
    'speciesname=s'      => \$speciesname,
    'outdir=s'           => \$outdir,
    't|threads=i'        => \$threads,
    'tmpdir=s'           => \$tmpdir,
    'keep_tmp!'          => \$keep_tmp,
    'allow_ambiguous!'   => \$allow_ambiguous,
    'minNS=i'            => \$min_ns,
    'phased!'            => \$phased,
    'legacy_format!'     => \$legacy_format,
    'force!'             => \$force,
    'help'               => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

for my $arg ([singmaf => $singmaf], [refsize => $refsize], [refname => $refname]) {
    die "ERROR: --$arg->[0] is required.\n\n" . usage() unless defined $arg->[1] && $arg->[1] ne '';
}

die "ERROR: --singmaf directory does not exist: $singmaf\n" unless -d $singmaf;
die "ERROR: --refsize file does not exist: $refsize\n" unless -s $refsize;

$speciesname    ||= 'Undefined';
$outdir         ||= 'SNPVCF_fast';
$threads        ||= 4;
$min_ns         ||= 1;
$phased = 1     unless defined $phased;
$legacy_format ||= 0;
$allow_ambiguous ||= 0;
$tmpdir         ||= File::Spec->catdir($outdir, 'tmp.snp_fast');

make_path($outdir) unless -d $outdir;
if (-d $tmpdir && $force) {
    remove_tree($tmpdir);
}
make_path($tmpdir) unless -d $tmpdir;

my ($chroms_ref, $chr_len_ref) = read_refsize($refsize);
my @chromosomes = @{$chroms_ref};
my %chr_len = %{$chr_len_ref};
my %chr_alias = build_chr_aliases($refname, \@chromosomes);

my @queries = find_query_mafs($singmaf, $refname, $query_list);
die "ERROR: No MAF files found in $singmaf. Expected files like $refname.<query>.sing.maf[.gz]\n" unless @queries;

print STDERR "[", timestamp(), "] Reference chromosomes: ", scalar(@chromosomes), "\n";
print STDERR "[", timestamp(), "] Query genomes: ", scalar(@queries), "\n";
print STDERR "[", timestamp(), "] Temporary directory: $tmpdir\n";

parse_all_queries(\@queries, \%chr_alias, $tmpdir, $threads);
write_final_vcf(\@queries, \@chromosomes, \%chr_len, $tmpdir, $outdir, $refname, $speciesname, $min_ns, $phased, $legacy_format);

unless ($keep_tmp) {
    remove_tree($tmpdir);
}

print STDERR "[", timestamp(), "] Done. Final VCF: $outdir/All.sorted.SNP.vcf\n";
exit 0;

################################################################################
# Main steps
################################################################################

sub parse_all_queries {
    my ($queries_ref, $chr_alias_ref, $tmpdir, $threads) = @_;

    $threads = 1 unless defined $threads && $threads =~ /^\d+$/ && $threads > 0;

    if ($threads <= 1) {
        for my $q (@{$queries_ref}) {
            parse_one_query($q, $chr_alias_ref, $tmpdir);
        }
        return;
    }

    # Use Perl built-in fork/wait instead of Parallel::ForkManager.
    # This avoids an extra CPAN dependency on clusters where users cannot install modules.
    my %children;
    my $had_error = 0;

    for my $q (@{$queries_ref}) {
        while (scalar(keys %children) >= $threads) {
            my $pid = wait();
            last if $pid == -1;
            my $status = $?;
            my $name = delete $children{$pid};
            if ($status != 0) {
                warn "ERROR: Child process failed while parsing query $name. Exit status: $status\n";
                $had_error = 1;
            }
        }

        my $pid = fork();
        die "ERROR: fork failed: $!\n" unless defined $pid;

        if ($pid == 0) {
            eval {
                parse_one_query($q, $chr_alias_ref, $tmpdir);
                1;
            } or do {
                my $err = $@ || 'unknown error';
                warn "ERROR: Failed to parse query $q->{query}: $err\n";
                exit 1;
            };
            exit 0;
        }

        $children{$pid} = $q->{query};
    }

    while (scalar(keys %children) > 0) {
        my $pid = wait();
        last if $pid == -1;
        my $status = $?;
        my $name = delete $children{$pid};
        if ($status != 0) {
            warn "ERROR: Child process failed while parsing query $name. Exit status: $status\n";
            $had_error = 1;
        }
    }

    die "ERROR: At least one query MAF parsing job failed. Please check the messages above.\n" if $had_error;
}

sub parse_one_query {
    my ($q, $chr_alias_ref, $tmpdir) = @_;
    my $query  = $q->{query};
    my $maf    = $q->{maf};
    my $qdir   = File::Spec->catdir($tmpdir, safe_name($query));
    make_path($qdir) unless -d $qdir;

    print STDERR "[", timestamp(), "] Parse MAF: $query\n";

    my (%snp_fh, %int_fh);
    my $fh = open_maybe_gzip($maf);
    my $block = '';
    while (my $line = <$fh>) {
        if ($line =~ /^a\s/ && $block ne '') {
            process_maf_block($block, $query, $qdir, $chr_alias_ref, \%snp_fh, \%int_fh);
            $block = $line;
        } else {
            $block .= $line;
        }
    }
    process_maf_block($block, $query, $qdir, $chr_alias_ref, \%snp_fh, \%int_fh) if $block ne '';
    close $fh;

    for my $h (values %snp_fh) { close $h; }
    for my $h (values %int_fh) { close $h; }

    print STDERR "[", timestamp(), "] Finished: $query\n";
}

sub process_maf_block {
    my ($block, $query, $qdir, $chr_alias_ref, $snp_fh_ref, $int_fh_ref) = @_;
    return if $block =~ /^\s*$/;
    return if $block =~ /^##maf/;

    my @s_lines;
    for my $line (split /\n/, $block) {
        next unless $line =~ /^s\s+/;
        my @f = split /\s+/, $line;
        next unless @f >= 7;
        push @s_lines, {
            src    => $f[1],
            start  => $f[2],
            size   => $f[3],
            strand => $f[4],
            srcsz  => $f[5],
            seq    => uc($f[6]),
        };
    }
    return unless @s_lines >= 2;

    my ($target_idx, $chr);
    for my $i (0 .. $#s_lines) {
        my $src = $s_lines[$i]{src};
        if (exists $chr_alias_ref->{$src}) {
            ($target_idx, $chr) = ($i, $chr_alias_ref->{$src});
            last;
        }
    }
    return unless defined $target_idx;

    my $query_idx;
    for my $i (0 .. $#s_lines) {
        next if $i == $target_idx;
        $query_idx = $i;
        last;
    }
    return unless defined $query_idx;

    my $t = $s_lines[$target_idx];
    my $q = $s_lines[$query_idx];
    return unless length($t->{seq}) == length($q->{seq});

    my $snp_fh = get_out_fh($snp_fh_ref, File::Spec->catfile($qdir, "$chr.snps.tsv"));
    my $int_fh = get_out_fh($int_fh_ref, File::Spec->catfile($qdir, "$chr.callable.bed"));

    my $tseq = $t->{seq};
    my $qseq = $q->{seq};
    my $aln_len = length($tseq);
    my $reverse = defined($t->{strand}) && $t->{strand} eq '-';
    my $ungapped_seen = 0;
    my ($run_start, $run_end, $prev_pos);

    for (my $i = 0; $i < $aln_len; $i++) {
        my $tb = substr($tseq, $i, 1);
        my $qb = substr($qseq, $i, 1);

        if ($tb ne '-') {
            $ungapped_seen++;
            my $pos = $reverse
                ? $t->{srcsz} - ($t->{start} + $ungapped_seen - 1)
                : $t->{start} + $ungapped_seen;
            if ($reverse) {
                $tb = revcomp_base($tb);
                $qb = revcomp_base($qb);
            }

            my $ref_callable = $allow_ambiguous
                ? $tb =~ /^[ACGTRYSWKMBDHVN]$/
                : $tb =~ /^[ACGT]$/;
            my $query_callable = $allow_ambiguous
                ? $qb =~ /^[ACGTRYSWKMBDHVN]$/
                : $qb =~ /^[ACGT]$/;
            if ($ref_callable && $query_callable) {
                if (!defined $run_start) {
                    ($run_start, $run_end, $prev_pos) = ($pos, $pos, $pos);
                } elsif (abs($pos - $prev_pos) == 1) {
                    $run_start = $pos if $pos < $run_start;
                    $run_end   = $pos if $pos > $run_end;
                    $prev_pos  = $pos;
                } else {
                    print {$int_fh} join("\t", $chr, $run_start, $run_end), "\n";
                    ($run_start, $run_end, $prev_pos) = ($pos, $pos, $pos);
                }

                if ($tb ne $qb) {
                    print {$snp_fh} join("\t", $pos, $tb, $qb), "\n";
                }
            } else {
                if (defined $run_start) {
                    print {$int_fh} join("\t", $chr, $run_start, $run_end), "\n";
                    undef $run_start;
                    undef $run_end;
                    undef $prev_pos;
                }
            }
        }
    }

    if (defined $run_start) {
        print {$int_fh} join("\t", $chr, $run_start, $run_end), "\n";
    }
}

sub write_final_vcf {
    my ($queries_ref, $chroms_ref, $chr_len_ref, $tmpdir, $outdir, $refname, $speciesname, $min_ns, $phased, $legacy_format) = @_;

    my $vcf_file = File::Spec->catfile($outdir, 'All.sorted.SNP.vcf');
    open my $vcf, '>', $vcf_file or die "ERROR: Cannot write $vcf_file: $!\n";

    print_vcf_header($vcf, $queries_ref, $chroms_ref, $chr_len_ref, $refname, $speciesname, $phased, $legacy_format);

    for my $chr (@{$chroms_ref}) {
        print STDERR "[", timestamp(), "] Merge chromosome: $chr\n";
        my ($sites_ref, $refbase_ref, $alt_ref) = collect_chr_sites($chr, $queries_ref, $tmpdir);
        my @positions = sort { $a <=> $b } keys %{$sites_ref};
        next unless @positions;

        my (%sample_snps, %sample_intervals);
        for my $q (@{$queries_ref}) {
            my $query = $q->{query};
            my $qdir = File::Spec->catdir($tmpdir, safe_name($query));
            $sample_snps{$query}      = load_sample_snps(File::Spec->catfile($qdir, "$chr.snps.tsv"));
            $sample_intervals{$query} = load_intervals(File::Spec->catfile($qdir, "$chr.callable.bed"));
        }

        for my $pos (@positions) {
            my $ref = $refbase_ref->{$pos};
            next unless defined $ref && valid_base($ref);

            my @alts = grep { exists $alt_ref->{$pos}{$_} && $_ ne $ref } qw(A C G T);
            next unless @alts;

            my %alt_index;
            for my $i (0 .. $#alts) {
                $alt_index{$alts[$i]} = $i + 1;
            }

            my (@gts, %ac);
            my $ns = 0;
            for my $q (@{$queries_ref}) {
                my $query = $q->{query};
                my $base;

                if (exists $sample_snps{$query}{$pos}) {
                    $base = $sample_snps{$query}{$pos};
                    $base = '.' unless valid_base($base);
                } elsif (is_pos_callable($sample_intervals{$query}, $pos)) {
                    $base = $ref;
                } else {
                    $base = '.';
                }

                if ($base eq '.') {
                    push @gts, missing_gt($phased, $legacy_format);
                } elsif ($base eq $ref) {
                    $ns++;
                    push @gts, hom_gt(0, $phased, $legacy_format);
                } elsif (exists $alt_index{$base}) {
                    $ns++;
                    $ac{$base} += 2;
                    push @gts, hom_gt($alt_index{$base}, $phased, $legacy_format);
                } else {
                    # A different clean base not retained in ALT due to inconsistent REF or filtering.
                    push @gts, missing_gt($phased, $legacy_format);
                }
            }
            next if $ns < $min_ns;

            my $an = $ns * 2;
            my @ac_values = map { $ac{$_} || 0 } @alts;
            my $total_ac = 0;
            $total_ac += $_ for @ac_values;
            next if $total_ac == 0;
            my @af_values = map { $an > 0 ? format_float(($_ || 0) / $an) : '0' } @ac_values;
            my $info = join(';',
                "NS=$ns",
                "AN=$an",
                'AC=' . join(',', @ac_values),
                'AF=' . join(',', @af_values),
            );
            my $format = $legacy_format ? 'GT:GQ:DP:HQ' : 'GT';
            print {$vcf} join("\t", $chr, $pos, '.', $ref, join(',', @alts), '.', 'PASS', $info, $format, @gts), "\n";
        }
    }
    close $vcf;
}

################################################################################
# File readers and parsers
################################################################################

sub read_refsize {
    my ($file) = @_;
    my (@chroms, %len);
    open my $fh, '<', $file or die "ERROR: Cannot read $file: $!\n";
    while (<$fh>) {
        chomp;
        next if /^\s*$/ || /^#/;
        my @f = split /\s+/;
        next unless @f >= 2;
        push @chroms, $f[0] unless exists $len{$f[0]};
        $len{$f[0]} = $f[1];
    }
    close $fh;
    die "ERROR: No chromosome sizes read from $file\n" unless @chroms;
    return (\@chroms, \%len);
}

sub build_chr_aliases {
    my ($refname, $chroms_ref) = @_;
    my %alias;
    for my $chr (@{$chroms_ref}) {
        $alias{$chr} = $chr;
        my $short = $chr;
        $short =~ s/^\Q$refname\E[\._-]?//;
        $alias{$short} = $chr if $short ne '';
        $alias{"$refname.$short"} = $chr if $short ne '';
        $alias{"${refname}_$short"} = $chr if $short ne '';
        $alias{"$refname-$short"} = $chr if $short ne '';
    }
    return %alias;
}

sub find_query_mafs {
    my ($dir, $refname, $query_list) = @_;
    my @queries;
    my %seen;

    if (defined $query_list && $query_list ne '') {
        open my $QL, '<', $query_list or die "ERROR: Cannot open --query_list/--glist $query_list: $!\n";
        while (my $q = <$QL>) {
            chomp $q;
            $q =~ s/^\s+|\s+$//g;
            next if $q eq '' || $q =~ /^#/;
            next if $q eq $refname;  # allow g.list to contain the reference
            my $file = "$refname.$q.sing.maf";
            my $path = File::Spec->catfile($dir, $file);
            if (!-s $path && -s "$path.gz") {
                $path = "$path.gz";
            }
            die "ERROR: MAF file for query '$q' does not exist or is empty: $path\n" unless -s $path;
            die "ERROR: Duplicated query in --query_list/--glist: $q\n" if $seen{$q}++;
            push @queries, { query => $q, maf => $path };
        }
        close $QL;
        return @queries;
    }

    opendir my $dh, $dir or die "ERROR: Cannot open $dir: $!\n";
    while (my $file = readdir $dh) {
        next if $file =~ /^\./;
        my $query;
        if ($file =~ /^\Q$refname\E\.(.+)\.sing\.maf(?:\.gz)?$/) {
            $query = $1;
        } else {
            next;
        }
        my $path = File::Spec->catfile($dir, $file);
        next unless -s $path;
        die "ERROR: Duplicated query inferred from MAF files: $query\n" if $seen{$query}++;
        push @queries, { query => $query, maf => $path };
    }
    closedir $dh;
    @queries = sort { $a->{query} cmp $b->{query} } @queries;
    return @queries;
}

sub collect_chr_sites {
    my ($chr, $queries_ref, $tmpdir) = @_;
    my (%sites, %refbase, %alts);
    for my $q (@{$queries_ref}) {
        my $query = $q->{query};
        my $file = File::Spec->catfile($tmpdir, safe_name($query), "$chr.snps.tsv");
        next unless -s $file;
        open my $fh, '<', $file or die "ERROR: Cannot read $file: $!\n";
        while (<$fh>) {
            chomp;
            next if /^\s*$/;
            my ($pos, $ref, $alt) = split /\t/;
            next unless defined $pos && $pos =~ /^\d+$/;
            next unless valid_base($ref) && valid_base($alt);
            next if $ref eq $alt;
            if (exists $refbase{$pos} && $refbase{$pos} ne $ref) {
                # Conflicting REF usually means duplicated/overlapping alignment or alias issue.
                # Keep the first REF and ignore the conflicting record for safety.
                next;
            }
            $sites{$pos} = 1;
            $refbase{$pos} = $ref;
            $alts{$pos}{$alt} = 1;
        }
        close $fh;
    }
    return (\%sites, \%refbase, \%alts);
}

sub load_sample_snps {
    my ($file) = @_;
    my %snp;
    return \%snp unless -s $file;
    open my $fh, '<', $file or die "ERROR: Cannot read $file: $!\n";
    while (<$fh>) {
        chomp;
        next if /^\s*$/;
        my ($pos, $ref, $alt) = split /\t/;
        next unless defined $pos && $pos =~ /^\d+$/;
        next unless valid_base($alt);
        if (exists $snp{$pos} && $snp{$pos} ne $alt) {
            $snp{$pos} = 'AMBIG';
        } else {
            $snp{$pos} = $alt;
        }
    }
    close $fh;
    return \%snp;
}

sub load_intervals {
    my ($file) = @_;
    my @ints;
    return \@ints unless -s $file;
    open my $fh, '<', $file or die "ERROR: Cannot read $file: $!\n";
    while (<$fh>) {
        chomp;
        next if /^\s*$/;
        my ($chr, $s, $e) = split /\t/;
        next unless defined $s && defined $e && $s =~ /^\d+$/ && $e =~ /^\d+$/;
        ($s, $e) = ($e, $s) if $s > $e;
        push @ints, [$s, $e];
    }
    close $fh;
    return merge_intervals(\@ints);
}

sub merge_intervals {
    my ($ints_ref) = @_;
    my @sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$ints_ref};
    my @merged;
    for my $int (@sorted) {
        if (!@merged || $int->[0] > $merged[-1][1] + 1) {
            push @merged, [@$int];
        } else {
            $merged[-1][1] = $int->[1] if $int->[1] > $merged[-1][1];
        }
    }
    return \@merged;
}

################################################################################
# Small utilities
################################################################################

sub get_out_fh {
    my ($href, $file) = @_;
    return $href->{$file} if exists $href->{$file};
    open my $fh, '>', $file or die "ERROR: Cannot write $file: $!\n";
    $href->{$file} = $fh;
    return $fh;
}

sub open_maybe_gzip {
    my ($file) = @_;
    if ($file =~ /\.gz$/) {
        open my $fh, '-|', 'gzip', '-dc', $file or die "ERROR: Cannot open gzip -dc $file: $!\n";
        return $fh;
    }
    open my $fh, '<', $file or die "ERROR: Cannot read $file: $!\n";
    return $fh;
}

sub maf_target_pos_1based {
    my ($start0, $ungapped_seen, $strand, $src_size) = @_;
    if (defined $strand && $strand eq '-') {
        return $src_size - ($start0 + $ungapped_seen - 1);
    }
    return $start0 + $ungapped_seen;
}

sub orient_bases_to_reference {
    my ($tb, $qb, $target_strand) = @_;
    if (defined $target_strand && $target_strand eq '-') {
        $tb = revcomp_base($tb);
        $qb = revcomp_base($qb);
    }
    return ($tb, $qb);
}

sub revcomp_base {
    my ($b) = @_;
    return '-' if !defined $b || $b eq '-';
    $b =~ tr/ACGTRYMKSWBDHVNacgtrymkswbdhvn/TGCAYRKMWSVHDBNtgcayrkmwsvhdbn/;
    return uc($b);
}

sub is_callable_base {
    my ($b) = @_;
    return 0 unless defined $b;
    return 1 if $b =~ /^[ACGT]$/;
    return 1 if $allow_ambiguous && $b =~ /^[ACGTRYSWKMBDHVN]$/;
    return 0;
}

sub valid_base {
    my ($b) = @_;
    return defined $b && $b =~ /^[ACGT]$/;
}

sub is_pos_callable {
    my ($ints_ref, $pos) = @_;
    return 0 unless $ints_ref && @{$ints_ref};
    my ($lo, $hi) = (0, scalar(@{$ints_ref}) - 1);
    while ($lo <= $hi) {
        my $mid = int(($lo + $hi) / 2);
        my ($s, $e) = @{$ints_ref->[$mid]};
        if ($pos < $s) {
            $hi = $mid - 1;
        } elsif ($pos > $e) {
            $lo = $mid + 1;
        } else {
            return 1;
        }
    }
    return 0;
}

sub hom_gt {
    my ($idx, $phased, $legacy_format) = @_;
    my $sep = $phased ? '|' : '/';
    my $gt = "$idx$sep$idx";
    return $legacy_format ? "$gt:48:10:51,51" : $gt;
}

sub missing_gt {
    my ($phased, $legacy_format) = @_;
    my $sep = $phased ? '|' : '/';
    my $gt = ".$sep.";
    return $legacy_format ? "$gt:.:.:." : $gt;
}

sub format_float {
    my ($x) = @_;
    my $s = sprintf('%.6f', $x);
    $s =~ s/0+$//;
    $s =~ s/\.$/.0/;
    return $s;
}

sub safe_name {
    my ($name) = @_;
    $name =~ s/[^A-Za-z0-9_.-]/_/g;
    return $name;
}

sub timestamp {
    return strftime('%F %T', localtime);
}

sub print_vcf_header {
    my ($fh, $queries_ref, $chroms_ref, $chr_len_ref, $refname, $speciesname, $phased, $legacy_format) = @_;
    my $date = strftime('%Y%m%d', localtime);
    print {$fh} "##fileformat=VCFv4.2\n";
    print {$fh} "##fileDate=$date\n";
    print {$fh} "##source=9_SNPgenotyping_fast.pl\n";
    print {$fh} "##reference=$refname\n";
    for my $chr (@{$chroms_ref}) {
        my $len = $chr_len_ref->{$chr};
        print {$fh} "##contig=<ID=$chr,length=$len,assembly=$refname,species=\"$speciesname\">\n";
    }
    print {$fh} "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with non-missing genotype\">\n";
    print {$fh} "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of called alleles in called genotypes\">\n";
    print {$fh} "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in called genotypes, for each ALT allele\">\n";
    print {$fh} "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for each ALT allele\">\n";
    print {$fh} "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    print {$fh} "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    if ($legacy_format) {
        print {$fh} "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Fixed legacy genotype quality placeholder\">\n";
        print {$fh} "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Fixed legacy depth placeholder\">\n";
        print {$fh} "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Fixed legacy haplotype quality placeholder\">\n";
    }
    print {$fh} "##phasing=", ($phased ? 'phased' : 'unphased'), "\n";
    print {$fh} join("\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', map { $_->{query} } @{$queries_ref}), "\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl 9_SNPgenotyping_fast.pl \
    --singmaf ./sing \
    --refsize ./Ref/S184.sizes \
    --refname S184 \
    --query_list query.list \
    --speciesname "Capsicum annuum" \
    --outdir SNPVCF_fast \
    --t 8

Required:
  --singmaf       Directory containing <refname>.<query>.sing.maf or .maf.gz files
  --refsize       Reference chromosome sizes file: <chromosome> <length>
  --refname       Reference name/prefix, e.g. S184, MH63, IRGSP1
  --query_list/--glist Optional query/sample list; VCF columns follow this order

Optional:
  --refseq        Kept for compatibility with old command lines; not used here
  --speciesname   Species name written into VCF contig header [Undefined]
  --outdir        Output directory [SNPVCF_fast]
  --tmpdir        Temporary directory [outdir/tmp.snp_fast]
  --t             Number of query MAF files parsed in parallel [4]
  --minNS         Keep SNPs called in at least this many samples [1]
  --phased        Output phased homozygous GT such as 0|0 and 1|1 [default]
  --nophased      Output unphased GT such as 0/0 and 1/1
  --legacy_format Output old-style GT:GQ:DP:HQ columns with placeholder values
  --allow_ambiguous
                  Treat IUPAC/N query bases as callable while parsing intervals.
                  SNP ALT alleles in the final VCF are still restricted to A/C/G/T.
  --keep_tmp      Keep temporary SNP/interval files
  --force         Remove existing tmpdir before running
  --help          Show this help message

Output:
  <outdir>/All.sorted.SNP.vcf

Notes:
  Coordinates are 1-based VCF coordinates converted from MAF target coordinates.
  Only columns where both reference and query bases are callable are used as SNP-genotyping sites.
  Query gaps, target gaps, and ambiguous bases are treated as missing by default.
USAGE
}
