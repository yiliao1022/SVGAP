#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
# ============================================================================
# SVGAP standardized script: 06_genotype_DEL.pl
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
use File::Basename qw(basename);
use File::Path qw(make_path remove_tree);
use POSIX qw(strftime WIFEXITED WEXITSTATUS);

# Genotype deletions from pairwise MAF files for multiple query assemblies in parallel.
#
# Input BED/TSV format with header, for example:
# #CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq
# S184.Chr01 9480 16324 GPE021490.Chr01_8674_9921_6844 DEL 6844 ...
#
# MAF file name format by default:
#   <target>.<query>.sing.maf
# Example:
#   S184.GPE046410.sing.maf
#
# MAF assumption:
#   Pairwise MAF.
#   First s line  = reference / target.
#   Second s line = query / sample.
#
# BED coordinate assumption:
#   POS and END are 0-based half-open coordinates.
#   For VCF output, POS is converted to POS + 1.
#
# Genotype meaning:
#   0 = query is reference-like, no deletion
#   1 = query carries the deletion relative to reference
#   . = uncertain

my %opt = (
    threads          => 4,
    target           => '',
    maf_suffix       => '.sing.maf',
    min_cov          => 0.8,
    min_support      => 0.8,
    bridge_del_ratio => 0.2,
    bridge_ref_low   => 0.8,
    bridge_ref_high  => 1.25,
    vcf_mode         => 'symbolic',     # symbolic or sequence
    diploid          => 1,
    keep_tmp         => 0,
    ref_fasta        => '',
);

GetOptions(
    'bed=s'              => \$opt{bed},
    'maf_dir=s'          => \$opt{maf_dir},
    'target|refname=s'           => \$opt{target},
    'ref_fasta|ref-fasta|ref=s'  => \$opt{ref_fasta},
    'query_list|query-list|glist=s'       => \$opt{query_list},
    'out=s'              => \$opt{out},
    'threads|t=i'          => \$opt{threads},
    'maf_suffix=s'       => \$opt{maf_suffix},
    'min_cov=f'          => \$opt{min_cov},
    'min_support=f'      => \$opt{min_support},
    'bridge_del_ratio=f' => \$opt{bridge_del_ratio},
    'bridge_ref_low=f'   => \$opt{bridge_ref_low},
    'bridge_ref_high=f'  => \$opt{bridge_ref_high},
    'vcf_mode=s'         => \$opt{vcf_mode},
    'haploid!'           => sub { $opt{diploid} = 0 },
    'diploid!'           => \$opt{diploid},
    'keep_tmp!'          => \$opt{keep_tmp},
    'help|h'             => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

die usage() unless $opt{bed} && $opt{maf_dir} && $opt{target} && $opt{out};
die "--threads must be >= 1\n" unless $opt{threads} >= 1;
die "--vcf_mode must be symbolic or sequence\n" unless $opt{vcf_mode} eq 'symbolic' || $opt{vcf_mode} eq 'sequence';

my ($variants_ref, $header_ref, $col_ref) = read_bed_variants($opt{bed});
my @variants = @$variants_ref;
die "No variants were loaded from $opt{bed}\n" unless @variants;
my $ref_db;
if ($opt{ref_fasta}) {
    die "Cannot read --ref_fasta $opt{ref_fasta}\n" unless -s $opt{ref_fasta};
    $ref_db = DELRefFAI->new($opt{ref_fasta});
}

my @jobs = discover_maf_jobs(\%opt);
die "No MAF files found. Check --maf_dir, --target, --query_list and --maf_suffix.\n" unless @jobs;

my $tmp_dir = "$opt{out}.tmp.DEL_genotype.$$";
make_path($tmp_dir) or die "Cannot create $tmp_dir: $!\n";

print STDERR "[", scalar localtime, "] Loaded ", scalar(@variants), " variants from $opt{bed}\n";
print STDERR "[", scalar localtime, "] Found ", scalar(@jobs), " query MAF files\n";
print STDERR "[", scalar localtime, "] Running with $opt{threads} parallel workers\n";

run_parallel_genotyping(\@jobs, $tmp_dir, \%opt);

my $gt_tsv = "$opt{out}.gt.tsv";
my $vcf    = "$opt{out}.vcf";
combine_gt_tsv(\@variants, $header_ref, \@jobs, $tmp_dir, $gt_tsv);
write_vcf(\@variants, \@jobs, $tmp_dir, $vcf, \%opt, $ref_db);

unless ($opt{keep_tmp}) {
    remove_tree($tmp_dir);
}

print STDERR "[", scalar localtime, "] Done.\n";
print STDERR "GT table: $gt_tsv\n";
print STDERR "VCF:      $vcf\n";

exit 0;

sub usage {
    return <<'USAGE';
Usage:
  perl genotype_DEL_parallel_from_maf.pl \
    --bed DEL.bed \
    --maf_dir sing_maf_dir \
    --target S184 \
    --out DEL.genotyped \
    --threads 8

Optional:
  --query_list query.list     or --glist g.list       one query name per line; sample name is the query name only
  --ref_fasta reference.fa      rebuild sequence-resolved left-anchored DEL alleles
  --maf_suffix .sing.maf        default: .sing.maf; MAF name is <target>.<query><maf_suffix>
  --min_cov 0.8                minimum confident coverage of the DEL interval
  --min_support 0.8            minimum DEL-like or REF-like support among confident bases
  --bridge_del_ratio 0.2       between-block query/ref gap ratio <= this is treated as DEL
  --bridge_ref_low 0.8         between-block query/ref gap ratio lower bound for REF
  --bridge_ref_high 1.25       between-block query/ref gap ratio upper bound for REF
  --vcf_mode symbolic          legacy compatibility option; when --ref_fasta is supplied, a reference-checked explicit DEL is preferred
  --haploid                    output haploid VCF GT values: 0, 1, .
  --diploid                    output diploid VCF GT values: 0/0, 1/1, ./. [default]
  --keep_tmp                   keep temporary per-query GT files

Input BED/TSV columns:
  #CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq

Output:
  <out>.gt.tsv  original BED columns plus one GT column per query
  <out>.vcf     multi-sample VCF with GT calls
USAGE
}

sub read_bed_variants {
    my ($bed_file) = @_;

    open my $BED, '<', $bed_file or die "Cannot open $bed_file: $!\n";

    my @header;
    my %col;
    my @variants;
    my $line_no = 0;

    while (my $line = <$BED>) {
        $line_no++;
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split_line($line);

        if ($line =~ /^#/) {
            @header = @f;
            $header[0] =~ s/^#//;
            for my $i (0 .. $#header) {
                $col{$header[$i]} = $i;
            }
            next;
        }

        if (!%col) {
            @header = qw(CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq);
            for my $i (0 .. $#header) {
                $col{$header[$i]} = $i;
            }
        }

        for my $required (qw(CHROM POS END ID SVTYPE)) {
            die "BED header lacks required column '$required'\n" unless exists $col{$required};
        }

        my $chr    = $f[$col{CHROM}];
        my $start  = $f[$col{POS}];
        my $end    = $f[$col{END}];
        my $id     = $f[$col{ID}];
        my $svtype = $f[$col{SVTYPE}];
        my $svlen  = exists $col{SVLEN} ? ($f[$col{SVLEN}] // '.') : '.';
        my $refseq = exists $col{REFseq} ? ($f[$col{REFseq}] // '.') : '.';
        my $altseq = exists $col{ALTseq} ? ($f[$col{ALTseq}] // '.') : '.';

        die "Bad BED coordinates at line $line_no: $line\n"
            unless defined $chr && defined $start && defined $end && $start =~ /^\d+$/ && $end =~ /^\d+$/ && $end > $start;

        push @variants, {
            line   => $line,
            fields => \@f,
            chr    => $chr,
            start  => $start + 0,
            end    => $end + 0,
            id     => defined($id) && $id ne '' ? $id : ($chr . '_' . $start . '_' . $end),
            svtype => $svtype,
            svlen  => defined($svlen) && $svlen ne '' ? $svlen : ($end - $start),
            refseq => defined($refseq) && $refseq ne '' ? $refseq : '.',
            altseq => defined($altseq) && $altseq ne '' ? $altseq : '.',
        };
    }

    close $BED;

    if (!@header) {
        @header = qw(CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq);
    }

    return (\@variants, \@header, \%col);
}

sub split_line {
    my ($line) = @_;
    my @f = split /\t/, $line, -1;
    @f = split /\s+/, $line if @f < 3;
    return @f;
}

sub discover_maf_jobs {
    my ($opt_ref) = @_;
    my @jobs;
    my %seen;

    my $maf_dir = $opt_ref->{maf_dir};
    $maf_dir =~ s/\/$//;
    my $target = $opt_ref->{target};
    my $suffix = $opt_ref->{maf_suffix};

    if ($opt_ref->{query_list}) {
        open my $QL, '<', $opt_ref->{query_list} or die "Cannot open $opt_ref->{query_list}: $!\n";
        while (my $q = <$QL>) {
            chomp $q;
            $q =~ s/^\s+|\s+$//g;
            next if $q eq '' || $q =~ /^#/;

            # Allow either query names such as S001, or full MAF basenames such as S184.S001.sing.maf.
            my ($query, $maf);
            if ($q =~ /^\Q$target\E\.(.+)\Q$suffix\E$/) {
                $query = $1;
                $maf = "$maf_dir/$q";
            } else {
                $query = $q;
                $maf = "$maf_dir/$target.$query$suffix";
            }

            die "MAF file for query '$query' does not exist or is empty: $maf\n" unless -s $maf;
            die "Duplicated query name: $query\n" if $seen{$query}++;
            push @jobs, { sample => $query, maf => $maf };
        }
        close $QL;
    } else {
        opendir my $DIR, $maf_dir or die "Cannot open --maf_dir $maf_dir: $!\n";
        my @files = sort readdir($DIR);
        closedir $DIR;

        foreach my $base (@files) {
            next if $base =~ /^\./;
            my $maf = "$maf_dir/$base";
            next unless -f $maf && -s $maf;

            my $query = parse_query_from_maf_name($base, $target, $suffix);
            next unless defined $query && $query ne '';

            die "Duplicated query name: $query\n" if $seen{$query}++;
            push @jobs, { sample => $query, maf => $maf };
        }

        if (!@jobs) {
            my @examples = grep { /\.maf(?:\.gz)?$/ || /sing/ } @files;
            @examples = @examples[0 .. 9] if @examples > 10;
            my $expected = "$target.<query>$suffix";
            my $seen_msg = @examples ? join(", ", @examples) : "no maf-like files in directory";
            die "No MAF files found in $maf_dir. Expected names like: $expected\n" .
                "Examples seen: $seen_msg\n" .
                "Please check --target and --maf_suffix. For files such as S184.S001.sing.maf, use --target S184 --maf_suffix .sing.maf\n";
        }
    }

    return @jobs;
}

sub parse_query_from_maf_name {
    my ($base, $target, $suffix) = @_;
    if ($base =~ /^\Q$target\E\.(.+)\Q$suffix\E$/) {
        return $1;
    }
    return undef;
}

sub run_parallel_genotyping {
    my ($jobs_ref, $tmp_dir, $opt_ref) = @_;

    my %children;
    my $active = 0;
    my $failed = 0;

    foreach my $job (@$jobs_ref) {
        while ($active >= $opt_ref->{threads}) {
            $failed += wait_one_child(\%children, \$active);
        }

        my $pid = fork();
        die "Cannot fork: $!\n" unless defined $pid;

        if ($pid == 0) {
            my $sample = $job->{sample};
            my $out_gt = "$tmp_dir/$sample.gt";
            eval {
                genotype_one_query(
                    $opt_ref->{bed},
                    $job->{maf},
                    $out_gt,
                    $opt_ref->{min_cov},
                    $opt_ref->{min_support},
                    $opt_ref->{bridge_del_ratio},
                    $opt_ref->{bridge_ref_low},
                    $opt_ref->{bridge_ref_high},
                );
            };
            if ($@) {
                print STDERR "[ERROR] sample=$sample maf=$job->{maf}\n$@\n";
                exit 255;
            }
            exit 0;
        } else {
            $children{$pid} = $job->{sample};
            $active++;
        }
    }

    while ($active > 0) {
        $failed += wait_one_child(\%children, \$active);
    }

    die "$failed child process(es) failed. Check STDERR above.\n" if $failed;
}

sub wait_one_child {
    my ($children_ref, $active_ref) = @_;
    my $pid = wait();
    return 0 if $pid < 0;
    $$active_ref--;
    my $sample = exists $children_ref->{$pid} ? $children_ref->{$pid} : 'unknown';
    delete $children_ref->{$pid};

    my $status = $?;
    if (WIFEXITED($status) && WEXITSTATUS($status) == 0) {
        print STDERR "[", scalar localtime, "] Finished $sample\n";
        return 0;
    } else {
        my $code = WIFEXITED($status) ? WEXITSTATUS($status) : $status;
        print STDERR "[", scalar localtime, "] FAILED $sample, status=$code\n";
        return 1;
    }
}

sub genotype_one_query {
    my ($bed_file, $maf_file, $out_gt, $min_cov, $min_support, $bridge_del_ratio, $bridge_ref_low, $bridge_ref_high) = @_;

    my %intervals;
    my %blocks_by_ref_chr;

    read_maf($maf_file, \%intervals, \%blocks_by_ref_chr);

    add_between_block_intervals(
        \%intervals,
        \%blocks_by_ref_chr,
        $bridge_del_ratio,
        $bridge_ref_low,
        $bridge_ref_high,
    );

    foreach my $chr (keys %intervals) {
        @{$intervals{$chr}} = sort {
            $a->[0] <=> $b->[0] ||
            $a->[1] <=> $b->[1]
        } @{$intervals{$chr}};
    }

    open my $BED, '<', $bed_file or die "Cannot open $bed_file: $!\n";
    open my $OUT, '>', $out_gt   or die "Cannot write $out_gt: $!\n";

    my %col;
    my $line_no = 0;

    while (my $line = <$BED>) {
        $line_no++;
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split_line($line);

        if ($line =~ /^#/) {
            $f[0] =~ s/^#//;
            for my $i (0 .. $#f) {
                $col{$f[$i]} = $i;
            }
            next;
        }

        if (!%col) {
            my @default = qw(CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq);
            for my $i (0 .. $#default) {
                $col{$default[$i]} = $i;
            }
        }

        my $chr    = $f[$col{CHROM}];
        my $start  = $f[$col{POS}];
        my $end    = $f[$col{END}];
        my $svtype = exists $col{SVTYPE} ? ($f[$col{SVTYPE}] // 'DEL') : 'DEL';

        my $gt = '.';
        if (defined $svtype && uc($svtype) eq 'DEL') {
            $gt = genotype_one_del(
                $chr,
                $start,
                $end,
                \%intervals,
                $min_cov,
                $min_support,
            );
        }

        print $OUT "$gt\n";
    }

    close $OUT;
    close $BED;
}

sub read_maf {
    my ($maf_file, $intervals_ref, $blocks_ref) = @_;

    open my $MAF, '<', $maf_file or die "Cannot open $maf_file: $!\n";

    my @block;

    while (my $line = <$MAF>) {
        chomp $line;

        if ($line =~ /^\s*$/) {
            process_maf_block(\@block, $intervals_ref, $blocks_ref) if @block;
            @block = ();
        } else {
            push @block, $line;
        }
    }

    process_maf_block(\@block, $intervals_ref, $blocks_ref) if @block;

    close $MAF;
}

sub process_maf_block {
    my ($block_ref, $intervals_ref, $blocks_ref) = @_;

    my @s_lines;

    foreach my $line (@$block_ref) {
        push @s_lines, $line if $line =~ /^s\s+/;
    }

    return unless @s_lines >= 2;

    my $ref = parse_s_line($s_lines[0]);
    my $qry = parse_s_line($s_lines[1]);

    my $ref_seq = $ref->{seq};
    my $qry_seq = $qry->{seq};

    my $aln_len = length($ref_seq);

    die "Alignment length differs in one MAF block\n"
        unless length($qry_seq) == $aln_len;

    my ($ref_fwd_start, $ref_fwd_end) = fwd_interval($ref);
    my ($qry_fwd_start, $qry_fwd_end) = fwd_interval($qry);

    push @{$blocks_ref->{$ref->{src}}}, {
        ref_chr    => $ref->{src},
        ref_start  => $ref_fwd_start,
        ref_end    => $ref_fwd_end,
        qry_chr    => $qry->{src},
        qry_start  => $qry_fwd_start,
        qry_end    => $qry_fwd_end,
        qry_strand => $qry->{strand},
    };

    my $ref_pos  = initial_pos($ref);
    my $ref_step = ($ref->{strand} eq '+') ? 1 : -1;

    my $run_status = '';
    my $run_first;
    my $run_last;

    for (my $i = 0; $i < $aln_len; $i++) {
        my $rb = substr($ref_seq, $i, 1);
        my $qb = substr($qry_seq, $i, 1);

        next if $rb eq '-';

        my $coord = $ref_pos;
        $ref_pos += $ref_step;

        my $status = ($qb eq '-') ? 'DEL' : 'REF';

        if ($run_status eq '') {
            $run_status = $status;
            $run_first  = $coord;
            $run_last   = $coord;
        }
        elsif ($status eq $run_status &&
               (($ref_step == 1  && $coord == $run_last + 1) ||
                ($ref_step == -1 && $coord == $run_last - 1))) {
            $run_last = $coord;
        }
        else {
            add_interval($intervals_ref, $ref->{src}, $run_first, $run_last, $run_status);

            $run_status = $status;
            $run_first  = $coord;
            $run_last   = $coord;
        }
    }

    if ($run_status ne '') {
        add_interval($intervals_ref, $ref->{src}, $run_first, $run_last, $run_status);
    }
}

sub add_between_block_intervals {
    my ($intervals_ref, $blocks_ref, $bridge_del_ratio, $bridge_ref_low, $bridge_ref_high) = @_;

    foreach my $chr (keys %$blocks_ref) {
        my @blocks = sort {
            $a->{ref_start} <=> $b->{ref_start} ||
            $a->{ref_end}   <=> $b->{ref_end}
        } @{$blocks_ref->{$chr}};

        for (my $i = 0; $i < @blocks - 1; $i++) {
            my $b1 = $blocks[$i];
            my $b2 = $blocks[$i + 1];

            my $ref_gap_start = $b1->{ref_end};
            my $ref_gap_end   = $b2->{ref_start};

            next if $ref_gap_end <= $ref_gap_start;

            my $ref_gap_len = $ref_gap_end - $ref_gap_start;
            my $status = 'UNC';

            if ($b1->{qry_chr} eq $b2->{qry_chr} &&
                $b1->{qry_strand} eq $b2->{qry_strand}) {

                my $qry_gap_len;

                if ($b1->{qry_strand} eq '+') {
                    $qry_gap_len = $b2->{qry_start} - $b1->{qry_end};
                } else {
                    $qry_gap_len = $b1->{qry_start} - $b2->{qry_end};
                }

                if (defined $qry_gap_len && $qry_gap_len >= 0) {
                    my $ratio = $qry_gap_len / $ref_gap_len;

                    if ($ratio <= $bridge_del_ratio) {
                        $status = 'DEL';
                    }
                    elsif ($ratio >= $bridge_ref_low && $ratio <= $bridge_ref_high) {
                        $status = 'REF';
                    }
                    else {
                        $status = 'UNC';
                    }
                }
            }

            push @{$intervals_ref->{$chr}}, [
                $ref_gap_start,
                $ref_gap_end,
                $status,
            ];
        }
    }
}

sub genotype_one_del {
    my ($chr, $start, $end, $intervals_ref, $min_cov, $min_support) = @_;

    return '.' unless defined $chr && defined $start && defined $end;
    return '.' unless $start =~ /^\d+$/ && $end =~ /^\d+$/;

    my $len = $end - $start;
    return '.' if $len <= 0;
    return '.' unless exists $intervals_ref->{$chr};

    my $ref_like  = 0;
    my $del_like  = 0;
    my $uncertain = 0;

    foreach my $iv (@{$intervals_ref->{$chr}}) {
        my ($iv_s, $iv_e, $status) = @$iv;

        last if $iv_s >= $end;
        next if $iv_e <= $start;

        my $ov_s = $iv_s > $start ? $iv_s : $start;
        my $ov_e = $iv_e < $end   ? $iv_e : $end;

        next if $ov_e <= $ov_s;

        my $ov_len = $ov_e - $ov_s;

        if ($status eq 'DEL') {
            $del_like += $ov_len;
        }
        elsif ($status eq 'REF') {
            $ref_like += $ov_len;
        }
        else {
            $uncertain += $ov_len;
        }
    }

    my $confident = $ref_like + $del_like;
    return '.' if $confident <= 0;
    return '.' if $confident / $len < $min_cov;

    if ($del_like / $confident >= $min_support) {
        return '1';
    }
    elsif ($ref_like / $confident >= $min_support) {
        return '0';
    }
    else {
        return '.';
    }
}

sub parse_s_line {
    my ($line) = @_;

    my @f = split /\s+/, $line;

    die "Bad MAF s line:\n$line\n" unless @f >= 7;

    return {
        src     => $f[1],
        start   => $f[2] + 0,
        size    => $f[3] + 0,
        strand  => $f[4],
        srcSize => $f[5] + 0,
        seq     => $f[6],
    };
}

sub fwd_interval {
    my ($rec) = @_;

    if ($rec->{strand} eq '+') {
        return ($rec->{start}, $rec->{start} + $rec->{size});
    } else {
        my $s = $rec->{srcSize} - $rec->{start} - $rec->{size};
        my $e = $rec->{srcSize} - $rec->{start};
        return ($s, $e);
    }
}

sub initial_pos {
    my ($rec) = @_;

    if ($rec->{strand} eq '+') {
        return $rec->{start};
    } else {
        return $rec->{srcSize} - $rec->{start} - 1;
    }
}

sub add_interval {
    my ($intervals_ref, $chr, $first, $last, $status) = @_;

    my $s = $first < $last ? $first : $last;
    my $e = $first > $last ? $first : $last;

    $e = $e + 1;

    push @{$intervals_ref->{$chr}}, [$s, $e, $status];
}

sub combine_gt_tsv {
    my ($variants_ref, $header_ref, $jobs_ref, $tmp_dir, $out_file) = @_;

    my @samples = map { $_->{sample} } @$jobs_ref;
    my @gt_arrays;

    foreach my $sample (@samples) {
        my $file = "$tmp_dir/$sample.gt";
        open my $GT, '<', $file or die "Cannot open $file: $!\n";
        my @gt = <$GT>;
        chomp @gt;
        close $GT;
        die "GT count mismatch for $sample: expected " . scalar(@$variants_ref) . ", got " . scalar(@gt) . "\n"
            unless @gt == @$variants_ref;
        push @gt_arrays, \@gt;
    }

    open my $OUT, '>', $out_file or die "Cannot write $out_file: $!\n";
    print $OUT '#', join("\t", @$header_ref, @samples), "\n";

    for (my $i = 0; $i < @$variants_ref; $i++) {
        my @row = split_line($variants_ref->[$i]{line});
        my @gts = map { $gt_arrays[$_]->[$i] } 0 .. $#gt_arrays;
        print $OUT join("\t", @row, @gts), "\n";
    }

    close $OUT;
}

sub write_vcf {
    my ($variants_ref, $jobs_ref, $tmp_dir, $out_file, $opt_ref, $ref_db) = @_;

    my @samples = map { $_->{sample} } @$jobs_ref;
    my @gt_arrays;

    foreach my $sample (@samples) {
        my $file = "$tmp_dir/$sample.gt";
        open my $GT, '<', $file or die "Cannot open $file: $!\n";
        my @gt = <$GT>;
        chomp @gt;
        close $GT;
        die "GT count mismatch for $sample: expected " . scalar(@$variants_ref) . ", got " . scalar(@gt) . "\n"
            unless @gt == @$variants_ref;
        push @gt_arrays, \@gt;
    }

    open my $VCF, '>', $out_file or die "Cannot write $out_file: $!\n";

    print $VCF "##fileformat=VCFv4.2\n";
    print $VCF "##fileDate=", strftime('%Y%m%d', localtime), "\n";
    print $VCF "##source=genotype_DEL_parallel_from_maf.pl\n";
    print $VCF "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate converted from BED 0-based half-open END\">\n";
    print $VCF "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    print $VCF "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n";
    print $VCF "##INFO=<ID=BED_START,Number=1,Type=Integer,Description=\"Original 0-based BED POS column\">\n";
    print $VCF "##INFO=<ID=ALLELE_SOURCE,Number=1,Type=String,Description=\"How REF/ALT were generated\">\n";
    print $VCF "##ALT=<ID=DEL,Description=\"Deletion\">\n";
    print $VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype: 0 reference-like, 1 deletion-like, . uncertain\">\n";
    print $VCF '#CHROM', "\t", 'POS', "\t", 'ID', "\t", 'REF', "\t", 'ALT', "\t", 'QUAL', "\t", 'FILTER', "\t", 'INFO', "\t", 'FORMAT';
    print $VCF "\t", join("\t", @samples) if @samples;
    print $VCF "\n";

    for (my $i = 0; $i < @$variants_ref; $i++) {
        my $v = $variants_ref->[$i];
        next unless uc($v->{svtype}) eq 'DEL';

        my $chrom = $v->{chr};
        my $id    = sanitize_vcf_field($v->{id});
        my $len   = $v->{end} - $v->{start};
        my ($ref,$alt,$pos,$end,$svlen,$allele_source);
        if ($ref_db) {
            ($ref,$alt,$pos,$end,$svlen,$allele_source) = canonical_del_alleles($v,$ref_db);
        }
        if (!defined $ref) {
            # Correct symbolic DEL convention: POS is the left anchor and
            # deleted sequence is (POS, END]. Thus END-POS=abs(SVLEN).
            $pos = $v->{start} > 0 ? $v->{start} : 1;
            $end = $v->{start} > 0 ? $v->{end} : $v->{end} + 1;
            $svlen = -abs($len);
            $ref = $ref_db ? $ref_db->get_seq($chrom,$pos-1,$pos) : 'N';
            $ref = 'N' unless defined($ref) && $ref =~ /^[ACGTN]$/i;
            $alt = '<DEL>';
            $allele_source = 'SYMBOLIC_FALLBACK';
        }

        my $info = join(';',
            'END=' . $end,
            'SVTYPE=DEL',
            'SVLEN=' . $svlen,
            'BED_START=' . $v->{start},
            'ALLELE_SOURCE=' . $allele_source,
        );

        my @vcf_gts = map { gt_to_vcf($gt_arrays[$_]->[$i], $opt_ref->{diploid}) } 0 .. $#gt_arrays;

        print $VCF join("\t", $chrom, $pos, $id, $ref, $alt, '.', 'PASS', $info, 'GT', @vcf_gts), "\n";
    }

    close $VCF;
}

sub canonical_del_alleles {
    my ($v,$ref_db)=@_;
    my ($chrom,$start0,$end0)=@{$v}{qw(chr start end)};
    return unless $ref_db->has_contig($chrom);
    return if $start0<0 || $end0<$start0 || $end0>$ref_db->length($chrom);
    my $event_ref=$ref_db->get_seq($chrom,$start0,$end0);
    return unless defined($event_ref) && $event_ref =~ /^[ACGTN]+$/i;
    my $raw_ref=clean_del_seq($v->{refseq});
    my $raw_alt=clean_del_seq($v->{altseq});
    # Legacy PairwiseSV may include one common right-padding base.
    if ($raw_ref ne '' && length($raw_ref)==length($event_ref)+1 && length($raw_alt)>=1
        && substr($raw_ref,0,length($event_ref)) eq $event_ref
        && substr($raw_ref,-1) eq substr($raw_alt,-1)) {
        $raw_ref=substr($raw_ref,0,-1);
        $raw_alt=substr($raw_alt,0,-1);
    }
    my $event_alt=$raw_alt;
    return unless $event_alt =~ /^[ACGTN]*$/i;
    return unless length($event_ref)>length($event_alt);
    my($pos,$ref,$alt,$source);
    if($start0>0){
        my$anchor=$ref_db->get_seq($chrom,$start0-1,$start0);return unless defined($anchor)&&$anchor=~/^[ACGTN]$/i;
        $pos=$start0;$ref=$anchor.$event_ref;$alt=$anchor.$event_alt;$source='REFERENCE_LEFT_ANCHORED';
    }else{
        return if $end0 >= $ref_db->length($chrom);
        my$anchor=$ref_db->get_seq($chrom,$end0,$end0+1);return unless defined($anchor)&&$anchor=~/^[ACGTN]$/i;
        $pos=1;$ref=$event_ref.$anchor;$alt=$event_alt.$anchor;$source='REFERENCE_RIGHT_ANCHORED';
    }
    my$end=$pos+length($ref)-1;my$svlen=length($alt)-length($ref);
    return unless $svlen<0;
    return($ref,$alt,$pos,$end,$svlen,$source);
}
sub clean_del_seq { my($s)=@_;$s='' unless defined$s;$s=uc($s);$s=~s/\s+//g;return '' if $s eq'.'||$s eq'-'||$s eq'NA';$s=~s/-//g;return $s; }

sub vcf_alleles {
    my ($v, $mode) = @_;

    if ($mode eq 'sequence') {
        my $ref = uc($v->{refseq} // '.');
        my $alt = uc($v->{altseq} // '.');
        if (valid_vcf_sequence($ref) && valid_vcf_sequence($alt)) {
            return ($ref, $alt);
        }
    }

    return ('N', '<DEL>');
}

sub valid_vcf_sequence {
    my ($seq) = @_;
    return 0 unless defined $seq;
    return 0 if $seq eq '' || $seq eq '.' || $seq eq '-';
    return $seq =~ /^[ACGTNacgtn]+$/ ? 1 : 0;
}

sub gt_to_vcf {
    my ($gt, $diploid) = @_;
    if (!defined $gt || $gt eq '.') {
        return $diploid ? './.' : '.';
    }
    elsif ($gt eq '0') {
        return $diploid ? '0/0' : '0';
    }
    elsif ($gt eq '1') {
        return $diploid ? '1/1' : '1';
    }
    else {
        return $diploid ? './.' : '.';
    }
}

sub sanitize_vcf_field {
    my ($x) = @_;
    $x = '.' unless defined $x && $x ne '';
    $x =~ s/[\s;=]/_/g;
    return $x;
}


package DELRefFAI;
use strict; use warnings; use Fcntl qw(SEEK_SET);
sub new{my($c,$f)=@_;die "Reference FASTA must be uncompressed: $f\n" if$f=~/\.gz$/i;my$s=bless{fasta=>$f,idx=>{}},$c;if(-s"$f.fai"){$s->_load("$f.fai")}else{$s->_build}open(my$h,'<',$f)or die"Cannot open $f: $!\n";binmode$h;$s->{fh}=$h;return$s;}
sub _load{my($s,$f)=@_;open(my$h,'<',$f)or die"Cannot read $f: $!\n";while(<$h>){chomp;my@x=split/\t/;next unless@x>=5;$s->{idx}{$x[0]}={len=>$x[1]+0,offset=>$x[2]+0,line_bases=>$x[3]+0,line_width=>$x[4]+0};}close$h;}
sub _build{my($s)=@_;open(my$h,'<',$s->{fasta})or die$!;binmode$h;my($n,$l,$o,$b,$w);while(1){my$p=tell$h;my$x=<$h>;last unless defined$x;if($x=~/^>(\S+)/){$s->{idx}{$n}={len=>$l,offset=>$o,line_bases=>$b,line_width=>$w}if defined$n;$n=$1;($l,$o,$b,$w)=(0,undef,undef,undef);}else{next unless defined$n;my$rw=length$x;$x=~s/[\r\n]+$//;my$z=length$x;next unless$z;$o=$p unless defined$o;$b=$z unless defined$b;$w=$rw unless defined$w;$l+=$z;}}$s->{idx}{$n}={len=>$l,offset=>$o,line_bases=>$b,line_width=>$w}if defined$n;close$h;}
sub has_contig{return exists$_[0]{idx}{$_[1]}}sub length{return$_[0]{idx}{$_[1]}{len}}
sub get_seq{my($s,$c,$a,$b)=@_;return undef unless$s->has_contig($c);my$x=$s->{idx}{$c};$a=0 if$a<0;$b=$x->{len}if$b>$x->{len};return''if$b<=$a;my$n=$b-$a;my$p=$x->{offset}+int($a/$x->{line_bases})*$x->{line_width}+($a%$x->{line_bases});seek($s->{fh},$p,SEEK_SET)or die"FASTA seek failed\n";my$q='';while(CORE::length($q)<$n){my$u='';my$r=$n-CORE::length($q);my$want=$r+int($r/$x->{line_bases}+2)*($x->{line_width}-$x->{line_bases});my$g=read($s->{fh},$u,$want);last unless$g;$u=~s/[\r\n]//g;$q.=$u;}return uc substr($q,0,$n)}
1;
