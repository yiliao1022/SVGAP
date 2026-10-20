#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
# ============================================================================
# SVGAP standardized script: 07_genotype_INS.pl
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
use Cwd qw(abs_path);
use POSIX qw(strftime WIFEXITED WEXITSTATUS);

# ------------------------------------------------------------
# genotype_INS_parallel_onefile.v1.pl
#
# One-file parallel insertion genotyper.
# It integrates the single-query logic from
# genotype_INS_from_maf_queryfa_run.v6_shortref_semiglobal.pl
# and the multi-query batch/VCF workflow into one script.
#
# MAF filename rule by default:
#   <target>.<query>.sing.maf
# Example:
#   S184.S001.sing.maf
#
# Query genome FASTA filename search order:
#   genome_dir/<query>
#   genome_dir/<query>.fa/.fasta/.fna/.fas
#   or one unique genome_dir/<query>*.(fa|fasta|fna|fas)
#
# Query list:
#   one query name per line, e.g. S001
#   also accepts full MAF basename, e.g. S184.S001.sing.maf
#
# Genotype meaning:
#   0 = query is reference-like, no insertion
#   1 = query carries the insertion relative to reference
#   . = uncertain
# ------------------------------------------------------------

my %opt = (
    threads          => 4,
    target           => '',
    maf_suffix       => '.sing.maf',
    pos_tol          => 100,
    seq_id           => 0.80,
    flank            => 100,
    len_ratio_min    => 0.80,
    len_ratio_max    => 1.25,
    len_slop         => 10,
    min_refcov       => 0.80,
    block_extra      => 20,
    short_ref_len    => 20,
    long_seq_len     => 1000,
    sample_len       => 1000,
    sample_seed_len  => 80,
    sample_min_hits  => 1,
    sample_windows   => 3,
    vcf_mode         => 'allele',  # allele or symbolic
    diploid          => 1,
    bed_one_based    => 0,
    alt_is_insseq    => 0,
    keep_tmp         => 0,
    ref_fasta        => '',       # reference FASTA used to build canonical explicit INS alleles
);

GetOptions(
    'bed=s'                => \$opt{bed},
    'maf_dir|maf-dir=s'    => \$opt{maf_dir},
    'genome_dir|genome-dir=s' => \$opt{genome_dir},
    'ref_fasta|ref-fasta|ref=s' => \$opt{ref_fasta},
    'target|refname=s'             => \$opt{target},
    'query_list|query-list=s' => \$opt{query_list},
    'out=s'                => \$opt{out},
    'threads|t=i'            => \$opt{threads},
    'maf_suffix|maf-suffix=s' => \$opt{maf_suffix},

    'pos_tol|pos-tol=i'    => \$opt{pos_tol},
    'seq_id|seq-id=f'      => \$opt{seq_id},
    'flank=i'              => \$opt{flank},
    'len_ratio_min|len-ratio-min=f' => \$opt{len_ratio_min},
    'len_ratio_max|len-ratio-max=f' => \$opt{len_ratio_max},
    'len_slop|len-slop=i'  => \$opt{len_slop},
    'min_refcov|min-refcov=f' => \$opt{min_refcov},
    'block_extra|block-extra=i' => \$opt{block_extra},
    'short_ref_len|short-ref-len=i' => \$opt{short_ref_len},

    'long_seq_len|long-seq-len=i' => \$opt{long_seq_len},
    'sample_len|sample-len=i' => \$opt{sample_len},
    'sample_seed_len|sample-seed-len=i' => \$opt{sample_seed_len},
    'sample_min_hits|sample-min-hits=i' => \$opt{sample_min_hits},
    'sample_windows|sample-windows=i' => \$opt{sample_windows},

    'vcf_mode|vcf-mode=s'  => \$opt{vcf_mode},
    'haploid!'             => sub { $opt{diploid} = 0 },
    'diploid!'             => \$opt{diploid},
    'bed_one_based|bed-one-based!' => \$opt{bed_one_based},
    'alt_is_insseq|alt-is-insseq!' => \$opt{alt_is_insseq},
    'keep_tmp|keep-tmp!'   => \$opt{keep_tmp},
    'help|h'               => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

die usage() unless $opt{bed} && $opt{maf_dir} && $opt{genome_dir} && $opt{target} && $opt{out};
die "--threads must be >= 1\n" unless $opt{threads} >= 1;
die "--vcf_mode must be allele or symbolic\n" unless $opt{vcf_mode} eq 'allele' || $opt{vcf_mode} eq 'symbolic';
die "Cannot find BED file: $opt{bed}\n" unless -s $opt{bed};
die "Cannot find MAF directory: $opt{maf_dir}\n" unless -d $opt{maf_dir};
die "Cannot find genome directory: $opt{genome_dir}\n" unless -d $opt{genome_dir};

my ($variants_ref, $header_ref) = read_bed_variants($opt{bed});
my @variants = @$variants_ref;
die "No variants were loaded from $opt{bed}\n" unless @variants;

# Resolve the reference FASTA independently of the query FASTA files.  It is
# required for a real VCF anchor base and for rebuilding REF from the reference
# assembly.  This prevents loss of the insertion sequence when PairwiseSV emits
# a pure inserted sequence with an empty reference-side interval.
$opt{ref_fasta} ||= find_query_fasta($opt{genome_dir}, $opt{target});
die "Cannot find reference FASTA for $opt{target}; provide --ref_fasta FILE\n"
    unless defined($opt{ref_fasta}) && -s $opt{ref_fasta};
my $ref_db = INSRefFAI->new($opt{ref_fasta});

my @jobs = discover_jobs(\%opt);
die "No valid jobs found. Check --maf_dir, --genome_dir, --target, --query_list and --maf_suffix.\n" unless @jobs;

my $tmp_dir = "$opt{out}.tmp.INS_genotype.$$";
make_path($tmp_dir) or die "Cannot create $tmp_dir: $!\n";

print STDERR "[", scalar localtime, "] Loaded ", scalar(@variants), " variants from $opt{bed}\n";
print STDERR "[", scalar localtime, "] Found ", scalar(@jobs), " query jobs\n";
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
  perl genotype_INS_parallel_onefile.v1.pl \
    --bed INS.bed \
    --maf_dir sing_maf_dir \
    --genome_dir query_genome_dir \
    --target S184 \
    --out INS.genotyped \
    --threads 8

Optional:
  --query_list query.list     or --glist g.list        one query name per line; sample name is the query name only
  --ref_fasta reference.fa       reference FASTA for canonical REF/ALT; auto-detected from --genome_dir/--target
  --maf_suffix .sing.maf         default: .sing.maf; MAF name is <target>.<query><maf_suffix>

Parameters from the original single-query INS genotyper:
  --pos_tol 100                  breakpoint tolerance around the BED insertion position
  --seq_id 0.80                  minimum sequence identity for insertion-sequence evidence
  --flank 100                    flank size used to call reference-like sites
  --len_ratio_min 0.80           minimum candidate/expected insertion length ratio
  --len_ratio_max 1.25           maximum candidate/expected insertion length ratio
  --len_slop 10                  absolute length tolerance in bp
  --min_refcov 0.80              minimum reference-like coverage to call GT=0
  --block_extra 20               query gap must exceed reference gap by this many bp for between-block INS
  --short_ref_len 20             allow short REF alleles to be treated as insertion-like

Fast long-insertion sequence check:
  --long_seq_len 1000            use seed-based check for expected sequence >= this length
  --sample_len 1000              sampled segment length for long sequences
  --sample_seed_len 80           exact seed length for long sequences
  --sample_min_hits 1            minimum seed hits required
  --sample_windows 3             number of deterministic sequence windows to sample

VCF options:
  --vcf_mode allele              allele: build reference-checked explicit REF/ALT [default]; symbolic: emit <INS> fallback
  --haploid                      output haploid GT values: 0, 1, .
  --diploid                      output diploid GT values: 0/0, 1/1, ./. [default]
  --bed_one_based                input POS/END are already 1-based; default input is 0-based half-open
  --alt_is_insseq                in allele mode, treat ALTseq/query sequence as pure inserted sequence and set ALT=REF+ALTseq
  --keep_tmp                     keep temporary per-query GT files

Accepted BED formats:
  1) Header-based SVGAP-like TSV:
     #CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq

  2) Original six-column INS BED used by genotype_INS_from_maf_queryfa_run.v6_shortref_semiglobal.pl:
     chr start end size ref_seq query_seq

Output:
  <out>.gt.tsv   original BED columns plus one GT column per query
  <out>.vcf      multi-sample VCF with GT calls
USAGE
}

# -----------------------------
# BED parsing
# -----------------------------
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
            for my $i (0 .. $#header) { $col{$header[$i]} = $i; }
            next;
        }

        my %v;
        if (%col && exists $col{CHROM} && exists $col{POS} && exists $col{END}) {
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
                unless defined $chr && defined $start && defined $end && $start =~ /^\d+$/ && $end =~ /^\d+$/;

            %v = (
                line   => $line,
                fields => \@f,
                chr    => $chr,
                start  => $start + 0,
                end    => $end + 0,
                id     => defined($id) && $id ne '' ? $id : ($chr . '_' . $start . '_' . $end),
                svtype => defined($svtype) && $svtype ne '' ? uc($svtype) : 'INS',
                svlen  => defined($svlen) && $svlen ne '' ? $svlen : '.',
                refseq => defined($refseq) && $refseq ne '' ? $refseq : '.',
                altseq => defined($altseq) && $altseq ne '' ? $altseq : '.',
                bedfmt => 'header',
            );
        } else {
            # Original six-column format: chr start end size ref_seq query_seq
            die "Bad six-column INS BED line $line_no: $line\n" unless @f >= 6;
            my ($chr,$start,$end,$size,$refseq,$altseq) = @f[0..5];
            die "Bad BED coordinates at line $line_no: $line\n"
                unless defined $chr && defined $start && defined $end && $start =~ /^\d+$/ && $end =~ /^\d+$/;
            if (!@header) {
                @header = qw(CHROM POS END SVLEN REFseq ALTseq);
            }
            %v = (
                line   => $line,
                fields => \@f,
                chr    => $chr,
                start  => $start + 0,
                end    => $end + 0,
                id     => $chr . '_' . $start . '_' . $end,
                svtype => 'INS',
                svlen  => $size,
                refseq => $refseq,
                altseq => $altseq,
                bedfmt => 'sixcol',
            );
        }

        # Retain all rows in the GT table, but non-INS rows become missing in the VCF/genotype.
        $v{expected} = derive_expected_insertion(\%v, \%opt);
        push @variants, \%v;
    }
    close $BED;

    @header = qw(CHROM POS END SVLEN REFseq ALTseq) unless @header;
    return (\@variants, \@header);
}

sub derive_expected_insertion {
    my ($v, $opt_ref) = @_;
    my $ref = clean_seq($v->{refseq});
    my $alt = clean_seq($v->{altseq});

    my $ref_len = length($ref);
    my $alt_len = length($alt);
    my $expected_len;

    if ($opt_ref->{alt_is_insseq}) {
        $expected_len = $alt_len;
    } else {
        $expected_len = $alt_len - $ref_len;
    }

    if (!defined($expected_len) || $expected_len <= 0) {
        if (defined $v->{svlen} && $v->{svlen} =~ /^-?\d+$/ && abs($v->{svlen}) > 0) {
            $expected_len = abs($v->{svlen});
        } elsif ($alt_len > 0) {
            $expected_len = $alt_len;
        } else {
            $expected_len = ($v->{end} >= $v->{start}) ? ($v->{end} - $v->{start}) : 0;
        }
    }

    my $expected_seq = $alt;
    if ($opt_ref->{alt_is_insseq}) {
        $expected_seq = $alt;
    }

    return {
        refseq      => $ref,
        altseq      => $alt,
        len         => $expected_len,
        seq         => $expected_seq,
        ref_len     => $ref_len,
        alt_len     => $alt_len,
        pure_point  => ($ref_len <= 1) ? 1 : 0,
        short_ref   => ($ref_len <= $opt_ref->{short_ref_len}) ? 1 : 0,
    };
}

sub split_line {
    my ($line) = @_;
    my @f = split /\t/, $line, -1;
    @f = split /\s+/, $line if @f < 3;
    return @f;
}

# -----------------------------
# Job discovery
# -----------------------------
sub discover_jobs {
    my ($opt_ref) = @_;
    my $maf_dir = $opt_ref->{maf_dir};
    my $genome_dir = $opt_ref->{genome_dir};
    $maf_dir =~ s/\/$//;
    $genome_dir =~ s/\/$//;
    my $target = $opt_ref->{target};
    my $suffix = $opt_ref->{maf_suffix};

    my @jobs;
    my %seen;

    if ($opt_ref->{query_list}) {
        open my $QL, '<', $opt_ref->{query_list} or die "Cannot open $opt_ref->{query_list}: $!\n";
        while (my $q = <$QL>) {
            chomp $q;
            $q =~ s/^\s+|\s+$//g;
            next if $q eq '' || $q =~ /^#/;

            my ($query, $maf);
            if ($q =~ /^\Q$target\E\.(.+)\Q$suffix\E$/) {
                $query = $1;
                $maf = "$maf_dir/$q";
            } elsif ($q =~ /^(.+?)\.([^\.]+)\.sing\.maf$/) {
                $query = $2;
                $maf = "$maf_dir/$q";
            } elsif ($q =~ /^(.+?)\.([^\.]+)\.maf$/) {
                $query = $2;
                $maf = "$maf_dir/$q";
            } else {
                $query = $q;
                $maf = "$maf_dir/$target.$query$suffix";
            }

            die "Duplicated query name in --query_list: $query\n" if $seen{$query}++;
            die "MAF file for query '$query' does not exist or is empty: $maf\n" unless -s $maf;
            my $fa = find_query_fasta($genome_dir, $query);
            die "Cannot find query genome FASTA for '$query' in $genome_dir\n" unless defined $fa && -s $fa;
            push @jobs, { sample => $query, maf => $maf, fasta => $fa };
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
            my $fa = find_query_fasta($genome_dir, $query);
            if (!defined $fa) {
                warn "Skip $base: cannot find query genome FASTA for $query in $genome_dir\n";
                next;
            }
            push @jobs, { sample => $query, maf => $maf, fasta => $fa };
        }
        if (!@jobs) {
            my @examples = grep { /\.maf(?:\.gz)?$/ || /sing/ } @files;
            @examples = @examples[0 .. 9] if @examples > 10;
            my $expected = "$target.<query>$suffix";
            my $seen_msg = @examples ? join(", ", @examples) : "no maf-like files in directory";
            die "No MAF files found in $maf_dir. Expected names like: $expected\nExamples seen: $seen_msg\n";
        }
    }

    return @jobs;
}

sub parse_query_from_maf_name {
    my ($base, $target, $suffix) = @_;
    return $1 if $base =~ /^\Q$target\E\.(.+)\Q$suffix\E$/;
    return undef;
}

sub find_query_fasta {
    my ($genome_dir, $query) = @_;

    my $no_suffix = "$genome_dir/$query";
    return $no_suffix if -e $no_suffix && -f $no_suffix;

    foreach my $suf (qw(.fa .fasta .fna .fas)) {
        my $file = "$genome_dir/$query$suf";
        return $file if -e $file && -f $file;
    }

    opendir my $DH, $genome_dir or die "Cannot open $genome_dir: $!\n";
    my @candidates = grep {
        /^\Q$query\E(?:\.|_|-)?/ && /\.(fa|fasta|fna|fas)$/i && -f "$genome_dir/$_"
    } readdir($DH);
    closedir $DH;

    if (@candidates == 1) {
        return "$genome_dir/$candidates[0]";
    }
    if (@candidates > 1) {
        warn "Multiple FASTA candidates found for $query in $genome_dir. Please keep only one or use exact name '$query.fa':\n";
        warn "  " . join("\n  ", @candidates) . "\n";
    }
    return undef;
}

# -----------------------------
# Parallel execution
# -----------------------------
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
                genotype_one_query($job, $out_gt, $opt_ref);
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
    my ($job, $out_gt, $opt_ref) = @_;

    my %query_seq_by_chr = read_fasta_to_hash($job->{fasta});

    my %ins_candidates;
    my %ref_coverage;
    my %blocks_by_ref_chr;

    read_maf($job->{maf}, \%ins_candidates, \%ref_coverage, \%blocks_by_ref_chr);
    add_between_block_insertions(\%ins_candidates, \%ref_coverage, \%blocks_by_ref_chr, $opt_ref->{block_extra});

    foreach my $chr (keys %ins_candidates) {
        @{$ins_candidates{$chr}} = sort {
            $a->{bp_start} <=> $b->{bp_start} ||
            $a->{bp_end}   <=> $b->{bp_end}   ||
            $a->{len}      <=> $b->{len}
        } @{$ins_candidates{$chr}};
    }
    foreach my $chr (keys %ref_coverage) {
        @{$ref_coverage{$chr}} = merge_intervals(@{$ref_coverage{$chr}});
    }

    open my $OUT, '>', $out_gt or die "Cannot write $out_gt: $!\n";

    foreach my $v (@variants) {
        my $gt = '.';
        if (uc($v->{svtype}) eq 'INS') {
            my @res = genotype_one_ins(
                $v,
                \%ins_candidates,
                \%ref_coverage,
                \%query_seq_by_chr,
                $opt_ref,
            );
            $gt = $res[0];
        }
        print $OUT "$gt\n";
    }

    close $OUT;
}

# -----------------------------
# MAF parsing and INS candidates
# -----------------------------
sub read_fasta_to_hash {
    my ($fa) = @_;
    open my $IN, '<', $fa or die "Cannot open FASTA $fa: $!\n";
    my %h;
    my $id = '';
    while (my $line = <$IN>) {
        chomp $line;
        if ($line =~ /^>(\S+)/) {
            $id = $1;
            $h{$id} = '';
            next;
        }
        next if $line =~ /^\s*$/;
        die "FASTA sequence found before header in $fa\n" if $id eq '';
        $line =~ s/\s+//g;
        $h{$id} .= uc($line);
    }
    close $IN;
    return %h;
}

sub read_maf {
    my ($maf_file, $ins_ref, $cov_ref, $blocks_ref) = @_;
    open my $MAF, '<', $maf_file or die "Cannot open MAF $maf_file: $!\n";
    my @block;
    while (my $line = <$MAF>) {
        chomp $line;
        if ($line =~ /^\s*$/) {
            process_maf_block($ins_ref, $cov_ref, $blocks_ref, @block) if @block;
            @block = ();
        } else {
            push @block, $line;
        }
    }
    process_maf_block($ins_ref, $cov_ref, $blocks_ref, @block) if @block;
    close $MAF;
}

sub process_maf_block {
    my ($ins_ref, $cov_ref, $blocks_ref, @lines) = @_;
    my @s = grep { /^s\s+/ } @lines;
    return unless @s >= 2;

    my $ref = parse_s_line($s[0]);
    my $qry = parse_s_line($s[1]);
    my ($rseq, $qseq) = ($ref->{seq}, $qry->{seq});
    my $aln_len = length($rseq);
    die "Alignment length differs in one MAF block\n" unless length($qseq) == $aln_len;

    my ($rs, $re) = fwd_interval($ref);
    my ($qs, $qe) = fwd_interval($qry);

    push @{$blocks_ref->{$ref->{src}}}, {
        ref_chr    => $ref->{src},
        ref_start  => $rs,
        ref_end    => $re,
        qry_chr    => $qry->{src},
        qry_start  => $qs,
        qry_end    => $qe,
        qry_strand => $qry->{strand},
    };
    push @{$cov_ref->{$ref->{src}}}, [$rs, $re];

    my $rpos  = initial_pos($ref);
    my $rstep = ($ref->{strand} eq '+') ? 1 : -1;
    my $qpos  = initial_pos($qry);
    my $qstep = ($qry->{strand} eq '+') ? 1 : -1;

    my $i = 0;
    while ($i < $aln_len) {
        my $rb = substr($rseq, $i, 1);
        my $qb = substr($qseq, $i, 1);

        if ($rb eq '-' && $qb ne '-') {
            my $bp = $rpos;
            my $ins_seq = '';
            my @qcoords;
            while ($i < $aln_len) {
                $rb = substr($rseq, $i, 1);
                $qb = substr($qseq, $i, 1);
                last unless ($rb eq '-' && $qb ne '-');
                $ins_seq .= $qb;
                push @qcoords, $qpos;
                $qpos += $qstep;
                $i++;
            }
            my ($q_s, $q_e) = coord_range_half_open(@qcoords);
            if ($qry->{strand} eq '-') {
                $ins_seq = reverse_complement($ins_seq);
            }
            push @{$ins_ref->{$ref->{src}}}, {
                bp_start   => $bp,
                bp_end     => $bp,
                len        => length($ins_seq),
                seq        => $ins_seq,
                source     => 'inside_block',
                qry_chr    => $qry->{src},
                qry_start  => $q_s,
                qry_end    => $q_e,
                qry_strand => $qry->{strand},
            };
            next;
        }

        $rpos += $rstep if $rb ne '-';
        $qpos += $qstep if $qb ne '-';
        $i++;
    }
}

sub add_between_block_insertions {
    my ($ins_ref, $cov_ref, $blocks_ref, $extra) = @_;
    foreach my $chr (keys %$blocks_ref) {
        my @b = sort {
            $a->{ref_start} <=> $b->{ref_start} ||
            $a->{ref_end}   <=> $b->{ref_end}
        } @{$blocks_ref->{$chr}};

        for (my $i = 0; $i < @b - 1; $i++) {
            my ($b1, $b2) = ($b[$i], $b[$i + 1]);
            next unless $b1->{qry_chr} eq $b2->{qry_chr};
            next unless $b1->{qry_strand} eq $b2->{qry_strand};

            my $ref_gap = $b2->{ref_start} - $b1->{ref_end};
            next if $ref_gap < 0;

            my ($qry_gap, $qgs, $qge);
            if ($b1->{qry_strand} eq '+') {
                $qry_gap = $b2->{qry_start} - $b1->{qry_end};
                $qgs = $b1->{qry_end};
                $qge = $b2->{qry_start};
            } else {
                $qry_gap = $b1->{qry_start} - $b2->{qry_end};
                $qgs = $b2->{qry_end};
                $qge = $b1->{qry_start};
            }

            next unless defined $qry_gap;
            next if $qry_gap < 0;

            if ($qry_gap > $ref_gap + $extra) {
                my $ins_len = $qry_gap - $ref_gap;
                push @{$ins_ref->{$chr}}, {
                    bp_start   => $b1->{ref_end},
                    bp_end     => $b2->{ref_start},
                    len        => $ins_len,
                    seq        => 'NA',
                    source     => 'between_blocks',
                    qry_chr    => $b1->{qry_chr},
                    qry_start  => $qgs,
                    qry_end    => $qge,
                    qry_strand => $b1->{qry_strand},
                };
            } elsif ($ref_gap > 0 && $qry_gap > 0) {
                my $ratio = $qry_gap / $ref_gap;
                push @{$cov_ref->{$chr}}, [$b1->{ref_end}, $b2->{ref_start}] if $ratio >= 0.80 && $ratio <= 1.25;
            }
        }
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
    my ($r) = @_;
    if ($r->{strand} eq '+') {
        return ($r->{start}, $r->{start} + $r->{size});
    } else {
        my $s = $r->{srcSize} - $r->{start} - $r->{size};
        my $e = $r->{srcSize} - $r->{start};
        return ($s, $e);
    }
}

sub initial_pos {
    my ($r) = @_;
    return ($r->{strand} eq '+') ? $r->{start} : $r->{srcSize} - $r->{start} - 1;
}

sub coord_range_half_open {
    my @c = @_;
    return ('NA', 'NA') unless @c;
    my ($min, $max) = ($c[0], $c[0]);
    foreach my $x (@c) {
        $min = $x if $x < $min;
        $max = $x if $x > $max;
    }
    return ($min, $max + 1);
}

# -----------------------------
# Genotyping logic
# -----------------------------
sub genotype_one_ins {
    my ($v, $ins_ref, $cov_ref, $seq_hash, $opt_ref) = @_;

    my $chr = $v->{chr};
    my $bp_start = $v->{start};
    my $bp_end   = $v->{end};
    ($bp_start, $bp_end) = sort { $a <=> $b } ($bp_start, $bp_end);

    my $expected = $v->{expected};
    my $expected_len = $expected->{len};
    my $expected_seq = $expected->{seq};
    my $short_ref_insertion_like = $expected->{short_ref};

    my $search_start = $bp_start - $opt_ref->{pos_tol};
    my $search_end   = $bp_end   + $opt_ref->{pos_tol};
    $search_start = 0 if $search_start < 0;

    my $best_id = 'NA';
    my $nearby_query_gap = 0;
    my $nearby_note = 'NA';
    my $nearby_ins;

    if (exists $ins_ref->{$chr}) {
        foreach my $ins (@{$ins_ref->{$chr}}) {
            next if $ins->{bp_end}   < $search_start;
            next if $ins->{bp_start} > $search_end;

            my $length_ok = length_match($ins->{len}, $expected_len,
                                         $opt_ref->{len_ratio_min}, $opt_ref->{len_ratio_max}, $opt_ref->{len_slop});

            my $large_between_block_gap = 0;
            if ($short_ref_insertion_like && $ins->{source} eq 'between_blocks' && defined $expected_len && $expected_len > 0) {
                if ($ins->{len} >= $expected_len * 0.50) {
                    $large_between_block_gap = 1;
                    $nearby_query_gap = 1;
                    $nearby_note = 'nearby_between_block_query_gap_not_reference_like';
                    $nearby_ins = $ins;
                }
            }

            next unless ($length_ok || $large_between_block_gap);

            my $cand_seq = $ins->{seq};
            if ((!defined $cand_seq || $cand_seq eq 'NA' || $cand_seq eq '') && %$seq_hash) {
                $cand_seq = fetch_candidate_seq($seq_hash, $ins);
            }

            my $seq_id = 'NA';
            if (defined $cand_seq && $cand_seq ne 'NA' && $cand_seq ne '' && defined $expected_seq && $expected_seq ne '' && $expected_seq ne 'NA') {
                $seq_id = best_sequence_identity($cand_seq, $expected_seq, $opt_ref);
                $best_id = $seq_id if $best_id eq 'NA' || $seq_id > $best_id;

                # If length is not OK, sequence similarity must support the rescue.
                next if $seq_id < $opt_ref->{seq_id};
            } else {
                # Without sequence evidence, do not call rescued large-gap candidates.
                next unless $length_ok;
            }

            return (
                '1',
                'INS_candidate',
                $ins->{len},
                $seq_id,
                $ins->{bp_start} . '-' . $ins->{bp_end},
                $ins->{source},
                $ins->{qry_chr},
                $ins->{qry_start},
                $ins->{qry_end},
                'PASS',
            );
        }
    }

    if ($nearby_query_gap) {
        my $ins = $nearby_ins;
        return (
            '.',
            'nearby_query_gap',
            defined($ins) ? $ins->{len} : 'NA',
            $best_id,
            defined($ins) ? ($ins->{bp_start} . '-' . $ins->{bp_end}) : 'NA',
            defined($ins) ? $ins->{source} : 'NA',
            defined($ins) ? $ins->{qry_chr} : 'NA',
            defined($ins) ? $ins->{qry_start} : 'NA',
            defined($ins) ? $ins->{qry_end} : 'NA',
            $nearby_note,
        );
    }

    my $ref_s = $bp_start - $opt_ref->{flank};
    my $ref_e = $bp_end   + $opt_ref->{flank};
    $ref_s = 0 if $ref_s < 0;
    $ref_e = $ref_s + 1 if $ref_e <= $ref_s;

    my $cov = coverage_fraction($chr, $ref_s, $ref_e, $cov_ref);
    if ($cov >= $opt_ref->{min_refcov}) {
        return ('0', 'reference_covered', 0, $best_id, 'NA', 'ref_coverage', 'NA', 'NA', 'NA', 'PASS');
    }

    return ('.', 'insufficient_evidence', 'NA', $best_id, 'NA', 'NA', 'NA', 'NA', 'NA', 'uncertain');
}

sub length_match {
    my ($obs, $exp, $min_ratio, $max_ratio, $slop) = @_;
    return 0 unless defined $obs && defined $exp;
    return 0 unless $obs =~ /^\d+(?:\.\d+)?$/ && $exp =~ /^\d+(?:\.\d+)?$/;
    return 0 if $obs <= 0 || $exp <= 0;
    return 1 if abs($obs - $exp) <= $slop;
    my $ratio = $obs / $exp;
    return ($ratio >= $min_ratio && $ratio <= $max_ratio) ? 1 : 0;
}

sub fetch_candidate_seq {
    my ($h, $ins) = @_;
    my $seq = fetch_seq_from_hash($h, $ins->{qry_chr}, $ins->{qry_start}, $ins->{qry_end});
    return 'NA' if !defined $seq || $seq eq 'NA' || $seq eq '';
    $seq = reverse_complement($seq) if defined($ins->{qry_strand}) && $ins->{qry_strand} eq '-';
    return $seq;
}

sub fetch_seq_from_hash {
    my ($h, $chr, $s, $e) = @_;
    return 'NA' unless defined $chr && exists $h->{$chr};
    return 'NA' if !defined $s || !defined $e || $s eq 'NA' || $e eq 'NA' || $e <= $s;
    $s = 0 if $s < 0;
    return 'NA' if $s >= length($h->{$chr});
    my $len = $e - $s;
    $len = length($h->{$chr}) - $s if $s + $len > length($h->{$chr});
    return substr($h->{$chr}, $s, $len);
}

sub coverage_fraction {
    my ($chr, $start, $end, $cov_ref) = @_;
    return 0 unless exists $cov_ref->{$chr};
    my $len = $end - $start;
    return 0 if $len <= 0;
    my $covered = 0;
    foreach my $iv (@{$cov_ref->{$chr}}) {
        my ($s, $e) = @$iv;
        last if $s >= $end;
        next if $e <= $start;
        my $os = $s > $start ? $s : $start;
        my $oe = $e < $end   ? $e : $end;
        $covered += $oe - $os if $oe > $os;
    }
    return $covered / $len;
}

sub merge_intervals {
    my @iv = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @_;
    my @m;
    foreach my $x (@iv) {
        next unless defined $x->[0] && defined $x->[1] && $x->[1] > $x->[0];
        if (!@m || $x->[0] > $m[-1][1]) {
            push @m, [$x->[0], $x->[1]];
        } else {
            $m[-1][1] = $x->[1] if $x->[1] > $m[-1][1];
        }
    }
    return @m;
}

# -----------------------------
# Sequence identity
# -----------------------------
sub best_sequence_identity {
    my ($cand, $expected, $opt_ref) = @_;
    $cand = clean_seq($cand);
    $expected = clean_seq($expected);
    return 0 if $cand eq '' || $expected eq '';

    if (length($expected) >= $opt_ref->{long_seq_len} || length($cand) >= $opt_ref->{long_seq_len}) {
        return seed_hit_identity($cand, $expected, $opt_ref);
    }

    my $id1 = semiglobal_identity($cand, $expected);
    my $rc = reverse_complement($cand);
    my $id2 = semiglobal_identity($rc, $expected);
    return $id1 > $id2 ? $id1 : $id2;
}

sub seed_hit_identity {
    my ($cand, $expected, $opt_ref) = @_;
    my $seed_len = $opt_ref->{sample_seed_len};
    $seed_len = length($expected) if $seed_len > length($expected);
    return semiglobal_identity($cand, $expected) if $seed_len < 10;

    my @starts = deterministic_window_starts(length($expected), $opt_ref->{sample_windows}, $opt_ref->{sample_len});
    my $hits = 0;
    foreach my $st (@starts) {
        my $seg_len = $opt_ref->{sample_len};
        $seg_len = length($expected) - $st if $st + $seg_len > length($expected);
        next if $seg_len < $seed_len;
        my $seg = substr($expected, $st, $seg_len);

        my @seed_starts = (0);
        push @seed_starts, int(($seg_len - $seed_len) / 2) if $seg_len - $seed_len > 10;
        push @seed_starts, $seg_len - $seed_len if $seg_len - $seed_len > 0;

        my %seen_seed;
        foreach my $ss (@seed_starts) {
            my $seed = substr($seg, $ss, $seed_len);
            next if $seen_seed{$seed}++;
            next if $seed =~ /^N+$/;
            if (index($cand, $seed) >= 0 || index(reverse_complement($cand), $seed) >= 0) {
                $hits++;
                return 1.0 if $hits >= $opt_ref->{sample_min_hits};
            }
        }
    }
    return 0.0;
}

sub deterministic_window_starts {
    my ($len, $n, $win) = @_;
    $n = 1 if $n < 1;
    return (0) if $len <= $win || $n == 1;
    my @s;
    for (my $i = 0; $i < $n; $i++) {
        my $st = int(($len - $win) * $i / ($n - 1));
        push @s, $st;
    }
    return @s;
}

sub semiglobal_identity {
    my ($text, $pat) = @_;
    $text = clean_seq($text);
    $pat  = clean_seq($pat);
    my $n = length($text);
    my $m = length($pat);
    return 0 if $n == 0 || $m == 0;

    # For very short sequences, simple containment/position identity is enough and faster.
    if ($m <= 30 || $n <= 30) {
        return simple_best_identity($text, $pat);
    }

    # Semi-global alignment: align the full pattern to any substring of text.
    # Score matches/mismatches/gaps and recover identity from the best end state.
    # To keep memory modest, store only two DP rows.
    my $match = 2;
    my $mismatch = -1;
    my $gap = -2;

    my @prev_score = (0) x ($m + 1);
    my @prev_match = (0) x ($m + 1);
    my @prev_len   = (0) x ($m + 1);

    for my $j (1 .. $m) {
        $prev_score[$j] = $prev_score[$j-1] + $gap;
        $prev_match[$j] = 0;
        $prev_len[$j] = $j;
    }

    my ($best_score, $best_match, $best_len) = (-10**9, 0, 1);

    for my $i (1 .. $n) {
        my @cur_score = (0) x ($m + 1);  # free gaps in text before alignment start
        my @cur_match = (0) x ($m + 1);
        my @cur_len   = (0) x ($m + 1);

        for my $j (1 .. $m) {
            my $a = substr($text, $i-1, 1);
            my $b = substr($pat,  $j-1, 1);
            my $is_match = ($a eq $b || $a eq 'N' || $b eq 'N') ? 1 : 0;

            my $diag_s = $prev_score[$j-1] + ($is_match ? $match : $mismatch);
            my $up_s   = $prev_score[$j] + $gap;
            my $left_s = $cur_score[$j-1] + $gap;

            if ($diag_s >= $up_s && $diag_s >= $left_s) {
                $cur_score[$j] = $diag_s;
                $cur_match[$j] = $prev_match[$j-1] + ($is_match ? 1 : 0);
                $cur_len[$j]   = $prev_len[$j-1] + 1;
            } elsif ($up_s >= $left_s) {
                $cur_score[$j] = $up_s;
                $cur_match[$j] = $prev_match[$j];
                $cur_len[$j]   = $prev_len[$j] + 1;
            } else {
                $cur_score[$j] = $left_s;
                $cur_match[$j] = $cur_match[$j-1];
                $cur_len[$j]   = $cur_len[$j-1] + 1;
            }
        }

        if ($cur_score[$m] > $best_score) {
            ($best_score, $best_match, $best_len) = ($cur_score[$m], $cur_match[$m], $cur_len[$m]);
        }

        @prev_score = @cur_score;
        @prev_match = @cur_match;
        @prev_len   = @cur_len;
    }

    return $best_len > 0 ? ($best_match / $best_len) : 0;
}

sub simple_best_identity {
    my ($text, $pat) = @_;
    my $n = length($text);
    my $m = length($pat);
    return 0 if $n == 0 || $m == 0;
    return 1.0 if index($text, $pat) >= 0;
    my $best = 0;
    my $min_shift = -$m + 1;
    my $max_shift = $n - 1;
    for (my $shift = $min_shift; $shift <= $max_shift; $shift++) {
        my ($same, $len) = (0, 0);
        for (my $j = 0; $j < $m; $j++) {
            my $i = $shift + $j;
            next if $i < 0 || $i >= $n;
            my $a = substr($text, $i, 1);
            my $b = substr($pat, $j, 1);
            $same++ if $a eq $b || $a eq 'N' || $b eq 'N';
            $len++;
        }
        next if $len == 0;
        my $den = $m > $n ? $m : $n;
        my $id = $same / $den;
        $best = $id if $id > $best;
    }
    return $best;
}

sub reverse_complement {
    my ($s) = @_;
    $s = reverse($s // '');
    $s =~ tr/ACGTNacgtn/TGCANtgcan/;
    return uc($s);
}

sub clean_seq {
    my ($seq) = @_;
    return '' unless defined $seq;
    $seq = uc($seq);
    $seq =~ s/-//g;
    $seq =~ s/[^ACGTN]//g;
    return $seq;
}

# -----------------------------
# Output
# -----------------------------
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
    print $VCF "##source=genotype_INS_parallel_onefile.v1.pl\n";
    print $VCF "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate derived from BED END\">\n";
    print $VCF "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    print $VCF "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n";
    print $VCF "##INFO=<ID=BED_START,Number=1,Type=Integer,Description=\"Original BED start/POS column\">\n";
    print $VCF "##INFO=<ID=ALLELE_SOURCE,Number=1,Type=String,Description=\"How REF/ALT were generated\">\n";
    print $VCF "##ALT=<ID=INS,Description=\"Insertion\">\n";
    print $VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype: 0 reference-like, 1 insertion-like, . uncertain\">\n";
    print $VCF '#CHROM', "\t", 'POS', "\t", 'ID', "\t", 'REF', "\t", 'ALT', "\t", 'QUAL', "\t", 'FILTER', "\t", 'INFO', "\t", 'FORMAT';
    print $VCF "\t", join("\t", @samples) if @samples;
    print $VCF "\n";

    for (my $i = 0; $i < @$variants_ref; $i++) {
        my $v = $variants_ref->[$i];
        next unless uc($v->{svtype}) eq 'INS';

        my $chrom = $v->{chr};
        my $id    = sanitize_vcf_field($v->{id});
        my ($ref, $alt, $pos, $end, $svlen, $allele_source);
        if ($opt_ref->{vcf_mode} eq 'allele') {
            ($ref, $alt, $pos, $end, $svlen, $allele_source) =
                canonical_ins_alleles($v, $opt_ref, $ref_db);
        }

        # Keep a valid symbolic fallback for unusual legacy rows, but retain a
        # clear flag so 13_prepare_vg_vcf.pl can report/skip it rather than
        # silently placing unknown insertion sequence into a graph.
        if (!defined($ref)) {
            my $start0 = $opt_ref->{bed_one_based} ? $v->{start} - 1 : $v->{start};
            $start0 = 0 if $start0 < 0;
            $pos = $start0 > 0 ? $start0 : 1;
            $ref = $ref_db->get_seq($chrom, $pos - 1, $pos);
            $ref = 'N' unless valid_vcf_sequence($ref);
            $alt = '<INS>';
            $end = $pos;
            $svlen = abs($v->{expected}{len} || 0);
            $allele_source = 'SYMBOLIC_FALLBACK';
        }
        my $info = join(';',
            'END=' . $end,
            'SVTYPE=INS',
            'SVLEN=' . $svlen,
            'BED_START=' . $v->{start},
            'ALLELE_SOURCE=' . $allele_source,
        );

        my @vcf_gts = map { gt_to_vcf($gt_arrays[$_]->[$i], $opt_ref->{diploid}) } 0 .. $#gt_arrays;
        print $VCF join("\t", $chrom, $pos, $id, $ref, $alt, '.', 'PASS', $info, 'GT', @vcf_gts), "\n";
    }

    close $VCF;
}

sub canonical_ins_alleles {
    my ($v, $opt_ref, $ref_db) = @_;
    my $chrom = $v->{chr};
    return unless $ref_db->has_contig($chrom);

    my $start0 = $opt_ref->{bed_one_based} ? $v->{start} - 1 : $v->{start};
    my $end0   = $opt_ref->{bed_one_based} ? $v->{end}       : $v->{end};
    $start0 = 0 if $start0 < 0;
    $end0 = $start0 if $end0 < $start0;
    my $chr_len = $ref_db->length($chrom);
    return if $start0 > $chr_len || $end0 > $chr_len;

    my $event_ref = $ref_db->get_seq($chrom, $start0, $end0);
    return unless defined($event_ref) && $event_ref =~ /^[ACGTN]*$/;

    my $raw_ref = clean_seq($v->{refseq});
    my $raw_alt = clean_seq($v->{altseq});
    return unless valid_vcf_sequence($raw_alt);

    # Old PairwiseSV sequence output added one right-padding base to both sides.
    # Strip that common base when it can be verified against the reference.
    if (valid_vcf_sequence($raw_ref)
        && length($raw_ref) == length($event_ref) + 1
        && length($raw_alt) >= 1
        && substr($raw_ref, 0, length($event_ref)) eq $event_ref
        && substr($raw_ref, -1) eq substr($raw_alt, -1)) {
        $raw_ref = substr($raw_ref, 0, -1);
        $raw_alt = substr($raw_alt, 0, -1);
    }

    my $event_alt;
    if ($opt_ref->{alt_is_insseq}) {
        $event_alt = $event_ref . $raw_alt;
    } else {
        # With --exact_seq_len, ALTseq is the complete query-side replacement
        # sequence.  For a zero-length reference interval it is naturally the
        # pure inserted sequence.
        $event_alt = $raw_alt;
    }
    return unless defined($event_alt) && $event_alt =~ /^[ACGTN]+$/;
    return unless length($event_alt) > length($event_ref);

    my ($pos, $ref, $alt, $source);
    if ($start0 > 0) {
        my $anchor = $ref_db->get_seq($chrom, $start0 - 1, $start0);
        return unless valid_vcf_sequence($anchor);
        $pos = $start0;             # 1-based coordinate of the left anchor
        $ref = $anchor . $event_ref;
        $alt = $anchor . $event_alt;
        $source = 'REFERENCE_LEFT_ANCHORED';
    } else {
        # No left base exists at the start of a contig.  Use the first base
        # following the event as a right anchor, as required by VCF padding.
        return if $end0 >= $chr_len;
        my $anchor = $ref_db->get_seq($chrom, $end0, $end0 + 1);
        return unless valid_vcf_sequence($anchor);
        $pos = 1;
        $ref = $event_ref . $anchor;
        $alt = $event_alt . $anchor;
        $source = 'REFERENCE_RIGHT_ANCHORED';
    }

    my $end = $pos + length($ref) - 1;
    my $svlen = length($alt) - length($ref);
    return unless $svlen > 0;
    return ($ref, $alt, $pos, $end, $svlen, $source);
}

sub valid_vcf_sequence {
    my ($seq) = @_;
    return 0 unless defined $seq;
    return 0 if $seq eq '' || $seq eq '.' || $seq eq '-';
    return $seq =~ /^[ACGTN]+$/ ? 1 : 0;
}

sub gt_to_vcf {
    my ($gt, $diploid) = @_;
    if (!defined $gt || $gt eq '.') {
        return $diploid ? './.' : '.';
    } elsif ($gt eq '0') {
        return $diploid ? '0/0' : '0';
    } elsif ($gt eq '1') {
        return $diploid ? '1/1' : '1';
    } else {
        return $diploid ? './.' : '.';
    }
}

sub sanitize_vcf_field {
    my ($x) = @_;
    $x = '.' unless defined $x && $x ne '';
    $x =~ s/[\s;=]/_/g;
    return $x;
}


################################################################################
# Lightweight reference FASTA random access (0-based, half-open get_seq).
################################################################################
package INSRefFAI;
use strict;
use warnings;
use Fcntl qw(SEEK_SET);

sub new {
    my ($class, $fasta) = @_;
    die "Reference FASTA must be uncompressed: $fasta\n" if $fasta =~ /\.gz$/i;
    my $self = bless { fasta=>$fasta, idx=>{}, order=>[] }, $class;
    if (-s "$fasta.fai") { $self->_load_fai("$fasta.fai"); }
    else { $self->_build_index(); }
    open(my $fh, '<', $fasta) or die "Cannot open reference FASTA $fasta: $!\n";
    binmode($fh);
    $self->{fh} = $fh;
    return $self;
}

sub _load_fai {
    my ($self, $fai) = @_;
    open(my $IN, '<', $fai) or die "Cannot read $fai: $!\n";
    while (<$IN>) {
        chomp;
        next if /^\s*$/;
        my @f = split(/\t/);
        next unless @f >= 5;
        $self->{idx}{$f[0]} = { len=>$f[1]+0, offset=>$f[2]+0, line_bases=>$f[3]+0, line_width=>$f[4]+0 };
        push @{$self->{order}}, $f[0];
    }
    close $IN;
}

sub _build_index {
    my ($self) = @_;
    my $fa = $self->{fasta};
    open(my $IN, '<', $fa) or die "Cannot read $fa: $!\n";
    binmode($IN);
    my ($name, $len, $offset, $line_bases, $line_width);
    while (1) {
        my $p = tell($IN);
        my $line = <$IN>;
        last unless defined $line;
        if ($line =~ /^>(\S+)/) {
            if (defined $name) {
                $self->{idx}{$name} = { len=>$len, offset=>$offset, line_bases=>$line_bases, line_width=>$line_width };
                push @{$self->{order}}, $name;
            }
            $name = $1;
            ($len, $offset, $line_bases, $line_width) = (0, undef, undef, undef);
        } else {
            next unless defined $name;
            my $raw_width = length($line);
            $line =~ s/[\r\n]+$//;
            my $bases = length($line);
            next if $bases == 0;
            $offset = $p unless defined $offset;
            $line_bases = $bases unless defined $line_bases;
            $line_width = $raw_width unless defined $line_width;
            $len += $bases;
        }
    }
    if (defined $name) {
        $self->{idx}{$name} = { len=>$len, offset=>$offset, line_bases=>$line_bases, line_width=>$line_width };
        push @{$self->{order}}, $name;
    }
    close $IN;
}

sub has_contig { return exists $_[0]{idx}{$_[1]}; }
sub length     { return $_[0]{idx}{$_[1]}{len}; }
sub get_seq {
    my ($self, $chr, $start, $end) = @_;
    return undef unless $self->has_contig($chr);
    my $x = $self->{idx}{$chr};
    $start = 0 if $start < 0;
    $end = $x->{len} if $end > $x->{len};
    return '' if $end <= $start;
    my $need = $end - $start;
    my $line = $x->{line_bases};
    my $width = $x->{line_width};
    my $byte = $x->{offset} + int($start / $line) * $width + ($start % $line);
    seek($self->{fh}, $byte, SEEK_SET) or die "Reference FASTA seek failed for $chr:$start-$end\n";
    my $seq = '';
    while (CORE::length($seq) < $need) {
        my $buf = '';
        my $remaining = $need - CORE::length($seq);
        my $want = $remaining + int($remaining / $line + 2) * ($width - $line);
        my $n = read($self->{fh}, $buf, $want);
        last unless $n;
        $buf =~ s/[\r\n]//g;
        $seq .= $buf;
    }
    return uc(substr($seq, 0, $need));
}
1;
