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
use Getopt::Long;
use File::Basename qw(basename);
use File::Path qw(make_path remove_tree);

# -----------------------------------------------------------------------------
# SVGAP small InDel genotyper, relaxed-coordinate split-BED version
# Original author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Updated date: 2026-06-13
#
# v6.1 update:
#   Supports split small-InDel inputs from 05_CombineSV.pl:
#     --del_bed CombinedSV/All.DELs.50bpsmall.bed.combined.sorted.txt
#     --ins_bed CombinedSV/All.INTs.50bpsmall.bed.combined.sorted.txt
#   A single mixed file can still be supplied with --bed.
#
# This version keeps the original within-block and block-junction logic, but adds
# relaxed genotyping for 1-50 bp InDels. Small InDels are often represented with
# a 1-5 bp shifted anchor in MAF/VCF, especially in homopolymers or tandem
# repeats. The new --relaxed_gt mode searches nearby MAF columns for a simple
# target-only gap (DEL) or query-only gap (INS) with the expected length before
# calling a site reference-like. This improves genotype recall while preserving
# conservative length/type matching.
# -----------------------------------------------------------------------------

my $maf_dir = '';
my $target  = '';
my $bed     = '';
my $del_bed = '';
my $ins_bed = '';
my $out     = 'small_indel_genotype';
my $max_len = 50;
my $min_len = 1;
my $maf_suffix = '.sing.maf';
my $other_gt = './.';
my $use_block_gap = 1;
my $max_junction_gap = 1000;
my $threads = 1;
my $keep_tmp = 0;
my $query_list = '';
my $gt_pos_slop = 5;
my $len_slop = 0;
my $relaxed_gt = 1;
my $ref_fasta = '';
# Coordinate system of input BED-like file.
# 0: MAF/BED-style 0-based start, half-open END. This is the default.
# 1: already 1-based POS/END.
my $bed_coord = 0;
my $help = 0;

GetOptions(
    'maf_dir=s'    => \$maf_dir,
    'target=s'     => \$target,
    'refname=s'    => \$target,
    'ref_fasta|ref-fasta|ref=s' => \$ref_fasta,
    'bed=s'        => \$bed,
    'variants=s'   => \$bed,
    'del_bed|del-bed|deletion_bed|deletion-bed=s' => \$del_bed,
    'ins_bed|ins-bed|insertion_bed|insertion-bed=s' => \$ins_bed,
    'out=s'        => \$out,
    'max_len=i'    => \$max_len,
    'min_len=i'    => \$min_len,
    'maf_suffix=s' => \$maf_suffix,
    'query_list|query-list|glist=s' => \$query_list,
    'other_gt=s'        => \$other_gt,
    'use_block_gap!'    => \$use_block_gap,
    'max_junction_gap=i'=> \$max_junction_gap,
    'threads|p=i'       => \$threads,
    'keep_tmp!'         => \$keep_tmp,
    'bed_coord=i'       => \$bed_coord,
    'gt_pos_slop|gt-pos-slop=i' => \$gt_pos_slop,
    'len_slop|len-slop=i' => \$len_slop,
    'relaxed_gt|relaxed-gt!' => \$relaxed_gt,
    'help'              => \$help,
) or usage();

usage() if $help;
usage("ERROR: --maf_dir is required") unless $maf_dir;
usage("ERROR: --target is required, for example --target S184") unless $target;
usage("ERROR: --ref_fasta is required to produce sequence-resolved vg-compatible alleles") unless $ref_fasta;
usage("ERROR: cannot read --ref_fasta $ref_fasta") unless -s $ref_fasta;
usage("ERROR: provide --bed, or --del_bed and/or --ins_bed") unless ($bed || $del_bed || $ins_bed);
usage("ERROR: --max_len must be larger than --min_len") unless $max_len > $min_len;
usage("ERROR: --threads must be >= 1") unless $threads >= 1;
usage("ERROR: --bed_coord must be 0 or 1") unless $bed_coord == 0 || $bed_coord == 1;
usage("ERROR: --gt_pos_slop must be >= 0") unless $gt_pos_slop >= 0;
usage("ERROR: --len_slop must be >= 0") unless $len_slop >= 0;

my $ref_db = SmallIndelRefFAI->new($ref_fasta);

my @input_beds;
push @input_beds, [$bed, 'MIXED'] if $bed;
push @input_beds, [$del_bed, 'DEL_BED'] if $del_bed;
push @input_beds, [$ins_bed, 'INS_BED'] if $ins_bed;

my ($variants_ref, $vars_by_chr_ref, $bed_header_ref) = read_bed_files(\@input_beds, $min_len, $max_len, $bed_coord);
my @variants = @$variants_ref;
my %vars_by_chr = %$vars_by_chr_ref;
my @bed_header = @$bed_header_ref;

if (!@variants) {
    die "No INS/DEL records with $min_len <= abs(SVLEN) < $max_len were found in the input BED file(s)\n";
}

for my $chr (keys %vars_by_chr) {
    my @idx = sort { $variants[$a]{pos} <=> $variants[$b]{pos} || $variants[$a]{end} <=> $variants[$b]{end} } @{$vars_by_chr{$chr}};
    $vars_by_chr{$chr} = \@idx;
}

my @jobs;
if ($query_list) {
    @jobs = jobs_from_query_list($query_list, $maf_dir, $target, $maf_suffix);
} else {
    opendir(my $dh, $maf_dir) or die "Cannot open --maf_dir $maf_dir: $!\n";
    my @maf_files = sort grep { /\Q$maf_suffix\E(?:\.gz)?$/ } readdir($dh);
    closedir($dh);

    if (!@maf_files) {
        die "No MAF files ending with '$maf_suffix' or '$maf_suffix.gz' were found in $maf_dir\n";
    }

    for my $mf (@maf_files) {
        my $sample = sample_from_filename($mf, $target, $maf_suffix);
        push @jobs, { maf => $mf, file => "$maf_dir/$mf", sample => $sample };
    }
}
my @samples = map { $_->{sample} } @jobs;

if ($threads == 1) {
    for my $job (@jobs) {
        warn "Processing $job->{maf} as sample $job->{sample} ...\n";
        process_maf($job->{file}, $job->{sample}, \@variants, \%vars_by_chr, $use_block_gap, $max_junction_gap);
    }
} else {
    run_parallel_jobs(\@jobs, \@variants, \%vars_by_chr, $threads, $use_block_gap, $max_junction_gap, $out, $keep_tmp);
}

write_matrix("$out.matrix.tsv", \@variants, \@samples, \@bed_header);
write_vcf("$out.vcf", \@variants, \@samples, $bed_coord, $ref_db);
write_stats("$out.sample_stats.tsv", \@variants, \@samples);

warn "Done. Output files:\n";
warn "  $out.matrix.tsv\n";
warn "  $out.vcf\n";
warn "  $out.sample_stats.tsv\n";

sub usage {
    my ($msg) = @_;
    print STDERR "$msg\n" if defined $msg;
    print STDERR <<'USAGE';

Usage:
  perl genotype_small_indels_from_sing_maf.v6.relaxed.pl \
    --maf_dir sing_maf_dir \
    --target S184 \
    --ref_fasta S184.fa \
    --del_bed CombinedSV/All.DELs.50bpsmall.bed.combined.sorted.txt \
    --ins_bed CombinedSV/All.INTs.50bpsmall.bed.combined.sorted.txt \
    --out small_indel_genotype \
    --threads 8

Important coordinate behavior:
  By default, the input BED-like file is treated as 0-based half-open because
  it is generated from MAF coordinates. Internal MAF genotyping coordinates are
  converted to 1-based. The final VCF is rebuilt independently from the
  reference FASTA with a real shared anchor: for an event starting at BED START
  > 0, VCF POS equals BED START (the 1-based position of the left anchor).
  END and SVLEN are derived from the actual explicit REF/ALT alleles. The matrix
  output keeps the original BED coordinates unchanged.

  If your input file is already 1-based, add:
    --bed_coord 1

Input BED format expected:
  #CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq
  S184.Chr01 5 6 GPE046410.Chr01_6291_6291_1 DEL 1 1 ... TG G

Main assumptions:
  1. Default: POS/END are MAF/BED-style 0-based half-open coordinates.
     The VCF writer uses --ref_fasta to add a real left anchor and derives
     END/SVLEN from the explicit allele lengths.
  2. REFseq and ALTseq may be exact unpadded event sequences (recommended with
     04_PairwiseSV.pl --exact_seq_len) or legacy right-padded alleles. Both are
     converted to canonical explicit VCF alleles.
  3. MAF files are pairwise single-coverage alignments named as:
       target.query.sing.maf
     for example:
       S184.GPE046410.sing.maf

Options:
  --maf_dir       Directory containing single-coverage MAF files
  --target/--refname  Target/reference sample prefix, for example S184
  --ref_fasta     Uncompressed reference FASTA used to rebuild REF/ALT [required]
  --bed           Optional single mixed BED file of small insertions/deletions
  --del_bed       Small deletion BED file, e.g. All.DELs.50bpsmall.bed.combined.sorted.txt
  --ins_bed       Small insertion BED file, e.g. All.INTs.50bpsmall.bed.combined.sorted.txt
  --out           Output prefix [small_indel_genotype]
  --max_len       Keep indels with abs(SVLEN) < max_len [50]
  --min_len       Keep indels with abs(SVLEN) >= min_len [1]
  --maf_suffix    MAF file suffix [.sing.maf]
  --query_list     Optional query list. If provided, VCF sample order follows this file.
  --other_gt            GT for non-reference/non-alternative alleles [./.]
  --use_block_gap       Also infer 1/1 calls from simple gaps between adjacent MAF blocks [on]
  --no-use_block_gap    Disable block-junction gap genotyping
  --max_junction_gap    Maximum target/query gap length considered at a block junction [1000]
  --threads, -p        Number of MAF files/query samples processed in parallel [1]
  --keep_tmp           Keep temporary per-sample genotype files [off]
  --bed_coord       Input coordinate system: 0=MAF/BED 0-based [default], 1=1-based
  --relaxed_gt      Search nearby MAF gap representation before calling reference [on]
  --no-relaxed_gt   Disable relaxed nearby-gap genotyping
  --gt_pos_slop     Position tolerance for relaxed small InDel genotyping [5]
  --len_slop        Allowed length difference for relaxed gap length matching [0]

Outputs:
  *.matrix.tsv          Original BED columns plus one GT column per query sample
  *.vcf                 VCF-like genotype file
  *.sample_stats.tsv    Numbers of 0/0, 1/1 and missing calls per sample

Parallel mode:
  Parallelization is by query sample/MAF file. Each child process reads one
  target.query.sing.maf file, genotypes all input variants for that query, and
  writes a temporary per-sample genotype file. The parent process then merges
  all samples into the final matrix and VCF. Memory usage is approximately
  multiplied by the number of active workers because each worker parses its own
  MAF file and keeps a copy of the variant table.

Block-junction gap genotyping:
  Some small indels, especially around 20-49 bp, may be represented as the
  unaligned interval between two adjacent MAF blocks rather than as '-' columns
  inside one block. When --use_block_gap is enabled, the script checks adjacent
  collinear MAF blocks from the same query chromosome and strand. It adds a 1/1
  genotype for simple junction events when the target-only gap length matches a
  DEL record, or the query-only gap length matches an INS record. This mode is
  length/coordinate based. For insertions, the inserted sequence itself is not
  present in the MAF gap, so exact ALTseq matching requires the query genome FASTA.

USAGE
    exit(defined $msg ? 1 : 0);
}

sub read_bed_files {
    my ($files_ref, $min_len, $max_len, $bed_coord) = @_;
    my (@all_variants, %vars_by_chr, @header);
    for my $entry (@$files_ref) {
        my ($file, $source_label) = @$entry;
        next unless defined $file && $file ne '';
        my ($v_ref, undef, $h_ref) = read_bed($file, $min_len, $max_len, $bed_coord, $source_label);
        @header = @$h_ref if !@header;
        push @all_variants, @$v_ref;
    }

    @all_variants = sort {
        $a->{chrom} cmp $b->{chrom}
        || $a->{pos} <=> $b->{pos}
        || $a->{end} <=> $b->{end}
        || $a->{type} cmp $b->{type}
        || $a->{id} cmp $b->{id}
    } @all_variants;

    for my $i (0 .. $#all_variants) {
        $all_variants[$i]{idx} = $i;
        push @{$vars_by_chr{$all_variants[$i]{chrom}}}, $i;
    }
    return (\@all_variants, \%vars_by_chr, \@header);
}

sub read_bed {
    my ($file, $min_len, $max_len, $bed_coord, $source_label) = @_;
    $source_label ||= 'BED';
    open(my $fh, '<', $file) or die "Cannot open BED $file: $!\n";

    my @header;
    my %hidx;
    my @variants;
    my %vars_by_chr;

    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        if ($line =~ /^#/) {
            $line =~ s/^#//;
            @header = split(/\s+/, $line);
            for my $i (0 .. $#header) {
                my $k = $header[$i];
                $k =~ s/^#//;
                $hidx{$k} = $i;
                $hidx{uc($k)} = $i;
            }
            next;
        }

        if (!@header) {
            @header = qw(CHROM POS END ID SVTYPE SVLEN MERGE_SAMPLES_NUM MERGE_SAMPLES QUERY Q_POS Q_END FIVE THREE REFseq ALTseq);
            for my $i (0 .. $#header) {
                $hidx{$header[$i]} = $i;
                $hidx{uc($header[$i])} = $i;
            }
        }

        my @f = split(/\s+/, $line);
        next if @f < 6;

        my $chrom  = get_field(\@f, \%hidx, 'CHROM', 0);
        my $pos_raw = get_field(\@f, \%hidx, 'POS', 1);
        my $end_raw = get_field(\@f, \%hidx, 'END', 2);
        my $id      = get_field(\@f, \%hidx, 'ID', 3);
        my $type    = uc(get_field(\@f, \%hidx, 'SVTYPE', 4));
        my $svlen   = get_field(\@f, \%hidx, 'SVLEN', 5);
        my $refseq  = get_field(\@f, \%hidx, 'REFseq', 13);
        my $altseq  = get_field(\@f, \%hidx, 'ALTseq', 14);

        next unless $type eq 'INS' || $type eq 'DEL';
        my $abs_len = abs($svlen + 0);
        next unless $abs_len >= $min_len && $abs_len < $max_len;

        $pos_raw += 0;
        $end_raw += 0;

        # Convert input coordinates to internal 1-based target/reference coordinates.
        # MAF uses 0-based starts. For this BED-like file, the default is therefore:
        #   internal POS = raw POS + 1
        #   internal END = raw END for deletion intervals [POS, END)
        # For insertions with raw POS=raw END, use the converted POS as the anchor.
        my ($pos, $end);
        if ($bed_coord == 0) {
            $pos = $pos_raw + 1;
            $end = $end_raw;
            $end = $pos if $type eq 'INS' && $end < $pos;
        } else {
            $pos = $pos_raw;
            $end = $end_raw;
            $end = $pos if $type eq 'INS' && $end < $pos;
        }

        if ($end < $pos && $type ne 'INS') {
            my $tmp = $pos;
            $pos = $end;
            $end = $tmp;
        }

        my $ref_norm = norm_allele($refseq);
        my $alt_norm = norm_allele($altseq);
        my $vcf_pos = $pos;
        my $vcf_end = vcf_end_from_alleles($vcf_pos, $type, $abs_len, $ref_norm, $alt_norm);

        my $idx = scalar(@variants);
        my %row;
        for my $i (0 .. $#header) {
            $row{$header[$i]} = defined $f[$i] ? $f[$i] : '.';
        }

        push @variants, {
            idx      => $idx,
            raw      => \@f,
            row      => \%row,
            chrom    => $chrom,
            raw_pos  => $pos_raw,
            raw_end  => $end_raw,
            pos      => $pos,
            end      => $end,
            vcf_pos  => $vcf_pos,
            vcf_end  => $vcf_end,
            id       => $id,
            type     => $type,
            svlen    => $abs_len,
            ref      => $ref_norm,
            alt      => $alt_norm,
            ref_raw  => $refseq,
            alt_raw  => $altseq,
            gt       => {},
            source_bed => $source_label,
        };
        push @{$vars_by_chr{$chrom}}, $idx;
    }
    close $fh;

    return (\@variants, \%vars_by_chr, \@header);
}

sub get_field {
    my ($f_ref, $hidx_ref, $name, $default_idx) = @_;
    my $idx;
    if (exists $hidx_ref->{$name}) {
        $idx = $hidx_ref->{$name};
    } elsif (exists $hidx_ref->{uc($name)}) {
        $idx = $hidx_ref->{uc($name)};
    } else {
        $idx = $default_idx;
    }
    return defined $f_ref->[$idx] ? $f_ref->[$idx] : '.';
}

sub norm_allele {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/\s+//g;
    return '' if $s eq '.' || $s eq '-' || uc($s) eq 'NA';
    $s =~ s/-//g;
    return uc($s);
}


sub vcf_end_from_alleles {
    my ($vcf_pos, $type, $svlen, $ref, $alt) = @_;

    # For explicit small indel alleles, the VCF record spans REF from POS to
    # POS + length(REF) - 1. For an insertion represented as one anchor base
    # plus inserted sequence, END therefore equals POS. For symbolic alleles,
    # fall back to SVLEN where possible.
    if (defined $ref && length($ref) > 0) {
        return $vcf_pos + length($ref) - 1;
    }
    if ($type eq 'DEL' && defined $svlen && $svlen > 0) {
        return $vcf_pos + $svlen - 1;
    }
    return $vcf_pos;
}

sub sample_from_filename {
    my ($file, $target, $suffix) = @_;
    my $b = basename($file);
    $b =~ s/\.gz$//;
    $b =~ s/\Q$suffix\E$//;
    $b =~ s/^\Q$target\E\.//;
    return $b;
}

sub jobs_from_query_list {
    my ($list, $maf_dir, $target, $suffix) = @_;
    open(my $fh, '<', $list) or die "Cannot open --query_list $list: $!\n";
    my @jobs;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '' || $line =~ /^#/;
        my ($sample) = split(/\s+/, $line);
        next if !defined $sample || $sample eq '';
        next if $sample eq $target;  # query.list may accidentally contain the reference

        my @candidates = (
            "$maf_dir/$target.$sample$suffix",
            "$maf_dir/$target.$sample$suffix.gz",
            "$maf_dir/$sample$suffix",
            "$maf_dir/$sample$suffix.gz",
        );
        my ($file) = grep { -s $_ } @candidates;
        die "Cannot find MAF for query '$sample'. Tried: @candidates\n" unless defined $file;
        my $maf = basename($file);
        push @jobs, { maf => $maf, file => $file, sample => $sample };
    }
    close $fh;
    die "No query jobs were loaded from --query_list $list\n" unless @jobs;
    return @jobs;
}


sub run_parallel_jobs {
    my ($jobs_ref, $variants_ref, $vars_by_chr_ref, $threads, $use_block_gap, $max_junction_gap, $out_prefix, $keep_tmp) = @_;

    my $tmp_dir = "$out_prefix.parallel_tmp.$$";
    make_path($tmp_dir) or die "Cannot create temporary directory $tmp_dir: $!\n";

    my %pid_to_job_idx;
    my @gt_files;

    for (my $i = 0; $i < @$jobs_ref; $i++) {
        while (scalar(keys %pid_to_job_idx) >= $threads) {
            wait_for_one_child(\%pid_to_job_idx);
        }

        my $job = $jobs_ref->[$i];
        my $gt_file = "$tmp_dir/sample_$i.gt.tsv";
        $gt_files[$i] = $gt_file;

        my $pid = fork();
        die "Cannot fork: $!\n" unless defined $pid;

        if ($pid == 0) {
            eval {
                warn "Processing $job->{maf} as sample $job->{sample} in child $$ ...\n";
                process_maf($job->{file}, $job->{sample}, $variants_ref, $vars_by_chr_ref, $use_block_gap, $max_junction_gap);
                write_one_sample_gt($gt_file, $variants_ref, $job->{sample});
            };
            if ($@) {
                warn "Child $$ failed for sample $job->{sample}: $@\n";
                exit 1;
            }
            exit 0;
        } else {
            $pid_to_job_idx{$pid} = $i;
        }
    }

    while (scalar(keys %pid_to_job_idx) > 0) {
        wait_for_one_child(\%pid_to_job_idx);
    }

    for (my $i = 0; $i < @$jobs_ref; $i++) {
        my $sample = $jobs_ref->[$i]{sample};
        read_one_sample_gt($gt_files[$i], $variants_ref, $sample);
    }

    if ($keep_tmp) {
        warn "Temporary files kept in $tmp_dir\n";
    } else {
        remove_tree($tmp_dir);
    }
}

sub wait_for_one_child {
    my ($pid_to_job_idx_ref) = @_;
    my $pid = wait();
    die "wait() failed: $!\n" if $pid < 0;
    my $status = $?;
    my $idx = delete $pid_to_job_idx_ref->{$pid};
    if ($status != 0) {
        my $exit_code = $status >> 8;
        die "A child process failed. PID=$pid, job_index=$idx, exit_code=$exit_code\n";
    }
}

sub write_one_sample_gt {
    my ($file, $variants_ref, $sample) = @_;
    open(my $fh, '>', $file) or die "Cannot write temporary genotype file $file: $!\n";
    print $fh join("\t", qw(IDX GT)), "\n";
    for my $v (@$variants_ref) {
        my $gt = exists $v->{gt}{$sample} ? $v->{gt}{$sample} : './.';
        print $fh join("\t", $v->{idx}, $gt), "\n";
    }
    close $fh;
}

sub read_one_sample_gt {
    my ($file, $variants_ref, $sample) = @_;
    open(my $fh, '<', $file) or die "Cannot read temporary genotype file $file: $!\n";
    my $header = <$fh>;
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my ($idx, $gt) = split(/\t/, $line);
        next unless defined $idx && defined $gt;
        $variants_ref->[$idx]{gt}{$sample} = $gt;
    }
    close $fh;
}

sub process_maf {
    my ($file, $sample, $variants_ref, $vars_by_chr_ref, $use_block_gap, $max_junction_gap) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "-|", "gzip -dc '$file'") or die "Cannot open gzipped MAF $file: $!\n";
    } else {
        open($fh, '<', $file) or die "Cannot open MAF $file: $!\n";
    }

    my @s_lines;
    my %block_meta_by_chr;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^s\s+/) {
            my @x = split(/\s+/, $line);
            push @s_lines, {
                src     => $x[1],
                start   => $x[2] + 0,
                size    => $x[3] + 0,
                strand  => $x[4],
                srcSize => $x[5] + 0,
                text    => $x[6],
            };
        } elsif ($line =~ /^\s*$/) {
            if (@s_lines) {
                my $meta = handle_block(\@s_lines, $sample, $variants_ref, $vars_by_chr_ref);
                push @{$block_meta_by_chr{$meta->{chr}}}, $meta if defined $meta;
            }
            @s_lines = ();
        }
    }
    if (@s_lines) {
        my $meta = handle_block(\@s_lines, $sample, $variants_ref, $vars_by_chr_ref);
        push @{$block_meta_by_chr{$meta->{chr}}}, $meta if defined $meta;
    }
    close $fh;

    if ($use_block_gap) {
        handle_block_junction_gaps(\%block_meta_by_chr, $sample, $variants_ref, $vars_by_chr_ref, $max_junction_gap);
    }
}

sub handle_block {
    my ($s_lines_ref, $sample, $variants_ref, $vars_by_chr_ref) = @_;

    my ($tline, $qline);
    for my $s (@$s_lines_ref) {
        if (exists $vars_by_chr_ref->{$s->{src}}) {
            $tline = $s;
            last;
        }
    }
    return undef unless $tline;

    for my $s (@$s_lines_ref) {
        next if $s == $tline;
        $qline = $s;
        last;
    }
    return undef unless $qline;

    my $chr = $tline->{src};
    my $block_start = $tline->{start} + 1;
    my $block_end   = $tline->{start} + $tline->{size};

    my $idx_list = $vars_by_chr_ref->{$chr} || return make_block_meta($tline, $qline);
    my $left = lower_bound_variant_pos($idx_list, $variants_ref, $block_start - $gt_pos_slop);

    my @t = split(//, uc($tline->{text}));
    my @q = split(//, uc($qline->{text}));
    return undef unless @t == @q;

    my %coord_to_idx;
    my $coord = $tline->{start} + 1;
    for my $i (0 .. $#t) {
        if ($t[$i] ne '-') {
            $coord_to_idx{$coord} = $i;
            $coord++;
        }
    }

    for (my $p = $left; $p < @$idx_list; $p++) {
        my $vidx = $idx_list->[$p];
        my $v = $variants_ref->[$vidx];
        last if $v->{pos} > $block_end + $gt_pos_slop;
        next if exists $v->{gt}{$sample} && $v->{gt}{$sample} ne './.';

        my $gt = genotype_one_variant($v, \@t, \@q, \%coord_to_idx, $other_gt);
        if (defined $gt) {
            $v->{gt}{$sample} = $gt;
            $v->{method}{$sample} = 'within_block' if $gt ne './.';
        }
    }

    return make_block_meta($tline, $qline);
}

sub make_block_meta {
    my ($tline, $qline) = @_;

    my $q_plus_start;
    my $q_plus_end;
    if ($qline->{strand} eq '+') {
        $q_plus_start = $qline->{start} + 1;
        $q_plus_end   = $qline->{start} + $qline->{size};
    } else {
        # Convert to coordinates on the query plus strand.
        $q_plus_start = $qline->{srcSize} - ($qline->{start} + $qline->{size}) + 1;
        $q_plus_end   = $qline->{srcSize} - $qline->{start};
    }

    return {
        chr          => $tline->{src},
        t_start      => $tline->{start} + 1,
        t_end        => $tline->{start} + $tline->{size},
        qsrc         => $qline->{src},
        qstrand      => $qline->{strand},
        q_plus_start => $q_plus_start,
        q_plus_end   => $q_plus_end,
    };
}

sub handle_block_junction_gaps {
    my ($blocks_by_chr_ref, $sample, $variants_ref, $vars_by_chr_ref, $max_junction_gap) = @_;

    for my $chr (keys %$blocks_by_chr_ref) {
        my @blocks = sort { $a->{t_start} <=> $b->{t_start} || $a->{t_end} <=> $b->{t_end} } @{$blocks_by_chr_ref->{$chr}};
        next unless @blocks > 1;
        my $idx_list = $vars_by_chr_ref->{$chr} || next;

        for (my $b = 0; $b < @blocks - 1; $b++) {
            my $left  = $blocks[$b];
            my $right = $blocks[$b + 1];

            next unless $left->{qsrc} eq $right->{qsrc};
            next unless $left->{qstrand} eq $right->{qstrand};

            my $tgap_len = $right->{t_start} - $left->{t_end} - 1;
            next if $tgap_len < 0;
            next if $tgap_len > $max_junction_gap;

            my $qgap_len;
            if ($left->{qstrand} eq '+') {
                next unless $right->{q_plus_start} > $left->{q_plus_end};
                $qgap_len = $right->{q_plus_start} - $left->{q_plus_end} - 1;
            } else {
                next unless $left->{q_plus_start} > $right->{q_plus_end};
                $qgap_len = $left->{q_plus_start} - $right->{q_plus_end} - 1;
            }
            next if !defined $qgap_len || $qgap_len < 0;
            next if $qgap_len > $max_junction_gap;

            # Only use simple junction gaps: target-only sequence suggests a DEL in query;
            # query-only sequence suggests an INS in query. Cases with gaps on both sides are
            # treated as complex and left unchanged.
            next if $tgap_len > 0 && $qgap_len > 0;
            next if $tgap_len == 0 && $qgap_len == 0;

            my $search_start = $left->{t_end} - 2;
            $search_start = 1 if $search_start < 1;
            my $search_end = $right->{t_start} + 2;

            my $p0 = lower_bound_variant_pos($idx_list, $variants_ref, $search_start - 60);
            for (my $p = $p0; $p < @$idx_list; $p++) {
                my $vidx = $idx_list->[$p];
                my $v = $variants_ref->[$vidx];
                last if $v->{pos} > $search_end + 60;

                next if exists $v->{gt}{$sample} && $v->{gt}{$sample} ne './.';

                my $gt = genotype_one_variant_from_junction_gap($v, $left, $right, $tgap_len, $qgap_len);
                if (defined $gt) {
                    $v->{gt}{$sample} = $gt;
                    $v->{method}{$sample} = 'block_gap';
                }
            }
        }
    }
}

sub genotype_one_variant_from_junction_gap {
    my ($v, $left, $right, $tgap_len, $qgap_len) = @_;

    my $type = $v->{type};
    my $len_from_alleles = abs(length($v->{ref}) - length($v->{alt}));
    my $expected_len = $len_from_alleles > 0 ? $len_from_alleles : $v->{svlen};
    return undef unless $expected_len > 0;

    if ($type eq 'DEL') {
        return undef unless $tgap_len > 0 && $qgap_len == 0;
        return undef unless abs($tgap_len - $expected_len) <= $len_slop;

        my $gap_start = $left->{t_end} + 1;
        my $gap_end   = $right->{t_start} - 1;

        # Allow one padded anchor base on either side.
        return undef unless $v->{pos} <= $gap_start + $gt_pos_slop;
        return undef unless $v->{end} >= $gap_end - $gt_pos_slop;
        return '1/1';
    }
    elsif ($type eq 'INS') {
        return undef unless $qgap_len > 0 && $tgap_len == 0;
        return undef unless abs($qgap_len - $expected_len) <= $len_slop;

        # Insertion is between left target end and right target start.
        # BED records may use the left anchor, the right anchor, or POS=END at the breakpoint.
        my $break_left  = $left->{t_end};
        my $break_right = $right->{t_start};
        return undef unless ($v->{pos} <= $break_right + $gt_pos_slop && $v->{end} >= $break_left - $gt_pos_slop);
        return '1/1';
    }

    return undef;
}

sub lower_bound_variant_pos {
    my ($idx_list, $variants_ref, $value) = @_;
    my ($lo, $hi) = (0, scalar(@$idx_list));
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        my $vpos = $variants_ref->[$idx_list->[$mid]]{pos};
        if ($vpos < $value) {
            $lo = $mid + 1;
        } else {
            $hi = $mid;
        }
    }
    return $lo;
}

sub genotype_one_variant {
    my ($v, $t_ref, $q_ref, $coord_to_idx_ref, $other_gt) = @_;

    my $type = $v->{type};
    my $pos  = $v->{pos};
    my $end  = $v->{end};
    my $ref  = $v->{ref};
    my $alt  = $v->{alt};

    return undef if $alt eq '' && $ref eq '';

    my ($target_allele, $query_allele);

    if ($type eq 'DEL') {
        my $region_start = $pos;
        my $region_end   = $end;

        # If REFseq length is inconsistent with POS-END, trust REFseq length.
        # This is useful for VCF-like padded alleles.
        if (length($ref) > 0 && length($ref) != ($region_end - $region_start + 1)) {
            $region_end = $region_start + length($ref) - 1;
        }

        unless (exists $coord_to_idx_ref->{$region_start} && exists $coord_to_idx_ref->{$region_end}) {
            if ($relaxed_gt) {
                my $relaxed = relaxed_gap_genotype($v, $t_ref, $q_ref, $coord_to_idx_ref);
                return $relaxed if defined $relaxed;
            }
            return undef;
        }
        my $i1 = $coord_to_idx_ref->{$region_start};
        my $i2 = $coord_to_idx_ref->{$region_end};
        if ($i2 < $i1) {
            if ($relaxed_gt) {
                my $relaxed = relaxed_gap_genotype($v, $t_ref, $q_ref, $coord_to_idx_ref);
                return $relaxed if defined $relaxed;
            }
            return undef;
        }

        ($target_allele, $query_allele) = alleles_from_columns($t_ref, $q_ref, $i1, $i2);
    }
    elsif ($type eq 'INS') {
        if (length($ref) > 0) {
            my $anchor_start = $pos;
            my $anchor_end   = $pos + length($ref) - 1;
            unless (exists $coord_to_idx_ref->{$anchor_start} && exists $coord_to_idx_ref->{$anchor_end}) {
                if ($relaxed_gt) {
                    my $relaxed = relaxed_gap_genotype($v, $t_ref, $q_ref, $coord_to_idx_ref);
                    return $relaxed if defined $relaxed;
                }
                return undef;
            }
            my $i1 = $coord_to_idx_ref->{$anchor_start};
            my $i2 = $coord_to_idx_ref->{$anchor_end};
            if ($i2 < $i1) {
                if ($relaxed_gt) {
                    my $relaxed = relaxed_gap_genotype($v, $t_ref, $q_ref, $coord_to_idx_ref);
                    return $relaxed if defined $relaxed;
                }
                return undef;
            }

            # Include query-insertion columns immediately after the reference anchor.
            # If the anchor is at the end of this MAF block and no insertion columns
            # are present, do not call 0/0 here; the insertion may be represented as
            # a query-only gap between this block and the next block.
            my $orig_i2 = $i2;
            my $j = $i2 + 1;
            while ($j < @$t_ref && $t_ref->[$j] eq '-') {
                $i2 = $j;
                $j++;
            }
            return undef if $i2 == $orig_i2 && $j >= @$t_ref;
            ($target_allele, $query_allele) = alleles_from_columns($t_ref, $q_ref, $i1, $i2);
        } else {
            # Unpadded insertion: genotype insertion sequence immediately after POS.
            unless (exists $coord_to_idx_ref->{$pos}) {
                if ($relaxed_gt) {
                    my $relaxed = relaxed_gap_genotype($v, $t_ref, $q_ref, $coord_to_idx_ref);
                    return $relaxed if defined $relaxed;
                }
                return undef;
            }
            my $anchor_idx = $coord_to_idx_ref->{$pos};
            my $j = $anchor_idx + 1;
            my $query_ins = '';
            while ($j < @$t_ref && $t_ref->[$j] eq '-') {
                $query_ins .= $q_ref->[$j] if $q_ref->[$j] ne '-';
                $j++;
            }
            return undef if $query_ins eq '' && $j >= @$t_ref;
            $target_allele = '';
            $query_allele = $query_ins;
        }
    }
    else {
        return undef;
    }

    $target_allele = norm_allele($target_allele);
    $query_allele  = norm_allele($query_allele);

    # Exact ALT call first.
    return '1/1' if $query_allele eq $alt;

    # Small InDels often have shifted anchors in repetitive sequence.  Before
    # rejecting the block or calling 0/0, search a small nearby window for a
    # simple gap of the expected length in the MAF columns.
    if ($relaxed_gt) {
        my $relaxed = relaxed_gap_genotype($v, $t_ref, $q_ref, $coord_to_idx_ref);
        return $relaxed if defined $relaxed;
    }

    # When REFseq is available, require that the target side of the MAF block matches it.
    # If it does not match, this block is probably not suitable for genotyping this record.
    if (length($ref) > 0 && $target_allele ne $ref) {
        return undef;
    }

    return '0/0' if $query_allele eq $ref;

    # If the ALT allele is the inserted-only sequence, also accept anchor+inserted sequence.
    if ($type eq 'INS' && length($ref) > 0 && length($alt) > length($ref)) {
        return '1/1' if $query_allele eq $alt;
    }

    return $other_gt;
}

sub relaxed_gap_genotype {
    my ($v, $t_ref, $q_ref, $coord_to_idx_ref) = @_;

    my $type = $v->{type};
    my $expected_len = abs(length($v->{ref}) - length($v->{alt}));
    $expected_len = $v->{svlen} if !$expected_len;
    return undef unless $expected_len > 0;

    my $center = $v->{pos};
    for my $shift (-$gt_pos_slop .. $gt_pos_slop) {
        my $p = $center + $shift;
        next if $p < 1;

        if ($type eq 'DEL') {
            my $s = $p;
            my $e = $p + $expected_len - 1;
            next unless exists $coord_to_idx_ref->{$s} && exists $coord_to_idx_ref->{$e};
            my $i1 = $coord_to_idx_ref->{$s};
            my $i2 = $coord_to_idx_ref->{$e};
            next if $i2 < $i1;

            my ($target_allele, $query_allele) = alleles_from_columns($t_ref, $q_ref, $i1, $i2);
            $target_allele = norm_allele($target_allele);
            $query_allele  = norm_allele($query_allele);
            next unless length($target_allele) >= $expected_len - $len_slop;
            next unless length($target_allele) <= $expected_len + $len_slop;
            return '1/1' if length($query_allele) <= $len_slop;
        }
        elsif ($type eq 'INS') {
            next unless exists $coord_to_idx_ref->{$p};
            my $anchor_idx = $coord_to_idx_ref->{$p};
            my $j = $anchor_idx + 1;
            my $query_ins = '';
            while ($j < @$t_ref && $t_ref->[$j] eq '-') {
                $query_ins .= $q_ref->[$j] if $q_ref->[$j] ne '-';
                $j++;
            }
            my $ins_len = length(norm_allele($query_ins));
            return '1/1' if abs($ins_len - $expected_len) <= $len_slop && $ins_len > 0;
        }
    }
    return undef;
}

sub alleles_from_columns {
    my ($t_ref, $q_ref, $i1, $i2) = @_;
    my ($ta, $qa) = ('', '');
    for my $i ($i1 .. $i2) {
        $ta .= $t_ref->[$i];
        $qa .= $q_ref->[$i];
    }
    $ta =~ s/-//g;
    $qa =~ s/-//g;
    return ($ta, $qa);
}

sub write_matrix {
    my ($file, $variants_ref, $samples_ref, $header_ref) = @_;
    open(my $outfh, '>', $file) or die "Cannot write $file: $!\n";
    print $outfh join("\t", @$header_ref, @$samples_ref), "\n";
    for my $v (@$variants_ref) {
        my @raw = @{$v->{raw}};
        # Keep the same number/order of fields as the BED header.
        if (@raw < @$header_ref) {
            push @raw, ('.') x (@$header_ref - @raw);
        }
        @raw = @raw[0 .. $#$header_ref];
        my @gts = map { exists $v->{gt}{$_} ? $v->{gt}{$_} : './.' } @$samples_ref;
        print $outfh join("\t", @raw, @gts), "\n";
    }
    close $outfh;
}

sub write_vcf {
    my ($file, $variants_ref, $samples_ref, $bed_coord, $ref_db) = @_;
    open(my $outfh, '>', $file) or die "Cannot write $file: $!\n";
    print $outfh "##fileformat=VCFv4.2\n";
    print $outfh "##source=genotype_small_indels_from_sing_maf.vg_ready.pl\n";
    print $outfh "##input_bed_coordinate_system=" . ($bed_coord == 0 ? "0-based_MAF/BED" : "1-based") . "\n";
    print $outfh "##reference=$ref_fasta\n";
    print $outfh "##relaxed_gt=$relaxed_gt;gt_pos_slop=$gt_pos_slop;len_slop=$len_slop\n";
    print $outfh "##input_bed=$bed\n" if $bed;
    print $outfh "##input_del_bed=$del_bed\n" if $del_bed;
    print $outfh "##input_ins_bed=$ins_bed\n" if $ins_bed;
    print $outfh "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type\">\n";
    print $outfh "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"ALT length minus REF length for explicit allele\">\n";
    print $outfh "##INFO=<ID=END,Number=1,Type=Integer,Description=\"1-based end coordinate based on REF allele span\">\n";
    print $outfh "##INFO=<ID=RAWPOS,Number=1,Type=Integer,Description=\"Original POS value from the input BED-like file\">\n";
    print $outfh "##INFO=<ID=RAWEND,Number=1,Type=Integer,Description=\"Original END value from the input BED-like file\">\n";
    print $outfh "##INFO=<ID=ALLELE_SOURCE,Number=1,Type=String,Description=\"How REF/ALT were generated\">\n";
    print $outfh "##ALT=<ID=INS,Description=\"Unresolved insertion sequence\">\n";
    print $outfh "##ALT=<ID=DEL,Description=\"Deletion\">\n";
    print $outfh "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $outfh join("\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', @$samples_ref), "\n";

    for my $v (@$variants_ref) {
        my ($ref,$alt,$vcf_pos,$vcf_end,$svlen,$source) = canonical_small_indel_alleles($v,$bed_coord,$ref_db);
        if (!defined $ref) {
            my $start0 = $bed_coord == 0 ? $v->{raw_pos} : $v->{raw_pos}-1;
            $start0 = 0 if $start0 < 0;
            $vcf_pos = $start0 > 0 ? $start0 : 1;
            $ref = $ref_db->get_seq($v->{chrom},$vcf_pos-1,$vcf_pos);
            $ref = 'N' unless defined($ref) && $ref =~ /^[ACGTN]$/i;
            $alt = $v->{type} eq 'INS' ? '<INS>' : '<DEL>';
            $vcf_end = $v->{type} eq 'INS' ? $vcf_pos : ($bed_coord == 0 ? $v->{raw_end} : $v->{raw_end});
            $svlen = $v->{type} eq 'INS' ? abs($v->{svlen}) : -abs($v->{svlen});
            $source = 'SYMBOLIC_FALLBACK';
        }
        my $info = join(';', 'SVTYPE=' . $v->{type}, 'SVLEN=' . $svlen, 'END=' . $vcf_end,
            'RAWPOS=' . $v->{raw_pos}, 'RAWEND=' . $v->{raw_end}, 'ALLELE_SOURCE=' . $source);
        my @gts = map { exists $v->{gt}{$_} ? $v->{gt}{$_} : './.' } @$samples_ref;
        print $outfh join("\t", $v->{chrom}, $vcf_pos, $v->{id}, $ref, $alt, '.', 'PASS', $info, 'GT', @gts), "\n";
    }
    close $outfh;
}

sub canonical_small_indel_alleles {
    my ($v,$bed_coord,$ref_db) = @_;
    my $chrom = $v->{chrom};
    return unless $ref_db->has_contig($chrom);
    my $start0 = $bed_coord == 0 ? $v->{raw_pos} : $v->{raw_pos}-1;
    my $end0   = $bed_coord == 0 ? $v->{raw_end} : $v->{raw_end};
    $start0=0 if $start0<0;
    $end0=$start0 if $end0<$start0;
    my $chr_len=$ref_db->length($chrom);
    return if $start0>$chr_len || $end0>$chr_len;
    my $event_ref=$ref_db->get_seq($chrom,$start0,$end0);
    return unless defined($event_ref) && $event_ref =~ /^[ACGTN]*$/i;
    my $raw_ref=norm_allele($v->{ref_raw});
    my $raw_alt=norm_allele($v->{alt_raw});

    # Compatibility with old PairwiseSV output that included one right anchor.
    if ($raw_ref ne '' && length($raw_ref)==length($event_ref)+1 && length($raw_alt)>=1
        && substr($raw_ref,0,length($event_ref)) eq $event_ref
        && substr($raw_ref,-1) eq substr($raw_alt,-1)) {
        $raw_ref=substr($raw_ref,0,-1);
        $raw_alt=substr($raw_alt,0,-1);
    }
    my $event_alt=$raw_alt;
    return unless $event_alt =~ /^[ACGTN]*$/i;
    if ($v->{type} eq 'INS') { return unless length($event_alt)>length($event_ref); }
    else                     { return unless length($event_ref)>length($event_alt); }

    my ($pos,$ref,$alt,$source);
    if ($start0>0) {
        my $anchor=$ref_db->get_seq($chrom,$start0-1,$start0);
        return unless defined($anchor) && $anchor =~ /^[ACGTN]$/i;
        $pos=$start0;
        $ref=$anchor.$event_ref;
        $alt=$anchor.$event_alt;
        $source='REFERENCE_LEFT_ANCHORED';
    } else {
        return if $end0 >= $chr_len;
        my $anchor=$ref_db->get_seq($chrom,$end0,$end0+1);
        return unless defined($anchor) && $anchor =~ /^[ACGTN]$/i;
        $pos=1;
        $ref=$event_ref.$anchor;
        $alt=$event_alt.$anchor;
        $source='REFERENCE_RIGHT_ANCHORED';
    }
    my $end=$pos+length($ref)-1;
    my $svlen=length($alt)-length($ref);
    return if $svlen==0;
    return ($ref,$alt,$pos,$end,$svlen,$source);
}

sub write_stats {
    my ($file, $variants_ref, $samples_ref) = @_;
    open(my $outfh, '>', $file) or die "Cannot write $file: $!\n";
    print $outfh join("\t", qw(Sample Total Ref_0_0 Alt_1_1 Missing Other)), "\n";
    for my $s (@$samples_ref) {
        my ($ref, $alt, $miss, $other) = (0, 0, 0, 0);
        for my $v (@$variants_ref) {
            my $gt = exists $v->{gt}{$s} ? $v->{gt}{$s} : './.';
            if ($gt eq '0/0') { $ref++; }
            elsif ($gt eq '1/1') { $alt++; }
            elsif ($gt eq './.') { $miss++; }
            else { $other++; }
        }
        print $outfh join("\t", $s, scalar(@$variants_ref), $ref, $alt, $miss, $other), "\n";
    }
    close $outfh;
}


################################################################################
# Lightweight reference FASTA random access (0-based, half-open get_seq).
################################################################################
package SmallIndelRefFAI;
use strict;
use warnings;
use Fcntl qw(SEEK_SET);
sub new {
    my ($class,$fasta)=@_;
    die "Reference FASTA must be uncompressed: $fasta\n" if $fasta =~ /\.gz$/i;
    my $self=bless {fasta=>$fasta,idx=>{}},$class;
    if (-s "$fasta.fai") { $self->_load_fai("$fasta.fai"); } else { $self->_build_index(); }
    open(my $fh,'<',$fasta) or die "Cannot open reference FASTA $fasta: $!\n"; binmode($fh); $self->{fh}=$fh;
    return $self;
}
sub _load_fai { my ($s,$f)=@_; open(my $I,'<',$f) or die "Cannot read $f: $!\n"; while(<$I>){chomp;next if /^\s*$/;my @x=split(/\t/);next unless @x>=5;$s->{idx}{$x[0]}={len=>$x[1]+0,offset=>$x[2]+0,line_bases=>$x[3]+0,line_width=>$x[4]+0};}close $I; }
sub _build_index {
    my ($s)=@_; open(my $I,'<',$s->{fasta}) or die "Cannot read $s->{fasta}: $!\n";binmode($I);
    my($n,$l,$o,$lb,$lw);while(1){my $p=tell($I);my $line=<$I>;last unless defined $line;if($line=~/^>(\S+)/){$s->{idx}{$n}={len=>$l,offset=>$o,line_bases=>$lb,line_width=>$lw} if defined $n;$n=$1;($l,$o,$lb,$lw)=(0,undef,undef,undef);}else{next unless defined $n;my $rw=length($line);$line=~s/[\r\n]+$//;my $b=length($line);next unless $b;$o=$p unless defined $o;$lb=$b unless defined $lb;$lw=$rw unless defined $lw;$l+=$b;}}$s->{idx}{$n}={len=>$l,offset=>$o,line_bases=>$lb,line_width=>$lw} if defined $n;close $I;
}
sub has_contig{return exists $_[0]{idx}{$_[1]};} sub length{return $_[0]{idx}{$_[1]}{len};}
sub get_seq { my($s,$c,$a,$b)=@_;return undef unless $s->has_contig($c);my$x=$s->{idx}{$c};$a=0 if $a<0;$b=$x->{len} if $b>$x->{len};return '' if $b<=$a;my$n=$b-$a;my$byte=$x->{offset}+int($a/$x->{line_bases})*$x->{line_width}+($a%$x->{line_bases});seek($s->{fh},$byte,SEEK_SET) or die "FASTA seek failed $c:$a-$b\n";my$q='';while(CORE::length($q)<$n){my$buf='';my$r=$n-CORE::length($q);my$want=$r+int($r/$x->{line_bases}+2)*($x->{line_width}-$x->{line_bases});my$got=read($s->{fh},$buf,$want);last unless $got;$buf=~s/[\r\n]//g;$q.=$buf;}return uc(substr($q,0,$n)); }
1;
