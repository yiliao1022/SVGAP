#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
# ============================================================================
# SVGAP standardized script: 04_PairwiseSV.pl
# Original author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Standardized date: 2026-06-12
# Update: v11 uses core Perl fork()/waitpid() for parallel processing, so
#         Parallel::ForkManager is no longer required for --threads > 1.
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
use Getopt::Long;
use File::Basename qw(basename);
use File::Path qw(make_path);

################################################################################
# 4_PairwiseSV_fast_seq.ratio_indel.split.pl
#
# Faster SV calling from Chain/Net/SynNet files, with an optional --seq mode.
# This version uses stricter DEL/INS classification: simple DEL/INS are
# reported only when one side of the gap contains the main SV sequence and
# the opposite side is short relative to the main side. By default, the
# opposite-side sequence must be <= 10% of the main-side sequence
# (--max_opposite_ratio 0.10). Otherwise the event is reported as
# CMPLX/CMPLX_PAV, not as DEL/INS.
#
# Compared with the original 4_PairwiseSV.pl:
#   - Parses SynNet files directly instead of constructing a Tree::DAG_Node tree.
#   - Keeps only one top-level fill branch in memory at a time.
#   - Does not read genome FASTA files unless --seq is specified.
#   - When --seq is specified, FASTA sequences are extracted by random access
#     using a Perl implementation of .fai indexing, rather than loading whole
#     genomes into memory.
#
# Output coordinates are 0-based, half-open, consistent with BED style.
################################################################################

my ($synnet, $syndir, $outdir, $genomes_dir, $threads, $help);
my ($ref_fa, $qry_fa, $do_seq, $exact_seq_len);
my ($invlen, $cnvlen, $dellen, $inslen, $trllen, $toplen, $secondlen, $ali_rate);
my ($max_level, $keep_tmp, $quiet, $seq_context, $max_full_seq, $edge_seq_len);
my ($max_opposite_len, $max_opposite_ratio);

GetOptions(
    "synnet|syn:s"     => \$synnet,
    "syndir:s"         => \$syndir,
    "outdir:s"         => \$outdir,
    "gd:s"             => \$genomes_dir,
    "ref_fa:s"         => \$ref_fa,
    "qry_fa|query_fa:s"=> \$qry_fa,
    "seq"              => \$do_seq,
    "exact_seq_len"    => \$exact_seq_len,
    "t|threads:i"      => \$threads,
    "invlen:i"         => \$invlen,
    "cnvlen:i"         => \$cnvlen,
    "dellen:i"         => \$dellen,
    "inslen:i"         => \$inslen,
    "trllen:i"         => \$trllen,
    "filltoplen:i"     => \$toplen,
    "fill2len:i"       => \$secondlen,
    "rate:f"           => \$ali_rate,
    "max_level:i"      => \$max_level,
    "seq_context:i"    => \$seq_context,
    "max_full_seq:i"   => \$max_full_seq,
    "edge_seq_len:i"   => \$edge_seq_len,
    "max_opposite_len:i" => \$max_opposite_len,  # deprecated; kept for backward compatibility
    "max_opposite_ratio:f" => \$max_opposite_ratio,
    "keep_tmp"         => \$keep_tmp,
    "quiet"            => \$quiet,
    "help"             => \$help,
) or die usage();

if ($help || (!defined $synnet && !defined $syndir)) {
    print usage();
    exit;
}

$outdir       ||= ".";
$threads      ||= 4;
$invlen       ||= 1000;
$cnvlen       ||= 500;
$dellen       ||= 20000;
$inslen       ||= $dellen;
$trllen       ||= 5000;
$ali_rate     ||= 0.25;
$toplen       ||= 20000;
$secondlen    ||= 10000;
$max_level    ||= 3;
$seq_context  ||= 20;
$max_full_seq ||= 1000000;
$edge_seq_len ||= 500;
$max_opposite_ratio = 0.10 unless defined $max_opposite_ratio;
$max_opposite_len = undef unless defined $max_opposite_len;  # optional absolute cap; not used unless explicitly set

make_path($outdir) unless -d $outdir;

my @jobs;
if (defined $synnet) {
    my ($rname, $qname) = infer_pair_from_file($synnet);
    my $rfa = $ref_fa;
    my $qfa = $qry_fa;
    if ($do_seq) {
        $rfa ||= find_genome_path($genomes_dir, $rname) if defined $genomes_dir && defined $rname;
        $qfa ||= find_genome_path($genomes_dir, $qname) if defined $genomes_dir && defined $qname;
        die "ERROR: --seq requires --ref_fa and --qry_fa, or --gd with genome files matching the SynNet file name.\n"
            unless defined $rfa && defined $qfa;
    }
    push @jobs, [$synnet, output_dir_for_file($synnet, $outdir), $rfa, $qfa];
} else {
    opendir(my $dh, $syndir) or die "ERROR: cannot open directory $syndir: $!\n";
    while (my $f = readdir($dh)) {
        next if $f =~ /^\./;
        next unless $f =~ /\.synnet$/ || $f =~ /chain\.filter\.(?:t|q)net\.synnet(?:\.refiltered)?\.synnet$/;
        my ($rname, $qname) = infer_pair_from_file($f);
        my ($rfa, $qfa) = ($ref_fa, $qry_fa);
        if ($do_seq) {
            $rfa ||= find_genome_path($genomes_dir, $rname) if defined $genomes_dir && defined $rname;
            $qfa ||= find_genome_path($genomes_dir, $qname) if defined $genomes_dir && defined $qname;
            die "ERROR: --seq cannot find FASTA for $f. Need --gd, or --ref_fa/--qry_fa for a single pair.\n"
                unless defined $rfa && defined $qfa;
        }
        push @jobs, ["$syndir/$f", output_dir_for_file($f, $outdir), $rfa, $qfa];
    }
    closedir($dh);
}

die "ERROR: no synnet files found.\n" unless @jobs;

if ($threads > 1 && @jobs > 1) {
    run_jobs_with_core_fork(\@jobs, $threads);
} else {
    call_sv_fast($_->[0], $_->[1], $_->[2], $_->[3]) for @jobs;
}

exit;

sub run_jobs_with_core_fork {
    my ($jobs, $max_threads) = @_;
    $max_threads ||= 1;
    $max_threads = 1 if $max_threads < 1;

    # Avoid requiring Parallel::ForkManager.  fork()/waitpid() are core Perl
    # functions and are available on normal Linux/HPC systems.
    my %children;
    my $failed = 0;

    JOB: for my $job (@$jobs) {
        while (keys(%children) >= $max_threads) {
            my $pid = wait();
            last if $pid == -1;
            my $status = $?;
            delete $children{$pid};
            $failed++ if $status != 0;
        }
        die "ERROR: previous PairwiseSV child process failed; stop launching new jobs.\n" if $failed;

        my $pid = fork();
        die "ERROR: cannot fork PairwiseSV worker: $!\n" unless defined $pid;
        if ($pid == 0) {
            eval {
                call_sv_fast($job->[0], $job->[1], $job->[2], $job->[3]);
                1;
            } or do {
                my $err = $@ || 'unknown error';
                warn $err;
                exit 1;
            };
            exit 0;
        } else {
            $children{$pid} = 1;
        }
    }

    while (keys %children) {
        my $pid = wait();
        last if $pid == -1;
        my $status = $?;
        delete $children{$pid};
        $failed++ if $status != 0;
    }

    die "ERROR: one or more PairwiseSV worker processes failed.\n" if $failed;
}

sub usage {
    return <<"END";
Usage:
  perl 4_PairwiseSV_fast_seq.pl --synnet B73vsTzi8.chain.filter.tnet.synnet.refiltered.synnet --outdir ./CleanSV
  perl 4_PairwiseSV_fast_seq.pl --syndir /path/to/CleanSynnet --outdir ./CleanSV --t 8

With sequence extraction:
  perl 4_PairwiseSV_fast_seq.pl --syndir /path/to/CleanSynnet --gd /path/to/genomes --outdir ./CleanSV --t 8 --seq

Required:
  --synnet <file>       one Chain/Net/SynNet file
  OR
  --syndir <dir>        folder containing Chain/Net/SynNet files

Optional:
  --outdir <dir>        output folder [default: .]
  --t <int>             number of parallel files to process [default: 4]
  --gd <dir>            folder containing genome FASTA files. Used by --seq.
  --ref_fa <file>       reference FASTA for --seq, useful for one --synnet file
  --qry_fa <file>       query FASTA for --seq, useful for one --synnet file
  --seq                 output real flanking/SV sequences instead of '.'
  --exact_seq_len       extract exactly Ref_end-Ref_beg and Sample_end-Sample_beg.
                        Without this option, SV sequences mimic the original script
                        and extract one extra base at the right side.
  --seq_context <int>   flanking sequence length [default: 20]
  --max_full_seq <int>  max SV length for full sequence extraction [default: 1000000]
  --edge_seq_len <int>  for very large SVs, output left/right edges of this length [default: 500]
  --invlen <int>        minimum inversion size [default: 1000]
  --cnvlen <int>        minimum CNV/tandem-duplication-like event size [default: 500]
  --dellen <int>        maximum deletion size to keep as simple DEL [default: 20000]
  --inslen <int>        maximum insertion size to keep as simple INS [default: same as dellen]
  --max_opposite_ratio <float>
                        maximum opposite-side / main-side length ratio allowed
                        for simple DEL/INS [default: 0.10].
                        DEL requires Query_gap / Ref_gap <= this value;
                        INS requires Ref_gap / Query_gap <= this value.
                        Example: with 0.10, Ref_gap=100 bp allows Query_gap <=10 bp.
  --max_opposite_len <int>
                        optional absolute cap for the opposite-side length.
                        If set together with --max_opposite_ratio, both criteria
                        must be satisfied. This option is kept for compatibility
                        and is not used by default.
  --trllen <int>        minimum translocation size [default: 5000]
  --filltoplen <int>    minimum top-fill aligned bp to parse [default: 20000]
  --fill2len <int>      minimum secondary-fill aligned bp to parse [default: 10000]
  --rate <float>        minimum alignable proportion for local CNV [default: 0.25]
  --max_level <int>     nested fill levels to parse [default: 3]
  --keep_tmp            keep unsorted temporary files
  --quiet               suppress progress messages
  --help                print this help

Notes:
  1. Coordinates in output are 0-based, half-open, consistent with BED style.
  2. Without --seq, the script is fastest and writes '.' for sequence columns.
  3. With --seq, FASTA files must be uncompressed. If .fai files are absent,
     this script builds an in-memory FASTA index by scanning the FASTA once.
  4. DEL and INS are written separately:
       *.del.sort for deletions and *.ins.sort for insertions.
  5. --seq uses random-access extraction, so it does not load whole maize genomes
     into RAM. Runtime will still increase because many sequence intervals are queried.
END
}

sub infer_pair_from_file {
    my ($file) = @_;
    my $base = basename($file);
    if ($base =~ /^(.+?)vs(.+?)\.chain\.filter\.tnet\.synnet(?:\.refiltered\.synnet)?$/) {
        return ($1, $2);
    } elsif ($base =~ /^(.+?)vs(.+?)\.chain\.filter\.qnet\.synnet(?:\.refiltered\.synnet)?$/) {
        return ($2, $1);
    } elsif ($base =~ /^(.+?)vs(.+?)\./) {
        return ($1, $2);
    }
    return (undef, undef);
}

sub output_dir_for_file {
    my ($file, $outdir) = @_;
    my $base = basename($file);
    if ($base =~ /^(.+?)vs(.+?)\.chain\.filter\.tnet\.synnet(?:\.refiltered\.synnet)?$/) {
        return "$outdir/$1.$2.rawSV";
    } elsif ($base =~ /^(.+?)vs(.+?)\.chain\.filter\.qnet\.synnet(?:\.refiltered\.synnet)?$/) {
        return "$outdir/$2.$1.rawSV";
    } elsif ($base =~ /^(.+?)vs(.+?)\./) {
        return "$outdir/$1.$2.rawSV";
    } else {
        $base =~ s/\.synnet$//;
        return "$outdir/$base.rawSV";
    }
}

sub find_genome_path {
    my ($dir, $name) = @_;
    return undef unless defined $dir && defined $name;
    my @candidates = (
        "$dir/$name",
        "$dir/$name.fa",
        "$dir/$name.fasta",
        "$dir/$name.fna",
        "$dir/$name.FA",
        "$dir/$name.FASTA",
        "$dir/$name.FNA",
    );
    for my $p (@candidates) {
        return $p if -s $p;
    }
    return undef;
}

sub call_sv_fast {
    my ($synnet_file, $output, $ref_fasta, $qry_fasta) = @_;
    $output =~ s/\/$//;
    make_path($output) unless -d $output;

    my $base = basename($synnet_file);
    warn "[INFO] Processing $base\n" unless $quiet;

    my $SEQ;
    if ($do_seq) {
        warn "[INFO] Building/loading FASTA index for reference: $ref_fasta\n" unless $quiet;
        my $ref_db = FastFAI->new($ref_fasta);
        warn "[INFO] Building/loading FASTA index for query: $qry_fasta\n" unless $quiet;
        my $qry_db = FastFAI->new($qry_fasta);
        $SEQ = { ref => $ref_db, qry => $qry_db };
    }

    my $del_tmp  = "$output/$base.del.tmp";
    my $ins_tmp  = "$output/$base.ins.tmp";
    my $inv_tmp  = "$output/$base.inv.tmp";
    my $cxsv_tmp = "$output/$base.cxsv.tmp";
    my $cnv_tmp  = "$output/$base.cnv.tmp";
    my $trl_tmp  = "$output/$base.trl.tmp";

    open(my $DEL, ">", $del_tmp)  or die "ERROR: cannot write $del_tmp: $!\n";
    open(my $INS, ">", $ins_tmp)  or die "ERROR: cannot write $ins_tmp: $!\n";
    open(my $INV, ">", $inv_tmp)  or die "ERROR: cannot write $inv_tmp: $!\n";
    open(my $CXS, ">", $cxsv_tmp) or die "ERROR: cannot write $cxsv_tmp: $!\n";
    open(my $CNV, ">", $cnv_tmp)  or die "ERROR: cannot write $cnv_tmp: $!\n";
    open(my $TRL, ">", $trl_tmp)  or die "ERROR: cannot write $trl_tmp: $!\n";

    my %FH = (DEL=>$DEL, INS=>$INS, INV=>$INV, CXS=>$CXS, CNV=>$CNV, TRL=>$TRL);

    open(my $IN, "<", $synnet_file) or die "ERROR: cannot read $synnet_file: $!\n";

    my $current_ref_chr;
    my $root = { kind => 'root', indent => -1 };
    my $top_fill;
    my @stack = ($root);
    my $n_top = 0;

    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/ || $line =~ /^#/;

        if ($line =~ /^net\s+(\S+)\s+(\d+)/) {
            process_top_fill($top_fill, \%FH, $SEQ) if defined $top_fill;
            undef $top_fill;
            $current_ref_chr = $1;
            $root = { kind => 'root', indent => -1, ref_chr => $current_ref_chr };
            @stack = ($root);
            next;
        }

        next unless $line =~ /^(\s*)(fill|gap)\s+/;
        my $indent = length($1);
        my $kind = $2;

        if ($kind eq 'fill' && $indent <= 1) {
            process_top_fill($top_fill, \%FH, $SEQ) if defined $top_fill;
            $top_fill = parse_node($line, $indent, $current_ref_chr);
            @stack = ($root, $top_fill);
            $n_top++;
            next;
        }

        next unless defined $top_fill;
        my $node = parse_node($line, $indent, $current_ref_chr);
        while (@stack && $stack[-1]{indent} >= $indent) { pop @stack; }
        my $parent = $stack[-1] || $top_fill;
        push @{ $parent->{children} }, $node;
        push @stack, $node;
    }
    close $IN;

    process_top_fill($top_fill, \%FH, $SEQ) if defined $top_fill;

    close $DEL; close $INS; close $INV; close $CXS; close $CNV; close $TRL;

    my $head1 = join("\t", split(/ /, "#Reference Ref_beg Ref_end SVtypes Sample Sample_beg Sample_end SVlength Strand 5_flanking 3_flanking Ref_5_flanking_20bp Ref_3_flanking_20bp Sample_5_flanking_20bp Sample_3_flanking_20bp Ref_SV_sequence Sample_SV_sequence")) . "\n";
    my $head2 = join("\t", split(/ /, "#Reference Ref_beg Ref_end Sample Sample_beg Sample_end Strand Ref_len Sample_len TRLtypes Level")) . "\n";

    sort_with_header($del_tmp,  "$output/$base.del.sort",  $head1, 1);
    sort_with_header($ins_tmp,  "$output/$base.ins.sort",  $head1, 1);
    sort_with_header($inv_tmp,  "$output/$base.inv.sort",  $head1, 1);
    sort_with_header($cxsv_tmp, "$output/$base.cxsv.sort", $head1, 1);
    sort_with_header($cnv_tmp,  "$output/$base.cnv.sort",  $head1, 1);
    sort_with_header($trl_tmp,  "$output/$base.trl.sort",  $head2, 1);

    unlink $del_tmp, $ins_tmp, $inv_tmp, $cxsv_tmp, $cnv_tmp, $trl_tmp unless $keep_tmp;
    warn "[INFO] Finished $base; parsed $n_top top fills\n" unless $quiet;
}

sub parse_node {
    my ($line, $indent, $ref_chr) = @_;
    $line =~ s/^\s+//;
    my @f = split /\s+/, $line;
    my $kind = shift @f;

    my %n = (
        kind     => $kind,
        indent   => $indent,
        ref_chr  => $ref_chr,
        children => [],
    );

    if ($kind eq 'fill' || $kind eq 'gap') {
        @n{qw(t_beg t_len q_chr strand q_beg q_len)} = @f[0..5];
        $n{t_beg} += 0; $n{t_len} += 0; $n{q_beg} += 0; $n{q_len} += 0;
        $n{t_end} = $n{t_beg} + $n{t_len};
        $n{q_end} = $n{q_beg} + $n{q_len};
        if ($kind eq 'fill') {
            my %attr;
            for (my $i=6; $i < @f-1; $i += 2) {
                $attr{$f[$i]} = $f[$i+1];
            }
            $n{ali}   = exists $attr{ali}   ? $attr{ali}+0   : $n{t_len};
            $n{type}  = exists $attr{type}  ? $attr{type}    : 'NA';
            $n{qFar}  = exists $attr{qFar}  ? $attr{qFar}+0  : 0;
            $n{qOver} = exists $attr{qOver} ? $attr{qOver}+0 : 0;
            $n{score} = exists $attr{score} ? $attr{score}+0 : 0;
            $n{qDup}  = exists $attr{qDup}  ? $attr{qDup}+0  : 0;
        }
    }
    return \%n;
}

sub process_top_fill {
    my ($fill, $FH, $SEQ) = @_;
    return unless defined $fill;
    return if ($fill->{ali} || 0) < $toplen;
    process_fill_branch($fill, $fill->{ref_chr}, 1, $FH, $SEQ);
}

sub process_fill_branch {
    my ($fill, $ref_chr, $level, $FH, $SEQ) = @_;
    return unless defined $fill && $fill->{kind} eq 'fill';

    my $fill_beg = $fill->{t_beg};
    my $fill_end = $fill->{t_end};
    my $tlen     = $fill->{t_len};
    my $qlen     = $fill->{q_len};
    my $qbeg     = $fill->{q_beg};
    my $qend     = $fill->{q_end};
    my $qchr     = $fill->{q_chr};
    my $strand   = $fill->{strand};
    my $type     = $fill->{type} || 'NA';

    # Potential translocations / non-syntenic fills.
    if ($type =~ /top/) {
        if ($tlen > 10 * $trllen) {
            print_trl($FH->{TRL}, $ref_chr, $fill_beg, $fill_end, $qchr, $qbeg, $qend, $strand, $tlen, $qlen, 'ChrEND', $type);
        }
    } elsif ($type =~ /syn|inv/) {
        if ($tlen > $trllen && ($fill->{qFar} || 0) > 2 * $tlen) {
            print_trl($FH->{TRL}, $ref_chr, $fill_beg, $fill_end, $qchr, $qbeg, $qend, $strand, $tlen, $qlen, 'IntraChr', $type);
        }
    } elsif ($type =~ /nonSyn/) {
        if ($tlen > $trllen) {
            print_trl($FH->{TRL}, $ref_chr, $fill_beg, $fill_end, $qchr, $qbeg, $qend, $strand, $tlen, $qlen, 'InterChr', $type);
        }
    }

    my @gaps = grep { $_->{kind} eq 'gap' } @{ $fill->{children} || [] };
    for (my $i=0; $i<@gaps; $i++) {
        my $prev = ($i == 0)       ? { t_beg => $fill_beg, t_len => 0, t_end => $fill_beg } : $gaps[$i-1];
        my $next = ($i == $#gaps)  ? { t_beg => $fill_end, t_len => 0, t_end => $fill_end } : $gaps[$i+1];
        my @downfills = grep { $_->{kind} eq 'fill' } @{ $gaps[$i]->{children} || [] };
        classify_gap($gaps[$i], $ref_chr, $prev, $next, \@downfills, $FH, $SEQ);
    }

    # Recursively parse secondary/tertiary fill branches.
    return if $level >= $max_level;
    for my $gap (@gaps) {
        for my $child (grep { $_->{kind} eq 'fill' } @{ $gap->{children} || [] }) {
            next if ($child->{ali} || 0) < $secondlen;
            process_fill_branch($child, $ref_chr, $level + 1, $FH, $SEQ);
        }
    }
}

sub classify_gap {
    my ($gap, $ref_chr, $prev, $next, $downfills, $FH, $SEQ) = @_;

    my $t_beg = $gap->{t_beg};
    my $t_end = $gap->{t_end};
    my $t_len = $gap->{t_len};
    my $q_chr = $gap->{q_chr};
    my $q_beg = $gap->{q_beg};
    my $q_end = $gap->{q_end};
    my $q_len = $gap->{q_len};
    my $strand = $gap->{strand};

    my $forward  = $t_beg - $prev->{t_beg} - $prev->{t_len};
    my $backward = $next->{t_beg} - $t_beg - $t_len;

    my @seqcols = get_seq_columns($SEQ, $ref_chr, $t_beg, $t_end, $q_chr, $q_beg, $q_end, $strand);

    # Ratio-based simple-INDEL definition.
    # DEL: reference gap has the major sequence; query gap is absent or only
    #      a short remnant relative to the reference gap.
    # INS: query gap has the major sequence; reference gap is absent or only
    #      a short remnant relative to the query gap.
    # If the opposite-side remnant exceeds --max_opposite_ratio, the event is
    # considered replacement-like/complex rather than a simple DEL/INS.
    my $max_len = ($t_len > $q_len) ? $t_len : $q_len;
    my $is_simple_del = 0;
    my $is_simple_ins = 0;
    if ($t_len > 0) {
        $is_simple_del = 1 if (($q_len / $t_len) <= $max_opposite_ratio);
        $is_simple_del = 0 if (defined $max_opposite_len && $q_len > $max_opposite_len);
    }
    if ($q_len > 0) {
        $is_simple_ins = 1 if (($t_len / $q_len) <= $max_opposite_ratio);
        $is_simple_ins = 0 if (defined $max_opposite_len && $t_len > $max_opposite_len);
    }

    if (!@$downfills) {
        if ($is_simple_del) {
            my $svt = ($t_len > 10000) ? 'DEL_PAV' : 'DEL';
            print_event($FH->{DEL}, $ref_chr, $t_beg, $t_end, $svt, $q_chr, $q_beg, $q_end, $t_len, $strand, $forward, $backward, @seqcols);
        } elsif ($is_simple_ins) {
            my $svt = ($q_len > 10000) ? 'INS_PAV' : 'INS';
            print_event($FH->{INS}, $ref_chr, $t_beg, $t_end, $svt, $q_chr, $q_beg, $q_end, $q_len, $strand, $forward, $backward, @seqcols);
        } else {
            print_event($FH->{CXS}, $ref_chr, $t_beg, $t_end, 'CMPLX_PAV', $q_chr, $q_beg, $q_end, $max_len, $strand, $forward, $backward, @seqcols);
        }
        return;
    }

    my ($ali, $inv_len, $total_cov) = (0, 0, 0);
    my (@q_coords, @t_coords);

    for my $f (@$downfills) {
        $total_cov += $f->{t_len};
        my $local_limit = 2 * $max_len;
        if (($f->{type} || '') =~ /syn|inv/ && ($f->{qFar} || 0) < $local_limit) {
            my $a = $f->{ali} || 0;
            $ali += $a;
            $inv_len += $a if ($f->{type} || '') =~ /inv/;
            push @q_coords, $f->{q_beg}, $f->{q_end};
            push @t_coords, $f->{t_beg}, $f->{t_end};
        }

        if (($f->{type} || '') =~ /inv/ && ($f->{qFar} || 0) < 2 * (($t_len > 0) ? $t_len : 1)
            && (($f->{qOver} || 0) >= $invlen || ($f->{t_len} > 2*$invlen && $f->{q_len} > 2*$invlen))
            && ($f->{ali} || 0) >= $invlen) {
            print_event($FH->{INV}, $ref_chr, $f->{t_beg}, $f->{t_end}, 'INV', $f->{q_chr}, $f->{q_beg}, $f->{q_end}, $f->{t_len}, $gap->{strand}, $forward, $backward, @seqcols);
        }
    }

    @q_coords = sort { $a <=> $b } @q_coords;
    @t_coords = sort { $a <=> $b } @t_coords;
    my ($q_min, $q_max) = @q_coords ? ($q_coords[0], $q_coords[-1]) : ($q_beg, $q_end);

    my $cnv_like = 0;
    my $inv_cnv_like = 0;
    if ($max_len > 0) {
        $inv_cnv_like = 1 if ($inv_len > ($q_len + $t_len)/4 && $inv_len > $cnvlen);
        # Use the larger local gap as the denominator. This avoids calling
        # replacement-like events as DEL_CNV or INS_CNV in the DEL/INS files.
        $cnv_like = 1 if (!$inv_cnv_like && $ali/$max_len > $ali_rate && $ali > $cnvlen);
    }

    if ($inv_cnv_like) {
        print_event($FH->{CNV}, $ref_chr, $t_beg, $t_end, 'InvCNV', $q_chr, $q_min, $q_max, $max_len, $strand, $forward, $backward, @seqcols);
    } elsif ($cnv_like) {
        print_event($FH->{CNV}, $ref_chr, $t_beg, $t_end, 'CNV', $q_chr, $q_min, $q_max, $max_len, $strand, $forward, $backward, @seqcols);
    }

    # Only ratio-qualified one-sided gaps are allowed into the DEL/INS files.
    # When the opposite-side remnant is too large relative to the main side,
    # the event is output as complex even if one side is much longer.
    if ($is_simple_del && !$cnv_like && !$inv_cnv_like) {
        if ($t_len <= $dellen) {
            my $svt = ($t_len > 10000 && $total_cov/($t_len || 1) < 0.1) ? 'DEL_PAV' : 'DEL';
            print_event($FH->{DEL}, $ref_chr, $t_beg, $t_end, $svt, $q_chr, $q_beg, $q_end, $t_len, $strand, $forward, $backward, @seqcols);
        } else {
            my $svt = ($t_len > 10000 && $total_cov/($t_len || 1) < 0.1) ? 'CMPLX_PAV' : 'CMPLX';
            print_event($FH->{CXS}, $ref_chr, $t_beg, $t_end, $svt, $q_chr, $q_beg, $q_end, $t_len, $strand, $forward, $backward, @seqcols);
        }
    } elsif ($is_simple_ins && !$cnv_like && !$inv_cnv_like) {
        if ($q_len <= $inslen) {
            my $svt = ($q_len > 10000 && $total_cov/($q_len || 1) < 0.1) ? 'INS_PAV' : 'INS';
            print_event($FH->{INS}, $ref_chr, $t_beg, $t_end, $svt, $q_chr, $q_beg, $q_end, $q_len, $strand, $forward, $backward, @seqcols);
        } else {
            my $svt = ($q_len > 10000 && $total_cov/($q_len || 1) < 0.1) ? 'CMPLX_PAV' : 'CMPLX';
            print_event($FH->{CXS}, $ref_chr, $t_beg, $t_end, $svt, $q_chr, $q_beg, $q_end, $q_len, $strand, $forward, $backward, @seqcols);
        }
    } else {
        my $den = $max_len || 1;
        my $svt = ($max_len > 10000 && $total_cov/$den < 0.1) ? 'CMPLX_PAV' : 'CMPLX';
        print_event($FH->{CXS}, $ref_chr, $t_beg, $t_end, $svt, $q_chr, $q_beg, $q_end, $max_len, $strand, $forward, $backward, @seqcols);
    }
}

sub get_seq_columns {
    my ($SEQ, $ref_chr, $rb, $re, $q_chr, $qb, $qe, $strand) = @_;
    return ('.', '.', '.', '.', '.', '.') unless $do_seq && defined $SEQ;

    my $extra = $exact_seq_len ? 0 : 1;  # original script used t_len+1/q_len+1

    my $ref5 = $SEQ->{ref}->get_seq($ref_chr, $rb - $seq_context, $rb);
    my $ref3 = $SEQ->{ref}->get_seq($ref_chr, $re, $re + $seq_context);
    my $qry5 = $SEQ->{qry}->get_seq($q_chr, $qb - $seq_context, $qb);
    my $qry3 = $SEQ->{qry}->get_seq($q_chr, $qe, $qe + $seq_context);

    my $ref_sv = get_sv_sequence($SEQ->{ref}, $ref_chr, $rb, $re + $extra);
    my $qry_sv = get_sv_sequence($SEQ->{qry}, $q_chr, $qb, $qe + $extra);

    if (defined $strand && $strand eq '-') {
        $qry_sv = revcom($qry_sv) if defined $qry_sv && $qry_sv ne '.';
        $qry5   = revcom($qry5)   if defined $qry5   && $qry5 ne '.';
        $qry3   = revcom($qry3)   if defined $qry3   && $qry3 ne '.';
    }

    for ($ref5, $ref3, $qry5, $qry3, $ref_sv, $qry_sv) {
        $_ = '.' if !defined($_) || $_ eq '';
    }
    return ($ref5, $ref3, $qry5, $qry3, $ref_sv, $qry_sv);
}

sub get_sv_sequence {
    my ($db, $chr, $start, $end) = @_;
    my $len = $end - $start;
    return '.' if $len <= 0;
    if ($len <= $max_full_seq) {
        return $db->get_seq($chr, $start, $end);
    }
    my $left  = $db->get_seq($chr, $start, $start + $edge_seq_len);
    my $right = $db->get_seq($chr, $end - $edge_seq_len, $end);
    $left  = '.' if !defined($left) || $left eq '';
    $right = '.' if !defined($right) || $right eq '';
    return join("NNNNNNNNNN", $left, $right);
}

sub print_event {
    my ($fh, $ref, $rb, $re, $svt, $sample, $qb, $qe, $len, $strand, $five, $three, @seqcols) = @_;
    $rb = 0 if !defined($rb) || $rb < 0;
    $re = $rb if !defined($re) || $re < $rb;
    $qb = 0 if !defined($qb) || $qb < 0;
    $qe = $qb if !defined($qe) || $qe < $qb;
    $len = 0 if !defined $len;
    @seqcols = ('.', '.', '.', '.', '.', '.') unless @seqcols == 6;
    print $fh join("\t", $ref, $rb, $re, $svt, $sample, $qb, $qe, $len, $strand, $five, $three, @seqcols), "\n";
}

sub print_trl {
    my ($fh, $ref, $rb, $re, $sample, $qb, $qe, $strand, $rlen, $qlen, $trltype, $level) = @_;
    print $fh join("\t", $ref, $rb, $re, $sample, $qb, $qe, $strand, $rlen, $qlen, $trltype, $level), "\n";
}

sub sort_with_header {
    my ($tmp, $final, $header, $sort_by_ref_start) = @_;
    open(my $OUT, ">", $final) or die "ERROR: cannot write $final: $!\n";
    print $OUT $header;
    if (-s $tmp) {
        if ($sort_by_ref_start) {
            open(my $SORT, "-|", "sort", "-k1,1", "-k2,2n", $tmp) or die "ERROR: cannot run sort: $!\n";
            while (my $line = <$SORT>) { print $OUT $line; }
            close $SORT;
        } else {
            open(my $IN, "<", $tmp) or die "ERROR: cannot read $tmp: $!\n";
            while (my $line = <$IN>) { print $OUT $line; }
            close $IN;
        }
    }
    close $OUT;
}

sub revcom {
    my $seq = shift;
    return $seq unless defined $seq && $seq ne '.';
    my $revcomp = reverse $seq;
    $revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
    return $revcomp;
}

################################################################################
# Lightweight FASTA random access. It can read existing .fai files, or build an
# in-memory index by scanning the FASTA once. The FASTA must be uncompressed.
################################################################################
package FastFAI;
use strict;
use warnings;
use Fcntl qw(SEEK_SET);

sub new {
    my ($class, $fasta) = @_;
    die "ERROR: FASTA path is not defined.\n" unless defined $fasta;
    die "ERROR: FASTA does not exist or is empty: $fasta\n" unless -s $fasta;
    die "ERROR: --seq currently requires uncompressed FASTA, but got: $fasta\n" if $fasta =~ /\.gz$/;

    my $self = bless { fasta => $fasta, idx => {}, cache => {}, warned => {} }, $class;
    if (-s "$fasta.fai") {
        $self->_load_fai("$fasta.fai");
    } else {
        $self->_build_index();
    }
    open(my $fh, "<", $fasta) or die "ERROR: cannot open FASTA $fasta: $!\n";
    binmode($fh);
    $self->{fh} = $fh;
    return $self;
}

sub _load_fai {
    my ($self, $fai) = @_;
    open(my $IN, "<", $fai) or die "ERROR: cannot read $fai: $!\n";
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my @f = split /\t/, $line;
        next unless @f >= 5;
        $self->{idx}{$f[0]} = {
            len        => $f[1] + 0,
            offset     => $f[2] + 0,
            line_bases => $f[3] + 0,
            line_width => $f[4] + 0,
        };
    }
    close $IN;
}

sub _build_index {
    my ($self) = @_;
    my $fasta = $self->{fasta};
    open(my $IN, "<", $fasta) or die "ERROR: cannot read $fasta: $!\n";
    binmode($IN);

    my ($name, $len, $offset, $line_bases, $line_width);
    while (1) {
        my $pos = tell($IN);
        my $line = <$IN>;
        last unless defined $line;
        if ($line =~ /^>(\S+)/) {
            if (defined $name) {
                $self->{idx}{$name} = {
                    len        => $len,
                    offset     => $offset,
                    line_bases => $line_bases,
                    line_width => $line_width,
                };
            }
            $name = $1;
            ($len, $offset, $line_bases, $line_width) = (0, undef, undef, undef);
        } else {
            next unless defined $name;
            my $raw_width = length($line);
            $line =~ s/[\r\n]+$//;
            my $bases = length($line);
            next if $bases == 0;
            if (!defined $offset) {
                $offset = $pos;
                $line_bases = $bases;
                $line_width = $raw_width;
            }
            $len += $bases;
        }
    }
    if (defined $name) {
        $self->{idx}{$name} = {
            len        => $len,
            offset     => $offset,
            line_bases => $line_bases,
            line_width => $line_width,
        };
    }
    close $IN;
}

sub get_seq {
    my ($self, $chr, $start, $end) = @_;
    return '.' unless defined $chr && defined $start && defined $end;
    my $real_chr = $self->_resolve_chr($chr);
    return '.' unless defined $real_chr;
    my $idx = $self->{idx}{$real_chr};
    return '.' unless defined $idx;

    $start = int($start);
    $end   = int($end);
    return '' if $end <= $start;

    my $left_pad = 0;
    my $right_pad = 0;
    if ($start < 0) {
        $left_pad = -$start;
        $start = 0;
    }
    if ($end > $idx->{len}) {
        $right_pad = $end - $idx->{len};
        $end = $idx->{len};
    }
    my $seq = '';
    $seq .= 'N' x $left_pad if $left_pad > 0;
    $seq .= $self->_read_interval($idx, $start, $end) if $end > $start;
    $seq .= 'N' x $right_pad if $right_pad > 0;
    return $seq;
}

sub _resolve_chr {
    my ($self, $chr) = @_;
    return $self->{cache}{$chr} if exists $self->{cache}{$chr};
    if (exists $self->{idx}{$chr}) {
        return $self->{cache}{$chr} = $chr;
    }

    my @try;
    if ($chr =~ /^[^.]+\.(.+)$/) {
        push @try, $1;          # B73.chr01 -> chr01
    }
    if ($chr =~ /(chr\S+)$/) {
        push @try, $1;          # any.prefix.chr01 -> chr01
    }
    if ($chr =~ /^chr0([1-9])$/) {
        push @try, "chr$1";     # chr01 -> chr1
    }
    if ($chr =~ /^chr([1-9])$/) {
        push @try, "chr0$1";    # chr1 -> chr01
    }

    for my $c (@try) {
        if (exists $self->{idx}{$c}) {
            return $self->{cache}{$chr} = $c;
        }
    }

    # As a last resort, match by suffix. This is useful when FASTA headers keep a sample prefix.
    for my $c (keys %{ $self->{idx} }) {
        if ($c =~ /\Q$chr\E$/ || $chr =~ /\Q$c\E$/) {
            return $self->{cache}{$chr} = $c;
        }
    }

    warn "[WARN] chromosome '$chr' not found in $self->{fasta}; sequence columns will be '.'.\n"
        unless $self->{warned}{$chr}++;
    $self->{cache}{$chr} = undef;
    return undef;
}

sub _read_interval {
    my ($self, $idx, $start, $end) = @_;
    return '' if $end <= $start;
    my $lb = $idx->{line_bases};
    my $lw = $idx->{line_width};
    my $off = $idx->{offset};
    my $fh = $self->{fh};

    my $pos = $start;
    my $remain = $end - $start;
    my $seq = '';
    while ($remain > 0) {
        my $line_no = int($pos / $lb);
        my $in_line = $pos % $lb;
        my $take = $lb - $in_line;
        $take = $remain if $take > $remain;
        my $byte_pos = $off + $line_no * $lw + $in_line;
        seek($fh, $byte_pos, SEEK_SET) or die "ERROR: seek failed in $self->{fasta}: $!\n";
        my $buf = '';
        my $got = read($fh, $buf, $take);
        die "ERROR: read failed in $self->{fasta}: $!\n" unless defined $got;
        $seq .= $buf;
        $pos += $take;
        $remain -= $take;
    }
    $seq =~ s/\s+//g;
    return $seq;
}

1;
