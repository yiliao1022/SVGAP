#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/pub";
use lib "$Bin/../lib";
use lib "$Bin/../pub";
use lib "$Bin/../../lib";
use lib "$Bin/../../pub";
use File::Path qw(make_path);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

####################################################################################
# Name: 01_Convert2Axt.pl
#
# Purpose:
#   1) Convert whole-genome alignment outputs from LAST, LASTZ, MUMmer,
#      minimap2 and AnchorWave to AXT format.
#   2) Optionally split FASTA sequences into overlapping subcontigs for faster
#      or more parallel alignment.
#   3) For alignments generated against split query subcontigs, reformat MAF
#      coordinates back to the original chromosome coordinates before mafToAxt.
#
# Integrated from:
#   - 1_Convert2Axt.pl
#   - SplitFa.pl
#   - Mafreformat.pl
#
# Original author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Standardized date: 2026-06-14
####################################################################################

my %opt = (
    ali        => '',
    input      => '.',
    wk         => '.',
    tname      => 'Ref',
    qname      => 'Que',
    mode       => 'auto',   # auto|chr|ctg. auto detects .chr/.ctg suffix; fallback is chr.
    t          => 4,
    query      => 0,        # produce query-as-reference AXT as in the old --query option
    keep_tmp   => 0,
    emit_variants => 0,     # optional paftools call output; unused by the SVGAP pipeline

    # External tool directories or executable paths.
    mum320     => "$Bin/pub/MUMMER/",
    ucsc       => "$Bin/pub/UCSC/",
    tba        => "$Bin/pub/TBA/",
    k8         => "$Bin/pub/MINIMAP2/",
    paftools   => "$Bin/pub/MINIMAP2/",

    # FASTA splitting options.
    splitfa    => [],       # repeatable: --splitfa genome.fa --splitfa another.fa
    bin        => 10_000_000,
    slide      => 10_000,
    split_out  => '',

    # Size files for coordinate reformatting.  Usually produced by --splitfa.
    sizes      => [],       # repeatable or comma-separated

    help       => 0,
);

GetOptions(
    'ali=s'        => \$opt{ali},
    'input=s'      => \$opt{input},
    'wk=s'         => \$opt{wk},
    'tname=s'      => \$opt{tname},
    'qname=s'      => \$opt{qname},
    'mode=s'       => \$opt{mode},
    'coord=s'      => \$opt{mode},
    't=i'          => \$opt{t},
    'threads=i'    => \$opt{t},   # alias for compatibility with SVGAP_pipeline.pl
    'query!'       => \$opt{query},
    'keep_tmp!'    => \$opt{keep_tmp},
    'emit_variants|emit-variants!' => \$opt{emit_variants},

    'mum320=s'     => \$opt{mum320},
    'ucsc=s'       => \$opt{ucsc},
    'tba=s'        => \$opt{tba},
    'k8=s'         => \$opt{k8},
    'paftools=s'   => \$opt{paftools},

    'splitfa=s@'   => \$opt{splitfa},
    'fasta=s@'     => \$opt{splitfa},
    'bin=i'        => \$opt{bin},
    'slide=i'      => \$opt{slide},
    'split_out=s'  => \$opt{split_out},

    'sizes=s@'     => \$opt{sizes},
    'qsize=s@'     => \$opt{sizes},
    'size=s@'      => \$opt{sizes},

    'help|h'       => \$opt{help},
) or die usage();

$opt{ali} = lc($opt{ali} || '');
$opt{mode} = lc($opt{mode} || 'auto');

if ($opt{help}) {
    print usage();
    exit 0;
}

validate_options(\%opt);
make_path($opt{wk}) unless -d $opt{wk};
$opt{split_out} ||= $opt{wk};
make_path($opt{split_out}) unless -d $opt{split_out};

# Step 1: optionally split FASTA files.  If --ali is not provided, the script
# behaves as a FASTA-splitting utility and exits after this step.
if (@{ $opt{splitfa} }) {
    for my $fa (@{ $opt{splitfa} }) {
        split_fasta($fa, $opt{split_out}, $opt{bin}, $opt{slide});
    }
    if (!$opt{ali}) {
        print STDERR "Done. FASTA splitting finished.\n";
        exit 0;
    }
}

# Step 2: convert alignment outputs to AXT.
if (!$opt{ali}) {
    die "ERROR: --ali is required for alignment-to-AXT conversion. Use --help for examples.\n";
}

my $target_t_dir = "$opt{wk}/Target_$opt{tname}";
my $target_q_dir = "$opt{wk}/Target_$opt{qname}";
make_path($target_t_dir) unless -d $target_t_dir;
make_path($target_q_dir) if $opt{query} && !-d $target_q_dir;

my %sizes = collect_size_files(\%opt);
my @jobs  = collect_alignment_jobs(\%opt);
die "ERROR: no alignment files found in $opt{input} for --ali $opt{ali}.\n" unless @jobs;

print STDERR "Found ".scalar(@jobs)." alignment file(s).\n";
run_jobs(\@jobs, $opt{t}, sub {
    my ($job) = @_;
    convert_one_alignment($job, \%opt, \%sizes, $target_t_dir, $target_q_dir);
});

print STDERR "Done. AXT files are in $target_t_dir" . ($opt{query} ? " and $target_q_dir" : "") . ".\n";

####################################################################################
# Option validation and usage
####################################################################################
sub validate_options {
    my ($opt) = @_;
    die "ERROR: --mode/--coord must be auto, chr or ctg.\n"
        unless $opt->{mode} =~ /^(auto|chr|ctg)$/;
    die "ERROR: --t must be >= 1.\n" if $opt->{t} < 1;
    die "ERROR: --bin must be > 0.\n" if $opt->{bin} <= 0;
    die "ERROR: --slide must be >= 0.\n" if $opt->{slide} < 0;
    die "ERROR: --slide must be smaller than --bin.\n" if $opt->{slide} >= $opt->{bin};

    if ($opt->{ali}) {
        die "ERROR: --ali must be one of last, lastz, maf, mummer/mummer4, minimap2 or anchorwave.\n"
            unless $opt->{ali} =~ /^(last|lastz|maf|mummer|mummer4|minimap2|anchorwave)$/;
    }
}

sub usage {
    return <<'USAGE';

Usage 1. Direct chromosome-level alignment conversion, default mode:
  perl 1_Convert2Axt_integrated.pl \
      --ali minimap2 \
      --input ./alignments \
      --wk ./Convert2Axt_out \
      --tname Ref \
      --qname Que \
      --t 8

  Input filenames can be like:
      RefvsQue.chr.paf       for minimap2 chr mode
      RefvsQue.chr.maf       for LAST chr mode
      RefvsQue.maf           for generic MAF/LASTZ/AnchorWave chr mode
      RefvsQue.delta         for MUMmer chr mode

Usage 2. Split a FASTA into overlapping subcontigs before alignment:
  perl 1_Convert2Axt_integrated.pl \
      --splitfa Que.fa \
      --wk ./split_fa \
      --bin 10000000 \
      --slide 10000

  Outputs:
      ./split_fa/Que.subcontigs.fa
      ./split_fa/Que.sizes
      ./split_fa/Que.subcontigs.map.tsv

  The split FASTA headers are:
      >Chr01_1_10010000
      >Chr01_9990001_20000000

Usage 3. Convert alignments generated from split query subcontigs:
  perl 1_Convert2Axt_integrated.pl \
      --ali minimap2 \
      --input ./alignments \
      --wk ./Convert2Axt_out \
      --mode ctg \
      --sizes ./split_fa/Que.sizes \
      --tname Ref \
      --qname Que \
      --t 8

  Input filenames can be like:
      RefvsQue.ctg.paf       for minimap2 ctg mode
      RefvsQue.ctg.maf       for LAST ctg mode

Main options:
  --ali        [last|lastz|maf|mummer|mummer4|minimap2|anchorwave]
  --input      folder containing alignment files
  --wk         output / working folder
  --mode       auto|chr|ctg. Default: auto. If no suffix is detected, chr is used.
  --query      also produce query-as-reference AXT files
  --emit_variants
               also run paftools.js call and keep its diagnostic variant table.
               This is disabled by default because SVGAP does not consume it.
  --tname      target assembly name for output folder name. Default: Ref
  --qname      query assembly name for output folder name. Default: Que
  --t          number of parallel processes. Default: 4

FASTA splitting options:
  --splitfa    FASTA to split. Can be used multiple times.
  --bin        main window size. Default: 10000000
  --slide      overlap extension around boundaries. Default: 10000
  --split_out  folder for split FASTA and size/map files. Default: --wk

Coordinate reformatting options for ctg mode:
  --sizes      original chromosome size file produced by --splitfa.
               Can be used multiple times or as comma-separated paths.
               If omitted, the script also searches --wk and --input for *.sizes.

External tool paths:
  --ucsc       directory containing mafToAxt, or path to mafToAxt
  --tba        directory containing maf_sort/maf_order, or path prefix directory
  --mum320     directory containing delta2maf, or path to delta2maf
  --k8         directory containing k8, or path to k8
  --paftools   directory containing paftools.js, or path to paftools.js

Notes:
  - chr mode is the default behavior for unsplit chromosome-level alignments.
  - ctg mode is for alignments against split query subcontigs. The script converts
    split-subcontig MAF coordinates back to original chromosome coordinates before
    running mafToAxt.
  - Parallel conversion uses Perl's built-in fork/wait and needs no CPAN module.

USAGE
}

####################################################################################
# FASTA splitting
####################################################################################
sub split_fasta {
    my ($fasta, $outdir, $bin, $slide) = @_;
    die "ERROR: cannot read FASTA $fasta\n" unless -s $fasta;
    make_path($outdir) unless -d $outdir;

    my $base = basename($fasta);
    $base =~ s/\.(fa|fasta|fna)$//i;

    my $out_fa  = "$outdir/$base.subcontigs.fa";
    my $size_f  = "$outdir/$base.sizes";
    my $map_f   = "$outdir/$base.subcontigs.map.tsv";

    open(my $IN,  '<', $fasta)  or die "ERROR: cannot read $fasta: $!\n";
    open(my $OUT, '>', $out_fa) or die "ERROR: cannot write $out_fa: $!\n";
    open(my $SZ,  '>', $size_f) or die "ERROR: cannot write $size_f: $!\n";
    open(my $MAP, '>', $map_f)  or die "ERROR: cannot write $map_f: $!\n";
    print $MAP join("\t", qw/SUBCONTIG ORIGINAL_CHR START END LENGTH/), "\n";

    local $/ = "\n>";
    while (my $rec = <$IN>) {
        $rec =~ s/^>//;
        $rec =~ s/>$//;
        next if $rec =~ /^\s*$/;
        chomp $rec;

        my ($head, @seq_lines) = split /\n/, $rec;
        my ($id) = $head =~ /^(\S+)/;
        die "ERROR: failed to parse FASTA header in $fasta\n" unless defined $id && $id ne '';
        my $seq = uc(join('', @seq_lines));
        $seq =~ s/\s+//g;
        my $len = length($seq);
        next if $len == 0;
        print $SZ "$id\t$len\n";

        for (my $i = 0; ; $i++) {
            my $beg0 = ($i == 0) ? 0 : $i * $bin - $slide;
            $beg0 = 0 if $beg0 < 0;
            last if $beg0 >= $len;

            my $readlen = $bin + $slide;
            my $end0 = $beg0 + $readlen;
            $end0 = $len if $end0 > $len;
            my $subseq = substr($seq, $beg0, $end0 - $beg0);
            last if $subseq eq '';

            my $start1 = $beg0 + 1;
            my $end1   = $end0;
            my $sub_id = $id . "_" . $start1 . "_" . $end1;
            print $OUT ">$sub_id\n";
            print_wrapped($OUT, $subseq, 80);
            print $MAP join("\t", $sub_id, $id, $start1, $end1, length($subseq)), "\n";

            last if $end0 >= $len;
        }
    }
    close $IN;
    close $OUT;
    close $SZ;
    close $MAP;

    print STDERR "Split FASTA written: $out_fa\n";
    print STDERR "Original sizes written: $size_f\n";
    print STDERR "Subcontig map written: $map_f\n";
}

sub print_wrapped {
    my ($fh, $seq, $width) = @_;
    $width ||= 80;
    for (my $i = 0; $i < length($seq); $i += $width) {
        print $fh substr($seq, $i, $width), "\n";
    }
}

####################################################################################
# Alignment job discovery and execution
####################################################################################
sub collect_alignment_jobs {
    my ($opt) = @_;
    opendir(my $DIR, $opt->{input}) or die "ERROR: cannot open input directory $opt->{input}: $!\n";
    my @jobs;
    while (my $file = readdir($DIR)) {
        next if $file =~ /^\./;
        my ($t, $q, $suffix_mode);

        if ($opt->{ali} eq 'minimap2') {
            if ($file =~ /^(.+?)vs(.+?)\.(chr|ctg)\.paf$/) {
                ($t, $q, $suffix_mode) = ($1, $2, $3);
            } elsif ($file =~ /^(.+?)vs(.+?)\.paf$/) {
                ($t, $q, $suffix_mode) = ($1, $2, undef);
            } else { next; }
        } elsif ($opt->{ali} eq 'last') {
            if ($file =~ /^(.+?)vs(.+?)\.(chr|ctg)\.maf$/) {
                ($t, $q, $suffix_mode) = ($1, $2, $3);
            } elsif ($file =~ /^(.+?)vs(.+?)\.maf$/) {
                ($t, $q, $suffix_mode) = ($1, $2, undef);
            } else { next; }
        } elsif ($opt->{ali} eq 'lastz' || $opt->{ali} eq 'anchorwave' || $opt->{ali} eq 'maf') {
            if ($file =~ /^(.+?)vs(.+?)(?:\.(chr|ctg))?\.maf$/) {
                ($t, $q, $suffix_mode) = ($1, $2, $3);
            } else { next; }
        } elsif ($opt->{ali} eq 'mummer' || $opt->{ali} eq 'mummer4') {
            if ($file =~ /^(.+?)vs(.+?)(?:\.(chr|ctg))?\.delta$/) {
                ($t, $q, $suffix_mode) = ($1, $2, $3);
            } else { next; }
        }

        my $mode = $opt->{mode} eq 'auto' ? ($suffix_mode || 'chr') : $opt->{mode};
        push @jobs, {
            file => $file,
            path => "$opt->{input}/$file",
            t    => $t,
            q    => $q,
            mode => $mode,
            ali  => $opt->{ali},
        };
    }
    closedir $DIR;

    @jobs = sort { $a->{file} cmp $b->{file} } @jobs;
    return @jobs;
}

sub run_jobs {
    my ($jobs, $nproc, $code) = @_;
    $nproc ||= 1;
    if ($nproc > 1 && @$jobs > 1) {
        my %children;
        my $failed = 0;
        for my $job (@$jobs) {
            while (keys(%children) >= $nproc) {
                my $pid = wait();
                last if $pid == -1;
                $failed++ if $? != 0;
                delete $children{$pid};
            }
            die "ERROR: a conversion worker failed.\n" if $failed;
            my $pid = fork();
            die "ERROR: cannot fork conversion worker: $!\n" unless defined $pid;
            if ($pid == 0) {
                eval { $code->($job); 1 } or do {
                    warn($@ || "unknown conversion error\n");
                    exit 1;
                };
                exit 0;
            }
            $children{$pid} = 1;
        }
        while (keys %children) {
            my $pid = wait();
            last if $pid == -1;
            $failed++ if $? != 0;
            delete $children{$pid};
        }
        die "ERROR: one or more conversion workers failed.\n" if $failed;
        return;
    }

    for my $job (@$jobs) {
        $code->($job);
    }
}

####################################################################################
# Conversion logic
####################################################################################
sub convert_one_alignment {
    my ($job, $opt, $sizes, $target_t_dir, $target_q_dir) = @_;
    my ($file, $t, $q, $mode, $ali) = @$job{qw/file t q mode ali/};
    my $prefix = "$t.$q";

    print STDERR "[$$] Converting $file ($ali, $mode mode) ...\n";

    if ($ali eq 'minimap2') {
        my $sorted_paf = "$opt->{wk}/$prefix.sorted.paf";
        my $raw_maf    = "$opt->{wk}/$prefix.sorted.maf";
        run_cmd("sort -k6,6 -k8,8n " . shq($job->{path}) . " > " . shq($sorted_paf));
        my $var_out = "$opt->{wk}/$prefix.sorted.var.txt";
        run_cmd(shq(exec_path($opt->{k8}, 'k8')) . " " . shq(exec_path($opt->{paftools}, 'paftools.js')) . " call " . shq($sorted_paf) . " > " . shq($var_out))
            if $opt->{emit_variants};
        run_cmd(shq(exec_path($opt->{k8}, 'k8')) . " " . shq(exec_path($opt->{paftools}, 'paftools.js')) . " view -f maf " . shq($sorted_paf) . " | sed 's/^a /a score=/g' > " . shq($raw_maf));

        my $maf_for_axt = prepare_maf_for_mode($raw_maf, $mode, $sizes);
        maf_to_axt_and_move($maf_for_axt, $t, $q, "$opt->{wk}/${t}vs${q}.axt", $target_t_dir, $opt);
        make_query_reference_axt($maf_for_axt, $t, $q, $target_q_dir, $opt) if $opt->{query};
        cleanup($opt, $sorted_paf, $raw_maf, ($opt->{emit_variants} ? $var_out : ()), ($maf_for_axt ne $raw_maf ? $maf_for_axt : ())) ;

    } elsif ($ali eq 'mummer' || $ali eq 'mummer4') {
        my $raw_maf = "$opt->{wk}/$prefix.sorted.maf";
        run_cmd(shq(exec_path($opt->{mum320}, 'delta2maf')) . " -r " . shq($job->{path}) . " > " . shq($raw_maf));
        my $maf_for_axt = prepare_maf_for_mode($raw_maf, $mode, $sizes);
        maf_to_axt_and_move($maf_for_axt, $t, $q, "$opt->{wk}/${t}vs${q}.axt", $target_t_dir, $opt);
        make_query_reference_axt($maf_for_axt, $t, $q, $target_q_dir, $opt) if $opt->{query};
        cleanup($opt, $raw_maf, ($maf_for_axt ne $raw_maf ? $maf_for_axt : ()));

    } elsif ($ali eq 'last') {
        my $header_maf = "$opt->{wk}/$prefix.header.maf";
        my $raw_maf    = "$opt->{wk}/$prefix.sorted.maf";
        ensure_maf_header($job->{path}, $header_maf);
        run_cmd(shq(exec_path($opt->{tba}, 'maf_sort')) . " " . shq($header_maf) . " " . shq($t) . " | sed 's/^a /a score=/g' > " . shq($raw_maf));
        my $maf_for_axt = prepare_maf_for_mode($raw_maf, $mode, $sizes);
        maf_to_axt_and_move($maf_for_axt, $t, $q, "$opt->{wk}/${t}vs${q}.axt", $target_t_dir, $opt);
        make_query_reference_axt($maf_for_axt, $t, $q, $target_q_dir, $opt) if $opt->{query};
        cleanup($opt, $header_maf, $raw_maf, ($maf_for_axt ne $raw_maf ? $maf_for_axt : ()));

    } elsif ($ali eq 'lastz' || $ali eq 'maf') {
        my $header_maf = "$opt->{wk}/$prefix.header.maf";
        my $raw_maf    = "$opt->{wk}/$prefix.sorted.maf";
        ensure_maf_header($job->{path}, $header_maf);
        run_cmd(shq(exec_path($opt->{tba}, 'maf_sort')) . " " . shq($header_maf) . " " . shq($t) . " > " . shq($raw_maf));
        my $maf_for_axt = prepare_maf_for_mode($raw_maf, $mode, $sizes);
        maf_to_axt_and_move($maf_for_axt, $t, $q, "$opt->{wk}/${t}vs${q}.axt", $target_t_dir, $opt);
        make_query_reference_axt($maf_for_axt, $t, $q, $target_q_dir, $opt) if $opt->{query};
        cleanup($opt, $header_maf, $raw_maf, ($maf_for_axt ne $raw_maf ? $maf_for_axt : ()));

    } elsif ($ali eq 'anchorwave') {
        my $trans_maf = "$opt->{wk}/$prefix.trans.maf";
        my $raw_maf   = "$opt->{wk}/$prefix.sorted.maf";

        # AnchorWave MAFs produced by genoAli/proali can have chromosome-only
        # sequence IDs (for example Chr01) and may not be accepted by TBA maf_sort.
        # Here we first add genome prefixes (ISO1.Chr01 / SIM001.Chr01) and then
        # sort MAF blocks internally by target sequence and target start, avoiding
        # the external maf_sort step for AnchorWave input.
        reformat_anchorwave_maf($job->{path}, $trans_maf, $t, $q);
        sort_maf_by_target($trans_maf, $raw_maf, $t);

        my $maf_for_axt = prepare_maf_for_mode($raw_maf, $mode, $sizes);
        maf_to_axt_and_move($maf_for_axt, $t, $q, "$opt->{wk}/${t}vs${q}.axt", $target_t_dir, $opt);
        make_query_reference_axt($maf_for_axt, $t, $q, $target_q_dir, $opt) if $opt->{query};
        cleanup($opt, $trans_maf, $raw_maf, ($maf_for_axt ne $raw_maf ? $maf_for_axt : ()));
    }

    print STDERR "[$$] Finished $file\n";
}

sub prepare_maf_for_mode {
    my ($maf, $mode, $sizes) = @_;
    return $maf if $mode eq 'chr';
    my $out = "$maf.reformat.maf";
    reformat_split_maf_to_chr($maf, $out, $sizes);
    return $out;
}

sub maf_to_axt_and_move {
    my ($maf, $t, $q, $axt, $target_dir, $opt) = @_;
    run_cmd(shq(exec_path($opt->{ucsc}, 'mafToAxt')) . " " . shq($maf) . " " . shq($t) . " " . shq($q) . " " . shq($axt));
    my $dest = "$target_dir/" . basename($axt);
    rename($axt, $dest) or die "ERROR: failed to move $axt to $dest: $!\n";
}

sub make_query_reference_axt {
    my ($maf, $t, $q, $target_q_dir, $opt) = @_;
    my $reordered = "$opt->{wk}/$q.$t.sorted.maf";
    my $axt       = "$opt->{wk}/${q}vs${t}.axt";
    run_cmd(shq(exec_path($opt->{tba}, 'maf_order')) . " " . shq($maf) . " " . shq($q) . " " . shq($t) . " > " . shq($reordered));
    run_cmd(shq(exec_path($opt->{ucsc}, 'mafToAxt')) . " " . shq($reordered) . " " . shq($q) . " " . shq($t) . " " . shq($axt));
    my $dest = "$target_q_dir/" . basename($axt);
    rename($axt, $dest) or die "ERROR: failed to move $axt to $dest: $!\n";
    cleanup($opt, $reordered);
}

####################################################################################
# MAF coordinate reformatting for split contigs
####################################################################################
sub collect_size_files {
    my ($opt) = @_;
    my @files;
    for my $entry (@{ $opt->{sizes} || [] }) {
        push @files, split /,/, $entry;
    }
    push @files, glob("$opt->{wk}/*.sizes");
    push @files, glob("$opt->{input}/*.sizes") if -d $opt->{input};
    push @files, glob("$opt->{split_out}/*.sizes") if $opt->{split_out} && -d $opt->{split_out};

    my (%seen, %sizes);
    for my $f (@files) {
        next unless defined $f && $f ne '' && -s $f;
        next if $seen{$f}++;
        open(my $IN, '<', $f) or die "ERROR: cannot read size file $f: $!\n";
        while (my $line = <$IN>) {
            chomp $line;
            next if $line =~ /^\s*$/;
            my @x = split /\t/, $line;
            next unless @x >= 2;
            $sizes{$x[0]} = $x[1];
        }
        close $IN;
    }
    return %sizes;
}

sub reformat_split_maf_to_chr {
    my ($in_maf, $out_maf, $sizes) = @_;
    die "ERROR: ctg mode requires a size file. Use --sizes Query.sizes or run --splitfa first.\n"
        unless %$sizes;

    open(my $IN,  '<', $in_maf)  or die "ERROR: cannot read $in_maf: $!\n";
    open(my $OUT, '>', $out_maf) or die "ERROR: cannot write $out_maf: $!\n";

    local $/ = "\n\n";
    my $converted = 0;
    while (my $block = <$IN>) {
        chomp $block;
        next if $block =~ /^\s*$/;
        my @lines = split /\n/, $block;
        my @headers = grep { /^\s*#/ } @lines;
        print $OUT join("\n", @headers), "\n\n" if @headers;
        @lines = grep { !/^\s*#/ } @lines;
        next unless @lines;
        for my $line (@lines) {
            if ($line =~ /^s\s+/) {
                my @u = split /\s+/, $line;
                # MAF s line: s src start size strand srcSize text
                if (@u >= 7) {
                    my ($orig, $sub_start, $sub_end) = parse_split_id($u[1]);
                    if (defined $orig) {
                        die "ERROR: no original chromosome size for '$orig'. Please provide the correct --sizes file.\n"
                            unless exists $sizes->{$orig};
                        my $orig_size = $sizes->{$orig};
                        my $new_start;
                        if ($u[4] eq '+') {
                            $new_start = $sub_start + $u[2] - 1;
                        } elsif ($u[4] eq '-') {
                            $new_start = $orig_size - $sub_end + $u[2];
                        } else {
                            die "ERROR: unknown MAF strand '$u[4]' in $in_maf\n";
                        }
                        $u[1] = $orig;
                        $u[2] = $new_start;
                        $u[5] = $orig_size;
                        $line = join("\t", @u);
                        $converted++;
                    }
                }
            }
        }
        print $OUT join("\n", @lines), "\n\n";
    }
    print $OUT "##eof maf\n";
    close $IN;
    close $OUT;

    die "ERROR: ctg mode was selected, but no split-contig coordinates were converted in $in_maf. Check subcontig names such as Chr01_1_100000 or Chr01.1.100000.\n"
        if $converted == 0;
}

sub parse_split_id {
    my ($src) = @_;
    # Compatible with both the original SplitFa.pl + sed convention:
    #   Chr01_1_100000
    #   Chr01.1.100000
    # Greedy parsing preserves original chromosome names containing dots/underscores.
    if ($src =~ /^(.+)[_\.](\d+)[_\.](\d+)$/) {
        my ($orig, $s, $e) = ($1, $2, $3);
        return ($orig, $s, $e) if $s >= 1 && $e >= $s;
    }
    return;
}

####################################################################################
# AnchorWave MAF reformatting
####################################################################################
sub reformat_anchorwave_maf {
    my ($in, $out, $t, $q) = @_;
    open(my $IN,  '<', $in)  or die "ERROR: cannot read $in: $!\n";
    open(my $OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";

    local $/ = "\n\n";
    while (my $block = <$IN>) {
        chomp $block;
        next if $block =~ /^\s*$/;
        my @lines = split /\n/, $block;
        my @headers = grep { /^\s*#/ } @lines;
        print $OUT join("\n", @headers), "\n\n" if @headers;
        @lines = grep { !/^\s*#/ } @lines;
        my @sidx = grep { $lines[$_] =~ /^s\s+/ } 0 .. $#lines;
        next unless @sidx >= 2;

        # Make sure each alignment block has a MAF score field. Some AnchorWave
        # outputs use a bare "a" line, which can break downstream tools.
        for my $i (0 .. $#lines) {
            if ($lines[$i] =~ /^a\s*$/) {
                $lines[$i] = 'a score=0';
            }
        }

        my @u1 = split /\s+/, $lines[$sidx[0]];
        my @u2 = split /\s+/, $lines[$sidx[1]];
        next unless @u1 >= 7 && @u2 >= 7;
        next if $u1[3] == 0 || $u2[3] == 0;

        # AnchorWave genoAli is often run after temporarily renaming both genomes
        # to the same chromosome IDs (Chr01, Chr02...).  Add genome prefixes back
        # so that mafToAxt can identify target and query species correctly.
        $u1[1] = "$t.$u1[1]" unless $u1[1] =~ /^\Q$t\E\./;
        $u2[1] = "$q.$u2[1]" unless $u2[1] =~ /^\Q$q\E\./;
        $lines[$sidx[0]] = join("\t", @u1);
        $lines[$sidx[1]] = join("\t", @u2);

        print $OUT join("\n", @lines), "\n\n";
    }
    close $IN;
    close $OUT;
}

sub sort_maf_by_target {
    my ($in, $out, $target) = @_;
    open(my $IN,  '<', $in)  or die "ERROR: cannot read $in: $!\n";

    my @headers;
    my @blocks;
    my $has_eof = 0;
    local $/ = "\n\n";
    while (my $block = <$IN>) {
        chomp $block;
        next if $block =~ /^\s*$/;
        my @lines = split /\n/, $block;
        for my $line (grep { /^\s*#/ } @lines) {
            if ($line =~ /^\s*##eof\b/i) {
                $has_eof = 1;
            } else {
                push @headers, $line;
            }
        }
        @lines = grep { !/^\s*#/ } @lines;
        next unless @lines;
        my ($src, $start) = ('', 0);
        for my $line (@lines) {
            next unless $line =~ /^s\s+/;
            my @x = split /\s+/, $line;
            next unless @x >= 4;
            if ($x[1] =~ /^\Q$target\E(?:\.|$)/) {
                ($src, $start) = ($x[1], $x[2]);
                last;
            }
        }

        # Keep only blocks that contain the target species. Blocks without target
        # cannot be converted by mafToAxt for this pair.
        next if $src eq '';
        push @blocks, { src => $src, start => $start + 0, block => $block };
    }
    close $IN;

    @blocks = sort {
        $a->{src} cmp $b->{src} || $a->{start} <=> $b->{start}
    } @blocks;

    open(my $OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";
    if (@headers) {
        my $has_maf_header = 0;
        for my $h (@headers) {
            $has_maf_header = 1 if $h =~ /^##maf/;
        }
        print $OUT "##maf version=1\n\n" unless $has_maf_header;
        print $OUT join("\n\n", @headers), "\n\n";
    } else {
        print $OUT "##maf version=1\n\n";
    }
    for my $b (@blocks) {
        print $OUT $b->{block}, "\n\n";
    }
    print $OUT "##eof maf\n" if $has_eof;
    close $OUT;
}

####################################################################################
# MAF header handling, shell utilities and cleanup
####################################################################################
sub ensure_maf_header {
    my ($in, $out) = @_;
    open(my $IN,  '<', $in)  or die "ERROR: cannot read $in: $!\n";
    open(my $OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";
    my $first = <$IN>;
    if (!defined $first) {
        close $IN; close $OUT;
        die "ERROR: empty MAF file: $in\n";
    }
    if ($first !~ /^##maf/) {
        print $OUT "##maf version=1\n";
    }
    print $OUT $first;
    while (my $line = <$IN>) {
        print $OUT $line;
    }
    close $IN;
    close $OUT;
}

sub exec_path {
    my ($path_or_dir, $exe) = @_;
    $path_or_dir ||= '';
    if ($path_or_dir eq '') {
        return $exe;
    }
    if (-d $path_or_dir) {
        $path_or_dir =~ s/\/$//;
        return "$path_or_dir/$exe";
    }
    if ($path_or_dir =~ /\/$/) {
        return "$path_or_dir$exe";
    }
    return $path_or_dir;
}

sub run_cmd {
    my ($cmd) = @_;
    print STDERR "CMD: $cmd\n";
    system($cmd);
    if ($? == -1) {
        die "ERROR: failed to execute command: $cmd\n$!\n";
    } elsif ($? & 127) {
        die sprintf("ERROR: command died with signal %d: %s\n", ($? & 127), $cmd);
    } elsif (($? >> 8) != 0) {
        die sprintf("ERROR: command exited with status %d: %s\n", ($? >> 8), $cmd);
    }
}

sub shq {
    my ($s) = @_;
    $s = '' unless defined $s;
    return "''" if $s eq '';
    $s =~ s/'/'"'"'/g;
    return "'$s'";
}

sub cleanup {
    my ($opt, @files) = @_;
    return if $opt->{keep_tmp};
    for my $f (@files) {
        next unless defined $f && $f ne '' && -e $f;
        unlink $f or warn "WARNING: failed to remove temporary file $f: $!\n";
    }
}
