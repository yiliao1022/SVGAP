#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
# ============================================================================
# SVGAP standardized script: 03_SynNetFilter.pl
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
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/pub";
use lib "$Bin/../lib";
use lib "$Bin/../pub";
use lib "$Bin/../../lib";
use lib "$Bin/../../pub";
use File::Path qw(make_path);
use File::Basename qw(basename);
use Cwd qw(abs_path getcwd);

# -----------------------------------------------------------------------------
# 3_SynNetFilter_standardized.pl
#
# Standardized step 3 for SVGAP-like workflows.
#
# Purpose:
#   1) Filter UCSC netSyntenic output (*.chain.filter.tnet/qnet.synnet) to remove
#      unreliable short fills inside gaps, especially short insertion-gap fragments.
#   2) Use the cleaned target-oriented synnet (tnet) to generate a reliable
#      single-coverage MAF file with netToAxt -> axtToMaf -> single_cov2.
#   3) Export syntenic gap BED files and homologous chromosome-pair summaries.
#
# Input expected from step 2:
#   - SYN_DIR/RefvsQuery.chain.filter.tnet.synnet
#   - SYN_DIR/RefvsQuery.chain.filter.qnet.synnet   (optional but processed)
#   - CHAIN_DIR/Ref.Query.chain.filter
#   - GENOME_DIR/Ref(.fa/.fasta/.fna or no suffix), Query(...)
#   - GENOME_DIR/Ref.sizes, Ref.2bit, Query.sizes, Query.2bit are generated if absent
#
# Authors: Yi Liao; standardized by ChatGPT
# -----------------------------------------------------------------------------

my %opt = (
    synnet      => undef,
    gd          => undef,
    chain       => undef,
    outdir      => 'SynNetFilter_out',

    ucsc        => "$Bin/pub/UCSC/",
    tba         => "$Bin/pub/TBA/",

    cleanSynnet => 'CleanSynNet',
    sing        => 'filtered_sing',
    axt         => 'axt',
    homopairs   => 'Hompairs',
    rawsv       => 'RawSV',
    tmp         => 'tmp',

    target      => undef,
    query       => undef,

    toplen      => 500000,
    nonsynlen   => 500000,
    synlen      => 500000,
    invlen      => 1000,
    alrate      => 0.10,
    qFarCutoff  => 5000000,
    adaptive_thresholds => 1,
    auto_fraction => 0.10,
    auto_min => 1000,
    minscore    => 1000,

    t           => 1,
    sing_only_tnet => 1,
    keep_temp   => 0,
    force       => 0,
    dry_run     => 0,
    help        => 0,
);

GetOptions(
    'synnet=s'      => \$opt{synnet},
    'gd=s'          => \$opt{gd},
    'chain=s'       => \$opt{chain},
    'outdir=s'      => \$opt{outdir},

    'ucsc=s'        => \$opt{ucsc},
    'tba=s'         => \$opt{tba},

    'cleanSynnet=s' => \$opt{cleanSynnet},
    'sing=s'        => \$opt{sing},
    'axt=s'         => \$opt{axt},
    'homopairs=s'   => \$opt{homopairs},
    'rawsv=s'       => \$opt{rawsv},
    'tmp=s'         => \$opt{tmp},

    'target|refname=s' => \$opt{target},
    'query=s'       => \$opt{query},

    'toplen=i'      => \$opt{toplen},
    'nonsynlen=i'   => \$opt{nonsynlen},
    'synlen=i'      => \$opt{synlen},
    'invlen=i'      => \$opt{invlen},
    'alrate=f'      => \$opt{alrate},
    'qFarCutoff=i'  => \$opt{qFarCutoff},
    'adaptive_thresholds|adaptive-thresholds!' => \$opt{adaptive_thresholds},
    'auto_fraction|auto-fraction=f' => \$opt{auto_fraction},
    'auto_min|auto-min=i' => \$opt{auto_min},
    'minscore=i'    => \$opt{minscore},

    't|threads=i'   => \$opt{t},
    'sing_only_tnet!' => \$opt{sing_only_tnet},
    'keep_temp!'    => \$opt{keep_temp},
    'force!'        => \$opt{force},
    'dry_run!'      => \$opt{dry_run},
    'help|h'        => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $required (qw/synnet gd chain/) {
    die "ERROR: --$required is required.\n" . usage() unless defined $opt{$required} && $opt{$required} ne '';
}

validate_options(\%opt);
normalize_paths(\%opt);

make_path($opt{outdir});
my $log_file = "$opt{outdir}/SynNetFilter.log";
open(my $LOG, '>', $log_file) or die "ERROR: cannot write log file $log_file: $!\n";
select((select($LOG), $| = 1)[0]);

for my $subdir (qw/cleanSynnet sing axt homopairs rawsv tmp/) {
    make_path($opt{outdir} . '/' . $opt{$subdir});
}

logmsg($LOG, "Input synnet directory: $opt{synnet}");
logmsg($LOG, "Input chain directory:  $opt{chain}");
logmsg($LOG, "Input genome directory: $opt{gd}");
logmsg($LOG, "Output directory:      $opt{outdir}");
logmsg($LOG, "Thresholds: toplen=$opt{toplen}; nonsynlen=$opt{nonsynlen}; synlen=$opt{synlen}; invlen=$opt{invlen}; alrate=$opt{alrate}; qFarCutoff=$opt{qFarCutoff}; adaptive=$opt{adaptive_thresholds}; auto_fraction=$opt{auto_fraction}; auto_min=$opt{auto_min}");

my %exe = find_required_programs(\%opt, $LOG);
my @jobs = discover_synnet_jobs(\%opt, $LOG);
die "ERROR: no *.chain.filter.tnet/qnet.synnet files found in $opt{synnet}\n" unless @jobs;

logmsg($LOG, "Discovered " . scalar(@jobs) . " synnet file(s) to process.");

# Prepare genome .sizes and .2bit files for all target/query IDs used by jobs.
my %genome_id;
for my $job (@jobs) {
    $genome_id{$job->{target}} = 1;
    $genome_id{$job->{query}}  = 1;
}
for my $id (sort keys %genome_id) {
    prepare_genome_files($id, \%opt, \%exe, $LOG);
}

run_jobs(\@jobs, \%opt, \%exe, $LOG);

logmsg($LOG, "Finished.");
close $LOG;

print STDERR "Done.\n";
print STDERR "Clean synnet: $opt{outdir}/$opt{cleanSynnet}/\n";
print STDERR "Filtered single-coverage MAF: $opt{outdir}/$opt{sing}/\n";
print STDERR "Raw syntenic gaps: $opt{outdir}/$opt{rawsv}/\n";
print STDERR "Log: $log_file\n";
exit 0;

# -----------------------------------------------------------------------------
# Main processing
# -----------------------------------------------------------------------------
sub run_jobs {
    my ($jobs, $opt, $exe, $LOG) = @_;

    if ($opt->{t} > 1 && @$jobs > 1) {
        my %children;
        my $failed = 0;
        for my $job (@$jobs) {
            while (keys(%children) >= $opt->{t}) {
                my $pid = wait();
                last if $pid == -1;
                $failed++ if $? != 0;
                delete $children{$pid};
            }
            die "ERROR: a SynNetFilter worker failed.\n" if $failed;
            my $pid = fork();
            die "ERROR: cannot fork SynNetFilter worker: $!\n" unless defined $pid;
            if ($pid == 0) {
                close $LOG;
                open(my $CHILD_LOG, '>>', "$opt->{outdir}/SynNetFilter.child.log") or die "ERROR: cannot write child log: $!\n";
                eval { process_one_job($job, $opt, $exe, $CHILD_LOG); 1 } or do {
                    warn($@ || "unknown SynNetFilter error\n");
                    close $CHILD_LOG;
                    exit 1;
                };
                close $CHILD_LOG;
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
        die "ERROR: one or more SynNetFilter workers failed.\n" if $failed;
    } else {
        for my $job (@$jobs) { process_one_job($job, $opt, $exe, $LOG); }
    }
}

sub process_one_job {
    my ($job, $opt, $exe, $LOG) = @_;

    my ($t, $q, $nettype, $file) = @{$job}{qw/target query nettype file/};
    my $in_synnet = "$opt->{synnet}/$file";
    my $prefix = "$t" . "vs" . "$q.chain.filter.$nettype.synnet";
    my $tmp_filtered   = "$opt->{outdir}/$opt->{tmp}/$prefix.refiltered.tmp.synnet";
    my $clean_filtered = "$opt->{outdir}/$opt->{cleanSynnet}/$prefix.refiltered.synnet";
    my $hom_file       = "$opt->{outdir}/$opt->{homopairs}/$prefix.pairedChr.tsv";
    my $raw_gap_bed    = "$opt->{outdir}/$opt->{rawsv}/$prefix.refiltered.rawsv.bed";

    logmsg($LOG, "Processing $file");

    my %pairs;
    filter_synnet($in_synnet, $tmp_filtered, \%pairs, $opt, $t, $q, $LOG);
    remove_empty_net_blocks($tmp_filtered, $clean_filtered);
    write_homologous_pairs($hom_file, \%pairs);
    write_gap_bed($clean_filtered, $raw_gap_bed);

    if (!$opt->{keep_temp}) {
        unlink $tmp_filtered if -e $tmp_filtered;
    }

    if ($nettype eq 'tnet' || !$opt->{sing_only_tnet}) {
        generate_single_cov_maf($t, $q, $nettype, $clean_filtered, $opt, $exe, $LOG);
    }
}

sub generate_single_cov_maf {
    my ($t, $q, $nettype, $clean_synnet, $opt, $exe, $LOG) = @_;

    # The reliable single-coverage MAF is usually produced from tnet, because it
    # reports single coverage relative to the target/reference genome.
    my $chain_file = find_chain_file($opt->{chain}, $t, $q);
    die "ERROR: cannot find chain.filter file for $t vs $q in $opt->{chain}\n" unless defined $chain_file;

    my $t_2bit  = "$opt->{gd}/$t.2bit";
    my $q_2bit  = "$opt->{gd}/$q.2bit";
    my $t_size  = "$opt->{gd}/$t.sizes";
    my $q_size  = "$opt->{gd}/$q.sizes";

    for my $f ($chain_file, $t_2bit, $q_2bit, $t_size, $q_size) {
        die "ERROR: required file not found or empty: $f\n" unless -s $f;
    }

    my $base = "$t" . "vs" . "$q.chain.filter.$nettype.synnet.refiltered";
    my $axt_out = "$opt->{outdir}/$opt->{axt}/$base.axt";
    my $maf_out = "$opt->{outdir}/$opt->{tmp}/$base.axt.maf";
    my $sing_out = "$opt->{outdir}/$opt->{sing}/$t.$q.sing.maf";

    run_cmd(join(' ', map { shell_quote($_) } ($exe->{netToAxt}, $clean_synnet, $chain_file, $t_2bit, $q_2bit, $axt_out)), $opt, $LOG);
    run_cmd(join(' ', map { shell_quote($_) } ($exe->{axtToMaf}, $axt_out, $t_size, $q_size, $maf_out)), $opt, $LOG);
    run_cmd(shell_quote($exe->{single_cov2}) . ' ' . shell_quote($maf_out) . ' > ' . shell_quote($sing_out), $opt, $LOG);

    if (!$opt->{keep_temp}) {
        unlink $maf_out if -e $maf_out;
    }
}

# -----------------------------------------------------------------------------
# Synnet filtering
# -----------------------------------------------------------------------------
sub filter_synnet {
    my ($input, $output, $pairs, $opt, $target_id, $query_id, $LOG) = @_;

    open(my $IN,  '<', $input)  or die "ERROR: cannot read $input: $!\n";
    open(my $OUT, '>', $output) or die "ERROR: cannot write $output: $!\n";

    my ($top_ok, $second_ok, $third_ok) = (0, 0, 0);
    my $target_chr = undef;
    my $target_size;

    while (my $line = <$IN>) {
        chomp $line;

        if ($line =~ /^\s*$/) {
            next;
        } elsif ($line =~ /^#/) {
            print $OUT "$line\n";
            next;
        } elsif ($line =~ /^net\s+(\S+)(?:\s+(\d+))?/) {
            my ($net_name, $net_size) = ($1, $2);
            print $OUT "$line\n";
            $target_chr = parse_seq_chr($net_name, $target_id);
            $target_size = defined($net_size) ? $net_size + 0 : undef;
            ($top_ok, $second_ok, $third_ok) = (0, 0, 0);
            next;
        }

        if ($line =~ /^\s{1}fill\b/) {
            my @f = split(/\s+/, trim($line));
            my $query_chr = parse_seq_chr(fill_qname(\@f) // '.', $query_id);
            my $len = fill_len(\@f);
            my $alrate = alignment_rate(\@f);
            my $threshold = effective_threshold($opt->{toplen}, $target_size, $opt);
            if ($len > $threshold && $alrate > $opt->{alrate}) {
                print $OUT "$line\n";
                record_pair($pairs, $target_chr, $query_chr);
                $top_ok = 1;
            } else {
                $top_ok = 0;
            }
            ($second_ok, $third_ok) = (0, 0);

        } elsif ($line =~ /^\s{2}gap\b/) {
            print $OUT "$line\n" if $top_ok;

        } elsif ($line =~ /^\s{3}fill\b/ && $line =~ /nonSyn/) {
            my @f = split(/\s+/, trim($line));
            my $query_chr = parse_seq_chr(fill_qname(\@f) // '.', $query_id);
            my $len = fill_len(\@f);
            my $alrate = alignment_rate(\@f);
            my $threshold = effective_threshold($opt->{nonsynlen}, $target_size, $opt);
            if ($top_ok && $len > $threshold && $alrate > $opt->{alrate}) {
                print $OUT "$line\n";
                record_pair($pairs, $target_chr, $query_chr);
                $second_ok = 1;
            } else {
                $second_ok = 0;
            }
            $third_ok = 0;

        } elsif ($line =~ /^\s{3}fill\b/ && $line =~ /(syn|inv)/) {
            my @f = split(/\s+/, trim($line));
            if ($top_ok && lower_fill_passes(\@f, $line, $opt, $target_size)) {
                print $OUT "$line\n";
                $second_ok = 1;
            } else {
                $second_ok = 0;
            }
            $third_ok = 0;

        } elsif ($line =~ /^\s{4}gap\b/) {
            print $OUT "$line\n" if $top_ok && $second_ok;

        } elsif ($line =~ /^\s{5}fill\b/ && $line =~ /nonSyn/) {
            my @f = split(/\s+/, trim($line));
            my $query_chr = parse_seq_chr(fill_qname(\@f) // '.', $query_id);
            my $len = fill_len(\@f);
            my $alrate = alignment_rate(\@f);
            my $threshold = effective_threshold($opt->{nonsynlen}, $target_size, $opt);
            if ($top_ok && $second_ok && $len > $threshold && $alrate > $opt->{alrate}) {
                print $OUT "$line\n";
                record_pair($pairs, $target_chr, $query_chr);
                $third_ok = 1;
            } else {
                $third_ok = 0;
            }

        } elsif ($line =~ /^\s{5}fill\b/ && $line =~ /(syn|inv)/) {
            my @f = split(/\s+/, trim($line));
            if ($top_ok && $second_ok && lower_fill_passes(\@f, $line, $opt, $target_size)) {
                print $OUT "$line\n";
                $third_ok = 1;
            } else {
                $third_ok = 0;
            }

        } elsif ($line =~ /^\s{6}gap\b/) {
            print $OUT "$line\n" if $top_ok && $second_ok && $third_ok;
        }
    }

    close $IN;
    close $OUT;
}

sub lower_fill_passes {
    my ($f, $line, $opt, $target_size) = @_;
    my $len = fill_len($f);
    my $alrate = alignment_rate($f);

    if ($line =~ /syn/) {
        my $threshold = effective_threshold($opt->{synlen}, $target_size, $opt);
        return ($len > $threshold && $alrate > $opt->{alrate}) ? 1 : 0;
    }

    if ($line =~ /inv/) {
        my $qfar  = value_after_key($f, 'qFar');
        my $qover = value_after_key($f, 'qOver');
        my $qsize = fill_qsize($f);

        return 1 if ($len > $opt->{invlen}
                  && $alrate > $opt->{alrate}
                  && defined($qover) && $qover > 0
                  && defined($qfar)  && $qfar < $opt->{qFarCutoff});

        return 1 if ($len > 5 * $opt->{invlen}
                  && $qsize > 5 * $opt->{invlen}
                  && defined($qfar) && $qfar < $len);
    }

    return 0;
}

sub effective_threshold {
    my ($configured, $sequence_size, $opt) = @_;
    return $configured unless $opt->{adaptive_thresholds};
    return $configured unless defined $sequence_size && $sequence_size > 0;
    my $scaled = int($sequence_size * $opt->{auto_fraction});
    $scaled = $opt->{auto_min} if $scaled < $opt->{auto_min};
    return $scaled < $configured ? $scaled : $configured;
}

sub alignment_rate {
    my ($f) = @_;
    my $len = fill_len($f);
    return 0 if !$len || $len !~ /^\d+(?:\.\d+)?$/;
    my $ali = value_after_key($f, 'ali');
    return 0 unless defined $ali && $ali =~ /^\d+(?:\.\d+)?$/;
    return $ali / $len;
}

# UCSC net fill lines look like:
#   fill tStart size qName qStart strand qSize chainId score ... ali N ...
# After trim()+split(), indexes are:
#   0=fill, 1=tStart, 2=size, 3=qName, 4=qStart, 5=strand, 6=qSize
# The original non-standard script used split(/\s/, line_with_leading_space),
# so its old index 3 happened to be size. In this standardized script we trim
# the line first, therefore size must be read from index 2.
sub fill_len {
    my ($f) = @_;
    my $v = $f->[2] // 0;
    return ($v =~ /^\d+(?:\.\d+)?$/) ? $v : 0;
}

sub fill_qname {
    my ($f) = @_;
    return $f->[3];
}

sub fill_qsize {
    my ($f) = @_;
    my $v = $f->[6] // 0;
    return ($v =~ /^\d+(?:\.\d+)?$/) ? $v : 0;
}

sub value_after_key {
    my ($f, $key) = @_;
    for (my $i = 0; $i < @$f - 1; $i++) {
        return $f->[$i + 1] if defined($f->[$i]) && $f->[$i] eq $key;
    }
    return undef;
}

sub remove_empty_net_blocks {
    my ($in, $out) = @_;
    open(my $IN,  '<', $in)  or die "ERROR: cannot read $in: $!\n";
    open(my $OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";

    my @block;
    my $has_fill = 0;

    while (my $line = <$IN>) {
        if ($line =~ /^#/) {
            flush_net_block($OUT, \@block, $has_fill);
            @block = ();
            $has_fill = 0;
            print $OUT $line;
            next;
        }

        if ($line =~ /^net\b/) {
            flush_net_block($OUT, \@block, $has_fill);
            @block = ($line);
            $has_fill = 0;
            next;
        }

        if (@block) {
            push @block, $line;
            $has_fill = 1 if $line =~ /^\s+fill\b/;
        }
    }
    flush_net_block($OUT, \@block, $has_fill);

    close $IN;
    close $OUT;
}

sub flush_net_block {
    my ($OUT, $block, $has_fill) = @_;
    return unless @$block;
    print $OUT @$block if $has_fill;
}

sub write_gap_bed {
    my ($clean_synnet, $out_bed) = @_;
    open(my $IN,  '<', $clean_synnet) or die "ERROR: cannot read $clean_synnet: $!\n";
    open(my $OUT, '>', $out_bed)     or die "ERROR: cannot write $out_bed: $!\n";

    my $chr;
    while (my $line = <$IN>) {
        chomp $line;
        $line =~ s/^\s+//;
        my @f = split(/\s+/, $line);
        if ($line =~ /^net\b/) {
            $chr = $f[1];
        } elsif ($line =~ /^gap\b/) {
            next unless defined $chr;
            my $start = $f[1];
            my $end = $f[1] + $f[2];
            print $OUT join("\t", $chr, $start, $end), "\n";
        }
    }

    close $IN;
    close $OUT;
}

sub write_homologous_pairs {
    my ($out, $pairs) = @_;
    open(my $OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";
    print $OUT join("\t", qw/TARGET_CHR QUERY_CHRS/), "\n";
    for my $chr (sort keys %$pairs) {
        my @q = sort keys %{ $pairs->{$chr} || {} };
        print $OUT join("\t", $chr, @q ? join(',', @q) : '.'), "\n";
    }
    close $OUT;
}

sub record_pair {
    my ($pairs, $target_chr, $query_chr) = @_;
    return unless defined $target_chr && defined $query_chr;
    $pairs->{$target_chr}{$query_chr} = 1;
}

sub parse_seq_chr {
    my ($name, $genome_id) = @_;
    return $name unless defined $name;
    if (defined $genome_id && $name =~ /^\Q$genome_id\E\.(.+)$/) {
        return $1;
    }
    # Historical SVGAP naming often used Genome.Chr; keep compatibility.
    my @p = split(/\./, $name);
    return $p[1] if @p == 2;
    return join('.', @p[1 .. $#p]) if @p > 2 && defined $genome_id && $p[0] eq $genome_id;
    return $name;
}

# -----------------------------------------------------------------------------
# Input discovery and genome preparation
# -----------------------------------------------------------------------------
sub discover_synnet_jobs {
    my ($opt, $LOG) = @_;
    opendir(my $DIR, $opt->{synnet}) or die "ERROR: cannot open synnet directory $opt->{synnet}: $!\n";

    my @jobs;
    while (my $file = readdir($DIR)) {
        next if $file =~ /^\./;
        next unless $file =~ /^(.+)vs(.+)\.chain\.filter\.(tnet|qnet)\.synnet$/;
        my ($t, $q, $nettype) = ($1, $2, $3);
        next if defined($opt->{target}) && $opt->{target} ne $t;
        next if defined($opt->{query})  && $opt->{query}  ne $q;
        push @jobs, { target => $t, query => $q, nettype => $nettype, file => $file };
    }
    closedir $DIR;

    @jobs = sort {
        $a->{target} cmp $b->{target}
        || $a->{query} cmp $b->{query}
        || $a->{nettype} cmp $b->{nettype}
    } @jobs;

    return @jobs;
}

sub prepare_genome_files {
    my ($id, $opt, $exe, $LOG) = @_;
    my $sizes = "$opt->{gd}/$id.sizes";
    my $twobit = "$opt->{gd}/$id.2bit";

    return if !$opt->{force} && -s $sizes && -s $twobit;

    my $fasta = find_genome_fasta($opt->{gd}, $id);
    die "ERROR: cannot find genome FASTA for ID '$id' in $opt->{gd}. Tried: $id, $id.fa, $id.fasta, $id.fna. A FASTA is needed because .sizes or .2bit is missing.\n"
        unless defined $fasta;

    if (!-s $sizes || $opt->{force}) {
        logmsg($LOG, "Generating sizes for $id");
        run_cmd(shell_quote($exe->{faSize}) . ' -detailed ' . shell_quote($fasta) . ' > ' . shell_quote($sizes), $opt, $LOG);
    }
    if (!-s $twobit || $opt->{force}) {
        logmsg($LOG, "Generating 2bit for $id");
        run_cmd(join(' ', map { shell_quote($_) } ($exe->{faToTwoBit}, $fasta, $twobit)), $opt, $LOG);
    }
}

sub find_genome_fasta {
    my ($gd, $id) = @_;
    for my $cand ("$gd/$id", "$gd/$id.fa", "$gd/$id.fasta", "$gd/$id.fna") {
        return $cand if -s $cand;
    }
    return undef;
}

sub find_chain_file {
    my ($chain_dir, $t, $q) = @_;
    for my $cand (
        "$chain_dir/$t.$q.chain.filter",
        "$chain_dir/$t" . "vs" . "$q.chain.filter",
        "$chain_dir/$t.$q.chain",
        "$chain_dir/$t" . "vs" . "$q.chain"
    ) {
        return $cand if -s $cand;
    }
    return undef;
}

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
sub validate_options {
    my ($opt) = @_;
    die "ERROR: --toplen must be >= 0\n"    if $opt->{toplen} < 0;
    die "ERROR: --nonsynlen must be >= 0\n" if $opt->{nonsynlen} < 0;
    die "ERROR: --synlen must be >= 0\n"     if $opt->{synlen} < 0;
    die "ERROR: --invlen must be >= 0\n"     if $opt->{invlen} < 0;
    die "ERROR: --alrate must be between 0 and 1\n" if $opt->{alrate} < 0 || $opt->{alrate} > 1;
    die "ERROR: --qFarCutoff must be >= 0\n" if $opt->{qFarCutoff} < 0;
    die "ERROR: --auto_fraction must be between 0 and 1\n"
        if $opt->{auto_fraction} < 0 || $opt->{auto_fraction} > 1;
    die "ERROR: --auto_min must be >= 0\n" if $opt->{auto_min} < 0;
    die "ERROR: --t must be >= 1\n" if $opt->{t} < 1;
    die "ERROR: --synnet is not a directory: $opt->{synnet}\n" unless -d $opt->{synnet};
    die "ERROR: --gd is not a directory: $opt->{gd}\n" unless -d $opt->{gd};
    die "ERROR: --chain is not a directory: $opt->{chain}\n" unless -d $opt->{chain};
}

sub normalize_paths {
    my ($opt) = @_;
    for my $k (qw/synnet gd chain/) {
        $opt->{$k} = abs_path($opt->{$k}) if defined $opt->{$k};
    }
    $opt->{outdir} = absolute_path_for_output($opt->{outdir});
    for my $k (qw/ucsc tba/) {
        $opt->{$k} =~ s{/+$}{};
        $opt->{$k} .= '/';
    }
}

# Avoid requiring File::Spec for old Perl installs.
sub absolute_path_for_output {
    my ($path) = @_;
    return $path if $path =~ m{^/};
    return getcwd() . '/' . $path;
}

sub find_required_programs {
    my ($opt, $LOG) = @_;
    my %exe;
    for my $name (qw/netToAxt axtToMaf faSize faToTwoBit/) {
        $exe{$name} = find_executable($name, $opt->{ucsc});
        die "ERROR: cannot find UCSC utility '$name'. Use --ucsc /path/to/UCSC/\n" unless defined $exe{$name};
    }
    $exe{single_cov2} = find_executable('single_cov2', $opt->{tba});
    die "ERROR: cannot find TBA program 'single_cov2'. Use --tba /path/to/TBA/\n" unless defined $exe{single_cov2};

    for my $k (sort keys %exe) {
        logmsg($LOG, "Program $k: $exe{$k}");
    }
    return %exe;
}

sub find_executable {
    my ($prog, $dir) = @_;
    my @candidates = (
        defined($dir) ? "$dir/$prog" : (),
        defined($dir) ? "$dir$prog"  : (),
        $prog,
    );
    for my $c (@candidates) {
        next unless defined $c;
        return $c if -x $c;
        my $path = `command -v $c 2>/dev/null`;
        chomp $path;
        return $path if $path && -x $path;
    }
    return undef;
}

sub run_cmd {
    my ($cmd, $opt, $LOG) = @_;
    logmsg($LOG, "CMD: $cmd");
    return if $opt->{dry_run};
    system($cmd);
    if ($? != 0) {
        my $exit = $? >> 8;
        die "ERROR: command failed with exit code $exit:\n$cmd\n";
    }
}

sub shell_quote {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/'/'\\''/g;
    return "'$s'";
}

sub trim {
    my ($s) = @_;
    $s =~ s/^\s+//;
    $s =~ s/\s+$//;
    return $s;
}

sub logmsg {
    my ($fh, $msg) = @_;
    my $t = scalar(localtime());
    print $fh "[$t] $msg\n";
}

sub usage {
    return <<'USAGE';

Usage:
  perl 3_SynNetFilter_standardized.pl \
      --synnet ChainNetSyn_out/synnet \
      --chain  ChainNetSyn_out/chain \
      --gd     ./genomes \
      --outdir SynNetFilter_out \
      --t 8

Required:
  --synnet       Directory containing *.chain.filter.tnet.synnet and/or qnet.synnet
  --chain        Directory containing *.chain.filter files from axtChain/chainPreNet
  --gd           Directory containing genome FASTA files and/or .sizes/.2bit files

Common optional parameters:
  --outdir       Output directory                           [SynNetFilter_out]
  --t            Number of parallel jobs                    [1]
  --ucsc         Directory of UCSC utilities                [script_dir/pub/UCSC/]
  --tba          Directory of TBA utilities                 [script_dir/pub/TBA/]

Filtering thresholds:
  --toplen       Minimum length for top syntenic fills       [500000]
  --nonsynlen    Minimum length for lower nonSyn fills       [500000]
  --synlen       Minimum length for lower syn fills          [500000]
  --invlen       Minimum length for inversion fills          [1000]
  --alrate       Minimum ali/fill_length proportion          [0.10]
  --qFarCutoff   Maximum qFar for accepting inv fills        [5000000]
  --adaptive_thresholds / --no-adaptive_thresholds
                 Scale the three fill thresholds down for short target contigs
                                                               [enabled]
  --auto_fraction
                 Per-contig adaptive threshold fraction       [0.10]
  --auto_min     Lower bound for adaptive thresholds          [1000]

Output subdirectories:
  --cleanSynnet  Cleaned/refiltered synnet directory         [CleanSynNet]
  --sing         Filtered single-coverage MAF directory      [filtered_sing]
  --axt          AXT generated from cleaned tnet directory   [axt]
  --homopairs    Homologous chromosome-pair summaries        [Hompairs]
  --rawsv        Raw syntenic-gap BED files                  [RawSV]

Pair selection:
  --target Ref   Only process synnet files with this target ID
  --query  Que   Only process synnet files with this query ID

Other:
  --keep_temp    Keep intermediate MAF files
  --dry_run      Print commands to log but do not execute them
  --help         Print this help message

Input file naming:
  RefvsQue.chain.filter.tnet.synnet
  RefvsQue.chain.filter.qnet.synnet
  Ref.Que.chain.filter

Main outputs:
  OUTDIR/CleanSynNet/RefvsQue.chain.filter.tnet.synnet.refiltered.synnet
  OUTDIR/filtered_sing/Ref.Que.sing.maf
  OUTDIR/RawSV/RefvsQue.chain.filter.tnet.synnet.refiltered.rawsv.bed
  OUTDIR/Hompairs/RefvsQue.chain.filter.tnet.synnet.pairedChr.tsv

Notes:
  - Single-coverage MAF is generated from cleaned tnet by default because it is
    single coverage relative to the target/reference genome.
  - The filtering removes short/unreliable fills nested inside syntenic gaps,
    which makes the final single_cov2 MAF more conservative for downstream SV calling.

USAGE
}
