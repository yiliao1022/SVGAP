#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
# ============================================================================
# SVGAP standardized script: 02_ChainNetSyn.pl
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
use Getopt::Long qw(GetOptions Configure);
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/pub";
use lib "$Bin/../lib";
use lib "$Bin/../pub";
use lib "$Bin/../../lib";
use lib "$Bin/../../pub";
use File::Path qw(make_path);
use Cwd qw(abs_path);

# -----------------------------------------------------------------------------
# 2_ChainNetSyn_standardized.pl
#
# Build UCSC chain/net/syntenic-net files from pairwise AXT alignments.
# Optionally generate target-oriented single-coverage MAF files used by SVGAP.
#
# Input AXT naming convention:
#   TargetIDvsQueryID.axt
#
# Genome naming convention:
#   The genome list contains genome IDs, one per line. For each ID, this script
#   looks for the FASTA in --gd as one of:
#       ID, ID.fa, ID.fasta, ID.fna
#   It generates/uses:
#       --gd/ID.sizes and --gd/ID.2bit
#
# Main UCSC workflow:
#   axtChain -> chainPreNet -> chainNet -> netSyntenic
#
# Optional single-coverage workflow:
#   netToAxt -> axtToMaf -> single_cov2
# -----------------------------------------------------------------------------

Configure('no_ignore_case');

my %opt = (
    gd              => undef,
    ad              => undef,
    lst             => undef,
    outdir          => 'ChainNetSyn_out',
    ucsc            => "$Bin/pub/UCSC/",
    tba             => "$Bin/pub/TBA/",
    linearGap       => "$Bin/pub/TBA/medium",
    minscore        => 1000,
    minSpace        => 1,
    t               => 1,
    sing            => 0,
    keep_raw_chain  => 0,
    keep_temp       => 0,
    force           => 0,
    dry_run         => 0,
    help            => 0,
);

GetOptions(
    'gd=s'                => \$opt{gd},
    'ad=s'                => \$opt{ad},
    'lst=s'               => \$opt{lst},
    'outdir=s'            => \$opt{outdir},
    'syn=s'               => \$opt{outdir},       # backward-compatible alias
    'ucsc=s'              => \$opt{ucsc},
    'tba=s'               => \$opt{tba},
    'linearGap=s'         => \$opt{linearGap},
    'lg=s'                => \$opt{linearGap},    # backward-compatible alias
    'minscore=i'          => \$opt{minscore},
    'minScore=i'          => \$opt{minscore},
    'minSpace=i'          => \$opt{minSpace},
    't|threads=i'         => \$opt{t},
    'sing!'               => \$opt{sing},
    'keep_raw_chain!'     => \$opt{keep_raw_chain},
    'keep_temp!'          => \$opt{keep_temp},
    'force!'              => \$opt{force},
    'dry_run!'            => \$opt{dry_run},
    'help|h'              => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

validate_options(\%opt);

# Normalize directories.
$opt{gd}     = abs_path($opt{gd})     || die "ERROR: cannot resolve --gd $opt{gd}\n";
$opt{ad}     = abs_path($opt{ad})     || die "ERROR: cannot resolve --ad $opt{ad}\n";
$opt{lst}    = abs_path($opt{lst})    || die "ERROR: cannot resolve --lst $opt{lst}\n";
$opt{outdir} = make_abs_dir($opt{outdir});

my $chain_dir  = "$opt{outdir}/chain";
my $net_dir    = "$opt{outdir}/net";
my $syn_dir    = "$opt{outdir}/synnet";
my $sing_dir   = "$opt{outdir}/singRaw";
my $tmp_dir    = "$opt{outdir}/tmp";
make_path($chain_dir, $net_dir, $syn_dir, $tmp_dir);
make_path($sing_dir) if $opt{sing};

my $log_file = "$opt{outdir}/ChainNetSyn.log";
open(my $LOG, '>', $log_file) or die "ERROR: cannot write $log_file: $!\n";
select((select($LOG), $| = 1)[0]);
logmsg($LOG, "Started 2_ChainNetSyn_standardized.pl");
logmsg($LOG, "Genome directory: $opt{gd}");
logmsg($LOG, "AXT directory: $opt{ad}");
logmsg($LOG, "Genome list: $opt{lst}");
logmsg($LOG, "Output directory: $opt{outdir}");
logmsg($LOG, "minScore=$opt{minscore}; minSpace=$opt{minSpace}; linearGap=$opt{linearGap}; threads=$opt{t}; sing=$opt{sing}");

# Resolve external tools.
my %tool;
for my $name (qw/faSize faToTwoBit axtChain chainPreNet chainNet netSyntenic/) {
    $tool{$name} = find_tool($opt{ucsc}, $name);
}
if ($opt{sing}) {
    for my $name (qw/netToAxt axtToMaf/) {
        $tool{$name} = find_tool($opt{ucsc}, $name);
    }
    $tool{single_cov2} = find_tool($opt{tba}, 'single_cov2');
}
check_linear_gap($opt{linearGap});

# Read genomes and prepare sizes/2bit.
my @genome_ids = read_list($opt{lst});
my %genome = prepare_genome_resources(\@genome_ids, \%opt, \%tool, $LOG);

# Discover AXT files.
my @jobs = discover_axt_jobs($opt{ad}, \%genome);
if (!@jobs) {
    die "ERROR: no AXT files found in $opt{ad}. Expected names like TargetIDvsQueryID.axt\n";
}
logmsg($LOG, "Discovered " . scalar(@jobs) . " AXT job(s)");

# Run independent genome-pair jobs with core Perl fork/wait.
if ($opt{t} > 1 && @jobs > 1) {
    my %children;
    my $failed = 0;
    for my $job (@jobs) {
        while (keys(%children) >= $opt{t}) {
            my $pid = wait();
            last if $pid == -1;
            $failed++ if $? != 0;
            delete $children{$pid};
        }
        die "ERROR: a ChainNetSyn worker failed.\n" if $failed;
        my $pid = fork();
        die "ERROR: cannot fork ChainNetSyn worker: $!\n" unless defined $pid;
        if ($pid == 0) {
            close $LOG;
            my $child_log_file = "$opt{outdir}/ChainNetSyn.child.log";
            open(my $CHILD_LOG, '>>', $child_log_file) or die "ERROR: cannot append $child_log_file: $!\n";
            eval {
                process_axt_job($job, \%opt, \%tool, \%genome, $chain_dir, $net_dir, $syn_dir, $sing_dir, $tmp_dir, $CHILD_LOG);
                1;
            } or do {
                warn($@ || "unknown ChainNetSyn error\n");
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
    die "ERROR: one or more ChainNetSyn workers failed.\n" if $failed;
} else {
    for my $job (@jobs) {
        process_axt_job($job, \%opt, \%tool, \%genome, $chain_dir, $net_dir, $syn_dir, $sing_dir, $tmp_dir, $LOG);
    }
}

logmsg($LOG, "Finished successfully");
logmsg($LOG, "Chain files: $chain_dir");
logmsg($LOG, "Net files: $net_dir");
logmsg($LOG, "Syntenic net files: $syn_dir");
logmsg($LOG, "Single-coverage MAF files: $sing_dir") if $opt{sing};
close $LOG;

print STDERR "Done.\n";
print STDERR "Output directory: $opt{outdir}\n";
print STDERR "Chain files: $chain_dir\n";
print STDERR "Net files: $net_dir\n";
print STDERR "Syntenic net files: $syn_dir\n";
print STDERR "Single-coverage MAF files: $sing_dir\n" if $opt{sing};
print STDERR "Log: $log_file\n";

# -----------------------------------------------------------------------------
# Main job
# -----------------------------------------------------------------------------
sub process_axt_job {
    my ($job, $opt, $tool, $genome, $chain_dir, $net_dir, $syn_dir, $sing_dir, $tmp_dir, $LOG) = @_;

    my ($t, $q, $axt) = @$job{qw/t q axt/};
    my $prefix = "$t.$q";

    die "ERROR: genome ID '$t' from AXT is not listed/found in --lst\n" unless exists $genome->{$t};
    die "ERROR: genome ID '$q' from AXT is not listed/found in --lst\n" unless exists $genome->{$q};

    my $raw_chain    = "$chain_dir/$prefix.chain";
    my $chain_filter = "$chain_dir/$prefix.chain.filter";
    my $tnet         = "$net_dir/$prefix.chain.filter.tnet";
    my $qnet         = "$net_dir/$prefix.chain.filter.qnet";
    my $tnet_syn     = "$syn_dir/${t}vs${q}.chain.filter.tnet.synnet";
    my $qnet_syn     = "$syn_dir/${t}vs${q}.chain.filter.qnet.synnet";

    if (!$opt->{force} && -s $tnet_syn && -s $qnet_syn && (!$opt->{sing} || -s "$sing_dir/$t.$q.sing.maf")) {
        logmsg($LOG, "SKIP $t vs $q: outputs already exist; use --force to overwrite");
        return;
    }

    logmsg($LOG, "Processing $t vs $q");

    run_cmd($LOG, $opt->{dry_run},
        $tool->{axtChain},
        "-linearGap=$opt->{linearGap}",
        "-minScore=$opt->{minscore}",
        $axt,
        $genome->{$t}{twobit},
        $genome->{$q}{twobit},
        $raw_chain,
    );

    run_cmd($LOG, $opt->{dry_run},
        $tool->{chainPreNet},
        $raw_chain,
        $genome->{$t}{sizes},
        $genome->{$q}{sizes},
        $chain_filter,
    );

    unlink $raw_chain if !$opt->{keep_raw_chain} && !$opt->{dry_run};

    run_cmd($LOG, $opt->{dry_run},
        $tool->{chainNet},
        "-minSpace=$opt->{minSpace}",
        "-minScore=$opt->{minscore}",
        $chain_filter,
        $genome->{$t}{sizes},
        $genome->{$q}{sizes},
        $tnet,
        $qnet,
    );

    run_cmd($LOG, $opt->{dry_run}, $tool->{netSyntenic}, $tnet, $tnet_syn);
    run_cmd($LOG, $opt->{dry_run}, $tool->{netSyntenic}, $qnet, $qnet_syn);

    if ($opt->{sing}) {
        my $sing_axt = "$tmp_dir/${t}vs${q}.chain.filter.tnet.synnet.axt";
        my $sing_maf = "$tmp_dir/${t}vs${q}.chain.filter.tnet.synnet.axt.maf";
        my $sing_out = "$sing_dir/$t.$q.sing.maf";

        run_cmd($LOG, $opt->{dry_run},
            $tool->{netToAxt},
            $tnet_syn,
            $chain_filter,
            $genome->{$t}{twobit},
            $genome->{$q}{twobit},
            $sing_axt,
        );

        run_cmd($LOG, $opt->{dry_run},
            $tool->{axtToMaf},
            $sing_axt,
            $genome->{$t}{sizes},
            $genome->{$q}{sizes},
            $sing_maf,
        );

        run_cmd_stdout($LOG, $opt->{dry_run}, $sing_out, $tool->{single_cov2}, $sing_maf);

        if (!$opt->{keep_temp} && !$opt->{dry_run}) {
            unlink $sing_axt if -e $sing_axt;
            unlink $sing_maf if -e $sing_maf;
        }
    }

    logmsg($LOG, "Finished $t vs $q");
}

# -----------------------------------------------------------------------------
# Genome preparation and AXT discovery
# -----------------------------------------------------------------------------
sub prepare_genome_resources {
    my ($ids, $opt, $tool, $LOG) = @_;
    my %genome;

    for my $id (@$ids) {
        next if $id =~ /\.sizes$/ || $id =~ /\.2bit$/;
        my $fa = find_genome_fasta($opt->{gd}, $id);
        my $sizes = "$opt->{gd}/$id.sizes";
        my $twobit = "$opt->{gd}/$id.2bit";

        if (!-s $sizes || $opt->{force}) {
            logmsg($LOG, "Generating size file for $id: $sizes");
            run_cmd_stdout($LOG, $opt->{dry_run}, $sizes, $tool->{faSize}, '-detailed', $fa);
        }
        if (!-s $twobit || $opt->{force}) {
            logmsg($LOG, "Generating 2bit file for $id: $twobit");
            run_cmd($LOG, $opt->{dry_run}, $tool->{faToTwoBit}, $fa, $twobit);
        }

        $genome{$id} = {
            fasta  => $fa,
            sizes  => $sizes,
            twobit => $twobit,
        };
    }

    return %genome;
}

sub discover_axt_jobs {
    my ($ad, $genome) = @_;
    opendir(my $DIR, $ad) or die "ERROR: cannot open AXT directory $ad: $!\n";
    my @jobs;
    while (my $file = readdir($DIR)) {
        next unless $file =~ /^(.+)vs(.+)\.axt$/;
        my ($t, $q) = ($1, $2);
        my $path = "$ad/$file";
        next unless -s $path;
        push @jobs, { t => $t, q => $q, axt => $path, file => $file };
    }
    closedir $DIR;

    @jobs = sort { $a->{t} cmp $b->{t} || $a->{q} cmp $b->{q} } @jobs;
    return @jobs;
}

sub find_genome_fasta {
    my ($dir, $id) = @_;
    my @candidates = (
        "$dir/$id",
        "$dir/$id.fa",
        "$dir/$id.fasta",
        "$dir/$id.fna",
    );
    for my $f (@candidates) {
        return $f if -s $f;
    }
    die "ERROR: cannot find genome FASTA for '$id' in $dir. Tried: " . join(', ', @candidates) . "\n";
}

sub read_list {
    my ($file) = @_;
    open(my $IN, '<', $file) or die "ERROR: cannot read $file: $!\n";
    my @ids;
    while (my $line = <$IN>) {
        chomp $line;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '' || $line =~ /^#/;
        my ($id) = split(/\s+/, $line);
        push @ids, $id;
    }
    close $IN;
    die "ERROR: no genome IDs found in $file\n" unless @ids;
    return @ids;
}

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
sub validate_options {
    my ($opt) = @_;
    die usage() if !$opt->{gd} || !$opt->{ad} || !$opt->{lst};
    die "ERROR: --gd does not exist: $opt->{gd}\n" unless -d $opt->{gd};
    die "ERROR: --ad does not exist: $opt->{ad}\n" unless -d $opt->{ad};
    die "ERROR: --lst does not exist: $opt->{lst}\n" unless -s $opt->{lst};
    die "ERROR: --t must be >= 1\n" if $opt->{t} < 1;
    die "ERROR: --minscore must be >= 0\n" if $opt->{minscore} < 0;
    die "ERROR: --minSpace must be >= 0\n" if $opt->{minSpace} < 0;
}

sub make_abs_dir {
    my ($dir) = @_;
    if (!-d $dir) {
        make_path($dir) or die "ERROR: cannot create directory $dir: $!\n";
    }
    return abs_path($dir) || die "ERROR: cannot resolve directory $dir\n";
}

sub normalize_tool_dir {
    my ($dir) = @_;
    $dir //= '';
    return '' if $dir eq '' || $dir eq 'PATH';
    $dir =~ s{/+$}{};
    return $dir;
}

sub find_tool {
    my ($dir, $name) = @_;
    $dir = normalize_tool_dir($dir);

    if ($dir ne '') {
        my $path = "$dir/$name";
        return $path if -x $path;
    }

    my $which = `command -v $name 2>/dev/null`;
    chomp $which;
    return $which if $which && -x $which;

    my $where = $dir ne '' ? "$dir/$name or PATH" : 'PATH';
    die "ERROR: cannot find executable '$name' in $where\n";
}

sub check_linear_gap {
    my ($linearGap) = @_;
    return if $linearGap !~ m{/};  # UCSC also accepts built-in names in some installations.
    die "ERROR: --linearGap file does not exist: $linearGap\n" unless -s $linearGap;
}

sub run_cmd {
    my ($LOG, $dry_run, @cmd) = @_;
    logmsg($LOG, 'CMD: ' . join(' ', map { shell_quote($_) } @cmd));
    return if $dry_run;
    system(@cmd) == 0 or die "ERROR: command failed with exit code " . ($? >> 8) . ": " . join(' ', @cmd) . "\n";
}

sub run_cmd_stdout {
    my ($LOG, $dry_run, $out, @cmd) = @_;
    logmsg($LOG, 'CMD: ' . join(' ', map { shell_quote($_) } @cmd) . ' > ' . shell_quote($out));
    return if $dry_run;

    open(my $OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";
    my $pid = open(my $IN, '-|', @cmd);
    die "ERROR: cannot run command: @cmd: $!\n" unless defined $pid;
    while (my $line = <$IN>) {
        print $OUT $line;
    }
    close $IN or die "ERROR: command failed with exit code " . ($? >> 8) . ": " . join(' ', @cmd) . "\n";
    close $OUT;
}

sub shell_quote {
    my ($s) = @_;
    return "''" if !defined($s) || $s eq '';
    return $s if $s =~ /^[A-Za-z0-9_\.\/\:\=\-]+$/;
    $s =~ s/'/'\\''/g;
    return "'$s'";
}

sub logmsg {
    my ($fh, $msg) = @_;
    my $t = scalar(localtime());
    print $fh "[$t] $msg\n";
    print STDERR "[$t] $msg\n";
}

sub usage {
    return <<'USAGE';

Usage:
  perl 2_ChainNetSyn_standardized.pl \
      --gd /path/to/genomes \
      --ad /path/to/axtfiles \
      --lst genome.lst \
      --outdir ChainNetSyn_out \
      --t 8 \
      --sing

Required:
  --gd        Directory containing genome FASTA files.
              Each genome ID in --lst should have one FASTA named ID, ID.fa,
              ID.fasta, or ID.fna.
  --ad        Directory containing pairwise AXT files named TargetIDvsQueryID.axt.
  --lst       Genome ID list, one genome ID per line.

Common options:
  --outdir    Output directory. Default: ChainNetSyn_out
  --syn       Backward-compatible alias of --outdir.
  --sing      Also generate target-oriented single-coverage MAF files with single_cov2.
  --t         Number of parallel jobs. Default: 1

UCSC/TBA utilities:
  --ucsc      Directory containing UCSC utilities, or PATH. Default: script/pub/UCSC/
  --tba       Directory containing single_cov2, or PATH. Default: script/pub/TBA/
  --linearGap ChainNet linearGap file/name. Default: script/pub/TBA/medium
  --lg        Backward-compatible alias of --linearGap.
  --minscore  minScore used by axtChain and chainNet. Default: 1000
  --minSpace  minSpace used by chainNet. Default: 1

Output structure:
  OUTDIR/chain/     *.chain.filter, and optionally raw *.chain
  OUTDIR/net/       *.chain.filter.tnet and *.chain.filter.qnet
  OUTDIR/synnet/    *.chain.filter.tnet.synnet and *.chain.filter.qnet.synnet
  OUTDIR/singRaw/   *.sing.maf, only when --sing is used
  OUTDIR/tmp/       temporary AXT/MAF files for --sing
  OUTDIR/ChainNetSyn.log

Notes:
  - Input AXT files should be named TargetIDvsQueryID.axt.
  - The TargetID and QueryID must match genome IDs in --lst.
  - The script automatically generates ID.sizes and ID.2bit in --gd when absent.
  - Raw chain files are removed by default after chainPreNet. Use --keep_raw_chain
    to keep them.
  - Temporary AXT/MAF files generated for --sing are removed by default. Use
    --keep_temp to keep them.
  - Use --force to overwrite existing outputs.
  - Use --dry_run to print commands without running them.

USAGE
}
