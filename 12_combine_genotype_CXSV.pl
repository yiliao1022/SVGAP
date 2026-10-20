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
use File::Find qw(find);
use File::Path qw(make_path);
use File::Basename qw(basename);
use POSIX qw(strftime);

################################################################################
# 12_combine_genotype_CXSV.pl
#
# Purpose:
#   Combine and genotype bidirectional SVGAP complex SV calls (*.cxsv.sort).
#
# Rationale:
#   PairwiseSV reports complex SVs when both the target/reference side and the
#   query/sample side contain unmatched sequence. For a reference-query pair,
#   useful evidence can exist in both directions:
#       REF.QUERY.rawSV/*.cxsv.sort
#       QUERY.REF.rawSV/*.cxsv.sort
#   This script projects both directions onto the original reference coordinate
#   system, merges redundant events, and outputs a population-level CXSV VCF.
#
# Input from 04_PairwiseSV.pl:
#   #Reference Ref_beg Ref_end SVtypes Sample Sample_beg Sample_end SVlength ...
#   ... Ref_SV_sequence Sample_SV_sequence
#   Coordinates are 0-based half-open, following SVGAP rawSV/BED convention.
#
# Output:
#   All.CXSV.records.projected.tsv
#   All.CXSV.combined.tsv
#   All.CXSV.genotype.matrix.tsv
#   All.CXSV.vcf                 legacy symbolic CXSV VCF
#   All.CXSV.vg.vcf              explicit sequence-resolved CXSV VCF when --ref_fasta is used
#   All.CXSV.vg.unresolved.tsv   clusters lacking a complete alternate sequence
#   combine_genotype_cxsv.log
#
# Genotype interpretation:
#   1|1: this sample has a projected raw CXSV record matching the merged site
#   0|0: no matching CXSV raw record was observed for this sample
#   ./. : optional missing call when --unknown_missing is used
#
# Author: Yi Liao; SVGAP module integrated with ChatGPT assistance
# Date: 2026-07-06
################################################################################

my %opt = (
    refname             => '',
    glist               => '',
    query_list          => '',
    rawsv_dir           => '.',
    outdir              => 'CXSV_genotype',
    break_tolerance     => 2000,
    min_ro              => 0.50,
    len_ratio           => 1.00,   # maximum relative difference for ref_len and alt_len when both are informative
    min_component_size  => 1,
    unknown_missing     => 0,
    keep_ref_sample     => 0,
    max_info_source_ids => 50,
    ref_fasta           => '',
    max_explicit_allele => 5_000_000,
    help                => 0,
);

GetOptions(
    'refname=s'             => \$opt{refname},
    'glist=s'               => \$opt{glist},
    'query_list|query-list=s' => \$opt{query_list},
    'rawsv_dir|rawsv-dir|svdir=s' => \$opt{rawsv_dir},
    'outdir=s'              => \$opt{outdir},
    'break_tolerance|break-tolerance=i' => \$opt{break_tolerance},
    'min_ro|min-ro=f'       => \$opt{min_ro},
    'len_ratio|len-ratio=f' => \$opt{len_ratio},
    'min_component_size|min-component-size=i' => \$opt{min_component_size},
    'unknown_missing|unknown-missing!' => \$opt{unknown_missing},
    'keep_ref_sample|keep-ref-sample!' => \$opt{keep_ref_sample},
    'max_info_source_ids|max-info-source-ids=i' => \$opt{max_info_source_ids},
    'ref_fasta|ref-fasta|ref=s' => \$opt{ref_fasta},
    'max_explicit_allele|max-explicit-allele=i' => \$opt{max_explicit_allele},
    'help|h'                => \$opt{help},
) or die usage();

if ($opt{help} || !$opt{refname} || (!$opt{glist} && !$opt{query_list})) {
    print usage();
    exit($opt{help} ? 0 : 1);
}

die "ERROR: --ref_fasta does not exist or is empty: $opt{ref_fasta}\n" if $opt{ref_fasta} && !-s $opt{ref_fasta};
die "ERROR: --max_explicit_allele must be >= 1\n" unless $opt{max_explicit_allele} >= 1;
make_path($opt{outdir}) unless -d $opt{outdir};
open(my $LOG, '>', "$opt{outdir}/combine_genotype_cxsv.log") or die "ERROR: cannot write log: $!\n";
logmsg($LOG, "CMD: $0 " . join(' ', map { "--$_=$opt{$_}" } sort keys %opt));

my $sample_list_file = $opt{query_list} || $opt{glist};
my @samples = read_sample_list($sample_list_file, $opt{refname}, $opt{keep_ref_sample});
die "ERROR: no query samples found in $sample_list_file after excluding $opt{refname}\n" unless @samples;
logmsg($LOG, "Reference: $opt{refname}");
logmsg($LOG, "Samples: " . join(',', @samples));

my @all_records;
for my $sample (@samples) {
    my @fwd_files = find_cxsv_files($opt{rawsv_dir}, "$opt{refname}.$sample.rawSV");
    my @rev_files = find_cxsv_files($opt{rawsv_dir}, "$sample.$opt{refname}.rawSV");
    logmsg($LOG, "$sample forward_files=" . scalar(@fwd_files) . " reverse_files=" . scalar(@rev_files));

    for my $file (@fwd_files) {
        my @r = parse_cxsv_file($file, $sample, 'forward', \%opt, $LOG);
        push @all_records, @r;
        logmsg($LOG, "  parsed forward $file records=" . scalar(@r));
    }
    for my $file (@rev_files) {
        my @r = parse_cxsv_file($file, $sample, 'reverse', \%opt, $LOG);
        push @all_records, @r;
        logmsg($LOG, "  parsed reverse $file records=" . scalar(@r));
    }
}

logmsg($LOG, "Projected CXSV records: " . scalar(@all_records));
write_projected_records("$opt{outdir}/All.CXSV.records.projected.tsv", \@all_records);

my @clusters = merge_records(\@all_records, \%opt);
assign_cluster_ids(\@clusters);
logmsg($LOG, "Combined CXSV clusters: " . scalar(@clusters));

write_combined_tsv("$opt{outdir}/All.CXSV.combined.tsv", \@clusters, $opt{max_info_source_ids});
write_genotype_matrix("$opt{outdir}/All.CXSV.genotype.matrix.tsv", \@clusters, \@samples, \%opt);
write_vcf("$opt{outdir}/All.CXSV.vcf", \@clusters, \@samples, \%opt);
if ($opt{ref_fasta}) {
    my $refdb = CXFastFAI->new($opt{ref_fasta});
    write_vg_vcf("$opt{outdir}/All.CXSV.vg.vcf", "$opt{outdir}/All.CXSV.vg.unresolved.tsv", \@clusters, \@samples, \%opt, $refdb);
}

logmsg($LOG, "Done");
close $LOG;
exit 0;

sub usage {
    return <<'USAGE';
Usage:
  perl 12_combine_genotype_CXSV.pl \
    --refname ISO1 \
    --glist genome.list \
    --rawsv_dir SVGAP_run/04_pairwise_rawSV \
    --outdir SVGAP_run/06_genotype/CXSV

Required:
  --refname STR        Original reference genome ID, e.g. ISO1
  --glist FILE         Genome list, one sample/genome per line. Reference is removed.
                       Alternatively use --query_list with query samples only.
  --rawsv_dir DIR      Directory containing pairwise *.rawSV folders from 04_PairwiseSV.pl

Main options:
  --break_tolerance N  Breakpoint tolerance for merging CXSVs [2000]
  --min_ro FLOAT       Minimum reciprocal overlap for interval merging [0.50]
  --len_ratio FLOAT    Maximum relative difference for ref_len and alt_len [1.00]
  --unknown_missing    Output ./. for non-carrier samples instead of 0|0
  --keep_ref_sample    Keep reference as a VCF sample column [off]
  --ref_fasta FILE     Reference FASTA; enables explicit sequence-resolved All.CXSV.vg.vcf
  --max_explicit_allele N  Maximum REF or ALT length written explicitly [5000000]

Input files searched recursively:
  REF.QUERY.rawSV/*.cxsv.sort
  QUERY.REF.rawSV/*.cxsv.sort
  also accepts .cxsv.sort.txt and .cxsv.sort.tsv

Outputs:
  All.CXSV.records.projected.tsv
  All.CXSV.combined.tsv
  All.CXSV.genotype.matrix.tsv
  All.CXSV.vcf                 legacy symbolic CXSV VCF
  All.CXSV.vg.vcf              explicit sequence-resolved CXSV VCF when --ref_fasta is used
  All.CXSV.vg.unresolved.tsv   clusters lacking a complete alternate sequence
USAGE
}

sub logmsg {
    my ($LOG, $msg) = @_;
    print $LOG '[' . strftime('%Y-%m-%d %H:%M:%S', localtime) . "] $msg\n";
}

sub read_sample_list {
    my ($file, $ref, $keep_ref) = @_;
    open(my $fh, '<', $file) or die "ERROR: cannot read $file: $!\n";
    my (@samples, %seen);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        $line =~ s/#.*$//;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '';
        my ($name) = split(/\s+/, $line);
        $name = basename($name);
        $name =~ s/\.gz$//i;
        $name =~ s/\.(fa|fasta|fna|fas)$//i;
        next if !$keep_ref && $name eq $ref;
        next if $seen{$name}++;
        push @samples, $name;
    }
    close $fh;
    return @samples;
}

sub find_cxsv_files {
    my ($root, $dirname) = @_;
    my @files;
    return @files unless -d $root;
    find({
        wanted => sub {
            return unless -f $_;
            my $path = $File::Find::name;
            return unless $path =~ /\Q$dirname\E(?:\/|\\)/;
            return unless basename($path) =~ /\.cxsv\.sort(?:\.(?:txt|tsv))?$/i;
            push @files, $path;
        },
        no_chdir => 1,
    }, $root);
    return sort @files;
}

sub parse_cxsv_file {
    my ($file, $sample, $direction, $opt, $LOG) = @_;
    open(my $fh, '<', $file) or do {
        logmsg($LOG, "WARNING: cannot open $file: $!");
        return ();
    };

    my (@records, %hidx, $line_no);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;

        if ($line =~ /^#/) {
            my $h = $line;
            $h =~ s/^#+//;
            $h =~ s/^\s+|\s+$//g;
            if ($h =~ /Reference|Ref_beg|Sample_beg|SVtypes|CHROM|POS|SVTYPE/i) {
                my @header = split(/\t|\s+/, $h);
                %hidx = build_header_index(@header);
            }
            next;
        }

        $line_no++;
        my @f = split(/\t/, $line, -1);
        @f = split(/\s+/, $line) if @f < 8;
        next if @f < 7;

        # SVGAP 04_PairwiseSV.pl header:
        # Reference Ref_beg Ref_end SVtypes Sample Sample_beg Sample_end SVlength Strand ... Ref_SV_sequence Sample_SV_sequence
        my $src_chrom = get_col(\@f, \%hidx, [('REFERENCE','CHROM','#CHROM','T_CHR','TARGET','TARGET_CHR','RNAME','REF_CHR')], 0);
        my $src_pos   = get_col(\@f, \%hidx, [qw(REF_BEG REF_START POS START T_POS TPOS TARGET_START)], 1);
        my $src_end   = get_col(\@f, \%hidx, [qw(REF_END END STOP T_END TEND TARGET_END)], 2);
        my $svtype    = get_col(\@f, \%hidx, [qw(SVTYPES SVTYPE TYPE)], 3);
        my $q_chrom   = get_col(\@f, \%hidx, [qw(SAMPLE QUERY Q_CHR QCHR QNAME QUERY_CHR ALT_CHR)], 4);
        my $q_pos     = get_col(\@f, \%hidx, [qw(SAMPLE_BEG SAMPLE_START Q_POS QPOS Q_START QSTART QUERY_START ALT_START)], 5);
        my $q_end     = get_col(\@f, \%hidx, [qw(SAMPLE_END Q_END QEND Q_STOP QSTOP QUERY_END ALT_END)], 6);
        my $svlen     = get_col(\@f, \%hidx, [qw(SVLENGTH SVLEN LEN LENGTH)], 7);
        my $strand    = get_col(\@f, \%hidx, [qw(STRAND)], 8);
        my $refseq    = get_col(\@f, \%hidx, [qw(REF_SV_SEQUENCE REFSEQ REFSEQ REF_SEQUENCE REF)], 15);
        my $altseq    = get_col(\@f, \%hidx, [qw(SAMPLE_SV_SEQUENCE ALTSEQ ALTSEQ ALT_SEQUENCE QUERYSEQ QUERY_SEQ)], 16);
        my $id        = get_col(\@f, \%hidx, [qw(ID NAME VARIANT_ID)], undef);

        my $t_pos = to_int($src_pos);
        my $t_end = to_int($src_end);
        my $s_pos = to_int($q_pos);
        my $s_end = to_int($q_end);
        next unless defined $src_chrom && defined $t_pos && defined $t_end;
        ($t_pos, $t_end) = sort_num_pair($t_pos, $t_end);
        ($s_pos, $s_end) = sort_num_pair($s_pos, $s_end) if defined $s_pos && defined $s_end;

        my ($ref_chrom, $ref_start, $ref_end, $alt_chrom, $alt_start, $alt_end);
        if ($direction eq 'forward') {
            $ref_chrom = $src_chrom;
            $ref_start = $t_pos;
            $ref_end   = $t_end;
            $alt_chrom = defined($q_chrom) ? $q_chrom : $sample;
            $alt_start = $s_pos;
            $alt_end   = $s_end;
        } else {
            # QUERY.REF.rawSV: the Sample/query columns are the original reference.
            if (defined $q_chrom && defined $s_pos && defined $s_end) {
                $ref_chrom = $q_chrom;
                $ref_start = $s_pos;
                $ref_end   = $s_end;
                $alt_chrom = $src_chrom;
                $alt_start = $t_pos;
                $alt_end   = $t_end;
            } else {
                logmsg($LOG, "WARNING: skip reverse record lacking Sample coordinates: $file line $line_no");
                next;
            }
        }

        ($ref_start, $ref_end) = sort_num_pair($ref_start, $ref_end);
        ($alt_start, $alt_end) = sort_num_pair($alt_start, $alt_end) if defined $alt_start && defined $alt_end;
        my $ref_len = $ref_end - $ref_start;
        my $alt_len = (defined $alt_start && defined $alt_end) ? ($alt_end - $alt_start) : abs_int($svlen);
        $alt_len = sequence_len($altseq) if (!defined($alt_len) || $alt_len <= 0) && defined($altseq);
        next if $ref_len < $opt->{min_component_size} && (!defined($alt_len) || $alt_len < $opt->{min_component_size});

        $svtype = 'CXSV' unless defined($svtype) && $svtype ne '';
        $id = basename($file) . ":$line_no" unless defined($id) && $id ne '';

        my ($projected_refseq, $projected_altseq) = ($refseq, $altseq);
        if ($direction eq 'reverse') {
            # Raw target is the non-reference assembly and raw Sample is the
            # original reference. Swap the sequence columns to the projected
            # reference/alternate orientation. PairwiseSV already oriented the
            # raw Sample sequence to the raw target; for a '-' relationship,
            # reverse-complement both alleles back to original-reference forward.
            ($projected_refseq, $projected_altseq) = ($altseq, $refseq);
            if (defined($strand) && $strand eq '-') {
                $projected_refseq = revcom_seq($projected_refseq);
                $projected_altseq = revcom_seq($projected_altseq);
            }
        }

        push @records, {
            sample      => $sample,
            direction   => $direction,
            ref_chrom   => $ref_chrom,
            ref_start   => $ref_start,
            ref_end     => $ref_end,
            ref_len     => $ref_len,
            alt_chrom   => defined($alt_chrom) ? $alt_chrom : '.',
            alt_start   => defined($alt_start) ? $alt_start : '.',
            alt_end     => defined($alt_end) ? $alt_end : '.',
            alt_len     => defined($alt_len) ? $alt_len : 0,
            strand      => defined($strand) ? $strand : '.',
            source_type => $svtype,
            source_id   => $id,
            source_file => $file,
            source_line => $line_no,
            refseq      => defined($projected_refseq) ? $projected_refseq : '.',
            altseq      => defined($projected_altseq) ? $projected_altseq : '.',
        };
    }
    close $fh;
    return @records;
}

sub build_header_index {
    my @h = @_;
    my %idx;
    for my $i (0..$#h) {
        my $k = norm_key($h[$i]);
        $idx{$k} = $i;
    }
    return %idx;
}

sub norm_key {
    my ($k) = @_;
    $k = '' unless defined $k;
    $k =~ s/^#//;
    $k = uc($k);
    $k =~ s/[^A-Z0-9_]+/_/g;
    return $k;
}

sub get_col {
    my ($fields, $hidx, $names, $fallback_index) = @_;
    for my $n (@$names) {
        my $k = norm_key($n);
        if (exists $hidx->{$k} && defined $fields->[$hidx->{$k}]) {
            return normalize_empty($fields->[$hidx->{$k}]);
        }
    }
    if (defined $fallback_index && defined $fields->[$fallback_index]) {
        return normalize_empty($fields->[$fallback_index]);
    }
    return undef;
}

sub normalize_empty {
    my ($v) = @_;
    return undef unless defined $v;
    $v =~ s/^\s+|\s+$//g;
    return undef if $v eq '' || $v eq '.' || uc($v) eq 'NA' || $v eq 'NaN';
    return $v;
}

sub to_int {
    my ($v) = @_;
    return undef unless defined $v;
    return undef unless $v =~ /^-?\d+(?:\.\d+)?$/;
    return int($v + 0);
}

sub abs_int {
    my ($v) = @_;
    return undef unless defined $v;
    return undef unless $v =~ /^-?\d+(?:\.\d+)?$/;
    return abs(int($v + 0));
}

sub sequence_len {
    my ($s) = @_;
    return undef unless defined $s && $s ne '.';
    $s =~ s/[^A-Za-z]//g;
    return length($s);
}

sub sort_num_pair {
    my ($a, $b) = @_;
    return ($a, $b) unless defined $a && defined $b;
    return $a <= $b ? ($a, $b) : ($b, $a);
}

sub write_projected_records {
    my ($file, $records) = @_;
    open(my $out, '>', $file) or die "ERROR: cannot write $file: $!\n";
    print $out join("\t", ('#CHROM', qw(POS END SAMPLE DIRECTION REF_LEN ALT_CHROM ALT_POS ALT_END ALT_LEN STRAND SOURCE_TYPE SOURCE_ID SOURCE_FILE SOURCE_LINE))) . "\n";
    for my $r (sort_records(@$records)) {
        print $out join("\t", @{$r}{qw(ref_chrom ref_start ref_end sample direction ref_len alt_chrom alt_start alt_end alt_len strand source_type source_id source_file source_line)}) . "\n";
    }
    close $out;
}

sub sort_records {
    return sort {
        $a->{ref_chrom} cmp $b->{ref_chrom} ||
        $a->{ref_start} <=> $b->{ref_start} ||
        $a->{ref_end}   <=> $b->{ref_end} ||
        $a->{sample}    cmp $b->{sample} ||
        $a->{direction} cmp $b->{direction}
    } @_;
}

sub merge_records {
    my ($records, $opt) = @_;
    my @sorted = sort_records(@$records);
    my @clusters;

    REC:
    for my $r (@sorted) {
        for (my $i = $#clusters; $i >= 0; $i--) {
            my $c = $clusters[$i];
            last if $c->{chrom} lt $r->{ref_chrom};
            next if $c->{chrom} ne $r->{ref_chrom};
            last if $c->{end_max} + $opt->{break_tolerance} < $r->{ref_start};
            if (compatible($c, $r, $opt)) {
                add_to_cluster($c, $r);
                next REC;
            }
        }
        my %c = (
            chrom     => $r->{ref_chrom},
            start_min => $r->{ref_start},
            start_max => $r->{ref_start},
            end_min   => $r->{ref_end},
            end_max   => $r->{ref_end},
            records   => [],
        );
        add_to_cluster(\%c, $r);
        push @clusters, \%c;
    }

    finalize_clusters(\@clusters);
    @clusters = sort { $a->{chrom} cmp $b->{chrom} || $a->{pos} <=> $b->{pos} || $a->{end} <=> $b->{end} } @clusters;
    return @clusters;
}

sub compatible {
    my ($c, $r, $opt) = @_;
    return 0 unless $c->{chrom} eq $r->{ref_chrom};

    my $bp_ok = (abs($c->{pos} - $r->{ref_start}) <= $opt->{break_tolerance} && abs($c->{end} - $r->{ref_end}) <= $opt->{break_tolerance});
    my $ov = interval_overlap($c->{pos}, $c->{end}, $r->{ref_start}, $r->{ref_end});
    my $len_c = max2(1, $c->{end} - $c->{pos});
    my $len_r = max2(1, $r->{ref_end} - $r->{ref_start});
    my $ro = $ov > 0 ? min2($ov / $len_c, $ov / $len_r) : 0;
    return 0 unless $bp_ok || $ro >= $opt->{min_ro};

    # Do not over-filter, but avoid merging very different complex alleles when lengths are informative.
    if ($opt->{len_ratio} >= 0) {
        return 0 unless length_compatible($c->{ref_len} || $len_c, $r->{ref_len}, $opt->{len_ratio});
        if (($c->{alt_len} || 0) > 0 && ($r->{alt_len} || 0) > 0) {
            return 0 unless length_compatible($c->{alt_len}, $r->{alt_len}, $opt->{len_ratio});
        }
    }
    return 1;
}

sub length_compatible {
    my ($a, $b, $max_diff) = @_;
    return 1 unless defined $a && defined $b && $a > 0 && $b > 0;
    my $min = min2($a, $b);
    my $max = max2($a, $b);
    return (($max - $min) / $max) <= $max_diff ? 1 : 0;
}

sub add_to_cluster {
    my ($c, $r) = @_;
    push @{$c->{records}}, $r;
    $c->{start_min} = min2($c->{start_min}, $r->{ref_start});
    $c->{start_max} = max2($c->{start_max}, $r->{ref_start});
    $c->{end_min}   = min2($c->{end_min},   $r->{ref_end});
    $c->{end_max}   = max2($c->{end_max},   $r->{ref_end});
    my @starts = map { $_->{ref_start} } @{$c->{records}};
    my @ends   = map { $_->{ref_end} } @{$c->{records}};
    my @rlens  = map { $_->{ref_len} } @{$c->{records}};
    my @alens  = map { $_->{alt_len} } @{$c->{records}};
    $c->{pos}     = median_int(@starts);
    $c->{end}     = median_int(@ends);
    $c->{ref_len} = median_int(@rlens);
    $c->{alt_len} = median_int(@alens);
}

sub finalize_clusters {
    my ($clusters) = @_;
    for my $c (@$clusters) {
        add_summary_fields($c);
    }
}

sub add_summary_fields {
    my ($c) = @_;
    my @starts = map { $_->{ref_start} } @{$c->{records}};
    my @ends   = map { $_->{ref_end} } @{$c->{records}};
    my @rlens  = map { $_->{ref_len} } @{$c->{records}};
    my @alens  = grep { defined $_ && $_ ne '.' } map { $_->{alt_len} } @{$c->{records}};
    $c->{pos} = median_int(@starts);
    $c->{end} = median_int(@ends);
    ($c->{pos}, $c->{end}) = sort_num_pair($c->{pos}, $c->{end});
    $c->{ref_len} = median_int(@rlens);
    $c->{alt_len} = @alens ? median_int(@alens) : 0;
    my (%samples, %dirs, %types, %strands);
    for my $r (@{$c->{records}}) {
        $samples{$r->{sample}} = 1;
        $dirs{$r->{direction}} = 1;
        $types{$r->{source_type}} = 1 if defined $r->{source_type};
        $strands{$r->{strand}} = 1 if defined $r->{strand};
    }
    $c->{supp_n}       = scalar(keys %samples);
    $c->{samples}      = [sort keys %samples];
    $c->{directions}   = [sort keys %dirs];
    $c->{source_types} = [sort keys %types];
    $c->{strands}      = [sort keys %strands];
}

sub assign_cluster_ids {
    my ($clusters) = @_;
    my $n = 1;
    for my $c (@$clusters) {
        $c->{id} = sprintf('CXSV%07d', $n++);
    }
}

sub write_combined_tsv {
    my ($file, $clusters, $max_ids) = @_;
    open(my $out, '>', $file) or die "ERROR: cannot write $file: $!\n";
    print $out join("\t", ('#CHROM', qw(POS END ID SVTYPE REFLEN ALTLEN SUPP_N SUPP_SAMPLES DIRS STRANDS SOURCE_TYPES N_RECORDS START_RANGE END_RANGE SOURCE_IDS))) . "\n";
    for my $c (@$clusters) {
        print $out join("\t",
            $c->{chrom}, $c->{pos}, $c->{end}, $c->{id}, 'CXSV', $c->{ref_len}, $c->{alt_len},
            $c->{supp_n}, join(',', @{$c->{samples}}), join(',', @{$c->{directions}}), join(',', @{$c->{strands}}), join(',', @{$c->{source_types}}),
            scalar(@{$c->{records}}), $c->{start_min} . '-' . $c->{start_max}, $c->{end_min} . '-' . $c->{end_max},
            join(',', source_id_list($c, $max_ids))
        ) . "\n";
    }
    close $out;
}

sub write_genotype_matrix {
    my ($file, $clusters, $samples, $opt) = @_;
    my @cols = @$samples;
    unshift @cols, $opt->{refname} if $opt->{keep_ref_sample};
    open(my $out, '>', $file) or die "ERROR: cannot write $file: $!\n";
    print $out join("\t", ('#CHROM', qw(POS END ID SVTYPE REFLEN ALTLEN SUPP_N SUPP_SAMPLES)), @cols) . "\n";
    for my $c (@$clusters) {
        my @gt = gt_for_cols($c, \@cols, $opt);
        print $out join("\t", $c->{chrom}, $c->{pos}, $c->{end}, $c->{id}, 'CXSV', $c->{ref_len}, $c->{alt_len}, $c->{supp_n}, join(',', @{$c->{samples}}), @gt) . "\n";
    }
    close $out;
}

sub write_vcf {
    my ($file, $clusters, $samples, $opt) = @_;
    my @cols = @$samples;
    unshift @cols, $opt->{refname} if $opt->{keep_ref_sample};
    open(my $out, '>', $file) or die "ERROR: cannot write $file: $!\n";
    print $out "##fileformat=VCFv4.2\n";
    print $out "##source=12_combine_genotype_CXSV.pl\n";
    print $out "##ALT=<ID=CXSV,Description=\"Complex structural variant with both reference-side and sample-side unmatched sequences\">\n";
    print $out "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
    print $out "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate on reference, converted from SVGAP 0-based half-open END to VCF-style END\">\n";
    print $out "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"ALTLEN minus REFLEN\">\n";
    print $out "##INFO=<ID=REFLEN,Number=1,Type=Integer,Description=\"Median reference-side unmatched length\">\n";
    print $out "##INFO=<ID=ALTLEN,Number=1,Type=Integer,Description=\"Median sample-side unmatched length\">\n";
    print $out "##INFO=<ID=SUPP_N,Number=1,Type=Integer,Description=\"Number of samples supporting the event\">\n";
    print $out "##INFO=<ID=SUPP_SAMPLES,Number=.,Type=String,Description=\"Samples supporting the event\">\n";
    print $out "##INFO=<ID=DIRS,Number=.,Type=String,Description=\"Calling directions supporting the event\">\n";
    print $out "##INFO=<ID=SOURCE_TYPES,Number=.,Type=String,Description=\"Original source SVtypes\">\n";
    print $out "##INFO=<ID=N_RECORDS,Number=1,Type=Integer,Description=\"Number of raw CXSV records merged\">\n";
    print $out "##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description=\"Limited sample:direction:source_id list\">\n";
    print $out "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $out join("\t", '#CHROM', qw(POS ID REF ALT QUAL FILTER INFO FORMAT), @cols) . "\n";
    for my $c (@$clusters) {
        my @gt = gt_for_cols($c, \@cols, $opt);
        my $svlen = ($c->{alt_len} || 0) - ($c->{ref_len} || ($c->{end} - $c->{pos}));
        my $info = join(';',
            'SVTYPE=CXSV', 'END=' . $c->{end}, 'SVLEN=' . $svlen,
            'REFLEN=' . $c->{ref_len}, 'ALTLEN=' . $c->{alt_len},
            'SUPP_N=' . $c->{supp_n}, 'SUPP_SAMPLES=' . escape_info(join(',', @{$c->{samples}})),
            'DIRS=' . escape_info(join(',', @{$c->{directions}})),
            'SOURCE_TYPES=' . escape_info(join(',', @{$c->{source_types}})),
            'N_RECORDS=' . scalar(@{$c->{records}}),
            'SOURCE_IDS=' . escape_info(join(',', source_id_list($c, $opt->{max_info_source_ids}))),
        );
        print $out join("\t", $c->{chrom}, $c->{pos} + 1, $c->{id}, 'N', '<CXSV>', '.', 'PASS', $info, 'GT', @gt) . "\n";
    }
    close $out;
}


sub write_vg_vcf {
    my ($file, $unresolved_file, $clusters, $samples, $opt, $refdb) = @_;
    my @cols = @$samples;
    unshift @cols, $opt->{refname} if $opt->{keep_ref_sample};
    open(my $out, '>', $file) or die "ERROR: cannot write $file: $!\n";
    open(my $bad, '>', $unresolved_file) or die "ERROR: cannot write $unresolved_file: $!\n";

    print $out "##fileformat=VCFv4.2\n";
    print $out "##source=12_combine_genotype_CXSV.pl_explicit_vg\n";
    print $out "##reference=$opt->{ref_fasta}\n";
    for my $chr ($refdb->contigs) {
        print $out "##contig=<ID=$chr,length=", $refdb->length($chr), ">\n";
    }
    print $out "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Sequence-resolved event class\">\n";
    print $out "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate derived from explicit REF span\">\n";
    print $out "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"ALT length minus REF length excluding shared anchor\">\n";
    print $out "##INFO=<ID=REFLEN,Number=1,Type=Integer,Description=\"Unanchored reference allele length\">\n";
    print $out "##INFO=<ID=ALTLEN,Number=1,Type=Integer,Description=\"Unanchored alternate allele length\">\n";
    print $out "##INFO=<ID=ORIG_SVTYPE,Number=.,Type=String,Description=\"Original SVGAP complex-SV classes\">\n";
    print $out "##INFO=<ID=REP_SAMPLE,Number=1,Type=String,Description=\"Sample supplying the explicit alternate sequence\">\n";
    print $out "##INFO=<ID=REP_DIRECTION,Number=1,Type=String,Description=\"Pairwise calling direction supplying the explicit alternate sequence\">\n";
    print $out "##INFO=<ID=SUPP_N,Number=1,Type=Integer,Description=\"Number of supporting samples\">\n";
    print $out "##INFO=<ID=SUPP_SAMPLES,Number=.,Type=String,Description=\"Supporting samples\">\n";
    print $out "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $out join("\t", '#CHROM', qw(POS ID REF ALT QUAL FILTER INFO FORMAT), @cols), "\n";
    print $bad join("\t", qw(ID CHROM CLUSTER_POS CLUSTER_END REASON N_RECORDS SOURCE_TYPES)), "\n";

    my ($written,$unresolved)=(0,0);
    for my $c (@$clusters) {
        my ($rep,$event_ref,$event_alt,$reason) = choose_explicit_record($c,$refdb,$opt);
        if (!$rep) {
            $unresolved++;
            print $bad join("\t", $c->{id},$c->{chrom},$c->{pos},$c->{end},$reason,scalar(@{$c->{records}}),join(',',@{$c->{source_types}})), "\n";
            next;
        }
        my ($pos,$ref,$alt,$end,$svlen,$err) = anchor_cxsv_alleles($refdb,$rep->{ref_chrom},$rep->{ref_start},$rep->{ref_end},$event_ref,$event_alt);
        if ($err) {
            $unresolved++;
            print $bad join("\t", $c->{id},$c->{chrom},$c->{pos},$c->{end},$err,scalar(@{$c->{records}}),join(',',@{$c->{source_types}})), "\n";
            next;
        }
        if (length($ref) > $opt->{max_explicit_allele} || length($alt) > $opt->{max_explicit_allele}) {
            $unresolved++;
            print $bad join("\t", $c->{id},$c->{chrom},$c->{pos},$c->{end},'explicit allele exceeds max_explicit_allele',scalar(@{$c->{records}}),join(',',@{$c->{source_types}})), "\n";
            next;
        }
        my $type = length($event_ref)==0 ? 'INS' : length($event_alt)==0 ? 'DEL' : 'CXSV';
        my @gt = gt_for_cols($c, \@cols, $opt);
        my $info = join(';',
            "SVTYPE=$type", "END=$end", "SVLEN=$svlen",
            'REFLEN=' . length($event_ref), 'ALTLEN=' . length($event_alt),
            'ORIG_SVTYPE=' . escape_info(join(',',@{$c->{source_types}})),
            'REP_SAMPLE=' . escape_info($rep->{sample}),
            'REP_DIRECTION=' . escape_info($rep->{direction}),
            'SUPP_N=' . $c->{supp_n},
            'SUPP_SAMPLES=' . escape_info(join(',',@{$c->{samples}})),
        );
        print $out join("\t",$rep->{ref_chrom},$pos,$c->{id},$ref,$alt,'.','PASS',$info,'GT',@gt),"\n";
        $written++;
    }
    close $out;
    close $bad;
    print STDERR "CXSV explicit vg VCF: written=$written unresolved=$unresolved file=$file\n";
}

sub choose_explicit_record {
    my ($c,$refdb,$opt)=@_;
    my @ranked;
    for my $r (@{$c->{records}}) {
        next unless $refdb->has_contig($r->{ref_chrom});
        next if $r->{ref_start} < 0 || $r->{ref_end} < $r->{ref_start} || $r->{ref_end} > $refdb->length($r->{ref_chrom});
        my $ref_len=$r->{ref_end}-$r->{ref_start};
        my $alt_len=defined($r->{alt_len}) ? $r->{alt_len}+0 : 0;
        next if $ref_len > $opt->{max_explicit_allele} || $alt_len > $opt->{max_explicit_allele};
        my $event_ref=$refdb->get_seq($r->{ref_chrom},$r->{ref_start},$r->{ref_end});
        my ($event_alt,$why)=recover_alt_event_sequence($r->{altseq},$alt_len);
        next unless defined $event_alt;
        my $score=abs($r->{ref_start}-$c->{pos})+abs($r->{ref_end}-$c->{end});
        $score += 100 if $r->{direction} ne 'forward';
        $score += 10 if $r->{sample} ne ($c->{samples}[0] // '');
        push @ranked,[$score,$r,$event_ref,$event_alt,$why];
    }
    return (undef,undef,undef,'no record has a complete alternate sequence matching its projected length') unless @ranked;
    @ranked=sort {$a->[0]<=>$b->[0] || $a->[1]{sample} cmp $b->[1]{sample}} @ranked;
    return ($ranked[0][1],$ranked[0][2],$ranked[0][3],$ranked[0][4]);
}

sub recover_alt_event_sequence {
    my ($seq,$expected)=@_;
    $expected=0 unless defined $expected && $expected>=0;
    return ('','empty alternate interval') if $expected==0 && (!defined($seq) || $seq eq '' || $seq eq '.');
    return (undef,'missing alternate sequence') unless defined $seq;
    $seq=uc($seq); $seq =~ s/\s+//g;
    return (undef,'alternate sequence contains non-IUPAC DNA characters') unless $seq =~ /^[ACGTN]+$/;
    if (length($seq)==$expected) { return ($seq,'exact sequence length'); }
    if (length($seq)==$expected+1) { return (substr($seq,0,$expected),'trimmed legacy +1 right anchor'); }
    if ($expected==0 && length($seq)==1) { return ('','trimmed legacy insertion anchor'); }
    return (undef,'alternate sequence length does not match projected ALT interval');
}

sub anchor_cxsv_alleles {
    my ($refdb,$chrom,$start0,$end0,$event_ref,$event_alt)=@_;
    return (undef,undef,undef,undef,undef,'reference event sequence mismatch') if uc($refdb->get_seq($chrom,$start0,$end0)) ne uc($event_ref);
    my ($pos,$ref,$alt);
    if ($start0>0) {
        my $anchor=$refdb->get_seq($chrom,$start0-1,$start0);
        return (undef,undef,undef,undef,undef,'cannot fetch left anchor') unless defined($anchor) && $anchor =~ /^[ACGTN]$/i;
        $pos=$start0; $ref=uc($anchor.$event_ref); $alt=uc($anchor.$event_alt);
    } else {
        return (undef,undef,undef,undef,undef,'cannot right-anchor event spanning entire contig') if $end0 >= $refdb->length($chrom);
        my $anchor=$refdb->get_seq($chrom,$end0,$end0+1);
        return (undef,undef,undef,undef,undef,'cannot fetch right anchor') unless defined($anchor) && $anchor =~ /^[ACGTN]$/i;
        $pos=1; $ref=uc($event_ref.$anchor); $alt=uc($event_alt.$anchor);
    }
    my $end=$pos+length($ref)-1;
    my $svlen=length($alt)-length($ref);
    return ($pos,$ref,$alt,$end,$svlen,undef);
}

sub revcom_seq {
    my ($seq)=@_;
    return $seq unless defined($seq) && $seq ne '.';
    my $x=reverse($seq); $x =~ tr/ACGTNacgtn/TGCANtgcan/; return $x;
}

sub gt_for_cols {
    my ($c, $cols, $opt) = @_;
    my %carrier = map { $_ => 1 } @{$c->{samples}};
    my @gt;
    for my $s (@$cols) {
        if ($s eq $opt->{refname}) {
            push @gt, '0|0';
        } elsif ($carrier{$s}) {
            push @gt, '1|1';
        } else {
            push @gt, $opt->{unknown_missing} ? './.' : '0|0';
        }
    }
    return @gt;
}

sub source_id_list {
    my ($c, $max_ids) = @_;
    my @src;
    my $n = 0;
    for my $r (@{$c->{records}}) {
        last if $n >= $max_ids;
        my $id = join(':', $r->{sample}, $r->{direction}, $r->{source_id});
        $id =~ s/[\s;,=]+/_/g;
        push @src, $id;
        $n++;
    }
    push @src, 'TRUNCATED' if scalar(@{$c->{records}}) > $max_ids;
    return @src;
}

sub escape_info {
    my ($s) = @_;
    $s = '.' unless defined $s && $s ne '';
    $s =~ s/\s+/_/g;
    $s =~ s/;/,/g;
    $s =~ s/=/-/g;
    return $s;
}

sub interval_overlap {
    my ($s1, $e1, $s2, $e2) = @_;
    my $s = max2($s1, $s2);
    my $e = min2($e1, $e2);
    return $e > $s ? $e - $s : 0;
}

sub median_int {
    my @v = sort { $a <=> $b } grep { defined $_ && $_ ne '' && $_ ne '.' } @_;
    return 0 unless @v;
    my $n = @v;
    return $n % 2 ? int($v[int($n/2)] + 0.5) : int(($v[$n/2 - 1] + $v[$n/2]) / 2 + 0.5);
}

sub min2 { $_[0] <= $_[1] ? $_[0] : $_[1] }
sub max2 { $_[0] >= $_[1] ? $_[0] : $_[1] }


################################################################################
# Lightweight FASTA random access for explicit CXSV alleles.
################################################################################
package CXFastFAI;
use strict;
use warnings;
use Fcntl qw(SEEK_SET);
sub new {
    my ($class,$fasta)=@_;
    die "ERROR: CXSV explicit output requires uncompressed FASTA: $fasta\n" if $fasta =~ /\.gz$/i;
    my $self=bless {fasta=>$fasta,idx=>{},order=>[]},$class;
    if (-s "$fasta.fai") { $self->_load_fai("$fasta.fai"); } else { $self->_build_index(); }
    open(my $fh,'<',$fasta) or die "ERROR: cannot open FASTA $fasta: $!\n"; binmode($fh); $self->{fh}=$fh; return $self;
}
sub _load_fai { my ($self,$fai)=@_; open(my $IN,'<',$fai) or die "ERROR: cannot read $fai: $!\n"; while(<$IN>){chomp;next if /^\s*$/;my @f=split(/\t/);next unless @f>=5;$self->{idx}{$f[0]}={len=>$f[1]+0,offset=>$f[2]+0,line_bases=>$f[3]+0,line_width=>$f[4]+0};push @{$self->{order}},$f[0];}close $IN; }
sub _build_index { my ($self)=@_;open(my $IN,'<',$self->{fasta}) or die $!;binmode($IN);my($name,$len,$offset,$lb,$lw);while(1){my$p=tell($IN);my$line=<$IN>;last unless defined$line;if($line=~/^>(\S+)/){if(defined$name){$self->{idx}{$name}={len=>$len,offset=>$offset,line_bases=>$lb,line_width=>$lw};push@{$self->{order}},$name;}$name=$1;($len,$offset,$lb,$lw)=(0,undef,undef,undef);}else{next unless defined$name;my$raw=CORE::length($line);$line=~s/[\r\n]+$//;my$b=CORE::length($line);next if$b==0;$offset=$p unless defined$offset;$lb=$b unless defined$lb;$lw=$raw unless defined$lw;$len+=$b;}}if(defined$name){$self->{idx}{$name}={len=>$len,offset=>$offset,line_bases=>$lb,line_width=>$lw};push@{$self->{order}},$name;}close$IN; }
sub contigs { return @{$_[0]{order}}; }
sub has_contig { return exists $_[0]{idx}{$_[1]}; }
sub length { return $_[0]{idx}{$_[1]}{len}; }
sub get_seq { my($self,$chr,$start,$end)=@_;return undef unless$self->has_contig($chr);my$x=$self->{idx}{$chr};$start=0 if$start<0;$end=$x->{len} if$end>$x->{len};return''if$end<=$start;my$need=$end-$start;my$line=$x->{line_bases};my$width=$x->{line_width};my$byte=$x->{offset}+int($start/$line)*$width+($start%$line);seek($self->{fh},$byte,SEEK_SET)or die"ERROR: FASTA seek failed\n";my$seq='';while(CORE::length($seq)<$need){my$buf='';my$n=read($self->{fh},$buf,$need-CORE::length($seq)+int(($need-CORE::length($seq))/$line+2)*($width-$line));last unless$n;$buf=~s/[\r\n]//g;$seq.=$buf;}return uc(substr($seq,0,$need)); }
1;
