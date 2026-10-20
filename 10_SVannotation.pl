#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
# ============================================================================
# SVGAP standardized script: 10_SVannotation.pl
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
use Getopt::Long;
use File::Path qw(make_path);
use Fcntl qw(SEEK_SET);

################################################################################
# 10_SVannotation_final.pl
#
# Fast SV annotation from VCF.
#
# Main annotations:
#   1. Tandem duplication / inverted tandem duplication-like events
#   2. TE-related SV sequence annotation by minimap2 against a TE library
#   3. CDS-related SV sequence annotation by minimap2 against a CDS FASTA
#
# Main features:
#   - Checks minimap2 before TE/CDS annotation starts.
#   - Supports --skip_tandem_dup to skip the slow tandem-duplication test.
#   - Uses FASTA random access instead of loading the whole reference genome.
#   - Supports symbolic DEL records such as REF=N ALT=<DEL> using POS/END/SVLEN.
#   - Supports sequence-resolved INS/DEL by trimming shared REF/ALT bases.
#   - Supports mixed INS/DEL VCF with --svtype AUTO.
#   - Parses PAF directly; no large external sort step is required.
#
# Coordinates:
#   - VCF POS and INFO/END are treated as 1-based inclusive coordinates.
#   - FASTA extraction internally uses 0-based half-open intervals.
################################################################################

my ($ref, $vcf, $te, $cds, $outfile, $svtype, $workdir, $help);
my ($minimap2, $threads, $keep_tmp, $quiet);
my ($kmer, $dup_seed_fraction, $skip_tandem_dup, $min_annot_pct, $min_aln_len, $min_mapq);
my ($max_dup_flank, $max_dup_seeds, $max_seq_len, $te_preset, $cds_preset);

GetOptions(
    "ref:s"               => \$ref,
    "vcf:s"               => \$vcf,
    "te:s"                => \$te,
    "cds:s"               => \$cds,
    "out:s"               => \$outfile,
    "svtype:s"            => \$svtype,
    "workdir:s"           => \$workdir,
    "threads|t:i"         => \$threads,
    "minimap2:s"          => \$minimap2,
    "kmer:i"              => \$kmer,
    "dup_seed_fraction:f" => \$dup_seed_fraction,
    "skip_tandem_dup"     => \$skip_tandem_dup,
    "min_annot_pct:f"     => \$min_annot_pct,
    "min_aln_len:i"       => \$min_aln_len,
    "min_mapq:i"          => \$min_mapq,
    "max_dup_flank:i"     => \$max_dup_flank,
    "max_dup_seeds:i"     => \$max_dup_seeds,
    "max_seq_len:i"       => \$max_seq_len,
    "te_preset:s"         => \$te_preset,
    "cds_preset:s"        => \$cds_preset,
    "keep_tmp"            => \$keep_tmp,
    "quiet"               => \$quiet,
    "help"                => \$help,
) or die usage();

if ($help || !defined $ref || !defined $vcf || !defined $outfile) {
    print usage();
    exit($help ? 0 : 1);
}

$threads           ||= 4;
$minimap2          ||= "minimap2";
$svtype            ||= "AUTO";
$svtype = uc($svtype);
$workdir           ||= "$outfile.work";
$kmer              ||= 10;
$dup_seed_fraction = 0.50 unless defined $dup_seed_fraction;
$skip_tandem_dup   ||= 0;
$min_annot_pct     = 1.0 unless defined $min_annot_pct;
$min_aln_len       ||= 20;
$min_mapq          = 0 unless defined $min_mapq;
$max_dup_flank     ||= 1000000;  # maximum flank length used for tandem-dup test
$max_dup_seeds     = 50000 unless defined $max_dup_seeds; # 0 means use all SV k-mers
$max_seq_len       ||= 0;        # 0 means no cap for writing SV sequence FASTA
$te_preset         ||= "";
$cds_preset        ||= "asm5";

validate_options();

my $need_minimap2 = 0;
$need_minimap2 = 1 if defined $te  && $te  ne '';
$need_minimap2 = 1 if defined $cds && $cds ne '';
if ($need_minimap2) {
    $minimap2 = check_executable($minimap2);
    logmsg("Using minimap2: $minimap2");
}

make_path($workdir) unless -d $workdir;

logmsg("Loading reference FASTA index: $ref");
my $refdb = FastFAI->new($ref);

my $sv_fa  = "$workdir/sv_sequences.fa";
my $meta_t = "$workdir/sv_metadata.tsv";
my $td_t   = "$workdir/tandem_duplication.tsv";

if ($skip_tandem_dup) {
    logmsg("Step 1. Extracting SV sequences; tandem duplication typing is skipped");
} else {
    logmsg("Step 1. Extracting SV sequences and typing tandem duplications");
}
my ($n_records, $n_with_seq) = prepare_sv_sequences($vcf, $sv_fa, $meta_t, $td_t, $refdb);
logmsg("  VCF records used: $n_records; records with SV sequence: $n_with_seq");

my (%te_ann, %cds_ann);
if (defined $te && -s $te && $n_with_seq > 0) {
    my $paf = "$workdir/TE.paf";
    logmsg("Step 2. Aligning SV sequences to TE library");
    run_minimap2($te, $sv_fa, $paf, $te_preset);
    logmsg("Step 3. Parsing TE alignments");
    %te_ann = parse_paf_annotation($paf, "TE");
} else {
    logmsg("Step 2. Skipping TE annotation because --te is not provided, file is empty, or no SV sequence is available");
}

if (defined $cds && -s $cds && $n_with_seq > 0) {
    my $paf = "$workdir/CDS.paf";
    logmsg("Step 4. Aligning SV sequences to CDS FASTA");
    run_minimap2($cds, $sv_fa, $paf, $cds_preset);
    logmsg("Step 5. Parsing CDS alignments");
    %cds_ann = parse_paf_annotation($paf, "CDS");
} else {
    logmsg("Step 4. Skipping CDS annotation because --cds is not provided, file is empty, or no SV sequence is available");
}

logmsg("Step 6. Writing annotated VCF");
write_annotated_vcf($vcf, $outfile, \%te_ann, \%cds_ann, $td_t);

unless ($keep_tmp) {
    # Keep sv_metadata.tsv and tandem_duplication.tsv for traceability.
    unlink $sv_fa if -e $sv_fa;
    unlink "$workdir/TE.paf"  if -e "$workdir/TE.paf";
    unlink "$workdir/CDS.paf" if -e "$workdir/CDS.paf";
}

logmsg("All jobs done. Output: $outfile");
exit;

sub usage {
    return <<"END";
Usage:
  perl 10_SVannotation_final.pl \
    --ref ref.fa \
    --vcf SV.vcf \
    --te TE.lib.fa \
    --cds CDS.fa \
    --out SV.annotated.vcf \
    --threads 8

Required:
  --ref <fa>                 reference genome FASTA; uncompressed FASTA required
  --vcf <vcf>                input SV VCF file
  --out <vcf>                output annotated VCF file

Optional annotation inputs:
  --te <fa>                  TE library FASTA, for TE annotation by minimap2
  --cds <fa>                 CDS FASTA, for CDS-related annotation by minimap2

Optional parameters:
  --svtype <AUTO|INS|DEL>    SV type. AUTO reads INFO/SVTYPE or infers from REF/ALT [AUTO]
  --workdir <dir>            intermediate directory [<out>.work]
  --threads <int>            minimap2 threads [4]
  --minimap2 <path>          minimap2 executable or full path [minimap2]

Tandem-duplication parameters:
  --skip_tandem_dup          skip tandem-dup typing; TANDUP=SKIP_TANDUP
  --kmer <int>               seed length for tandem-dup typing [10]
  --dup_seed_fraction <f>    minimum matched seed fraction for DUP typing [0.50]
  --max_dup_flank <int>      maximum flank length for DUP typing [1000000]
  --max_dup_seeds <int>      maximum SV k-mers tested per record; 0 uses all [50000]

TE/CDS minimap2 and PAF parsing parameters:
  --te_preset <str>          minimap2 preset for TE alignment [none]
  --cds_preset <str>         minimap2 preset for CDS alignment [asm5]
  --min_annot_pct <f>        minimum query coverage percentage to report TE/CDS [1.0]
  --min_aln_len <int>        minimum PAF alignment block length [20]
  --min_mapq <int>           minimum PAF mapping quality [0]

Other parameters:
  --max_seq_len <int>        maximum SV sequence length written to FASTA; 0 means no cap [0]
  --keep_tmp                 keep intermediate FASTA and PAF files
  --quiet                    suppress progress messages
  --help                     print this help message

Notes:
  1. For symbolic DEL records, for example REF=N ALT=<DEL>, the deleted sequence
     is extracted from the reference genome using POS and INFO/END. If END is
     absent, abs(SVLEN) is used when available.
  2. For sequence-resolved DEL/INS records, shared prefix and suffix are trimmed
     before extracting the SV allele sequence.
  3. For symbolic INS records, sequence annotation requires an inserted sequence
     in INFO fields such as SEQ, SVSEQ, INSSEQ, ALTSEQ, INSERTION_SEQUENCE or
     INSERTED_SEQUENCE.
END
}

sub validate_options {
    die "ERROR: --ref file does not exist or is empty: $ref\n" unless -s $ref;
    die "ERROR: --vcf file does not exist or is empty: $vcf\n" unless -s $vcf;
    die "ERROR: --threads must be >= 1\n" unless defined $threads && $threads >= 1;
    die "ERROR: --svtype must be AUTO, INS or DEL\n" unless $svtype =~ /^(AUTO|INS|DEL)$/;
    die "ERROR: --kmer must be >= 1\n" unless defined $kmer && $kmer >= 1;
    die "ERROR: --dup_seed_fraction must be > 0 and <= 1\n" unless $dup_seed_fraction > 0 && $dup_seed_fraction <= 1;
    die "ERROR: --max_dup_flank must be >= 1\n" unless defined $max_dup_flank && $max_dup_flank >= 1;
    die "ERROR: --max_dup_seeds must be >= 0\n" unless defined $max_dup_seeds && $max_dup_seeds >= 0;
    die "ERROR: --min_annot_pct must be >= 0\n" unless defined $min_annot_pct && $min_annot_pct >= 0;
    die "ERROR: --min_aln_len must be >= 0\n" unless defined $min_aln_len && $min_aln_len >= 0;
    die "ERROR: --min_mapq must be >= 0\n" unless defined $min_mapq && $min_mapq >= 0;
    die "ERROR: --max_seq_len must be >= 0\n" unless defined $max_seq_len && $max_seq_len >= 0;
    die "ERROR: --te file does not exist or is empty: $te\n" if defined $te && $te ne '' && !-s $te;
    die "ERROR: --cds file does not exist or is empty: $cds\n" if defined $cds && $cds ne '' && !-s $cds;
}

sub check_executable {
    my ($prog) = @_;
    die "ERROR: minimap2 executable is not defined.\n" unless defined $prog && $prog ne '';

    if ($prog =~ m{/}) {
        die "ERROR: minimap2 is not executable at '$prog'. Please check --minimap2.\n" unless -x $prog;
        return $prog;
    }

    for my $dir (split /:/, $ENV{PATH} // '') {
        next if $dir eq '';
        my $candidate = "$dir/$prog";
        return $candidate if -x $candidate;
    }

    die "ERROR: cannot find '$prog' in PATH. Please install minimap2 or specify its full path with --minimap2 /path/to/minimap2.\n";
}

sub prepare_sv_sequences {
    my ($vcf_file, $fasta_out, $meta_out, $td_out, $db) = @_;
    open(my $IN,   "<", $vcf_file)  or die "ERROR: cannot read VCF $vcf_file: $!\n";
    open(my $FA,   ">", $fasta_out) or die "ERROR: cannot write $fasta_out: $!\n";
    open(my $META, ">", $meta_out)  or die "ERROR: cannot write $meta_out: $!\n";
    open(my $TD,   ">", $td_out)    or die "ERROR: cannot write $td_out: $!\n";

    print $META join("\t", qw(QID CHROM POS ID SVTYPE SOURCE SV_START0 SV_END0 SVLEN HAS_SEQUENCE)), "\n";
    print $TD   join("\t", qw(QID TANDUP)), "\n";

    my ($idx_all, $idx_used, $with_seq) = (0, 0, 0);
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^#/;
        next if $line =~ /^\s*$/;
        $idx_all++;
        my @f = split /\t/, $line;
        next unless @f >= 8;
        my ($chrom, $pos, $id, $ref_allele, $alt_allele, $info) = @f[0,1,2,3,4,7];
        my %info = parse_info($info);
        my $type = infer_svtype($ref_allele, $alt_allele, \%info);
        next if $svtype ne 'AUTO' && $type ne $svtype;
        next unless $type eq 'INS' || $type eq 'DEL';

        $idx_used++;
        my $qid = qid_for_record($idx_all, $chrom, $pos, $id);
        my ($seq, $source, $s0, $e0) = sv_sequence_and_interval($db, $chrom, $pos, $ref_allele, $alt_allele, \%info, $type);
        my $has_seq = (defined $seq && $seq ne '' && $seq ne '.');
        my $svlen = $has_seq ? length($seq) : estimate_svlen($pos, $ref_allele, $alt_allele, \%info, $type);

        my $td = 'NOSEQ';
        if ($has_seq) {
            $with_seq++;
            my $seq_to_write = clean_seq($seq);
            if ($max_seq_len && length($seq_to_write) > $max_seq_len) {
                $seq_to_write = substr($seq_to_write, 0, $max_seq_len);
            }
            print $FA ">$qid\n";
            print_wrapped_fasta($FA, $seq_to_write, 80);

            if ($skip_tandem_dup) {
                $td = 'SKIP_TANDUP';
            } else {
                $td = tandem_dup_type($db, $chrom, $type, $seq, $s0, $e0);
            }
        } elsif ($skip_tandem_dup) {
            $td = 'SKIP_TANDUP';
        }

        print $TD "$qid\t$td\n";
        print $META join("\t", $qid, $chrom, $pos, $id, $type, ($source // 'NA'), ($s0 // 'NA'), ($e0 // 'NA'), ($svlen // 'NA'), ($has_seq ? 1 : 0)), "\n";
    }
    close $IN;
    close $FA;
    close $META;
    close $TD;
    return ($idx_used, $with_seq);
}

sub infer_svtype {
    my ($ref_allele, $alt_allele, $info) = @_;
    if (exists $info->{SVTYPE} && $info->{SVTYPE} ne '') {
        my $t = uc($info->{SVTYPE});
        return 'INS' if $t =~ /INS/;
        return 'DEL' if $t =~ /DEL/;
        return $t;
    }
    if (defined $alt_allele && $alt_allele =~ /^<([^>]+)>/) {
        my $t = uc($1);
        return 'INS' if $t =~ /INS/;
        return 'DEL' if $t =~ /DEL/;
        return $t;
    }
    my ($r, $a) = trim_shared_alleles($ref_allele, first_alt($alt_allele));
    return 'DEL' if length($r) > length($a);
    return 'INS' if length($a) > length($r);
    return 'UNK';
}

sub sv_sequence_and_interval {
    my ($db, $chrom, $pos, $ref_allele, $alt_allele, $info, $type) = @_;
    my $alt1 = first_alt($alt_allele);
    my $symbolic_alt = (defined $alt1 && $alt1 =~ /^</) ? 1 : 0;
    my $pos0 = $pos - 1;

    if ($type eq 'DEL') {
        if ($symbolic_alt || !defined($ref_allele) || $ref_allele eq 'N' || length($ref_allele) <= 1) {
            my $end1 = get_end1($pos, $ref_allele, $alt1, $info, $type);
            return (undef, 'NO_END_FOR_SYMBOLIC_DEL', $pos0, undef) unless defined $end1 && $end1 >= $pos;
            my $seq = $db->get_seq($chrom, $pos0, $end1);  # POS..END inclusive in VCF -> [POS-1, END)
            return (clean_seq($seq), 'REF_BY_POS_END', $pos0, $end1);
        } else {
            my ($r, $a, $prefix_len) = trim_shared_alleles($ref_allele, $alt1);
            my $start0 = $pos0 + $prefix_len;
            my $end0   = $start0 + length($r);
            if (length($r) == 0 && exists $info->{END}) {
                my $end1 = get_end1($pos, $ref_allele, $alt1, $info, $type);
                my $seq = $db->get_seq($chrom, $pos0, $end1);
                return (clean_seq($seq), 'REF_BY_POS_END_FALLBACK', $pos0, $end1);
            }
            return (clean_seq($r), 'VCF_REF_TRIMMED', $start0, $end0);
        }
    } elsif ($type eq 'INS') {
        if ($symbolic_alt) {
            for my $key (qw(SEQ SVSEQ INSSEQ ALTSEQ INSERTION_SEQUENCE INSERTED_SEQUENCE)) {
                if (exists $info->{$key} && $info->{$key} ne '') {
                    my $seq = clean_seq($info->{$key});
                    return ($seq, "INFO_$key", $pos, $pos); # inserted after POS
                }
            }
            return (undef, 'NO_SEQUENCE_FOR_SYMBOLIC_INS', $pos, $pos);
        } else {
            my ($r, $a, $prefix_len) = trim_shared_alleles($ref_allele, $alt1);
            my $ins_site0 = $pos0 + $prefix_len;
            return (clean_seq($a), 'VCF_ALT_TRIMMED', $ins_site0, $ins_site0);
        }
    }
    return (undef, 'UNSUPPORTED_SVTYPE', undef, undef);
}

sub get_end1 {
    my ($pos, $ref_allele, $alt_allele, $info, $type) = @_;
    if (exists $info->{END} && $info->{END} =~ /^-?\d+$/) {
        return $info->{END} + 0;
    }
    if (exists $info->{SVLEN} && $info->{SVLEN} =~ /^-?\d+/) {
        my $len = abs($info->{SVLEN});
        return $pos + $len - 1 if $len > 0;
    }
    if (defined $ref_allele && $ref_allele ne '' && $type eq 'DEL') {
        return $pos + length($ref_allele) - 1;
    }
    return undef;
}

sub estimate_svlen {
    my ($pos, $ref_allele, $alt_allele, $info, $type) = @_;
    if (exists $info->{SVLEN} && $info->{SVLEN} =~ /^-?\d+/) {
        return abs($info->{SVLEN});
    }
    if ($type eq 'DEL') {
        my $end1 = get_end1($pos, $ref_allele, first_alt($alt_allele), $info, $type);
        return $end1 - $pos + 1 if defined $end1;
    }
    my ($r, $a) = trim_shared_alleles($ref_allele, first_alt($alt_allele));
    return $type eq 'INS' ? length($a) : length($r);
}

sub trim_shared_alleles {
    my ($ref_allele, $alt_allele) = @_;
    $ref_allele = '' unless defined $ref_allele;
    $alt_allele = '' unless defined $alt_allele;
    return ($ref_allele, $alt_allele, 0) if $alt_allele =~ /^</;

    my ($r, $a) = ($ref_allele, $alt_allele);
    my $prefix = 0;
    while (length($r) > 0 && length($a) > 0 && substr($r, 0, 1) eq substr($a, 0, 1)) {
        substr($r, 0, 1) = '';
        substr($a, 0, 1) = '';
        $prefix++;
    }
    while (length($r) > 0 && length($a) > 0 && substr($r, -1, 1) eq substr($a, -1, 1)) {
        chop $r;
        chop $a;
    }
    return ($r, $a, $prefix);
}

sub tandem_dup_type {
    my ($db, $chrom, $type, $svseq, $sv_start0, $sv_end0) = @_;
    $svseq = clean_seq($svseq);
    return 'NOSEQ' unless defined $svseq && $svseq ne '' && $svseq ne '.';
    my $len = length($svseq);
    return 'SHORTSEQ' if $len < 2;

    my $flank_len = $len;
    $flank_len = $max_dup_flank if $flank_len > $max_dup_flank;

    my ($left_anchor, $right_anchor);
    if ($type eq 'DEL') {
        $left_anchor  = defined $sv_start0 ? $sv_start0 : 0;
        $right_anchor = defined $sv_end0   ? $sv_end0   : $left_anchor;
    } elsif ($type eq 'INS') {
        $left_anchor  = defined $sv_start0 ? $sv_start0 : 0;
        $right_anchor = $left_anchor;
    } else {
        return 'NONDUP';
    }

    my $seq5 = $db->get_seq($chrom, $left_anchor - $flank_len, $left_anchor);
    my $seq3 = $db->get_seq($chrom, $right_anchor, $right_anchor + $flank_len);
    $seq5 = '' if !defined($seq5) || $seq5 eq '.';
    $seq3 = '' if !defined($seq3) || $seq3 eq '.';
    $seq5 = clean_seq($seq5);
    $seq3 = clean_seq($seq3);

    my ($iden,    $orient)    = seed_identity_fast($svseq, $seq5, $seq3);
    my ($idenrev, $orientrev) = seed_identity_fast(revcom($svseq), $seq5, $seq3);

    if ($iden) {
        return $orient == 5 ? '5DUP' : $orient == 3 ? '3DUP' : 'DUP';
    } elsif ($idenrev) {
        return $orientrev == 5 ? '5INVDUP' : $orientrev == 3 ? '3INVDUP' : 'INVDUP';
    }
    return 'NONDUP';
}

sub seed_identity_fast {
    my ($seq, $flank5, $flank3) = @_;
    $seq    = clean_seq($seq);
    $flank5 = clean_seq($flank5);
    $flank3 = clean_seq($flank3);
    return (0, 0) if $seq eq '';

    my $len = length($seq);
    my $k = $kmer;
    $k = $len if $len < $k;
    return (0, 0) if $k <= 0;

    my $num = $len - $k + 1;
    return (0, 0) if $num <= 0;

    my %k5 = build_kmer_set($flank5, $k);
    my %k3 = build_kmer_set($flank3, $k);

    my @positions = seed_positions($num, $max_dup_seeds);
    my ($m5, $m3, $tested) = (0, 0, 0);
    for my $i (@positions) {
        my $seed = substr($seq, $i, $k);
        next if $seed =~ /N/i;
        $tested++;
        $m5++ if exists $k5{$seed};
        $m3++ if exists $k3{$seed};
    }
    return (0, 0) if $tested <= 0;

    my $p5 = $m5 / $tested;
    my $p3 = $m3 / $tested;
    if ($p5 >= $dup_seed_fraction || $p3 >= $dup_seed_fraction) {
        return $p5 >= $p3 ? (1, 5) : (1, 3);
    }
    return (0, 0);
}

sub build_kmer_set {
    my ($seq, $k) = @_;
    my %set;
    return %set unless defined $seq && length($seq) >= $k;
    my $n = length($seq) - $k + 1;
    for (my $i = 0; $i < $n; $i++) {
        my $kmer_seq = substr($seq, $i, $k);
        next if $kmer_seq =~ /N/i;
        $set{$kmer_seq} = 1;
    }
    return %set;
}

sub seed_positions {
    my ($num, $cap) = @_;
    return (0 .. $num - 1) if !defined($cap) || $cap <= 0 || $num <= $cap;

    my @pos;
    my $step = $num / $cap;
    my %seen;
    for (my $j = 0; $j < $cap; $j++) {
        my $p = int($j * $step);
        $p = $num - 1 if $p >= $num;
        next if $seen{$p}++;
        push @pos, $p;
    }
    return @pos;
}

sub run_minimap2 {
    my ($target_fa, $query_fa, $paf_out, $preset) = @_;
    die "ERROR: minimap2 executable is not available.\n" unless defined $minimap2 && -x $minimap2;
    open(my $OUT, ">", $paf_out) or die "ERROR: cannot write $paf_out: $!\n";
    my @cmd = ($minimap2, "-t", $threads, "-c");
    push @cmd, ("-x", $preset) if defined $preset && $preset ne '';
    push @cmd, ($target_fa, $query_fa);
    logmsg("  CMD: @cmd");
    open(my $MM, "-|", @cmd) or die "ERROR: cannot run minimap2. Command: @cmd\n";
    while (my $line = <$MM>) {
        print $OUT $line;
    }
    close $MM;
    my $status = $?;
    close $OUT;
    if ($status != 0) {
        my $exit_code = $status >> 8;
        my $signal    = $status & 127;
        die "ERROR: minimap2 failed. exit_code=$exit_code signal=$signal. Command: @cmd\n";
    }
}

sub parse_paf_annotation {
    my ($paf, $kind) = @_;
    my %hits;
    open(my $IN, "<", $paf) or die "ERROR: cannot read $paf: $!\n";
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my @f = split /\t/, $line;
        next unless @f >= 12;
        my ($qid, $qlen, $qs, $qe, $strand, $target, $tlen, $ts, $te, $matches, $alnlen, $mapq) = @f[0..11];
        next if defined $min_mapq && $mapq < $min_mapq;
        next if defined $min_aln_len && $alnlen < $min_aln_len;
        next if $qlen <= 0 || $qe <= $qs;
        my $label = annotation_label($target, $kind);
        push @{ $hits{$qid}{$label}{intervals} }, [$qs + 0, $qe + 0];
        $hits{$qid}{$label}{qlen} = $qlen + 0;
    }
    close $IN;

    my %ann;
    for my $qid (keys %hits) {
        my @items;
        for my $label (keys %{ $hits{$qid} }) {
            my $qlen = $hits{$qid}{$label}{qlen} || 1;
            my $cov = merged_interval_length($hits{$qid}{$label}{intervals});
            my $pct = 100 * $cov / $qlen;
            next if $pct < $min_annot_pct;
            my $pct_s = sprintf("%.1f", $pct);
            my $safe_label = sanitize_info_value($label);
            push @items, [$pct, "$safe_label:$pct_s%"];
        }
        if (@items) {
            @items = sort { $b->[0] <=> $a->[0] || $a->[1] cmp $b->[1] } @items;
            $ann{$qid} = "(" . join(";", map { $_->[1] } @items) . ")";
        }
    }
    return %ann;
}

sub annotation_label {
    my ($target, $kind) = @_;
    $target //= 'Unknown';
    return $target;
}

sub merged_interval_length {
    my ($ivals) = @_;
    return 0 unless defined $ivals && @$ivals;
    my @s = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$ivals;
    my ($cs, $ce) = @{ shift @s };
    my $total = 0;
    for my $iv (@s) {
        my ($s0, $e0) = @$iv;
        if ($s0 <= $ce) {
            $ce = $e0 if $e0 > $ce;
        } else {
            $total += $ce - $cs;
            ($cs, $ce) = ($s0, $e0);
        }
    }
    $total += $ce - $cs;
    return $total;
}

sub write_annotated_vcf {
    my ($vcf_file, $out_vcf, $te_hash, $cds_hash, $td_file) = @_;
    my %td;
    open(my $TD, "<", $td_file) or die "ERROR: cannot read $td_file: $!\n";
    while (my $line = <$TD>) {
        chomp $line;
        next if $line =~ /^QID\t/;
        my ($qid, $val) = split /\t/, $line;
        $td{$qid} = $val if defined $qid;
    }
    close $TD;

    open(my $IN,  "<", $vcf_file) or die "ERROR: cannot read VCF $vcf_file: $!\n";
    open(my $OUT, ">", $out_vcf)  or die "ERROR: cannot write $out_vcf: $!\n";
    my $idx_all = 0;
    while (my $line = <$IN>) {
        chomp $line;
        if ($line =~ /^##/) {
            print $OUT "$line\n";
            next;
        }
        if ($line =~ /^#CHROM/) {
            print $OUT "##INFO=<ID=TANDUP,Number=1,Type=String,Description=\"Tandem duplication annotation based on SV sequence similarity to 5-prime or 3-prime flanking reference sequences; SKIP_TANDUP indicates --skip_tandem_dup was used\">\n";
            print $OUT "##INFO=<ID=CDSANN,Number=1,Type=String,Description=\"CDS sequence annotation based on minimap2 alignment of the SV sequence to the CDS FASTA\">\n";
            print $OUT "##INFO=<ID=TEANN,Number=1,Type=String,Description=\"TE sequence annotation based on minimap2 alignment of the SV sequence to the TE library\">\n";
            print $OUT "$line\n";
            next;
        }
        if ($line =~ /^#/) {
            print $OUT "$line\n";
            next;
        }
        $idx_all++;
        my @f = split /\t/, $line;
        if (@f < 8) {
            print $OUT "$line\n";
            next;
        }
        my %info = parse_info($f[7]);
        my $type = infer_svtype($f[3], $f[4], \%info);
        if ($svtype ne 'AUTO' && $type ne $svtype) {
            print $OUT "$line\n";
            next;
        }
        unless ($type eq 'INS' || $type eq 'DEL') {
            print $OUT "$line\n";
            next;
        }

        my $qid = qid_for_record($idx_all, $f[0], $f[1], $f[2]);
        my $tdv  = exists $td{$qid}        ? $td{$qid}        : 'NF';
        my $cdsv = exists $cds_hash->{$qid} ? $cds_hash->{$qid} : 'NF';
        my $tev  = exists $te_hash->{$qid}  ? $te_hash->{$qid}  : 'NF';
        $f[7] = add_info_fields($f[7], TANDUP => $tdv, CDSANN => $cdsv, TEANN => $tev);
        print $OUT join("\t", @f), "\n";
    }
    close $IN;
    close $OUT;
}

sub add_info_fields {
    my ($info, %add) = @_;
    $info = '' if !defined($info) || $info eq '.';
    for my $k (sort keys %add) {
        my $v = sanitize_info_value($add{$k});
        $info .= ";" if $info ne '';
        $info .= "$k=$v";
    }
    return $info eq '' ? '.' : $info;
}

sub parse_info {
    my ($info) = @_;
    my %h;
    return %h if !defined $info || $info eq '.';
    for my $x (split /;/, $info) {
        next if $x eq '';
        my ($k, $v) = split /=/, $x, 2;
        $v = 1 unless defined $v;
        $h{$k} = $v;
    }
    return %h;
}

sub first_alt {
    my ($alt) = @_;
    return '' unless defined $alt;
    my ($a) = split /,/, $alt;
    return $a;
}

sub qid_for_record {
    my ($idx, $chrom, $pos, $id) = @_;
    $id = 'NA' if !defined($id) || $id eq '.' || $id eq '';
    my $qid = sprintf("SV%09d_%s_%s_%s", $idx, $chrom, $pos, $id);
    $qid =~ s/[^A-Za-z0-9_.:-]/_/g;
    return $qid;
}

sub sanitize_info_value {
    my ($v) = @_;
    $v = 'NF' unless defined $v && $v ne '';
    $v =~ s/[\s,]+/_/g;
    $v =~ s/;/,/g;
    $v =~ s/=/:/g;
    return $v;
}

sub clean_seq {
    my ($seq) = @_;
    return '' unless defined $seq;
    $seq = uc($seq);
    $seq =~ s/[^ACGTN]//g;
    return $seq;
}

sub print_wrapped_fasta {
    my ($fh, $seq, $width) = @_;
    $seq //= '';
    for (my $i = 0; $i < length($seq); $i += $width) {
        print $fh substr($seq, $i, $width), "\n";
    }
}

sub revcom {
    my ($seq) = @_;
    return '' unless defined $seq;
    my $rev = reverse $seq;
    $rev =~ tr/ACGTNacgtn/TGCANtgcan/;
    return $rev;
}

sub logmsg {
    my ($msg) = @_;
    return if $quiet;
    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
    my $t = sprintf('%04d-%02d-%02d %02d:%02d:%02d', $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
    print STDERR "[$t] $msg\n";
}

################################################################################
# Lightweight FASTA random access. It can read an existing .fai file, or build an
# in-memory index by scanning the FASTA once. FASTA must be uncompressed.
################################################################################
package FastFAI;
use strict;
use warnings;
use Fcntl qw(SEEK_SET);

sub new {
    my ($class, $fasta) = @_;
    die "ERROR: FASTA path is not defined.\n" unless defined $fasta;
    die "ERROR: FASTA does not exist or is empty: $fasta\n" unless -s $fasta;
    die "ERROR: this script requires uncompressed FASTA, but got: $fasta\n" if $fasta =~ /\.gz$/;

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

    my ($left_pad, $right_pad) = (0, 0);
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
    $seq =~ s/\s+//g;
    return uc($seq);
}

sub _resolve_chr {
    my ($self, $chr) = @_;
    return $self->{cache}{$chr} if exists $self->{cache}{$chr};
    if (exists $self->{idx}{$chr}) {
        return $self->{cache}{$chr} = $chr;
    }

    my @try;
    if ($chr =~ /^[^.]+\.(.+)$/) { push @try, $1; }
    if ($chr =~ /(chr\S+)$/i)    { push @try, $1; }
    if ($chr =~ /^chr0([1-9])$/i){ push @try, "chr$1"; push @try, "Chr$1"; }
    if ($chr =~ /^chr([1-9])$/i) { push @try, "chr0$1"; push @try, "Chr0$1"; }

    for my $c (@try) {
        if (exists $self->{idx}{$c}) {
            return $self->{cache}{$chr} = $c;
        }
    }

    for my $c (keys %{ $self->{idx} }) {
        if ($c =~ /\Q$chr\E$/ || $chr =~ /\Q$c\E$/) {
            return $self->{cache}{$chr} = $c;
        }
    }

    warn "[WARN] chromosome '$chr' not found in $self->{fasta}; sequence will be '.'.\n"
        unless $self->{warned}{$chr}++;
    $self->{cache}{$chr} = undef;
    return undef;
}

sub _read_interval {
    my ($self, $idx, $start, $end) = @_;
    return '' if $end <= $start;
    my $lb  = $idx->{line_bases};
    my $lw  = $idx->{line_width};
    my $off = $idx->{offset};
    my $fh  = $self->{fh};

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
