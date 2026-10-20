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
use File::Copy qw(copy);
use File::Basename qw(basename);
use Cwd qw(abs_path);
use POSIX qw(strftime);

# ============================================================================
# SVGAP standardized script: SVGAP_pipeline.pl
# Original author: Yi Liao
# Emails: yiliao@scau.edu.cn; yiliao1022@gmail.com
# Standardized date: 2026-06-14
# Update: v11 fixes qnet PairwiseSV input by copying synnet files instead of
#         creating relative symlinks that can break in nested output folders.
# Update: v12 adds 12_combine_genotype_CXSV.pl for bidirectional complex SV
#         merging and population genotyping.
# Update: v13 repairs skip-filter and empty-call paths, adds adaptive SynNet
#         thresholds, and reduces conversion/SNP overhead.
#
# Purpose:
#   One-command SVGAP workflow from one reference assembly and multiple query
#   assemblies to population SV/SNP/InDel VCF files.
# Release: SVGAP version2 (the previous repository workflow is version1).
#
# Major design in this version:
#   1. The pipeline can start from:
#        --start_from align : pairwise aligner outputs, converted by 01_Convert2Axt.pl
#        --start_from maf   : existing pairwise MAF files, converted to AXT
#        --start_from axt   : existing pairwise AXT files
#   2. A true independent reciprocal alignment is optional, but NOT required.
#      If only forward Ref-vs-Query alignment is available, SVGAP can use the
#      qnet generated from forward chainNet to recover query-side gaps.
#   3. If true Query-vs-Ref alignments are available, the pipeline can keep the
#      old reciprocal workflow and run chain/net/synnet on those independent
#      reciprocal AXT/MAF files.  This is controlled by --reciprocal_strategy.
#   4. The pipeline never generates fake reciprocal AXT by reversing forward AXT.
#      Such files are coordinate-reordered alignments, not independent
#      Query-vs-Ref alignments, and can be misleading for SV discovery.
#   5. VCF sample order is strictly controlled by g.list after removing refname.
# ============================================================================

my %opt = (
    refname => '',
    glist => '',
    genome_dir => '',

    start_from => 'align',          # align | maf | axt
    align_dir => '',
    aligner => 'minimap2',
    align_mode => 'auto',           # auto | chr | ctg
    align_sizes => [],              # original chromosome sizes for split-contig input
    maf_dir => '',
    maf_reciprocal_dir => '',       # optional, true independently generated Query-vs-Ref MAFs
    axt_dir => '',
    axt_forward_dir => '',
    axt_reciprocal_dir => '',       # optional, true independently generated QueryvsRef AXT files
    reciprocal_strategy => 'auto',  # auto | qnet | true | both

    outdir => 'SVGAP_run',
    threads => 8,
    speciesname => 'Unknown',
    ucsc => "$Bin/pub/UCSC/",
    tba => "$Bin/pub/TBA/",
    k8 => "$Bin/pub/MINIMAP2/",
    paftools => "$Bin/pub/MINIMAP2/",
    mum320 => "$Bin/pub/MUMMER/",
    linearGap => "$Bin/pub/TBA/medium",
    refsize => '',
    te => '',
    cds => '',
    minimap2 => 'minimap2',
    annotation => 1,

    # SynNet filtering. The fixed values remain upper bounds; adaptive mode
    # scales them down for short chromosomes/contigs instead of deleting every
    # top-level fill when a sequence is shorter than 500 kb.
    filter_toplen => 500_000,
    filter_nonsynlen => 500_000,
    filter_synlen => 500_000,
    # A 10 kb floor misses many well-aligned assembly inversions.  Keep 1 kb as
    # the conservative default and pass the same floor to both filtering and
    # PairwiseSV so the two stages cannot silently disagree.
    filter_invlen => 1_000,
    filter_alrate => 0.10,
    filter_qfar_cutoff => 5_000_000,
    filter_adaptive => 1,
    filter_auto_fraction => 0.10,
    filter_auto_min => 1_000,

    skip_convert => 0,
    skip_chainnet => 0,
    skip_filter => 0,
    skip_pairwise_sv => 0,
    skip_combine => 0,
    skip_genotype => 0,
    skip_cxsv => 0,
    skip_snp => 0,
    skip_annotation => 0,
    skip_vg_vcf => 0,
    max_vg_allele_len => 5_000_000,
    index_vg_vcf => 1,
    bgzip => 'bgzip',
    tabix => 'tabix',

    cxsv_break_tolerance => 2000,
    cxsv_min_ro => 0.50,
    cxsv_len_ratio => 1.00,
    cxsv_unknown_missing => 0,
    dry_run => 0,
    force => 0,
    version => 0,
    help => 0,
);

GetOptions(
    'refname|target=s' => \$opt{refname},
    'glist=s' => \$opt{glist},
    'genome_dir|gd=s' => \$opt{genome_dir},

    'start_from=s' => \$opt{start_from},
    'align_dir|input=s' => \$opt{align_dir},
    'aligner|ali=s' => \$opt{aligner},
    'align_mode|align-mode=s' => \$opt{align_mode},
    'align_sizes|align-sizes=s@' => \$opt{align_sizes},
    'maf_dir=s' => \$opt{maf_dir},
    'maf_reciprocal_dir=s' => \$opt{maf_reciprocal_dir},
    'axt_dir|ad=s' => \$opt{axt_dir},
    'axt_forward_dir=s' => \$opt{axt_forward_dir},
    'axt_reciprocal_dir=s' => \$opt{axt_reciprocal_dir},
    'reciprocal_strategy=s' => \$opt{reciprocal_strategy},

    'outdir=s' => \$opt{outdir},
    'threads|t=i' => \$opt{threads},
    'speciesname=s' => \$opt{speciesname},
    'ucsc=s' => \$opt{ucsc},
    'tba=s' => \$opt{tba},
    'k8=s' => \$opt{k8},
    'paftools=s' => \$opt{paftools},
    'mum320=s' => \$opt{mum320},
    'linearGap|lg=s' => \$opt{linearGap},
    'refsize=s' => \$opt{refsize},
    'te=s' => \$opt{te},
    'cds=s' => \$opt{cds},
    'minimap2=s' => \$opt{minimap2},
    'annotation!' => \$opt{annotation},

    'filter_toplen|filter-toplen=i' => \$opt{filter_toplen},
    'filter_nonsynlen|filter-nonsynlen=i' => \$opt{filter_nonsynlen},
    'filter_synlen|filter-synlen=i' => \$opt{filter_synlen},
    'filter_invlen|filter-invlen=i' => \$opt{filter_invlen},
    'filter_alrate|filter-alrate=f' => \$opt{filter_alrate},
    'filter_qfar_cutoff|filter-qfar-cutoff=i' => \$opt{filter_qfar_cutoff},
    'filter_adaptive|filter-adaptive!' => \$opt{filter_adaptive},
    'filter_auto_fraction|filter-auto-fraction=f' => \$opt{filter_auto_fraction},
    'filter_auto_min|filter-auto-min=i' => \$opt{filter_auto_min},

    'skip_convert!' => \$opt{skip_convert},
    'skip_chainnet!' => \$opt{skip_chainnet},
    'skip_filter!' => \$opt{skip_filter},
    'skip_pairwise_sv!' => \$opt{skip_pairwise_sv},
    'skip_combine!' => \$opt{skip_combine},
    'skip_genotype!' => \$opt{skip_genotype},
    'skip_cxsv!' => \$opt{skip_cxsv},
    'skip_snp!' => \$opt{skip_snp},
    'skip_annotation!' => \$opt{skip_annotation},
    'skip_vg_vcf|skip-vg-vcf!' => \$opt{skip_vg_vcf},
    'max_vg_allele_len|max-vg-allele-len=i' => \$opt{max_vg_allele_len},
    'index_vg_vcf|index-vg-vcf!' => \$opt{index_vg_vcf},
    'bgzip=s' => \$opt{bgzip},
    'tabix=s' => \$opt{tabix},

    'cxsv_break_tolerance|cxsv-break-tolerance=i' => \$opt{cxsv_break_tolerance},
    'cxsv_min_ro|cxsv-min-ro=f' => \$opt{cxsv_min_ro},
    'cxsv_len_ratio|cxsv-len-ratio=f' => \$opt{cxsv_len_ratio},
    'cxsv_unknown_missing|cxsv-unknown-missing!' => \$opt{cxsv_unknown_missing},
    'dry_run!' => \$opt{dry_run},
    'force!' => \$opt{force},
    'version' => \$opt{version},
    'help|h' => \$opt{help},
) or die usage();

$opt{aligner} = lc($opt{aligner} || '');
$opt{align_mode} = lc($opt{align_mode} || 'auto');

if ($opt{version}) {
    print "SVGAP version2\n";
    exit 0;
}

if ($opt{help} || !$opt{refname} || !$opt{glist} || !$opt{genome_dir}) {
    print usage();
    exit($opt{help} ? 0 : 1);
}

validate_options(\%opt);

$opt{outdir} =~ s{/+$}{};
make_path($opt{outdir});
my $log = "$opt{outdir}/SVGAP_pipeline.log";
open(my $LOG, '>', $log) or die "ERROR: cannot write $log: $!\n";
logmsg($LOG, "Started SVGAP pipeline");
logmsg($LOG, "start_from=$opt{start_from}");
logmsg($LOG, "requested reciprocal_strategy=$opt{reciprocal_strategy}");

my ($all_list, $query_list, $queries_ref) = prepare_lists(\%opt, $LOG);
my $ref_fasta = find_genome_fasta($opt{genome_dir}, $opt{refname});
$opt{refsize} ||= find_refsize($opt{genome_dir}, $opt{refname});
if (!$opt{refsize}) {
    $opt{refsize} = "$opt{genome_dir}/$opt{refname}.sizes";
}

my $axt_root  = "$opt{outdir}/01_axt";
my $chain_fwd = "$opt{outdir}/02_chain_forward";
my $chain_rec = "$opt{outdir}/02_chain_independent_reciprocal";
my $filter_fwd = "$opt{outdir}/03_filter_forward";
my $filter_rec = "$opt{outdir}/03_filter_independent_reciprocal";
my $pairsv = "$opt{outdir}/04_pairwise_sv";
my $combined = "$opt{outdir}/05_combined_sv";
my $geno = "$opt{outdir}/06_genotype";
my $cxsv = "$geno/CXSV";
my $snp = "$opt{outdir}/07_snp";
my $vcfdir = "$opt{outdir}/08_vcf";
my $ann = "$opt{outdir}/09_annotation";
make_path($axt_root,$chain_fwd,$chain_rec,$filter_fwd,$filter_rec,$pairsv,$combined,$geno,$cxsv,$snp,$vcfdir,$ann);

my ($fwd_axt, $independent_rec_axt) = prepare_or_resolve_axt(\%opt, $axt_root, $LOG);
my $reciprocal_strategy = choose_reciprocal_strategy(\%opt, $independent_rec_axt, $LOG);

# 2. Build chain/net/synnet for forward AXT. The forward chainNet step creates
#    both tnet and qnet. The qnet synnet is important: after 04_PairwiseSV.pl,
#    RefvsQuery.qnet is written as Query.Ref.rawSV, so reference insertions can
#    be recovered without a fake reciprocal AXT.
if (!$opt{skip_chainnet}) {
    run_cmd($LOG, $opt{dry_run}, "perl '$Bin/02_ChainNetSyn.pl' --gd '$opt{genome_dir}' --ad '$fwd_axt' --lst '$all_list' --outdir '$chain_fwd' --threads $opt{threads} --ucsc '$opt{ucsc}' --tba '$opt{tba}' --linearGap '$opt{linearGap}' --sing");
    if ($independent_rec_axt) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/02_ChainNetSyn.pl' --gd '$opt{genome_dir}' --ad '$independent_rec_axt' --lst '$all_list' --outdir '$chain_rec' --threads $opt{threads} --ucsc '$opt{ucsc}' --tba '$opt{tba}' --linearGap '$opt{linearGap}'");
    }
}

# 3. Filter synnet and create high-confidence target-oriented single-coverage MAF.
if (!$opt{skip_filter}) {
    my $filter_args = "--toplen $opt{filter_toplen} --nonsynlen $opt{filter_nonsynlen} --synlen $opt{filter_synlen} --invlen $opt{filter_invlen} --alrate $opt{filter_alrate} --qFarCutoff $opt{filter_qfar_cutoff} --auto_fraction $opt{filter_auto_fraction} --auto_min $opt{filter_auto_min}";
    $filter_args .= $opt{filter_adaptive} ? " --adaptive_thresholds" : " --no-adaptive_thresholds";
    run_cmd($LOG, $opt{dry_run}, "perl '$Bin/03_SynNetFilter.pl' --synnet '$chain_fwd/synnet' --chain '$chain_fwd/chain' --gd '$opt{genome_dir}' --outdir '$filter_fwd' --threads $opt{threads} --ucsc '$opt{ucsc}' --tba '$opt{tba}' $filter_args");
    if ($independent_rec_axt && (-d "$chain_rec/synnet" || $opt{dry_run})) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/03_SynNetFilter.pl' --synnet '$chain_rec/synnet' --chain '$chain_rec/chain' --gd '$opt{genome_dir}' --outdir '$filter_rec' --threads $opt{threads} --ucsc '$opt{ucsc}' --tba '$opt{tba}' $filter_args");
    }
} else {
    logmsg($LOG, "Skipping SynNet filtering; PairwiseSV will use raw synnet files and genotyping will use singRaw MAF files.");
}

my $pair_syn_fwd = $opt{skip_filter} ? "$chain_fwd/synnet" : "$filter_fwd/CleanSynNet";
my $pair_syn_rec = $opt{skip_filter} ? "$chain_rec/synnet" : "$filter_rec/CleanSynNet";
my $singdir = $opt{skip_filter} ? "$chain_fwd/singRaw" : "$filter_fwd/filtered_sing";

# 4. Pairwise SV calling.
#    reciprocal_strategy controls which CleanSynNet files enter PairwiseSV:
#      qnet : forward tnet + forward qnet only.
#      true : forward tnet + true reciprocal tnet only.
#      both : forward tnet/qnet + true reciprocal tnet/qnet.
#    auto selects true if independent reciprocal alignment was provided, otherwise qnet.
if (!$opt{skip_pairwise_sv}) {
    my $pair_input_root = "$opt{outdir}/04_pairwise_inputs";
    make_path($pair_input_root) unless -d $pair_input_root;

    if ($reciprocal_strategy eq 'qnet') {
        my $fwd_all = prepare_synnet_subset($pair_syn_fwd, "$pair_input_root/forward_tnet_qnet", 'all', $LOG, $opt{dry_run});
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/04_PairwiseSV.pl' --syndir '$fwd_all' --outdir '$pairsv' --gd '$opt{genome_dir}' --threads $opt{threads} --invlen $opt{filter_invlen} --seq --exact_seq_len");
    } elsif ($reciprocal_strategy eq 'true') {
        my $fwd_tnet = prepare_synnet_subset($pair_syn_fwd, "$pair_input_root/forward_tnet_only", 'tnet', $LOG, $opt{dry_run});
        my $rec_tnet = prepare_synnet_subset($pair_syn_rec, "$pair_input_root/reciprocal_tnet_only", 'tnet', $LOG, $opt{dry_run});
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/04_PairwiseSV.pl' --syndir '$fwd_tnet' --outdir '$pairsv' --gd '$opt{genome_dir}' --threads $opt{threads} --invlen $opt{filter_invlen} --seq --exact_seq_len");
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/04_PairwiseSV.pl' --syndir '$rec_tnet' --outdir '$pairsv' --gd '$opt{genome_dir}' --threads $opt{threads} --invlen $opt{filter_invlen} --seq --exact_seq_len");
    } elsif ($reciprocal_strategy eq 'both') {
        my $fwd_all = prepare_synnet_subset($pair_syn_fwd, "$pair_input_root/forward_tnet_qnet", 'all', $LOG, $opt{dry_run});
        my $rec_all = prepare_synnet_subset($pair_syn_rec, "$pair_input_root/reciprocal_tnet_qnet", 'all', $LOG, $opt{dry_run});
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/04_PairwiseSV.pl' --syndir '$fwd_all' --outdir '$pairsv' --gd '$opt{genome_dir}' --threads $opt{threads} --invlen $opt{filter_invlen} --seq --exact_seq_len");
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/04_PairwiseSV.pl' --syndir '$rec_all' --outdir '$pairsv' --gd '$opt{genome_dir}' --threads $opt{threads} --invlen $opt{filter_invlen} --seq --exact_seq_len");
    } else {
        die "ERROR: internal error: unknown reciprocal_strategy $reciprocal_strategy\n";
    }
}

# 5. Merge pairwise raw SVs into population catalogs.
if (!$opt{skip_combine}) {
    run_cmd($LOG, $opt{dry_run}, "perl '$Bin/05_CombineSV.pl' --sv '$pairsv' --refname '$opt{refname}' --outdir '$combined' --sort_threads $opt{threads}");
}

# 6. Genotype large DEL/INS, small InDels, and INV. VCF sample order follows query.list.
if (!$opt{skip_genotype}) {
    my $del_large = "$combined/All.DELs.50bplarge.bed.combined.sorted.txt";
    my $del_small = "$combined/All.DELs.50bpsmall.bed.combined.sorted.txt";
    my $ins_large = "$combined/All.INTs.50bplarge.bed.combined.sorted.txt";
    my $ins_small = "$combined/All.INTs.50bpsmall.bed.combined.sorted.txt";
    my $inv_file = "$combined/All.INVs.bed.combined.sorted.txt";

    # 6a. Combine and genotype complex SVs (CXSV) from bidirectional *.cxsv.sort rawSV files.
    if (!$opt{skip_cxsv}) {
        my $cxsv_cmd = "perl '$Bin/12_combine_genotype_CXSV.pl' --refname '$opt{refname}' --query_list '$query_list' --rawsv_dir '$pairsv' --outdir '$cxsv' --break_tolerance $opt{cxsv_break_tolerance} --min_ro $opt{cxsv_min_ro} --len_ratio $opt{cxsv_len_ratio} --ref_fasta '$ref_fasta' --max_explicit_allele $opt{max_vg_allele_len}";
        $cxsv_cmd .= " --unknown_missing" if $opt{cxsv_unknown_missing};
        run_cmd($LOG, $opt{dry_run}, $cxsv_cmd);
        if ($opt{dry_run} || -s "$cxsv/All.CXSV.vcf") {
            run_cmd($LOG, $opt{dry_run}, "cp '$cxsv/All.CXSV.vcf' '$vcfdir/All.CXSV.vcf'");
        } else {
            logmsg($LOG, "Skipping CXSV VCF copy because $cxsv/All.CXSV.vcf was not created or is empty");
        }
        if ($opt{dry_run} || -s "$cxsv/All.CXSV.vg.vcf") {
            run_cmd($LOG, $opt{dry_run}, "cp '$cxsv/All.CXSV.vg.vcf' '$vcfdir/All.CXSV.vg.vcf'");
            run_cmd($LOG, $opt{dry_run}, "cp '$cxsv/All.CXSV.vg.unresolved.tsv' '$vcfdir/All.CXSV.vg.unresolved.tsv'") if $opt{dry_run} || -s "$cxsv/All.CXSV.vg.unresolved.tsv";
        }
    } else {
        logmsg($LOG, "Skipping CXSV combine/genotyping because --skip_cxsv was set");
    }

    run_if_records($LOG, $opt{dry_run}, $del_large, "perl '$Bin/06_genotype_DEL.pl' --bed '$del_large' --maf_dir '$singdir' --ref_fasta '$ref_fasta' --refname '$opt{refname}' --query_list '$query_list' --out '$geno/DEL.large' --threads $opt{threads}");
    run_if_records($LOG, $opt{dry_run}, $ins_large, "perl '$Bin/07_genotype_INS.pl' --bed '$ins_large' --maf_dir '$singdir' --genome_dir '$opt{genome_dir}' --ref_fasta '$ref_fasta' --refname '$opt{refname}' --query_list '$query_list' --out '$geno/INS.large' --threads $opt{threads}");

    my $have_del_small = tabular_file_has_records($del_small);
    my $have_ins_small = tabular_file_has_records($ins_small);
    if (($have_del_small && $have_ins_small) || $opt{dry_run}) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/08_genotype_small_indels.pl' --maf_dir '$singdir' --ref_fasta '$ref_fasta' --refname '$opt{refname}' --query_list '$query_list' --del_bed '$del_small' --ins_bed '$ins_small' --out '$geno/smallIndel' --threads $opt{threads}");
    } elsif ($have_del_small) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/08_genotype_small_indels.pl' --maf_dir '$singdir' --ref_fasta '$ref_fasta' --refname '$opt{refname}' --query_list '$query_list' --del_bed '$del_small' --out '$geno/smallIndel' --threads $opt{threads}");
    } elsif ($have_ins_small) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/08_genotype_small_indels.pl' --maf_dir '$singdir' --ref_fasta '$ref_fasta' --refname '$opt{refname}' --query_list '$query_list' --ins_bed '$ins_small' --out '$geno/smallIndel' --threads $opt{threads}");
    } else {
        logmsg($LOG, "Skipping small InDel genotyping because neither $del_small nor $ins_small exists");
    }

    if ($opt{dry_run} || inv_file_has_records($inv_file)) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/11_INVgenotyping.pl' --inv '$inv_file' --maf_dir '$singdir' --refsize '$opt{refsize}' --refname '$opt{refname}' --speciesname '$opt{speciesname}' --query_list '$query_list' --out '$geno/INV.vcf' --threads $opt{threads}");
    } else {
        logmsg($LOG, "Skipping INV genotyping because no INV records were found in $inv_file");
    }

    my @large_sv_vcfs = grep { vcf_has_records($_) } ("$geno/DEL.large.vcf", "$geno/INS.large.vcf", "$geno/INV.vcf");
    if (@large_sv_vcfs || $opt{dry_run}) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/SVGAP_vcf_tools.pl' --mode merge --vcfs '" . join(',', @large_sv_vcfs) . "' --query_list '$query_list' --refname '$opt{refname}' --refsize '$opt{refsize}' --out '$vcfdir/All.largeSV.raw.vcf'");
    }
    if (vcf_has_records("$geno/smallIndel.vcf") || $opt{dry_run}) {
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/SVGAP_vcf_tools.pl' --mode merge --vcfs '$geno/smallIndel.vcf' --query_list '$query_list' --refname '$opt{refname}' --refsize '$opt{refsize}' --out '$vcfdir/All.smallIndel.vcf'");
    }
}

# 7. SNP genotyping.
if (!$opt{skip_snp}) {
    run_cmd($LOG, $opt{dry_run}, "perl '$Bin/09_SNPgenotyping.pl' --singmaf '$singdir' --refsize '$opt{refsize}' --refname '$opt{refname}' --query_list '$query_list' --speciesname '$opt{speciesname}' --outdir '$snp' --threads $opt{threads}");
    if (vcf_has_records("$snp/All.sorted.SNP.vcf")) {
        run_cmd($LOG, $opt{dry_run}, "cp '$snp/All.sorted.SNP.vcf' '$vcfdir/All.sorted.SNP.vcf'");
    } else {
        logmsg($LOG, "No SNP records were produced; leaving the header-only SNP VCF in $snp");
    }
}

# 8. Prepare coordinate-consistent, vg-ready VCFs.
#    DEL/INS/CXSV are explicit sequence alleles; INV remains symbolic <INV>.
if (!$opt{skip_vg_vcf}) {
    my @sv_candidates = (
        "$vcfdir/All.largeSV.raw.vcf",
        "$vcfdir/All.smallIndel.vcf",
        "$vcfdir/All.CXSV.vg.vcf",
    );
    my @sv_vcfs = $opt{dry_run} ? @sv_candidates : grep { vcf_has_records($_) } @sv_candidates;
    if (@sv_vcfs) {
        run_cmd($LOG, $opt{dry_run},
            "perl '$Bin/13_prepare_vg_vcf.pl' --ref '$ref_fasta' --vcfs '" . join(',',@sv_vcfs) .
            "' --query_list '$query_list' --refname '$opt{refname}' --refsize '$opt{refsize}' --max_explicit_allele $opt{max_vg_allele_len} --out '$vcfdir/All.SV.vg.vcf' --report '$vcfdir/All.SV.vg.report.tsv'");
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/14_validate_vg_vcf.pl' --vcf '$vcfdir/All.SV.vg.vcf' --ref '$ref_fasta' --out '$vcfdir/All.SV.vg.validation.tsv'");
        maybe_bgzip_tabix($LOG, $opt{dry_run}, "$vcfdir/All.SV.vg.vcf", \%opt);
    } else {
        logmsg($LOG, "Skipping All.SV.vg.vcf because no SV/InDel/CXSV component VCF exists");
    }

    my @all_candidates = (@sv_candidates, "$vcfdir/All.sorted.SNP.vcf");
    my @all_vcfs = $opt{dry_run} ? @all_candidates : grep { vcf_has_records($_) } @all_candidates;
    if (@all_vcfs) {
        run_cmd($LOG, $opt{dry_run},
            "perl '$Bin/13_prepare_vg_vcf.pl' --ref '$ref_fasta' --vcfs '" . join(',',@all_vcfs) .
            "' --query_list '$query_list' --refname '$opt{refname}' --refsize '$opt{refsize}' --max_explicit_allele $opt{max_vg_allele_len} --out '$vcfdir/All.SVGAP.vg.vcf' --report '$vcfdir/All.SVGAP.vg.report.tsv'");
        run_cmd($LOG, $opt{dry_run}, "perl '$Bin/14_validate_vg_vcf.pl' --vcf '$vcfdir/All.SVGAP.vg.vcf' --ref '$ref_fasta' --out '$vcfdir/All.SVGAP.vg.validation.tsv'");
        maybe_bgzip_tabix($LOG, $opt{dry_run}, "$vcfdir/All.SVGAP.vg.vcf", \%opt);
    }
}

# 9. Large-SV annotation.
if (!$opt{skip_annotation} && $opt{annotation} && vcf_has_records("$vcfdir/All.largeSV.raw.vcf")) {
    my $ann_cmd = "perl '$Bin/10_SVannotation.pl' --ref '$ref_fasta' --vcf '$vcfdir/All.largeSV.raw.vcf' --out '$vcfdir/All.largeSV.annotated.vcf' --workdir '$ann' --threads $opt{threads} --minimap2 '$opt{minimap2}'";
    $ann_cmd .= " --te '$opt{te}'" if $opt{te};
    $ann_cmd .= " --cds '$opt{cds}'" if $opt{cds};
    run_cmd($LOG, $opt{dry_run}, $ann_cmd);
}

logmsg($LOG, "Finished SVGAP pipeline");
print STDERR "Done. Main outputs are in $vcfdir\n";
close $LOG;

sub validate_options {
    my ($opt) = @_;
    die "ERROR: --start_from must be align, maf or axt.\n" . usage() unless $opt->{start_from} =~ /^(align|maf|axt)$/;
    die "ERROR: --aligner must be minimap2, last, lastz, mummer/mummer4, anchorwave, or maf.\n"
        unless $opt->{aligner} =~ /^(minimap2|last|lastz|mummer|mummer4|anchorwave|maf)$/;
    die "ERROR: --align_mode must be auto, chr, or ctg.\n"
        unless $opt->{align_mode} =~ /^(auto|chr|ctg)$/;
    die "ERROR: --threads must be >= 1.\n" unless $opt->{threads} >= 1;
    die "ERROR: --reciprocal_strategy must be auto, qnet, true or both.\n" . usage()
        unless $opt->{reciprocal_strategy} =~ /^(auto|qnet|true|both)$/;
    die "ERROR: --cxsv_break_tolerance must be >= 0.\n" unless $opt->{cxsv_break_tolerance} >= 0;
    die "ERROR: --cxsv_min_ro must be between 0 and 1.\n" unless $opt->{cxsv_min_ro} >= 0 && $opt->{cxsv_min_ro} <= 1;
    die "ERROR: --cxsv_len_ratio must be >= 0.\n" unless $opt->{cxsv_len_ratio} >= 0;
    die "ERROR: --max_vg_allele_len must be >= 1.\n" unless $opt->{max_vg_allele_len} >= 1;
    for my $name (qw/filter_toplen filter_nonsynlen filter_synlen filter_invlen filter_qfar_cutoff filter_auto_min/) {
        die "ERROR: --$name must be >= 0.\n" unless $opt->{$name} >= 0;
    }
    die "ERROR: --filter_alrate must be between 0 and 1.\n"
        unless $opt->{filter_alrate} >= 0 && $opt->{filter_alrate} <= 1;
    die "ERROR: --filter_auto_fraction must be between 0 and 1.\n"
        unless $opt->{filter_auto_fraction} >= 0 && $opt->{filter_auto_fraction} <= 1;

    return if $opt->{skip_convert};

    if ($opt->{start_from} eq 'align') {
        die "ERROR: --align_dir is required when --start_from align.\n" . usage() unless $opt->{align_dir};
    } elsif ($opt->{start_from} eq 'maf') {
        die "ERROR: --maf_dir is required when --start_from maf.\n" . usage() unless $opt->{maf_dir};
    } elsif ($opt->{start_from} eq 'axt') {
        die "ERROR: --axt_dir or --axt_forward_dir is required when --start_from axt.\n" . usage()
            unless $opt->{axt_dir} || $opt->{axt_forward_dir};
    }
}

sub prepare_or_resolve_axt {
    my ($opt, $axt_root, $LOG) = @_;
    my ($fwd, $rec);

    if ($opt->{skip_convert}) {
        ($fwd, $rec) = resolve_existing_axt_dirs($opt, $axt_root, $LOG);
        return ($fwd, $rec);
    }

    if ($opt->{start_from} eq 'align') {
        # Convert aligner output to forward AXT only. No fake reciprocal AXT is generated.
        my $convert_cmd = "perl '$Bin/01_Convert2Axt.pl' --ali '$opt->{aligner}' --input '$opt->{align_dir}' --wk '$axt_root' --tname '$opt->{refname}' --qname Query --mode '$opt->{align_mode}' --t $opt->{threads} --ucsc '$opt{ucsc}' --tba '$opt{tba}' --k8 '$opt{k8}' --paftools '$opt{paftools}' --mum320 '$opt{mum320}'";
        $convert_cmd .= join('', map { " --sizes '$_'" } @{ $opt->{align_sizes} || [] });
        run_cmd($LOG, $opt->{dry_run}, $convert_cmd);
        $fwd = "$axt_root/Target_$opt->{refname}";
    } elsif ($opt->{start_from} eq 'maf') {
        my $convert_cmd = "perl '$Bin/01_Convert2Axt.pl' --ali maf --input '$opt->{maf_dir}' --wk '$axt_root' --tname '$opt->{refname}' --qname Query --mode '$opt->{align_mode}' --t $opt->{threads} --ucsc '$opt{ucsc}' --tba '$opt{tba}'";
        $convert_cmd .= join('', map { " --sizes '$_'" } @{ $opt->{align_sizes} || [] });
        run_cmd($LOG, $opt->{dry_run}, $convert_cmd);
        $fwd = "$axt_root/Target_$opt->{refname}";
        if ($opt->{maf_reciprocal_dir}) {
            my $rec_root = "$axt_root/independent_reciprocal_from_maf";
            my $rec_convert_cmd = "perl '$Bin/01_Convert2Axt.pl' --ali maf --input '$opt->{maf_reciprocal_dir}' --wk '$rec_root' --tname Reciprocal --qname Ref --mode '$opt->{align_mode}' --t $opt->{threads} --ucsc '$opt{ucsc}' --tba '$opt{tba}'";
            $rec_convert_cmd .= join('', map { " --sizes '$_'" } @{ $opt->{align_sizes} || [] });
            run_cmd($LOG, $opt->{dry_run}, $rec_convert_cmd);
            $rec = "$rec_root/Target_Reciprocal";
        }
    } elsif ($opt->{start_from} eq 'axt') {
        ($fwd, $rec) = resolve_existing_axt_dirs($opt, $axt_root, $LOG);
    }

    die "ERROR: forward AXT directory does not exist: $fwd\n" unless $opt->{dry_run} || -d $fwd;
    if (!$opt->{dry_run}) {
        my @axt = glob("$fwd/*.axt");
        die "ERROR: no .axt files found in forward AXT directory: $fwd\n" unless @axt;
    }
    logmsg($LOG, "Forward AXT directory: $fwd");
    logmsg($LOG, "Independent reciprocal AXT directory: " . ($rec || 'NA'));
    logmsg($LOG, "Note: qnet generated from forward AXT will be used for reference insertion discovery.");
    return ($fwd, $rec);
}


sub choose_reciprocal_strategy {
    my ($opt, $independent_rec_axt, $LOG) = @_;
    my $strategy = $opt->{reciprocal_strategy};
    if ($strategy eq 'auto') {
        $strategy = $independent_rec_axt ? 'true' : 'qnet';
    }
    if (($strategy eq 'true' || $strategy eq 'both') && !$independent_rec_axt) {
        die "ERROR: --reciprocal_strategy $strategy requires true independent reciprocal input. Provide --axt_reciprocal_dir or --maf_reciprocal_dir.\n";
    }
    logmsg($LOG, "effective reciprocal_strategy=$strategy");
    if ($strategy eq 'qnet') {
        logmsg($LOG, "PairwiseSV will use forward tnet plus forward qnet. No independent reciprocal alignment is required.");
    } elsif ($strategy eq 'true') {
        logmsg($LOG, "PairwiseSV will use forward tnet plus true reciprocal tnet. Forward qnet will not be used for SV discovery to avoid redundant calls.");
    } elsif ($strategy eq 'both') {
        logmsg($LOG, "PairwiseSV will use both forward tnet/qnet and true reciprocal tnet/qnet. This is the most inclusive mode and may create redundant candidates before merge.");
    }
    return $strategy;
}

sub prepare_synnet_subset {
    my ($src, $dest, $which, $LOG, $dry) = @_;
    make_path($dest) unless -d $dest;
    if ($dry) {
        logmsg($LOG, "DRY-RUN subset CleanSynNet: source=$src dest=$dest type=$which");
        return $dest;
    }
    die "ERROR: CleanSynNet directory does not exist: $src\n" unless -d $src;
    opendir(my $DH, $src) or die "ERROR: cannot open CleanSynNet directory $src: $!\n";
    my @files;
    while (my $f = readdir($DH)) {
        next if $f =~ /^\./;
        next if $f =~ /\.out$/;
        next unless $f =~ /\.synnet$/;
        if ($which eq 'all') {
            push @files, $f if $f =~ /chain\.filter\.(?:t|q)net\.synnet/;
        } elsif ($which eq 'tnet') {
            push @files, $f if $f =~ /chain\.filter\.tnet\.synnet/;
        } elsif ($which eq 'qnet') {
            push @files, $f if $f =~ /chain\.filter\.qnet\.synnet/;
        } else {
            die "ERROR: unknown synnet subset type: $which\n";
        }
    }
    closedir($DH);
    die "ERROR: no $which synnet files found in $src\n" unless @files;

    for my $f (@files) {
        my $from = "$src/$f";
        my $to   = "$dest/$f";
        unlink $to if -e $to || -l $to;
        # A hard link avoids copying very large SynNet files and cannot suffer
        # from the broken-relative-symlink issue seen in older controllers.
        # Fall back to a real copy across filesystems or on unsupported storage.
        link($from, $to) || copy($from, $to)
            or die "ERROR: cannot link or copy $from to $to: $!\n";
    }
    logmsg($LOG, "Prepared PairwiseSV input subset: $dest ($which; n=" . scalar(@files) . ")");
    return $dest;
}

sub resolve_existing_axt_dirs {
    my ($opt, $default_axt_dir, $LOG) = @_;
    my ($fwd, $rec);
    $fwd = $opt->{axt_forward_dir} if $opt->{axt_forward_dir};
    $rec = $opt->{axt_reciprocal_dir} if $opt->{axt_reciprocal_dir};

    if (!$fwd && $opt->{axt_dir}) {
        my $ad = $opt->{axt_dir};
        $ad =~ s{/+$}{};
        if (-d "$ad/Target_$opt->{refname}") {
            $fwd = "$ad/Target_$opt->{refname}";
            $rec ||= "$ad/Target_Query" if -d "$ad/Target_Query";
            $rec ||= "$ad/independent_reciprocal" if -d "$ad/independent_reciprocal";
        } elsif (-d "$ad/forward") {
            $fwd = "$ad/forward";
            $rec ||= "$ad/reciprocal" if -d "$ad/reciprocal";
        } else {
            $fwd = $ad;  # flat directory containing RefvsQuery.axt
        }
    }

    die "ERROR: cannot resolve forward AXT directory. Use --axt_dir or --axt_forward_dir.\n" unless $fwd;
    die "ERROR: forward AXT directory does not exist: $fwd\n" unless $opt->{dry_run} || -d $fwd;
    if ($rec && !$opt->{dry_run} && !-d $rec) {
        die "ERROR: independent reciprocal AXT directory does not exist: $rec\n";
    }
    return ($fwd, $rec);
}

sub prepare_lists {
    my ($opt, $LOG) = @_;
    open(my $IN, '<', $opt->{glist}) or die "ERROR: cannot read --glist $opt->{glist}: $!\n";
    my (@all, @queries, %seen);
    while (<$IN>) {
        chomp;
        s/^\s+|\s+$//g;
        next if $_ eq '' || /^#/;
        die "ERROR: duplicated genome ID in --glist: $_\n" if $seen{$_}++;
        push @all, $_;
        push @queries, $_ unless $_ eq $opt->{refname};
    }
    close $IN;
    die "ERROR: --refname $opt->{refname} was not found in --glist\n" unless $seen{$opt->{refname}};
    die "ERROR: no query genomes remain after excluding --refname\n" unless @queries;
    my $cfg = "$opt->{outdir}/00_config";
    make_path($cfg);
    my $all_list = "$cfg/all.g.list";
    my $query_list = "$cfg/query.list";
    write_list($all_list, @all);
    write_list($query_list, @queries);
    logmsg($LOG, "Genome list: ".join(',', @all));
    logmsg($LOG, "Query order for VCF: ".join(',', @queries));
    return ($all_list, $query_list, \@queries);
}

sub write_list {
    my ($file, @items) = @_;
    open(my $OUT, '>', $file) or die "ERROR: cannot write $file: $!\n";
    print $OUT join("\n", @items), "\n";
    close $OUT;
}

sub find_genome_fasta {
    my ($gd, $id) = @_;
    for my $p ("$gd/$id", "$gd/$id.fa", "$gd/$id.fasta", "$gd/$id.fna", "$gd/$id.fas") {
        return $p if -s $p;
    }
    die "ERROR: cannot find genome FASTA for $id in $gd. Tried ID, ID.fa, ID.fasta, ID.fna and ID.fas\n";
}

sub find_refsize {
    my ($gd, $ref) = @_;
    for my $p ("$gd/$ref.sizes", "$gd/$ref.fa.sizes", "$gd/$ref.fasta.sizes") {
        return $p if -s $p;
    }
    return '';
}

sub run_if_records {
    my ($LOG, $dry, $file, $cmd) = @_;
    if ($dry || tabular_file_has_records($file)) { run_cmd($LOG, $dry, $cmd); }
    else { logmsg($LOG, "Skip because file has no data records: $file"); }
}

sub tabular_file_has_records {
    my ($file) = @_;
    return 0 unless defined $file && -s $file;
    open(my $IN, '<', $file) or return 0;
    while (my $line = <$IN>) {
        next if $line =~ /^\s*$/ || $line =~ /^#/;
        close $IN;
        return 1;
    }
    close $IN;
    return 0;
}

sub vcf_has_records {
    return tabular_file_has_records($_[0]);
}

# Return true only if the merged INV file contains at least one real INV record.
# Some SVGAP runs produce an INV catalog file with only a header, or with no valid
# INV rows. In that case INV genotyping should be skipped rather than treated as
# a fatal error.
sub inv_file_has_records {
    my ($file) = @_;
    return 0 unless defined $file && -s $file;
    open(my $IN, '<', $file) or return 0;
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;
        my @f = split(/	/, $line);
        next unless @f >= 6;
        my $svtype = uc($f[4] // '');
        next unless $svtype eq 'INV' || $svtype =~ /INV/;
        my ($start, $end) = @f[1,2];
        next unless defined $start && defined $end && $start =~ /^-?\d+(?:\.\d+)?$/ && $end =~ /^-?\d+(?:\.\d+)?$/;
        if ($end > $start) {
            close $IN;
            return 1;
        }
    }
    close $IN;
    return 0;
}

sub maybe_bgzip_tabix {
    my ($LOG,$dry,$vcf,$opt)=@_;
    return unless $opt->{index_vg_vcf};
    my $bgzip = $dry ? $opt->{bgzip} : resolve_executable($opt->{bgzip});
    my $tabix = $dry ? $opt->{tabix} : resolve_executable($opt->{tabix});
    if (!$bgzip || !$tabix) {
        logmsg($LOG, "WARNING: bgzip/tabix not found; kept plain $vcf. Install htslib or provide --bgzip/--tabix, then run bgzip -f $vcf && tabix -f -p vcf $vcf.gz");
        return;
    }
    run_cmd($LOG,$dry,"'$bgzip' -f -c '$vcf' > '$vcf.gz'");
    run_cmd($LOG,$dry,"'$tabix' -f -p vcf '$vcf.gz'");
}

sub resolve_executable {
    my ($cmd)=@_;
    return undef unless defined($cmd) && $cmd ne '';
    return abs_path($cmd) if -x $cmd;
    my $safe=$cmd; $safe =~ s/'/'\\''/g;
    my $found=`command -v '$safe' 2>/dev/null`;
    chomp $found;
    return $found if $found ne '' && -x $found;
    return undef;
}

sub run_cmd {
    my ($LOG, $dry, $cmd) = @_;
    logmsg($LOG, $cmd);
    if ($dry) { print STDERR "[DRY-RUN] $cmd\n"; return; }
    my $ret = system($cmd);
    if ($ret != 0) {
        die "ERROR: command failed with status $ret:\n$cmd\nSee $log for previous steps.\n";
    }
}

sub logmsg {
    my ($LOG, $msg) = @_;
    # Be tolerant if a scalar reference to the filehandle was accidentally passed.
    # Older pipeline versions passed \$LOG to prepare_synnet_subset(), which caused
    # 'Not a GLOB reference' when logmsg tried to print to it.
    $LOG = $$LOG if ref($LOG) eq 'SCALAR';
    my $t = strftime('%F %T', localtime);
    print $LOG "[$t] $msg\n" if defined $LOG;
    print STDERR "[$t] $msg\n";
}

sub usage {
    return <<'USAGE';
Usage:
  SVGAP version2

  # 1. Start from aligner outputs. 01_Convert2Axt.pl will convert them to AXT.
  perl SVGAP_pipeline.pl \
    --refname Ref \
    --glist g.list \
    --genome_dir ./genomes \
    --start_from align \
    --aligner minimap2 \
    --align_dir ./alignments \
    --outdir SVGAP_run \
    --threads 16

  # 2. Start directly from existing forward MAF files.
  #    File names should be RefvsQuery.maf, RefvsQuery.chr.maf, or RefvsQuery.ctg.maf.
  perl SVGAP_pipeline.pl \
    --refname Ref \
    --glist g.list \
    --genome_dir ./genomes \
    --start_from maf \
    --maf_dir ./MAF_forward \
    --outdir SVGAP_run_from_maf \
    --threads 16

  # 3. Start directly from existing forward AXT files.
  #    File names should be RefvsQuery.axt.
  perl SVGAP_pipeline.pl \
    --refname Ref \
    --glist g.list \
    --genome_dir ./genomes \
    --start_from axt \
    --axt_dir ./AXT_forward \
    --outdir SVGAP_run_from_axt \
    --threads 16

  # 4. Optional: provide true independently generated reciprocal Query-vs-Ref AXT.
  #    This preserves the old reciprocal-alignment workflow. Do NOT generate this
  #    by simply reversing forward AXT.
  perl SVGAP_pipeline.pl \
    --refname Ref \
    --glist g.list \
    --genome_dir ./genomes \
    --start_from axt \
    --axt_forward_dir ./AXT_forward \
    --axt_reciprocal_dir ./AXT_query_vs_ref_independent \
    --reciprocal_strategy true \
    --outdir SVGAP_run_with_true_reciprocal \
    --threads 16

Required:
  --refname       Reference genome ID. This ID must appear in g.list.
  --glist         Genome list. One genome ID per line. VCF sample order is exactly
                  g.list order after removing --refname.
  --genome_dir    Directory containing reference/query genome FASTA files.

Start options:
  --start_from    align, maf or axt [align]
  --align_dir     Pairwise aligner output directory. Required for --start_from align.
  --aligner       minimap2, last, lastz, mummer/mummer4, anchorwave, or maf [minimap2]
                  Names are case-insensitive.
  --align_mode    auto, chr, or ctg [auto]. Use ctg for split-query alignments.
  --align_sizes FILE      Original chromosome size file for ctg mode; repeatable.
  --maf_dir       Existing forward MAF directory. Required for --start_from maf.
  --axt_dir       Existing forward AXT directory, or a parent directory containing
                  Target_<refname>/.
  --axt_forward_dir       Directory with forward RefvsQuery.axt files.
  --axt_reciprocal_dir    Optional directory with true independent QueryvsRef.axt files.
  --maf_reciprocal_dir    Optional directory with true independent Query-vs-Ref MAFs.
  --reciprocal_strategy   auto, qnet, true, or both [auto]
                          auto: use true reciprocal if provided; otherwise use qnet
                          qnet: use forward tnet+qnet only
                          true: use forward tnet + true reciprocal tnet only
                          both: use forward tnet/qnet + reciprocal tnet/qnet

Important reciprocal-alignment logic:
  A true Query-vs-Ref alignment is optional. If only forward Ref-vs-Query AXT/MAF
  is available, SVGAP still runs. The forward chainNet step generates both:
      RefvsQuery.chain.filter.tnet.synnet
      RefvsQuery.chain.filter.qnet.synnet
  and qnet-derived calls can be used as Query.Ref.rawSV for reference insertions.

  If you independently ran Query-vs-Ref alignments, provide them with
  --axt_reciprocal_dir or --maf_reciprocal_dir and use --reciprocal_strategy true
  to preserve the old reciprocal workflow. In this mode PairwiseSV uses forward
  tnet plus reciprocal tnet, avoiding duplicated qnet-derived candidates.

  Do not create reciprocal AXT by simply reversing a forward AXT. That is only a
  coordinate reordering, not an independent Query-vs-Ref alignment.

SynNet filtering options:
  --skip_filter                Use raw chain/net synnet plus singRaw MAF directly.
  --filter_toplen N            Top-level fill threshold ceiling [500000].
  --filter_nonsynlen N         Lower non-syntenic fill threshold ceiling [500000].
  --filter_synlen N            Lower syntenic fill threshold ceiling [500000].
  --filter_invlen N            Minimum inversion fill/call length [1000].
  --filter_alrate FLOAT        Minimum aligned fraction [0.10].
  --filter_qfar_cutoff N       Maximum qFar for inversion fills [5000000].
  --filter_adaptive / --no-filter_adaptive
                               Scale fill thresholds down on short contigs [enabled].
  --filter_auto_fraction FLOAT Per-contig threshold fraction [0.10].
  --filter_auto_min N          Lower bound for adaptive thresholds [1000].

vg-ready VCF options:
  --skip_vg_vcf               Do not create final graph-construction VCFs.
  --max_vg_allele_len N       Maximum explicit REF/ALT allele length [5000000].
  --index_vg_vcf / --no-index_vg_vcf  bgzip and tabix final VCFs [enabled].
  --bgzip FILE                bgzip executable or command name [bgzip].
  --tabix FILE                tabix executable or command name [tabix].

CXSV options:
  --skip_cxsv                 Skip complex SV merging/genotyping.
  --cxsv_break_tolerance N    Breakpoint tolerance for merging CXSVs [2000].
  --cxsv_min_ro FLOAT         Minimum reciprocal overlap for CXSV merging [0.50].
  --cxsv_len_ratio FLOAT      Maximum relative length difference for CXSV merging [1.00].
  --cxsv_unknown_missing      Output ./. for non-carrier CXSV samples instead of 0|0.

Main outputs:
  OUTDIR/08_vcf/All.largeSV.raw.vcf
  OUTDIR/08_vcf/All.largeSV.annotated.vcf   if annotation is enabled
  OUTDIR/08_vcf/All.smallIndel.vcf
  OUTDIR/08_vcf/All.sorted.SNP.vcf
  OUTDIR/08_vcf/All.CXSV.vcf
  OUTDIR/08_vcf/All.CXSV.vg.vcf
  OUTDIR/08_vcf/All.SV.vg.vcf
  OUTDIR/08_vcf/All.SVGAP.vg.vcf
  OUTDIR/08_vcf/All.SVGAP.vg.vcf.gz + .tbi   when bgzip/tabix are available
  OUTDIR/08_vcf/*.vg.report.tsv
  OUTDIR/08_vcf/*.vg.validation.tsv
USAGE
}
