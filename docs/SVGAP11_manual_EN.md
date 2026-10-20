---
title: "SVGAP11-CXSV Clean Package Manual"
author: "Yi Liao @ scau.edu.cn"
date: "2026-07-07"
geometry: margin=0.75in
mainfont: "Noto Serif"
monofont: "DejaVu Sans Mono"
fontsize: 10pt
colorlinks: true
---

# SVGAP version2 Manual

Version naming: the previous repository workflow is **version1**; the repaired
and optimized workflow documented here is **version2**.

**Author tag**: every Perl script in this package now contains `Author: Yi Liao @ scau.edu.cn`.

**Purpose**: run assembly-to-assembly SV discovery, merging, complex-SV handling, SNP/small-InDel genotyping and optional annotation from one reference assembly and multiple query assemblies.
**Example dataset**: `ISO1` is the reference genome; `A3`, `A4` and `A5` are query genomes.

## 1. What was fixed and cleaned

### 1.1 Parallel execution without an extra module

The recurring warning was:

```text
WARNING: Parallel::ForkManager is not installed; falling back to single-thread mode
```

Parallel steps now use Perl's built-in `fork()` and `wait()`. They no longer
depend on `Parallel::ForkManager` and do not fall back to serial execution because
of CPAN module path differences. `--threads` directly controls the number of
concurrent sample or file jobs.

All `.pl` scripts still include package-local library paths through `FindBin`, for example:

```perl
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/pub";
use lib "$Bin/../lib";
use lib "$Bin/../pub";
use lib "$Bin/../../lib";
use lib "$Bin/../../pub";
```

Dependency check:

```bash
cd SVGAP11_clean_fixed
perl util/Test_dependency.pl
```

### 1.2 Author line added

Every Perl script now contains:

```perl
# Author: Yi Liao @ scau.edu.cn
```

### 1.3 Non-production files moved to trash

The following files were moved to `trash/`, not deleted:

```text
trash/benchmark_tools/
  simulate_population_variants_benchmark_shared.v4.pl
  13_compare_truth_SVGAP.pl
  14_diagnose_small_indel_benchmark.pl

trash/integration_notes/
  SVGAP_pipeline.add_12_CXSV.patch
  README_CXSV_integration.txt

trash/legacy_utilities/
  FaFileChecl.pl
  INS2TDP.pl
  Mafreformat.pl
  SplitFa.pl
```

They are not required for running the production ISO1/A3/A4/A5 workflow. They can be restored if simulation, benchmarking or legacy utility behavior is needed.

## 2. Clean package structure

```text
SVGAP11_clean_fixed/
  SVGAP_pipeline.pl                 # Main controller for the full workflow
  01_Convert2Axt.pl                 # Step 01: convert alignments to AXT
  02_ChainNetSyn.pl                 # Step 02: AXT -> chain/net/synnet
  03_SynNetFilter.pl                # Step 03: synnet filtering + single-coverage MAF
  04_PairwiseSV.pl                  # Step 04: pairwise SV calling
  05_CombineSV.pl                   # Step 05: merge SV candidates
  06_genotype_DEL.pl                # Step 06: deletion genotyping
  07_genotype_INS.pl                # Step 07: insertion genotyping
  08_genotype_small_indels.pl       # Step 08: small-InDel genotyping
  09_SNPgenotyping.pl               # Step 09: SNP genotyping
  10_SVannotation.pl                # Step 10: SV annotation
  11_INVgenotyping.pl               # Step 11: inversion genotyping
  12_combine_genotype_CXSV.pl       # Step 12: complex-SV merging and genotyping
  SVGAP_vcf_tools.pl                # VCF merging and sample-order utility
  lib/                              # Bundled Perl dependencies, including Parallel::ForkManager
  pub/                              # UCSC/TBA/MUMMER/Stretcher helpers
  util/Test_dependency.pl           # Dependency test script
  docs/                             # Manuals
  trash/                            # Archived non-production files
```

## 3. Input preparation: ISO1 reference + A3/A4/A5 queries

### 3.1 Recommended directory layout

```text
project/
  SVGAP11_clean_fixed/
  genomes/
    ISO1.fa
    A3.fa
    A4.fa
    A5.fa
  g.list
  01_axt/Target/                 # If starting from AXT
    ISO1vsA3.axt
    ISO1vsA4.axt
    ISO1vsA5.axt
```

### 3.2 g.list

`g.list` contains one genome ID per line. The reference ID must be included. The pipeline removes the reference from the VCF sample columns, so the final sample order will be `A3, A4, A5`.

```text
ISO1
A3
A4
A5
```

### 3.3 Genome FASTA names

For each genome ID, the scripts search for one of these names:

```text
ID
ID.fa
ID.fasta
ID.fna
```

For this example, use:

```text
genomes/ISO1.fa
genomes/A3.fa
genomes/A4.fa
genomes/A5.fa
```

Recommended FASTA headers:

```text
>ISO1.Chr01
>A3.Chr01
>A4.Chr01
>A5.Chr01
```

Avoid sequence IDs beginning with numbers.

### 3.4 AXT file names

If starting directly from AXT, forward AXT files must be named:

```text
ISO1vsA3.axt
ISO1vsA4.axt
ISO1vsA5.axt
```

If true independently generated reciprocal alignments are available, prepare:

```text
A3vsISO1.axt
A4vsISO1.axt
A5vsISO1.axt
```

Do not generate reciprocal AXT by simply reversing forward AXT. That only reorders coordinates; it is not an independent query-vs-reference alignment.

## 4. Running the full pipeline

### 4.1 Recommended mode: start from existing AXT

```bash
cd project/SVGAP11_clean_fixed

perl SVGAP_pipeline.pl \
  --refname ISO1 \
  --glist ../g.list \
  --genome_dir ../genomes \
  --start_from axt \
  --axt_dir ../01_axt/Target \
  --outdir ../SVGAP_run \
  --threads 12 \
  --speciesname ISO1_A3_A4_A5 \
  --skip_annotation
```

`--skip_annotation` is recommended for the first training run because annotation may require TE/CDS inputs and a working minimap2 installation. The core SV/SNP/InDel pipeline can be validated first.

### 4.2 Start from MAF

Forward MAF files should be named like:

```text
ISO1vsA3.maf
ISO1vsA4.maf
ISO1vsA5.maf
```

Command:

```bash
perl SVGAP_pipeline.pl \
  --refname ISO1 \
  --glist ../g.list \
  --genome_dir ../genomes \
  --start_from maf \
  --maf_dir ../MAF_forward \
  --outdir ../SVGAP_run_from_maf \
  --threads 12 \
  --speciesname ISO1_A3_A4_A5 \
  --skip_annotation
```

### 4.3 Start from aligner outputs

If the inputs are minimap2/LAST/LASTZ/MUMmer/AnchorWave outputs, `01_Convert2Axt.pl` converts them to AXT.

```bash
perl SVGAP_pipeline.pl \
  --refname ISO1 \
  --glist ../g.list \
  --genome_dir ../genomes \
  --start_from align \
  --aligner minimap2 \
  --align_dir ../alignments \
  --outdir ../SVGAP_run_from_align \
  --threads 12 \
  --speciesname ISO1_A3_A4_A5 \
  --skip_annotation
```

### 4.4 Use true reciprocal alignments

```bash
perl SVGAP_pipeline.pl \
  --refname ISO1 \
  --glist ../g.list \
  --genome_dir ../genomes \
  --start_from axt \
  --axt_forward_dir ../AXT_forward \
  --axt_reciprocal_dir ../AXT_query_vs_ref_independent \
  --reciprocal_strategy true \
  --outdir ../SVGAP_run_true_reciprocal \
  --threads 12 \
  --speciesname ISO1_A3_A4_A5 \
  --skip_annotation
```

If no true reciprocal alignments are available, the default `auto` strategy uses forward tnet/qnet and is appropriate.

## 5. Main output directories

The pipeline creates the following directories under `--outdir`:

```text
SVGAP_run/
  01_axt/                         # Prepared or converted AXT files
  02_chain_forward/               # Forward chain/net/synnet
  02_chain_independent_reciprocal/ # Optional reciprocal chain/net/synnet
  03_filter_forward/              # Filtered synnet + single-coverage MAF
  03_filter_independent_reciprocal/
  04_pairwise_inputs/             # Synnet subsets for PairwiseSV
  04_pairwise_sv/                 # Pairwise raw SVs
  05_combined_sv/                 # Merged DEL/INS/INV/CNV/TRL candidates
  06_genotype/                    # DEL/INS/small-InDel/INV/CXSV genotyping
  07_snp/                         # SNP calling/genotyping
  08_vcf/                         # Final VCF files
  09_annotation/                  # Optional annotation workdir
  SVGAP_pipeline.log              # Main log file
```

Key final VCFs:

```text
SVGAP_run/08_vcf/All.largeSV.raw.vcf
SVGAP_run/08_vcf/All.largeSV.annotated.vcf     # if annotation is enabled
SVGAP_run/08_vcf/All.smallIndel.vcf
SVGAP_run/08_vcf/All.sorted.SNP.vcf
SVGAP_run/08_vcf/All.CXSV.vcf
```

## 6. Detailed workflow steps and parameters

### Step 01: 01_Convert2Axt.pl

**Purpose**: convert LAST, LASTZ, MUMmer, minimap2, AnchorWave or existing MAF alignments to AXT. It can also split query FASTA into overlapping subcontigs and restore split-contig coordinates to original chromosome coordinates.

```bash
perl 01_Convert2Axt.pl \
  --ali minimap2 \
  --input ../alignments \
  --wk ../01_axt \
  --tname ISO1 \
  --qname A3 \
  --t 8
```

| Option | Meaning |
|---|---|
| `--ali` | Alignment type: `last`, `lastz`, `maf`, `mummer`, `minimap2`, `anchorwave`. |
| `--input` | Input alignment directory. |
| `--wk` | Working/output directory. |
| `--tname` | Target/reference name. Default: `Ref`. |
| `--qname` | Query name. Default: `Que`. |
| `--mode` / `--coord` | `auto`, `chr` or `ctg`. Use `ctg` for split query subcontig alignments. |
| `--t` / `--threads` | Number of parallel workers. Default: 4. |
| `--query` | Also generate query-as-reference AXT files. |
| `--keep_tmp` | Keep temporary files. |
| `--mum320` | Directory/path for MUMmer `delta2maf`. |
| `--ucsc` | Directory/path for UCSC utilities such as `mafToAxt`. |
| `--tba` | Directory/path for TBA utilities. |
| `--k8` | Path to `k8`. |
| `--paftools` | Path to `paftools.js`. |
| `--splitfa` / `--fasta` | FASTA file to split; can be repeated. |
| `--bin` | Main split window size. Default: 10,000,000 bp. |
| `--slide` | Overlap between windows. Default: 10,000 bp. |
| `--split_out` | Output directory for split FASTA. |
| `--sizes` / `--qsize` / `--size` | Original chromosome size file for coordinate restoration in ctg mode. |

### Step 02: 02_ChainNetSyn.pl

**Purpose**: build UCSC chain/net/syntenic-net from AXT and optionally generate target-oriented single-coverage MAF.

```bash
perl 02_ChainNetSyn.pl \
  --gd ../genomes \
  --ad ../01_axt/Target \
  --lst ../g.list \
  --outdir ../SVGAP_run/02_chain_forward \
  --threads 12 \
  --sing
```

| Option | Meaning |
|---|---|
| `--gd` | Genome FASTA directory. The script searches for `ID`, `ID.fa`, `ID.fasta`, or `ID.fna`. |
| `--ad` | AXT directory. Input names must look like `TargetIDvsQueryID.axt`. |
| `--lst` | Genome ID list. |
| `--outdir` / `--syn` | Output directory. |
| `--ucsc` | UCSC utility directory. |
| `--tba` | TBA utility directory. |
| `--linearGap` / `--lg` | ChainNet linearGap file. Default: `pub/TBA/medium`. |
| `--minscore` / `--minScore` | Minimum score for axtChain/chainNet. Default: 1000. |
| `--minSpace` | minSpace for chainNet. Default: 1. |
| `--t` / `--threads` | Number of parallel jobs. Default: 1. |
| `--sing` | Generate single-coverage MAF. |
| `--keep_raw_chain` | Keep raw chain files. |
| `--keep_temp` | Keep temporary AXT/MAF files. |
| `--force` | Overwrite existing outputs. |
| `--dry_run` | Print commands without executing them. |

### Step 03: 03_SynNetFilter.pl

**Purpose**: filter unreliable short fills in synnet, generate CleanSynNet, raw syntenic-gap files and filtered single-coverage MAF.

```bash
perl 03_SynNetFilter.pl \
  --synnet ../SVGAP_run/02_chain_forward/synnet \
  --chain ../SVGAP_run/02_chain_forward/chain \
  --gd ../genomes \
  --outdir ../SVGAP_run/03_filter_forward \
  --threads 12
```

| Option | Meaning |
|---|---|
| `--synnet` | Directory containing synnet files from Step 02. |
| `--chain` | Chain directory from Step 02. |
| `--gd` | Genome FASTA/.sizes/.2bit directory. |
| `--outdir` | Output directory. |
| `--ucsc` | UCSC utility directory. |
| `--tba` | TBA utility directory. |
| `--cleanSynnet` | Output subdirectory for cleaned synnet. Default: `CleanSynNet`. |
| `--sing` | Output subdirectory for filtered single MAF. Default: `filtered_sing`. |
| `--axt` | Output subdirectory for cleaned-synnet AXT. |
| `--homopairs` | Output subdirectory for homologous pair summaries. |
| `--rawsv` | Output subdirectory for raw SV gaps. |
| `--target` / `--refname` | Process only this target/reference. |
| `--query` | Process only this query. |
| `--toplen` | Minimum length for top syntenic fills. Default: 500000. |
| `--nonsynlen` | Minimum length for lower nonSyn fills. Default: 500000. |
| `--synlen` | Minimum length for lower syn fills. Default: 500000. |
| `--invlen` | Minimum length for inversion fills. Default: 10000. |
| `--alrate` | Minimum alignment/fill-length proportion. Default: 0.10. |
| `--qFarCutoff` | Maximum qFar for inversion fills. Default: 5000000. |
| `--minscore` | Minimum score threshold. |
| `--t` / `--threads` | Number of parallel jobs. |
| `--sing_only_tnet` | Generate single MAF from tnet only. |
| `--keep_temp` | Keep temporary files. |
| `--force` | Overwrite existing outputs. |
| `--dry_run` | Print commands only. |

### Step 04: 04_PairwiseSV.pl

**Purpose**: call pairwise SVs from cleaned synnet files, including DEL, INS, CNV, INV, TRL and complex PAV/CMPLX events. This version uses core Perl `fork()/waitpid()` for parallel execution.

```bash
perl 04_PairwiseSV.pl \
  --syndir ../SVGAP_run/03_filter_forward/CleanSynNet \
  --outdir ../SVGAP_run/04_pairwise_sv \
  --gd ../genomes \
  --threads 12 \
  --seq
```

| Option | Meaning |
|---|---|
| `--synnet` / `--syn` | Single cleaned synnet file. |
| `--syndir` | Directory of cleaned synnet files. |
| `--outdir` | Output directory. |
| `--gd` | Genome FASTA directory. |
| `--ref_fa` | Reference FASTA for single-file mode. |
| `--qry_fa` / `--query_fa` | Query FASTA. |
| `--seq` | Write REF/ALT sequence into output. |
| `--exact_seq_len` | Try to output full sequences. |
| `--t` / `--threads` | Number of parallel jobs. |
| `--invlen` | Minimum inversion length. |
| `--cnvlen` | Minimum CNV length. |
| `--dellen` | Minimum deletion length. |
| `--inslen` | Minimum insertion length. |
| `--trllen` | Minimum translocation length. |
| `--filltoplen` | Top-fill length threshold. |
| `--fill2len` | Second-level fill length threshold. |
| `--rate` | Alignment-rate threshold. |
| `--max_level` | Maximum synnet parsing level. |
| `--seq_context` | Flanking context length for sequence output. |
| `--max_full_seq` | Maximum length for full sequence output. |
| `--edge_seq_len` | Edge sequence length when long sequences are truncated. |
| `--max_opposite_len` | Deprecated compatibility option. |
| `--max_opposite_ratio` | Maximum opposite-side gap ratio for simple DEL/INS classification. |
| `--keep_tmp` | Keep temporary files. |
| `--quiet` | Reduce screen output. |

### Step 05: 05_CombineSV.pl

**Purpose**: merge pairwise raw SVs into population-level candidate events and split them by SV type.

```bash
perl 05_CombineSV.pl \
  --sv ../SVGAP_run/04_pairwise_sv \
  --refname ISO1 \
  --outdir ../SVGAP_run/05_combined_sv \
  --sort_threads 12
```

| Option | Meaning |
|---|---|
| `--sv` / `--svdir` / `--pairsv` | Pairwise raw SV directory. |
| `--refname` | Reference ID. |
| `--outdir` | Output directory. |
| `--ins_tol` | Insertion merging tolerance. |
| `--del_ov` | Minimum deletion overlap for merging. |
| `--len_ratio` | Allowed length-ratio difference. |
| `--cant` / `--candidate` | Minimum candidate support. |
| `--min_large_len` | Minimum large-SV length. |
| `--max_indel_len` | Maximum small-InDel length. |
| `--legacy_indel` | Use legacy InDel output logic. |
| `--write_ins_reference` | Write reference insertion table. |
| `--cnv_ov` | CNV overlap threshold. |
| `--cnv_len_ratio` | CNV length-ratio threshold. |
| `--cnv_bp_tol` | CNV breakpoint tolerance in bp. |
| `--inv_ov` | Inversion overlap threshold. |
| `--inv_len_ratio` | Inversion length-ratio threshold. |
| `--inv_bp_tol` | Inversion breakpoint tolerance in bp. |
| `--trl_ov` | TRL overlap threshold. |
| `--trl_len_ratio` | TRL length-ratio threshold. |
| `--trl_bp_tol` | TRL breakpoint tolerance in bp. |
| `--trl_require_qchr` | Require query chromosome match for TRL merging. |
| `--sort_mem` | Memory for sort, for example `2G`. |
| `--sort_threads` | Number of sort threads. |
| `--sort_bin` | Path to sort binary. |
| `--keep_temp` | Keep temporary files. |
| `--verbose` | More verbose logging. |

### Step 06: 06_genotype_DEL.pl

**Purpose**: genotype merged deletions using single-coverage MAF files.

```bash
perl 06_genotype_DEL.pl \
  --bed ../SVGAP_run/05_combined_sv/All.DELs.50bplarge.bed.combined.sorted.txt \
  --maf_dir ../SVGAP_run/03_filter_forward/filtered_sing \
  --refname ISO1 \
  --query_list ../SVGAP_run/query.list \
  --out ../SVGAP_run/06_genotype/DEL.large \
  --threads 12
```

| Option | Meaning |
|---|---|
| `--bed` | Merged DEL BED file. |
| `--maf_dir` | Single-coverage MAF directory. |
| `--target` / `--refname` | Reference ID. |
| `--query_list` / `--glist` | Query list; controls VCF sample order. |
| `--out` | Output prefix. |
| `--threads` / `--t` | Number of parallel query jobs. |
| `--maf_suffix` | MAF suffix, usually `.sing.maf`. |
| `--min_cov` | Minimum reference coverage. |
| `--min_support` | Minimum alternative support. |
| `--bridge_del_ratio` | Ratio threshold for bridge-based deletion calls. |
| `--bridge_ref_low` | Lower bound of bridge reference coverage. |
| `--bridge_ref_high` | Upper bound of bridge reference coverage. |
| `--vcf_mode` | VCF output mode. |
| `--haploid` / `--diploid` | Genotype ploidy. |
| `--keep_tmp` | Keep temporary per-sample files. |

### Step 07: 07_genotype_INS.pl

**Purpose**: genotype merged insertions using MAF context, genome FASTA and sequence/flank matching.

```bash
perl 07_genotype_INS.pl \
  --bed ../SVGAP_run/05_combined_sv/All.INTs.50bplarge.bed.combined.sorted.txt \
  --maf_dir ../SVGAP_run/03_filter_forward/filtered_sing \
  --genome_dir ../genomes \
  --refname ISO1 \
  --query_list ../SVGAP_run/query.list \
  --out ../SVGAP_run/06_genotype/INS.large \
  --threads 12
```

| Option | Meaning |
|---|---|
| `--bed` | Merged INS BED file. |
| `--maf_dir` | Single-coverage MAF directory. |
| `--genome_dir` | FASTA directory for sequence extraction/verification. |
| `--target` / `--refname` | Reference ID. |
| `--query_list` | Query list. |
| `--out` | Output prefix. |
| `--threads` / `--t` | Number of parallel query jobs. |
| `--maf_suffix` | MAF file suffix. |
| `--pos_tol` | Position tolerance in bp. |
| `--seq_id` | Sequence identity threshold. |
| `--flank` | Flank length for matching. |
| `--len_ratio_min` | Lower bound of ALT length ratio. |
| `--len_ratio_max` | Upper bound of ALT length ratio. |
| `--len_slop` | Absolute length tolerance in bp. |
| `--min_refcov` | Reference-coverage threshold. |
| `--block_extra` | Extra length for MAF block search. |
| `--short_ref_len` | Short reference-anchor threshold. |
| `--long_seq_len` | Long insertion threshold. |
| `--sample_len` | Window length for long-sequence sampling. |
| `--sample_seed_len` | Seed length for sampling. |
| `--sample_min_hits` | Minimum seed hits. |
| `--sample_windows` | Number of sampling windows. |
| `--vcf_mode` | VCF output mode. |
| `--haploid` / `--diploid` | Genotype ploidy. |
| `--bed_one_based` | Treat input BED as 1-based. |
| `--alt_is_insseq` | Treat ALT as the inserted sequence itself. |
| `--keep_tmp` | Keep temporary files. |

### Step 08: 08_genotype_small_indels.pl

**Purpose**: genotype 1-50 bp small insertions and deletions. It can read small DEL and small INS files together.

```bash
perl 08_genotype_small_indels.pl \
  --maf_dir ../SVGAP_run/03_filter_forward/filtered_sing \
  --refname ISO1 \
  --query_list ../SVGAP_run/query.list \
  --del_bed ../SVGAP_run/05_combined_sv/All.DELs.50bpsmall.bed.combined.sorted.txt \
  --ins_bed ../SVGAP_run/05_combined_sv/All.INTs.50bpsmall.bed.combined.sorted.txt \
  --out ../SVGAP_run/06_genotype/smallIndel \
  --threads 12
```

| Option | Meaning |
|---|---|
| `--maf_dir` | Single-coverage MAF directory. |
| `--target` / `--refname` | Reference ID. |
| `--bed` / `--variants` | Mixed small-InDel BED file. |
| `--del_bed` | Small deletion BED file. |
| `--ins_bed` | Small insertion BED file. |
| `--out` | Output prefix. Default: `small_indel_genotype`. |
| `--max_len` | Maximum InDel length. Default: 50. |
| `--min_len` | Minimum InDel length. Default: 1. |
| `--maf_suffix` | MAF suffix. Default: `.sing.maf`. |
| `--query_list` / `--glist` | Query list. |
| `--other_gt` | GT for non-reference/non-alternative alleles. Default: `./.`. |
| `--use_block_gap` / `--no-use_block_gap` | Infer genotypes from adjacent MAF block gaps. Default: on. |
| `--max_junction_gap` | Maximum junction-gap length. Default: 1000. |
| `--threads` / `--p` | Number of parallel query jobs. |
| `--keep_tmp` | Keep temporary files. |
| `--bed_coord` | Coordinate system: 0 for 0-based MAF/BED, 1 for 1-based. |
| `--gt_pos_slop` | Position tolerance for relaxed genotyping. Default: 5. |
| `--len_slop` | Length tolerance for relaxed genotyping. Default: 0. |
| `--relaxed_gt` / `--no-relaxed_gt` | Search nearby MAF gaps before calling reference. Default: on. |

### Step 09: 09_SNPgenotyping.pl

**Purpose**: call SNPs from single-coverage MAF and merge all query SNP sites into one population VCF.

```bash
perl 09_SNPgenotyping.pl \
  --singmaf ../SVGAP_run/03_filter_forward/filtered_sing \
  --refsize ../genomes/ISO1.sizes \
  --refname ISO1 \
  --query_list ../SVGAP_run/query.list \
  --speciesname ISO1_A3_A4_A5 \
  --outdir ../SVGAP_run/07_snp \
  --threads 12
```

| Option | Meaning |
|---|---|
| `--singmaf` | Single-coverage MAF directory or file. |
| `--refsize` | Reference chromosome size file. |
| `--refseq` | Kept for compatibility; usually not required. |
| `--refname` / `--target` | Reference ID. |
| `--query_list` / `--glist` | Query list controlling VCF sample order. |
| `--speciesname` | Project/species name for VCF header. |
| `--outdir` | Output directory. |
| `--t` / `--threads` | Number of parallel jobs. |
| `--tmpdir` | Temporary directory. |
| `--keep_tmp` | Keep temporary files. |
| `--allow_ambiguous` | Allow ambiguous bases. |
| `--minNS` | N/ambiguous filter threshold. |
| `--phased` | Output phased genotypes. |
| `--legacy_format` | Use legacy format. |
| `--force` | Overwrite existing outputs. |

### Step 10: 10_SVannotation.pl

**Purpose**: annotate SV VCF records, including tandem duplication, TE overlap and CDS overlap. For symbolic DEL alleles, the deleted sequence can be extracted from reference coordinates.

```bash
perl 10_SVannotation.pl \
  --ref ../genomes/ISO1.fa \
  --vcf ../SVGAP_run/08_vcf/All.largeSV.raw.vcf \
  --te ../annotation/TE.bed \
  --cds ../annotation/CDS.bed \
  --out ../SVGAP_run/08_vcf/All.largeSV.annotated.vcf \
  --workdir ../SVGAP_run/09_annotation \
  --threads 12 \
  --minimap2 minimap2
```

If minimap2 is unavailable, skip annotation in the main pipeline:

```bash
--skip_annotation
```

or skip tandem-duplication detection when running annotation alone:

```bash
--skip_tandem_dup
```

| Option | Meaning |
|---|---|
| `--ref` | Reference FASTA. |
| `--vcf` | Input VCF. |
| `--te` | TE annotation file. |
| `--cds` | CDS annotation file. |
| `--out` | Output VCF. |
| `--svtype` | Annotate only a selected SV type. |
| `--workdir` | Working directory. |
| `--threads` / `--t` | Number of threads. |
| `--minimap2` | Path to minimap2. |
| `--kmer` | K-mer size for tandem-duplication seeds. |
| `--dup_seed_fraction` | Seed hit fraction for tandem duplication. |
| `--skip_tandem_dup` | Skip tandem-duplication detection. |
| `--min_annot_pct` | Minimum TE/CDS overlap percentage. |
| `--min_aln_len` | Minimum minimap2 alignment length. |
| `--min_mapq` | Minimum minimap2 MAPQ. |
| `--max_dup_flank` | Maximum flank length for tandem-duplication checking. |
| `--max_dup_seeds` | Maximum number of duplication seeds. |
| `--max_seq_len` | Maximum sequence length used for annotation. |
| `--te_preset` | TE input preset. |
| `--cds_preset` | CDS input preset. |
| `--keep_tmp` | Keep temporary files. |
| `--quiet` | Reduce output. |

### Step 11: 11_INVgenotyping.pl

**Purpose**: genotype merged inversion candidates using MAF body coverage and flank coverage.

```bash
perl 11_INVgenotyping.pl \
  --inv ../SVGAP_run/05_combined_sv/All.INVs.bed.combined.sorted.txt \
  --maf_dir ../SVGAP_run/03_filter_forward/filtered_sing \
  --refsize ../genomes/ISO1.sizes \
  --refname ISO1 \
  --speciesname ISO1_A3_A4_A5 \
  --query_list ../SVGAP_run/query.list \
  --out ../SVGAP_run/06_genotype/INV.vcf \
  --threads 12
```

| Option | Meaning |
|---|---|
| `--inv` | Inversion BED file. |
| `--maf_dir` | Single-coverage MAF directory. |
| `--sing` / `--singmaf` | Single MAF input or directory. |
| `--refsize` | Reference size file. |
| `--refname` / `--target` | Reference ID. |
| `--speciesname` | Project name for VCF header. |
| `--query_list` / `--glist` | Query list. |
| `--out` | Output VCF. |
| `--outdir` | Output directory. |
| `--maf_suffix` | MAF suffix. |
| `--threads` / `--t` | Number of parallel query jobs. |
| `--bed_coord` | BED coordinate system. |
| `--flank_ratio` | Flank length as a fraction of inversion length. |
| `--min_flank` | Minimum flank length. |
| `--max_flank` | Maximum flank length. |
| `--min_body_cov` | Minimum body coverage. |
| `--min_flank_cov` | Minimum flank coverage. |
| `--dominant_frac` | Minimum dominant-orientation fraction. |
| `--min_body_bp` | Minimum body coverage in bp. |
| `--trust_carriers` | Trust candidate carrier samples. |
| `--body_only` | Use body evidence only. |
| `--phased` | Output phased genotypes. |
| `--keep_tmp` | Keep temporary files. |
| `--debug` | Output debug information. |

### Step 12: 12_combine_genotype_CXSV.pl

**Purpose**: merge and genotype complex SVs where reference and query gaps coexist and cannot be represented well as simple DEL/INS.

```bash
perl 12_combine_genotype_CXSV.pl \
  --refname ISO1 \
  --query_list ../SVGAP_run/query.list \
  --rawsv_dir ../SVGAP_run/04_pairwise_sv \
  --outdir ../SVGAP_run/06_genotype/CXSV \
  --break_tolerance 2000 \
  --min_ro 0.50 \
  --len_ratio 1.00
```

| Option | Meaning |
|---|---|
| `--refname` | Reference ID. |
| `--glist` | Genome list. |
| `--query_list` | Query list; preferred for sample order. |
| `--rawsv_dir` / `--svdir` | Step 04 raw SV directory. |
| `--outdir` | Output directory. |
| `--break_tolerance` | Breakpoint merging tolerance. Default: 2000. |
| `--min_ro` | Minimum reciprocal overlap. Default: 0.50. |
| `--len_ratio` | Maximum relative length difference. Default: 1.00. |
| `--min_component_size` | Minimum number of component fragments. |
| `--unknown_missing` | Output `./.` instead of `0|0` for non-carriers. |
| `--keep_ref_sample` | Keep reference sample in VCF. |
| `--max_info_source_ids` | Maximum number of source IDs kept in INFO. |

### SVGAP_pipeline.pl controller options

| Option | Meaning |
|---|---|
| `--refname` / `--target` | Reference genome ID; must appear in `g.list`. |
| `--glist` | Genome list; VCF sample order is the list after removing the reference. |
| `--genome_dir` / `--gd` | Genome FASTA directory. |
| `--start_from` | Start point: `align`, `maf`, or `axt`. Default: `align`. |
| `--align_dir` / `--input` | Aligner output directory. |
| `--aligner` / `--ali` | Aligner type: `minimap2`, `last`, `lastz`, `mummer`, `anchorwave`, `maf`. |
| `--maf_dir` | Forward MAF directory. |
| `--maf_reciprocal_dir` | True independent reciprocal MAF directory. |
| `--axt_dir` / `--ad` | Forward AXT directory. |
| `--axt_forward_dir` | Forward AXT directory. |
| `--axt_reciprocal_dir` | True independent reciprocal AXT directory. |
| `--reciprocal_strategy` | `auto`, `qnet`, `true`, or `both`. Default: `auto`. |
| `--outdir` | Main output directory. Default: `SVGAP_run`. |
| `--threads` / `--t` | Number of threads/jobs. |
| `--speciesname` | Project/species name in VCF headers. |
| `--ucsc` | UCSC utility path. |
| `--tba` | TBA utility path. |
| `--k8` | k8 path. |
| `--paftools` | paftools.js path. |
| `--mum320` | MUMmer/delta2maf path. |
| `--linearGap` / `--lg` | ChainNet linearGap file. |
| `--refsize` | Reference size file. |
| `--te` | TE annotation file. |
| `--cds` | CDS annotation file. |
| `--minimap2` | minimap2 path. |
| `--annotation` / `--no-annotation` | Enable/disable annotation. |
| `--skip_convert` | Skip Step 01. |
| `--skip_chainnet` | Skip Step 02. |
| `--skip_filter` | Skip Step 03. |
| `--skip_pairwise_sv` | Skip Step 04. |
| `--skip_combine` | Skip Step 05. |
| `--skip_genotype` | Skip DEL/INS/small-InDel/INV genotyping. |
| `--skip_cxsv` | Skip Step 12. |
| `--skip_snp` | Skip Step 09. |
| `--skip_annotation` | Skip Step 10. |
| `--cxsv_break_tolerance` | CXSV breakpoint merging tolerance. |
| `--cxsv_min_ro` | CXSV minimum reciprocal overlap. |
| `--cxsv_len_ratio` | CXSV length-ratio threshold. |
| `--cxsv_unknown_missing` | Output missing genotypes for CXSV non-carriers. |
| `--dry_run` | Print commands without running. |
| `--force` | Overwrite existing outputs. |

### SVGAP_vcf_tools.pl

**Purpose**: merge multiple VCFs and enforce sample order according to `query.list` or `g.list`.

```bash
perl SVGAP_vcf_tools.pl \
  --mode merge \
  --vcfs DEL.vcf,INS.vcf,INV.vcf \
  --query_list query.list \
  --refname ISO1 \
  --refsize ../genomes/ISO1.sizes \
  --out All.largeSV.raw.vcf
```

| Option | Meaning |
|---|---|
| `--mode` | Current main mode is `merge`. |
| `--vcfs` | Comma-separated VCF paths. |
| `--out` | Output VCF. |
| `--query_list` / `--glist` | File controlling sample order. |
| `--refname` / `--target` | Reference ID. |
| `--refsize` | Chromosome size file for sorting. |
| `--missing_gt` | Missing genotype. Default: `./.`. |

## 7. Output sanity checks

After the run, inspect:

```bash
ls -lh SVGAP_run/08_vcf/
less SVGAP_run/SVGAP_pipeline.log
```

Check VCF sample columns:

```bash
grep '^#CHROM' SVGAP_run/08_vcf/All.largeSV.raw.vcf
```

Expected order:

```text
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT A3 A4 A5
```

Count records:

```bash
grep -v '^#' SVGAP_run/08_vcf/All.largeSV.raw.vcf | wc -l
grep -v '^#' SVGAP_run/08_vcf/All.smallIndel.vcf | wc -l
grep -v '^#' SVGAP_run/08_vcf/All.sorted.SNP.vcf | wc -l
grep -v '^#' SVGAP_run/08_vcf/All.CXSV.vcf | wc -l
```

## 8. Troubleshooting

### 8.1 Parallel::ForkManager warning still appears

First confirm you are running the clean package, not an old copy:

```bash
pwd
head -5 SVGAP_pipeline.pl
```

The script should contain:

```perl
# Author: Yi Liao @ scau.edu.cn
```

The current production scripts no longer load `Parallel::ForkManager`. Run:

```bash
perl util/Test_dependency.pl
```

If the pipeline still prints that warning, the command is calling an older script path.

### 8.2 AXT files are not detected

AXT files must be non-empty and named:

```text
ISO1vsA3.axt
ISO1vsA4.axt
ISO1vsA5.axt
```

Zero-byte files are ignored.

### 8.3 Genome FASTA not found

The ID in `g.list` must match the FASTA prefix. For `A3`, one of the following must exist:

```text
A3
A3.fa
A3.fasta
A3.fna
```

### 8.4 minimap2 is missing

The core SV/SNP/InDel workflow does not always require system minimap2, but annotation can require it for tandem-duplication checks. If minimap2 is unavailable, use:

```bash
--skip_annotation
```

or, for annotation alone:

```bash
--skip_tandem_dup
```

### 8.5 Memory and threads

Higher `--threads` can speed up query-level parallelism but also increases memory use. For three queries, `--threads 3` to `--threads 12` is usually reasonable. For large genomes, start with `--threads 4` or `--threads 6`.

### 8.6 qnet versus reciprocal alignments

If no true query-vs-reference alignment exists, keep the default `auto` behavior. SVGAP will use forward chainNet tnet/qnet to recover reference insertion information. If true independent `A3vsISO1.axt`, `A4vsISO1.axt` and `A5vsISO1.axt` files exist, use `--reciprocal_strategy true`.

## 9. Recommended training workflow

For a training session:

1. Show `g.list` and genome FASTA naming.
2. Show AXT naming: `ISO1vsA3.axt`.
3. Run the controller with `--dry_run` so trainees can see each command.
4. Run the full pipeline first with `--skip_annotation`.
5. Focus on the four outputs in `08_vcf/`: large SV, small InDel, SNP and CXSV.
6. Explain how `query.list` controls VCF sample order.
7. Demonstrate annotation and benchmark tools separately.

Dry-run example:

```bash
perl SVGAP_pipeline.pl \
  --refname ISO1 \
  --glist ../g.list \
  --genome_dir ../genomes \
  --start_from axt \
  --axt_dir ../01_axt/Target \
  --outdir ../SVGAP_run_dryrun \
  --threads 12 \
  --skip_annotation \
  --dry_run
```
