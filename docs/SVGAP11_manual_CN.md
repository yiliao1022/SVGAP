---
title: "SVGAP11-CXSV 清洁版中文使用手册"
author: "Yi Liao @ scau.edu.cn"
date: "2026-07-07"
geometry: margin=0.75in
mainfont: "Noto Serif CJK SC"
monofont: "Noto Sans Mono CJK SC"
fontsize: 10pt
colorlinks: true
---

# SVGAP version2 中文使用手册

版本约定：此前仓库中的旧流程称为 **version1**；本手册描述的修复优化流程称为 **version2**。

**作者标注**：本包中每个 Perl 脚本均已加入 `Author: Yi Liao @ scau.edu.cn`。

**适用场景**：以一个参考基因组和多个 query 基因组为输入，完成 assembly-to-assembly 的 SV、复杂 SV、SNP 和 small InDel 检测、合并、genotyping 与可选注释。
**示例数据设定**：`ISO1` 为 reference，`A3`、`A4`、`A5` 为 query genomes。

## 1. 本次清理和修复内容

### 1.1 无额外模块的并行执行

你之前反复遇到的警告是：

```text
WARNING: Parallel::ForkManager is not installed; falling back to single-thread mode
```

新版的并行步骤统一使用 Perl 内置的 `fork()` 和 `wait()`，不再依赖
`Parallel::ForkManager`，也不会因为 CPAN 模块搜索路径不同而回退到单线程。
`--threads` 可直接控制同时运行的样本或文件任务数。

所有 `.pl` 脚本仍统一加入包内库路径，以便加载 SVGAP 自带模块：

```perl
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/pub";
use lib "$Bin/../lib";
use lib "$Bin/../pub";
use lib "$Bin/../../lib";
use lib "$Bin/../../pub";
```

可用下面的命令检查依赖和内置工具：

```bash
cd SVGAP11_clean_fixed
perl util/Test_dependency.pl
```

### 1.2 脚本作者信息

每个 Perl 脚本头部均加入：

```perl
# Author: Yi Liao @ scau.edu.cn
```

### 1.3 文件包清理

为保持生产流程简洁，下面文件没有删除，而是移动到 `trash/`：

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

这些文件不是正式运行 ISO1/A3/A4/A5 流程所必需的。如果后续要做模拟 benchmark 或恢复旧工具，可以从 `trash/` 取回。

## 2. 清洁版包结构

```text
SVGAP11_clean_fixed/
  SVGAP_pipeline.pl                 # 主控脚本，一键运行全流程
  01_Convert2Axt.pl                 # Step 01: alignment 转 AXT
  02_ChainNetSyn.pl                 # Step 02: AXT -> chain/net/synnet
  03_SynNetFilter.pl                # Step 03: synnet 过滤 + single-coverage MAF
  04_PairwiseSV.pl                  # Step 04: pairwise SV calling
  05_CombineSV.pl                   # Step 05: SV 合并
  06_genotype_DEL.pl                # Step 06: deletion genotyping
  07_genotype_INS.pl                # Step 07: insertion genotyping
  08_genotype_small_indels.pl       # Step 08: small InDel genotyping
  09_SNPgenotyping.pl               # Step 09: SNP genotyping
  10_SVannotation.pl                # Step 10: SV annotation
  11_INVgenotyping.pl               # Step 11: inversion genotyping
  12_combine_genotype_CXSV.pl       # Step 12: complex SV 合并和 genotyping
  SVGAP_vcf_tools.pl                # VCF 合并/样本顺序整理工具
  lib/                              # Perl 模块依赖，包括 Parallel::ForkManager
  pub/                              # UCSC/TBA/MUMMER/Stretcher 等工具或辅助脚本
  util/Test_dependency.pl           # 依赖测试脚本
  docs/                             # 本手册
  trash/                            # 非生产流程必需文件的归档目录
```

## 3. 输入文件准备：ISO1 reference + A3/A4/A5 queries

### 3.1 推荐目录

```text
project/
  SVGAP11_clean_fixed/
  genomes/
    ISO1.fa
    A3.fa
    A4.fa
    A5.fa
  g.list
  01_axt/Target/                 # 如果从 AXT 开始
    ISO1vsA3.axt
    ISO1vsA4.axt
    ISO1vsA5.axt
```

### 3.2 g.list 文件

`g.list` 中每行一个 genome ID。reference 也要写进去，主流程会自动把 reference 从 VCF 样本列中排除。因此最终 VCF 样本顺序为 `A3, A4, A5`。

```text
ISO1
A3
A4
A5
```

### 3.3 genome FASTA 命名

每个 ID 可以对应以下任一文件名：

```text
ID
ID.fa
ID.fasta
ID.fna
```

在本例中推荐使用：

```text
genomes/ISO1.fa
genomes/A3.fa
genomes/A4.fa
genomes/A5.fa
```

FASTA header 推荐带样本前缀，便于后续结果解释：

```text
>ISO1.Chr01
>A3.Chr01
>A4.Chr01
>A5.Chr01
```

不建议 header 以数字开头。

### 3.4 AXT 命名

如果直接从 AXT 开始，forward AXT 文件名必须是：

```text
ISO1vsA3.axt
ISO1vsA4.axt
ISO1vsA5.axt
```

如果有真正独立生成的 reciprocal alignment，可以准备：

```text
A3vsISO1.axt
A4vsISO1.axt
A5vsISO1.axt
```

注意：不要简单反转 forward AXT 来伪造 reciprocal AXT。那只是坐标重排，不是真正的 query-vs-reference alignment。

## 4. 一键运行主流程

### 4.1 推荐方式：从已有 AXT 开始

如果你已经有 `ISO1vsA3.axt`、`ISO1vsA4.axt`、`ISO1vsA5.axt`，建议这样运行：

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

这里加 `--skip_annotation` 是因为注释步骤需要 TE/CDS 输入或 minimap2 环境。第一次培训演示时建议先跳过 annotation，确保主体检测流程完成。

### 4.2 从 MAF 开始

如果已有 forward MAF，文件名推荐为：

```text
ISO1vsA3.maf
ISO1vsA4.maf
ISO1vsA5.maf
```

运行：

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

### 4.3 从 aligner 原始输出开始

如果从 minimap2/LAST/LASTZ/MUMmer/AnchorWave 的输出开始，由 `01_Convert2Axt.pl` 统一转换为 AXT。

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

### 4.4 使用真正 reciprocal alignment

如果你独立运行了 query-vs-reference alignment，可用：

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

如果没有 reciprocal alignment，默认使用 `qnet` 策略即可。

## 5. 主流程输出目录解释

主控脚本会在 `--outdir` 中生成：

```text
SVGAP_run/
  01_axt/                         # AXT 准备或转换结果
  02_chain_forward/               # forward chain/net/synnet
  02_chain_independent_reciprocal/ # 可选 reciprocal chain/net/synnet
  03_filter_forward/              # filtered synnet + single-coverage MAF
  03_filter_independent_reciprocal/
  04_pairwise_inputs/             # PairwiseSV 输入 synnet 子集
  04_pairwise_sv/                 # 每个 pair 的 raw SV
  05_combined_sv/                 # 合并后的 DEL/INS/INV/CNV/TRL 候选
  06_genotype/                    # DEL/INS/small InDel/INV/CXSV genotyping
  07_snp/                         # SNP calling/genotyping
  08_vcf/                         # 最终 VCF 输出
  09_annotation/                  # 可选 annotation 工作目录
  SVGAP_pipeline.log              # 主流程日志
```

最重要的最终文件通常在：

```text
SVGAP_run/08_vcf/All.largeSV.raw.vcf
SVGAP_run/08_vcf/All.largeSV.annotated.vcf     # 如果 annotation 开启
SVGAP_run/08_vcf/All.smallIndel.vcf
SVGAP_run/08_vcf/All.sorted.SNP.vcf
SVGAP_run/08_vcf/All.CXSV.vcf
```

## 6. 流程步骤详细说明

### Step 01：01_Convert2Axt.pl

**功能**：把 LAST、LASTZ、MUMmer、minimap2、AnchorWave 或已有 MAF 转换为 AXT。也可把 query FASTA 切分为 overlapping subcontigs，并把 split contig 坐标恢复到原始 chromosome 坐标。

**常用命令**：

```bash
perl 01_Convert2Axt.pl \
  --ali minimap2 \
  --input ../alignments \
  --wk ../01_axt \
  --tname ISO1 \
  --qname A3 \
  --t 8
```

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--ali` | 输入 alignment 类型：`last`、`lastz`、`maf`、`mummer`、`minimap2`、`anchorwave`。|
| `--input` | alignment 文件所在目录。|
| `--wk` | 工作目录/输出目录。|
| `--tname` | target/reference 名称，默认 `Ref`。|
| `--qname` | query 名称，默认 `Que`。|
| `--mode` / `--coord` | `auto`、`chr` 或 `ctg`。`ctg` 用于 split query subcontig 坐标恢复。|
| `--t` / `--threads` | 并行进程数，默认 4；使用 Perl 内置 `fork()/wait()`。|
| `--query` | 同时生成 query-as-reference AXT。|
| `--keep_tmp` | 保留中间文件。|
| `--mum320` | MUMmer `delta2maf` 所在目录或路径。|
| `--ucsc` | UCSC 工具目录，例如 `mafToAxt`。|
| `--tba` | TBA 工具目录。|
| `--k8` | `k8` 路径，用于 paftools。|
| `--paftools` | `paftools.js` 路径。|
| `--splitfa` / `--fasta` | 需要切分的 FASTA，可重复使用。|
| `--bin` | split 主窗口大小，默认 10,000,000 bp。|
| `--slide` | split 窗口重叠大小，默认 10,000 bp。|
| `--split_out` | split FASTA 输出目录。|
| `--sizes` / `--qsize` / `--size` | 原始 chromosome size 文件，用于 ctg 模式坐标恢复。|

### Step 02：02_ChainNetSyn.pl

**功能**：从 AXT 构建 UCSC chain/net/syntenic-net，并可生成 target-oriented single-coverage MAF。

**常用命令**：

```bash
perl 02_ChainNetSyn.pl \
  --gd ../genomes \
  --ad ../01_axt/Target \
  --lst ../g.list \
  --outdir ../SVGAP_run/02_chain_forward \
  --threads 12 \
  --sing
```

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--gd` | genome FASTA 目录。脚本会寻找 `ID`、`ID.fa`、`ID.fasta`、`ID.fna`。|
| `--ad` | AXT 文件目录，文件名必须类似 `TargetIDvsQueryID.axt`。|
| `--lst` | genome ID list。|
| `--outdir` / `--syn` | 输出目录。|
| `--ucsc` | UCSC 工具目录。|
| `--tba` | TBA 工具目录。|
| `--linearGap` / `--lg` | chainNet linearGap 文件，默认 `pub/TBA/medium`。|
| `--minscore` / `--minScore` | axtChain/chainNet 最低分数，默认 1000。|
| `--minSpace` | chainNet minSpace，默认 1。|
| `--t` / `--threads` | 并行任务数，默认 1。|
| `--sing` | 生成 single-coverage MAF。|
| `--keep_raw_chain` | 保留 raw chain 文件。|
| `--keep_temp` | 保留临时 AXT/MAF 文件。|
| `--force` | 覆盖已有结果。|
| `--dry_run` | 只打印命令，不真正执行。|

**主要输出**：

```text
chain/*.chain.filter
net/*.tnet, *.qnet
synnet/*.tnet.synnet, *.qnet.synnet
singRaw/*.sing.maf
```

### Step 03：03_SynNetFilter.pl

**功能**：过滤 synnet 中不可靠的短 fill，生成更保守的 CleanSynNet、RawSV gap 文件和 filtered single-coverage MAF。

**常用命令**：

```bash
perl 03_SynNetFilter.pl \
  --synnet ../SVGAP_run/02_chain_forward/synnet \
  --chain ../SVGAP_run/02_chain_forward/chain \
  --gd ../genomes \
  --outdir ../SVGAP_run/03_filter_forward \
  --threads 12
```

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--synnet` | Step 02 生成的 synnet 目录。|
| `--chain` | Step 02 生成的 chain 目录。|
| `--gd` | genome FASTA/.sizes/.2bit 目录。|
| `--outdir` | 输出目录。|
| `--ucsc` | UCSC 工具目录。|
| `--tba` | TBA 工具目录。|
| `--cleanSynnet` | CleanSynNet 输出子目录名，默认 `CleanSynNet`。|
| `--sing` | filtered single MAF 输出子目录名，默认 `filtered_sing`。|
| `--axt` | cleaned synnet 转 AXT 的输出子目录名。|
| `--homopairs` | homologous pair summary 输出子目录名。|
| `--rawsv` | raw SV gap 输出子目录名。|
| `--target` / `--refname` | 只处理指定 target/reference。|
| `--query` | 只处理指定 query。|
| `--toplen` | top syntenic fill 最小长度，默认 500000。|
| `--nonsynlen` | 下层 nonSyn fill 最小长度，默认 500000。|
| `--synlen` | 下层 syn fill 最小长度，默认 500000。|
| `--invlen` | inversion fill 最小长度，默认 1000。|
| `--alrate` | 最小 alignment/fill length 比例，默认 0.10。|
| `--qFarCutoff` | inversion fill 的最大 qFar，默认 5000000。|
| `--minscore` | 最低分数阈值。|
| `--t` / `--threads` | 并行任务数。|
| `--sing_only_tnet` | 只从 tnet 生成 single MAF。|
| `--keep_temp` | 保留临时文件。|
| `--force` | 覆盖已有输出。|
| `--dry_run` | 只打印命令。|

### Step 04：04_PairwiseSV.pl

**功能**：从 cleaned synnet 直接识别 pairwise SV，包括 DEL、INS、CNV、INV、TRL 和复杂 PAV/CMPLX 事件。新版使用 Perl 内置 `fork()/waitpid()`，不依赖 `Parallel::ForkManager`。

**常用命令**：

```bash
perl 04_PairwiseSV.pl \
  --syndir ../SVGAP_run/03_filter_forward/CleanSynNet \
  --outdir ../SVGAP_run/04_pairwise_sv \
  --gd ../genomes \
  --threads 12 \
  --seq
```

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--synnet` / `--syn` | 单个 cleaned synnet 文件。|
| `--syndir` | cleaned synnet 文件目录，适合批量运行。|
| `--outdir` | 输出目录。|
| `--gd` | genome FASTA 目录。|
| `--ref_fa` | reference FASTA，单文件运行时可显式指定。|
| `--qry_fa` / `--query_fa` | query FASTA。|
| `--seq` | 在输出中写入 REF/ALT 序列。|
| `--exact_seq_len` | 尽量输出完整序列。|
| `--t` / `--threads` | 并行任务数。|
| `--invlen` | INV 最小长度。|
| `--cnvlen` | CNV 最小长度。|
| `--dellen` | DEL 最小长度。|
| `--inslen` | INS 最小长度。|
| `--trllen` | translocation/TRL 最小长度。|
| `--filltoplen` | top fill 长度阈值。|
| `--fill2len` | second-level fill 长度阈值。|
| `--rate` | alignment rate 阈值。|
| `--max_level` | synnet 解析最大层级。|
| `--seq_context` | 输出序列时保留的 flank context。|
| `--max_full_seq` | 完整输出序列的最大长度。|
| `--edge_seq_len` | 长序列截断时两端保留长度。|
| `--max_opposite_len` | 旧参数，保留兼容。|
| `--max_opposite_ratio` | DEL/INS 判定中另一侧 gap 最大比例，默认建议 0.10。|
| `--keep_tmp` | 保留临时文件。|
| `--quiet` | 减少屏幕输出。|

### Step 05：05_CombineSV.pl

**功能**：把多个 pairwise raw SV 合并为 population-level 候选事件，并按 SV 类型输出。

**常用命令**：

```bash
perl 05_CombineSV.pl \
  --sv ../SVGAP_run/04_pairwise_sv \
  --refname ISO1 \
  --outdir ../SVGAP_run/05_combined_sv \
  --sort_threads 12
```

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--sv` / `--svdir` / `--pairsv` | Step 04 raw SV 目录。|
| `--refname` | reference ID。|
| `--outdir` | 输出目录。|
| `--ins_tol` | insertion 合并位置容忍度/比例。|
| `--del_ov` | deletion 合并最小 overlap。|
| `--len_ratio` | 合并时允许的长度比例差异。|
| `--cant` / `--candidate` | 最小候选支持数。|
| `--min_large_len` | large SV 最小长度。|
| `--max_indel_len` | small InDel 最大长度。|
| `--legacy_indel` | 使用旧的 indel 输出逻辑。|
| `--write_ins_reference` | 输出 reference insertion 表。|
| `--cnv_ov` | CNV 合并 overlap。|
| `--cnv_len_ratio` | CNV 长度比例阈值。|
| `--cnv_bp_tol` | CNV breakpoint 容忍 bp。|
| `--inv_ov` | INV 合并 overlap。|
| `--inv_len_ratio` | INV 长度比例阈值。|
| `--inv_bp_tol` | INV breakpoint 容忍 bp。|
| `--trl_ov` | TRL overlap 阈值。|
| `--trl_len_ratio` | TRL 长度比例阈值。|
| `--trl_bp_tol` | TRL breakpoint 容忍 bp。|
| `--trl_require_qchr` | TRL 合并时要求 query chromosome 匹配。|
| `--sort_mem` | sort 内存，例如 `2G`。|
| `--sort_threads` | sort 使用线程数。|
| `--sort_bin` | sort 程序路径。|
| `--keep_temp` | 保留临时文件。|
| `--verbose` | 输出更详细日志。|

### Step 06：06_genotype_DEL.pl

**功能**：基于 single-coverage MAF 对 merged deletions 进行 population genotyping，输出 matrix 和 VCF。

**常用命令**：

```bash
perl 06_genotype_DEL.pl \
  --bed ../SVGAP_run/05_combined_sv/All.DELs.50bplarge.bed.combined.sorted.txt \
  --maf_dir ../SVGAP_run/03_filter_forward/filtered_sing \
  --refname ISO1 \
  --query_list ../SVGAP_run/query.list \
  --out ../SVGAP_run/06_genotype/DEL.large \
  --threads 12
```

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--bed` | merged DEL BED 文件。|
| `--maf_dir` | single-coverage MAF 目录。|
| `--target` / `--refname` | reference ID。|
| `--query_list` / `--glist` | query 列表；VCF 样本顺序按该文件。|
| `--out` | 输出前缀。|
| `--threads` / `--t` | 并行 query 数。|
| `--maf_suffix` | MAF 后缀，默认常为 `.sing.maf`。|
| `--min_cov` | 参考覆盖判定的最小覆盖比例。|
| `--min_support` | ALT 支持最小比例。|
| `--bridge_del_ratio` | bridge alignment 判定 deletion 的比例阈值。|
| `--bridge_ref_low` | bridge reference coverage 下界。|
| `--bridge_ref_high` | bridge reference coverage 上界。|
| `--vcf_mode` | VCF 输出模式。|
| `--haploid` / `--diploid` | 输出单倍体或二倍体 genotype。|
| `--keep_tmp` | 保留中间 per-sample 文件。|

### Step 07：07_genotype_INS.pl

**功能**：对 merged insertions 进行 genotyping。支持长插入序列抽样匹配、flank 匹配和 VCF 输出。

**常用命令**：

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

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--bed` | merged INS BED 文件。|
| `--maf_dir` | single-coverage MAF 目录。|
| `--genome_dir` | query/reference FASTA 目录，用于插入序列提取/验证。|
| `--target` / `--refname` | reference ID。|
| `--query_list` | query 列表。|
| `--out` | 输出前缀。|
| `--threads` / `--t` | 并行 query 数。|
| `--maf_suffix` | MAF 文件后缀。|
| `--pos_tol` | 插入位置容忍 bp。|
| `--seq_id` | 序列 identity 阈值。|
| `--flank` | flank 匹配长度。|
| `--len_ratio_min` | ALT 长度比例下限。|
| `--len_ratio_max` | ALT 长度比例上限。|
| `--len_slop` | 长度绝对容忍 bp。|
| `--min_refcov` | reference coverage 判定阈值。|
| `--block_extra` | MAF block 扩展搜索长度。|
| `--short_ref_len` | 短 reference anchor 长度阈值。|
| `--long_seq_len` | 长插入序列阈值。|
| `--sample_len` | 长序列抽样窗口长度。|
| `--sample_seed_len` | 抽样 seed 长度。|
| `--sample_min_hits` | seed 最小命中数。|
| `--sample_windows` | 抽样窗口数量。|
| `--vcf_mode` | VCF 输出模式。|
| `--haploid` / `--diploid` | genotype ploidy。|
| `--bed_one_based` | 输入 BED 是否为 1-based。|
| `--alt_is_insseq` | ALT 是否为插入序列本身。|
| `--keep_tmp` | 保留临时文件。|

### Step 08：08_genotype_small_indels.pl

**功能**：对 1-50 bp small insertions/deletions 进行 genotyping。该脚本可同时读取 small DEL 和 small INS 文件。

**常用命令**：

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

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--maf_dir` | single-coverage MAF 目录。|
| `--target` / `--refname` | reference ID。|
| `--bed` / `--variants` | mixed small InDel BED 文件。|
| `--del_bed` | small deletion BED 文件。|
| `--ins_bed` | small insertion BED 文件。|
| `--out` | 输出前缀，默认 `small_indel_genotype`。|
| `--max_len` | 最大 indel 长度，默认 50。|
| `--min_len` | 最小 indel 长度，默认 1。|
| `--maf_suffix` | MAF 后缀，默认 `.sing.maf`。|
| `--query_list` / `--glist` | query 列表。|
| `--other_gt` | 非 reference/alternative allele 的 GT，默认 `./.`。|
| `--use_block_gap` / `--no-use_block_gap` | 是否用相邻 MAF block gap 推断 genotype，默认开启。|
| `--max_junction_gap` | junction gap 最大长度，默认 1000。|
| `--threads` / `--p` | 并行 query 数。|
| `--keep_tmp` | 保留临时 per-sample 文件。|
| `--bed_coord` | 输入坐标系统：0 表示 MAF/BED 0-based，1 表示 1-based。|
| `--gt_pos_slop` | relaxed genotyping 的位置容忍度，默认 5。|
| `--len_slop` | relaxed 模式长度容忍度，默认 0。|
| `--relaxed_gt` / `--no-relaxed_gt` | 是否启用 nearby-gap 搜索，默认开启。|

### Step 09：09_SNPgenotyping.pl

**功能**：从 single-coverage MAF 中检测 SNP，并把所有 query 的 SNP 位点整合成统一 population VCF。

**常用命令**：

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

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--singmaf` | single-coverage MAF 目录或文件。|
| `--refsize` | reference chromosome sizes 文件。|
| `--refseq` | 兼容旧参数，新版通常不需要。|
| `--refname` / `--target` | reference ID。|
| `--query_list` / `--glist` | query 列表，控制 VCF 样本顺序。|
| `--speciesname` | VCF header 中的物种/项目名。|
| `--outdir` | 输出目录。|
| `--t` / `--threads` | 并行任务数。|
| `--tmpdir` | 临时目录。|
| `--keep_tmp` | 保留临时文件。|
| `--allow_ambiguous` | 是否允许 ambiguous bases。|
| `--minNS` | N/ambiguous 过滤阈值。|
| `--phased` | 输出 phased genotype。|
| `--legacy_format` | 使用旧格式。|
| `--force` | 覆盖已有输出。|

### Step 10：10_SVannotation.pl

**功能**：对 SV VCF 做注释，包括 tandem duplication、TE overlap、CDS overlap 等。DEL 为 symbolic allele 时，会根据 reference 坐标提取序列。

**常用命令**：

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

如果系统没有 `minimap2`，可以先跳过 annotation：

```bash
--skip_annotation
```

或者跳过 tandem duplication 检测：

```bash
--skip_tandem_dup
```

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--ref` | reference FASTA。|
| `--vcf` | 待注释 VCF。|
| `--te` | TE annotation BED/GFF。|
| `--cds` | CDS annotation BED/GFF。|
| `--out` | 输出 VCF。|
| `--svtype` | 只注释某类 SV。|
| `--workdir` | 工作目录。|
| `--threads` / `--t` | 并行线程数。|
| `--minimap2` | minimap2 程序路径。|
| `--kmer` | tandem duplication seed k-mer 大小。|
| `--dup_seed_fraction` | tandem duplication seed 命中比例阈值。|
| `--skip_tandem_dup` | 跳过 tandem duplication 判断。|
| `--min_annot_pct` | TE/CDS overlap 最小百分比。|
| `--min_aln_len` | minimap2 比对最小长度。|
| `--min_mapq` | minimap2 最小 MAPQ。|
| `--max_dup_flank` | tandem duplication flank 最大长度。|
| `--max_dup_seeds` | tandem duplication seed 最大数量。|
| `--max_seq_len` | 用于注释的最大序列长度。|
| `--te_preset` | TE 注释输入 preset。|
| `--cds_preset` | CDS 注释输入 preset。|
| `--keep_tmp` | 保留临时文件。|
| `--quiet` | 减少输出。|

### Step 11：11_INVgenotyping.pl

**功能**：对 merged inversion 候选进行 genotyping，结合 MAF body coverage 和 flank coverage 判断样本是否携带 inversion。

**常用命令**：

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

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--inv` | inversion BED 文件。|
| `--maf_dir` | single-coverage MAF 目录。|
| `--sing` / `--singmaf` | 单个或目录形式的 sing MAF 输入。|
| `--refsize` | reference sizes 文件。|
| `--refname` / `--target` | reference ID。|
| `--speciesname` | VCF header 项目名。|
| `--query_list` / `--glist` | query 列表。|
| `--out` | 输出 VCF 文件。|
| `--outdir` | 输出目录。|
| `--maf_suffix` | MAF 后缀。|
| `--threads` / `--t` | 并行 query 数。|
| `--bed_coord` | BED 坐标系统。|
| `--flank_ratio` | flank 长度相对 inversion 长度比例。|
| `--min_flank` | 最小 flank 长度。|
| `--max_flank` | 最大 flank 长度。|
| `--min_body_cov` | inversion body 最小覆盖比例。|
| `--min_flank_cov` | flank 最小覆盖比例。|
| `--dominant_frac` | dominant orientation 最小比例。|
| `--min_body_bp` | body 最小覆盖 bp。|
| `--trust_carriers` | 信任候选携带样本。|
| `--body_only` | 只根据 body 进行判定。|
| `--phased` | 输出 phased genotype。|
| `--keep_tmp` | 保留临时文件。|
| `--debug` | 输出 debug 信息。|

### Step 12：12_combine_genotype_CXSV.pl

**功能**：合并和 genotyping complex SV。适合处理 reference 与 query gap 同时存在、简单 DEL/INS 无法准确表达的复杂变异。

**常用命令**：

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

**主要参数**：

| 参数 | 说明 |
|---|---|
| `--refname` | reference ID。|
| `--glist` | genome list。|
| `--query_list` | query list，优先控制样本顺序。|
| `--rawsv_dir` / `--svdir` | Step 04 raw SV 目录。|
| `--outdir` | 输出目录。|
| `--break_tolerance` | 断点合并容忍度，默认 2000。|
| `--min_ro` | 最小 reciprocal overlap，默认 0.50。|
| `--len_ratio` | 最大相对长度差异，默认 1.00。|
| `--min_component_size` | complex event 最小组成片段数。|
| `--unknown_missing` | 非携带样本输出 `./.` 而不是 `0|0`。|
| `--keep_ref_sample` | VCF 中保留 reference sample。|
| `--max_info_source_ids` | INFO 中最多保留的 source ID 数。|

### SVGAP_pipeline.pl 主控脚本参数

| 参数 | 说明 |
|---|---|
| `--refname` / `--target` | reference genome ID，必须出现在 `g.list` 中。|
| `--glist` | genome list；VCF 样本顺序等于去掉 reference 后的顺序。|
| `--genome_dir` / `--gd` | genome FASTA 目录。|
| `--start_from` | 起点：`align`、`maf`、`axt`，默认 `align`。|
| `--align_dir` / `--input` | aligner 原始输出目录。|
| `--aligner` / `--ali` | aligner 类型：`minimap2`、`last`、`lastz`、`mummer`、`anchorwave`、`maf`。|
| `--maf_dir` | forward MAF 目录。|
| `--maf_reciprocal_dir` | 独立 reciprocal MAF 目录。|
| `--axt_dir` / `--ad` | forward AXT 目录。|
| `--axt_forward_dir` | forward AXT 目录。|
| `--axt_reciprocal_dir` | 独立 reciprocal AXT 目录。|
| `--reciprocal_strategy` | `auto`、`qnet`、`true`、`both`。默认 `auto`。|
| `--outdir` | 主输出目录，默认 `SVGAP_run`。|
| `--threads` / `--t` | 线程/并行任务数。|
| `--speciesname` | VCF header 项目名。|
| `--ucsc` | UCSC 工具路径。|
| `--tba` | TBA 工具路径。|
| `--k8` | k8 路径。|
| `--paftools` | paftools.js 路径。|
| `--mum320` | MUMmer/delta2maf 路径。|
| `--linearGap` / `--lg` | chainNet linearGap 文件。|
| `--refsize` | reference sizes 文件。|
| `--te` | TE annotation 文件。|
| `--cds` | CDS annotation 文件。|
| `--minimap2` | minimap2 路径。|
| `--annotation` / `--no-annotation` | 是否进行 annotation。|
| `--skip_convert` | 跳过 01 转换。|
| `--skip_chainnet` | 跳过 02。|
| `--skip_filter` | 跳过 03。|
| `--skip_pairwise_sv` | 跳过 04。|
| `--skip_combine` | 跳过 05。|
| `--skip_genotype` | 跳过 06/07/08/11。|
| `--skip_cxsv` | 跳过 12。|
| `--skip_snp` | 跳过 09。|
| `--skip_annotation` | 跳过 10。|
| `--cxsv_break_tolerance` | CXSV 断点合并容忍度。|
| `--cxsv_min_ro` | CXSV 最小 reciprocal overlap。|
| `--cxsv_len_ratio` | CXSV 长度比例阈值。|
| `--cxsv_unknown_missing` | CXSV 非携带样本输出 missing。|
| `--dry_run` | 只打印命令。|
| `--force` | 覆盖已有输出。|

### SVGAP_vcf_tools.pl

**功能**：合并多个 VCF，并强制样本列顺序与 `query.list` 或 `g.list` 一致。

```bash
perl SVGAP_vcf_tools.pl \
  --mode merge \
  --vcfs DEL.vcf,INS.vcf,INV.vcf \
  --query_list query.list \
  --refname ISO1 \
  --refsize ../genomes/ISO1.sizes \
  --out All.largeSV.raw.vcf
```

| 参数 | 说明 |
|---|---|
| `--mode` | 当前主要为 `merge`。|
| `--vcfs` | 逗号分隔的 VCF 文件。|
| `--out` | 输出 VCF。|
| `--query_list` / `--glist` | 样本顺序文件。|
| `--refname` / `--target` | reference ID。|
| `--refsize` | chromosome size 文件，用于排序。|
| `--missing_gt` | 缺失 genotype，默认 `./.`。|

## 7. 结果检查建议

运行完成后，建议检查：

```bash
ls -lh SVGAP_run/08_vcf/
less SVGAP_run/SVGAP_pipeline.log
```

检查 VCF 样本列：

```bash
grep '^#CHROM' SVGAP_run/08_vcf/All.largeSV.raw.vcf
```

应该看到：

```text
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT A3 A4 A5
```

统计每类变异数量：

```bash
grep -v '^#' SVGAP_run/08_vcf/All.largeSV.raw.vcf | wc -l
grep -v '^#' SVGAP_run/08_vcf/All.smallIndel.vcf | wc -l
grep -v '^#' SVGAP_run/08_vcf/All.sorted.SNP.vcf | wc -l
grep -v '^#' SVGAP_run/08_vcf/All.CXSV.vcf | wc -l
```

## 8. 常见问题

### 8.1 仍然出现 Parallel::ForkManager warning

先确认你运行的是清洁版包中的脚本，而不是旧目录中的脚本：

```bash
which perl
pwd
head -5 SVGAP_pipeline.pl
```

脚本头部应出现：

```perl
# Author: Yi Liao @ scau.edu.cn
```

新版生产脚本不再加载 `Parallel::ForkManager`。运行下面的检查：

```bash
perl util/Test_dependency.pl
```

如果流程仍出现该 warning，说明命令调用的是旧脚本路径。

### 8.2 AXT 文件找不到

AXT 必须非空，并且命名为：

```text
ISO1vsA3.axt
ISO1vsA4.axt
ISO1vsA5.axt
```

大小为 0 的空文件会被跳过。

### 8.3 genome FASTA 找不到

确认 `g.list` 中 ID 与 FASTA 文件名前缀一致。例如 `A3` 必须对应：

```text
A3
A3.fa
A3.fasta
A3.fna
```

### 8.4 minimap2 缺失

主流程检测 SV/SNP/InDel 不一定需要系统有 `minimap2`，但 annotation 中 tandem duplication 检查会用到。若没有 minimap2，可先加：

```bash
--skip_annotation
```

或在单独运行 `10_SVannotation.pl` 时加：

```bash
--skip_tandem_dup
```

### 8.5 内存和线程

`--threads` 越大，并行 query 越多，速度越快，但内存也会增加。对于 3 个 query，如 A3/A4/A5，`--threads 3` 到 `--threads 12` 都可以；如果单个 genome 很大，建议先用 `--threads 4` 或 `--threads 6` 测试。

### 8.6 qnet 和 reciprocal 怎么选

没有真正独立的 query-vs-reference alignment 时，使用默认 `auto` 即可，流程会使用 forward chainNet 的 `qnet` 恢复 reference insertion 信息。如果你确实独立跑了 `A3vsISO1.axt`、`A4vsISO1.axt`、`A5vsISO1.axt`，可以使用 `--reciprocal_strategy true`。

## 9. 培训演示建议流程

培训时推荐按以下顺序讲解：

1. 展示 `g.list` 和 genome FASTA 命名。
2. 展示 AXT 命名规则 `ISO1vsA3.axt`。
3. 用 `--dry_run` 运行主控脚本，让学员看到每一步命令。
4. 正式运行时先加 `--skip_annotation`。
5. 运行完成后重点讲 `08_vcf/` 中的四类输出：large SV、small InDel、SNP、CXSV。
6. 解释 `query.list` 如何控制 VCF 样本列顺序。
7. 最后再单独演示 annotation 和 benchmark 工具。

`--dry_run` 示例：

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
