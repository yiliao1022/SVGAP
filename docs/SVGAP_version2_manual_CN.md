# SVGAP version2 中文运行手册

## 1. 版本说明

- **version1**：此前 Git 仓库中的原始流程。
- **version2**：当前替换后的流程，包含 CXSV、vg-ready VCF、短 contig
  自适应过滤、空变异兼容、`--skip_filter` 修复和并行提速。

检查版本：

```bash
perl SVGAP_pipeline.pl --version
```

应输出：

```text
SVGAP version2
```

## 2. version2 的主要变化

1. 用 `SVGAP_pipeline.pl` 统一运行从比对结果到群体 VCF 的完整流程。
2. 支持大 DEL、INS、INV、复杂结构变异、small InDel 和 SNP。
3. 输出经过参考序列校验的 `All.SV.vg.vcf` 和 `All.SVGAP.vg.vcf`。
4. 修复原始 qnet 的方向识别，保留 `Query.Ref` 插入证据。
5. `--skip_filter` 会直接使用原始 SynNet 和 `singRaw` MAF。
6. 没有变异时会跳过相应 genotyper，不再因只有表头的文件而退出。
7. SynNet 阈值会根据 contig 长度自动缩放，默认也适用于短 contig。
8. 各样本任务使用 Perl 内置 `fork()/wait()` 并行，无需额外 CPAN 模块。
9. Minimap2 PAF 转换默认不再执行流程未使用的 `paftools.js call`。

## 3. 环境检查

进入仓库目录后运行：

```bash
perl util/Test_dependency.pl
```

正常输出应包含：

```text
SVGAP dependency check: OK
Parallel workers: core Perl fork/wait (no CPAN module required)
```

外部的全基因组比对软件需根据输入方式准备。例如从 Minimap2 PAF 开始时，
先安装 `minimap2`；从已经生成的 MAF 或 AXT 开始时不需要再次运行 Minimap2。

## 4. 输入文件

### 4.1 基因组目录

每个组装使用一个 FASTA 文件。文件可以命名为 genome ID 本身，或使用
`.fa`、`.fasta`、`.fna` 后缀。例如：

```text
genomes/
  Ref
  Sample01
  Sample02
```

FASTA 序列名推荐使用：

```text
GenomeID.ChromosomeID
```

例如：

```text
>Ref.chr01
>Sample01.chr01
```

同一条染色体在不同组装中应保持一致的 chromosome ID。

### 4.2 genome list

`g.list` 每行一个 genome ID，包含参考基因组和所有 query：

```text
Ref
Sample01
Sample02
```

最终 VCF 的样本顺序按照该文件去掉 `--refname` 后的顺序确定。

### 4.3 比对输入命名

Minimap2 PAF：

```text
RefvsSample01.chr.paf
RefvsSample02.chr.paf
```

MAF：

```text
RefvsSample01.maf
RefvsSample02.maf
```

AXT：

```text
RefvsSample01.axt
RefvsSample02.axt
```

version2 可以利用 forward alignment 产生的 tnet 和 qnet，不要求把 forward
比对简单交换坐标来制造 reciprocal alignment。若已有真正独立运行的
query-vs-reference 比对，可通过 reciprocal 参数提供。

## 5. 推荐运行方法

### 5.1 从 Minimap2 PAF 开始

先对每个 query 运行比对：

```bash
minimap2 -t 16 -x asm5 -c --cs=long --mask-level 1 \
  genomes/Ref genomes/Sample01 \
  > alignments/RefvsSample01.chr.paf
```

`--mask-level 1` 会保留被同一主链覆盖区域中的局部反向链。默认的 `0.5`
在相邻变异较多时可能压掉真实倒位，从而使倒位事件和携带者同时漏检。
它会稍微增加 PAF 大小和后续 chain 时间；如果只追求最快运行且不关注倒位，
可省略该参数。

再运行完整流程：

```bash
perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from align \
  --aligner minimap2 \
  --align_dir alignments \
  --outdir SVGAP_version2_run \
  --threads 16
```

### 5.2 从 MAF 开始

```bash
perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from maf \
  --maf_dir maf_forward \
  --outdir SVGAP_version2_run \
  --threads 16
```

### 5.3 从 AXT 开始

```bash
perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from axt \
  --axt_forward_dir axt_forward \
  --outdir SVGAP_version2_run \
  --threads 16
```

## 6. 输出目录

```text
SVGAP_version2_run/
  00_config/          样本列表和本次运行配置
  01_axt/             转换后的 AXT
  02_chain_forward/   chain、net、SynNet 和 singRaw MAF
  03_filter_forward/  清理后的 SynNet 和 filtered_sing MAF
  04_pairwise_sv/     每个样本的原始 SV
  05_combined_sv/     群体合并后的候选集合
  06_genotype/        DEL、INS、small InDel、INV 和 CXSV 分型结果
  07_snp/             SNP VCF
  08_vcf/             最终 VCF
  09_annotation/      可选注释中间结果
```

主要结果：

| 文件 | 内容 |
|---|---|
| `All.largeSV.raw.vcf` | 大 DEL、INS 和 INV |
| `All.smallIndel.vcf` | small insertion/deletion |
| `All.CXSV.vcf` | 复杂结构变异 |
| `All.sorted.SNP.vcf` | SNP |
| `All.SV.vg.vcf` | 适合图基因组构建的 SV、small InDel、CXSV |
| `All.SVGAP.vg.vcf` | 上一文件再加入 SNP 的综合 VCF |

`All.SV.vg.validation.tsv` 和 `All.SVGAP.vg.validation.tsv` 记录参考等位基因、
坐标、`END`、`SVLEN`、排序和 ID 唯一性检查结果。

## 7. 常用参数

| 参数 | 说明 |
|---|---|
| `--threads N` | 独立样本或文件任务的最大并行数 |
| `--speciesname NAME` | 写入 VCF 的物种名 |
| `--skip_annotation` | 跳过可选的大 SV 注释 |
| `--skip_filter` | 使用原始 SynNet 和 `singRaw` MAF |
| `--skip_snp` | 不进行 SNP 分型 |
| `--skip_cxsv` | 不合并和分型复杂 SV |
| `--skip_vg_vcf` | 不生成最终 vg-ready VCF |
| `--no-index-vg-vcf` | 不运行 bgzip/tabix |
| `--reciprocal_strategy auto` | 自动选择 qnet 或真正的 reciprocal alignment |

默认建议保留 SynNet 过滤。短 contig 的阈值会自动调整。需要复现实验中的固定阈值时可用：

```text
--no-filter-adaptive
--filter_toplen N
--filter_nonsynlen N
--filter_synlen N
--filter_invlen N
```

`--filter_invlen` 同时控制 SynNet 过滤和 PairwiseSV 倒位调用，version2 默认
为 `1000` bp。旧的 10 kb 门槛会直接漏掉 1–10 kb 的倒位。

自适应参数：

```text
--filter_auto_fraction 0.10
--filter_auto_min 1000
```

实际阈值是固定阈值上限与按 contig 长度计算值中的较小者。

## 8. 断点续跑

若前面步骤已经成功，可配合以下参数跳过已有阶段：

```text
--skip_convert
--skip_chainnet
--skip_filter
--skip_pairwise_sv
--skip_combine
--skip_genotype
--skip_snp
--skip_annotation
--skip_vg_vcf
```

断点续跑前应确认 `--outdir` 中所需输入完整。`--skip_filter` 表示主动使用
原始结果；它与“复用已有 filtered 结果”含义不同。

## 9. 结果验证

单独验证最终 VCF：

```bash
perl 14_validate_vg_vcf.pl \
  --vcf SVGAP_version2_run/08_vcf/All.SVGAP.vg.vcf \
  --ref genomes/Ref \
  --out SVGAP_version2_run/08_vcf/All.SVGAP.vg.validation.tsv
```

报告末行应类似：

```text
#SUMMARY records=125;errors=0;warnings=0
```

出现 error 时不要直接把该 VCF 用于图构建，应先检查报告中对应记录的
参考碱基、坐标和等位基因长度。

## 10. 性能建议

1. `--threads` 不宜超过 query 数量太多；通常设为可用 CPU 数与 query 数的较小值。
2. FASTA、PAF、MAF 和输出目录放在本地 SSD 可减少大量随机 I/O。
3. 已有高质量 AXT 时从 `--start_from axt` 开始可省去格式转换。
4. 生产流程默认不生成未使用的 paftools 诊断变异表。单独运行
   `01_Convert2Axt.pl --emit_variants` 可恢复该输出。
5. PairwiseSV staging 在同一文件系统内使用硬链接，跨文件系统时自动复制。
6. INS 分型通常是 version2 后半段最慢的步骤；增加 query 级并行比提高单个
   样本内部线程更有效。
7. 倒位敏感模式使用 Minimap2 `--mask-level 1` 和默认的
   `--filter_invlen 1000`。若速度优先且只分析大于 10 kb 的倒位，可在上游
   省略 `--mask-level 1`，并设置 `--filter_invlen 10000`。

## 11. 常见问题

### 流程提示找不到 genome FASTA

确认 `g.list` 中的 ID 与 FASTA 文件名完全一致，并检查每个 FASTA 至少包含
一条非空序列。

### 找不到 PAF、MAF 或 AXT

检查文件名是否符合 `RefvsQuery` 形式，并确认 `Ref` 与 `--refname` 一致。

### 没有产生某一类 VCF

若对应类别没有真实数据记录，version2 会正常跳过该 genotyper。查看
`SVGAP_pipeline.log` 中的 `Skip because file has no data records` 信息。

### 最终 VCF 没有 `.gz` 和 `.tbi`

安装 htslib 的 `bgzip` 和 `tabix`，或通过 `--bgzip`、`--tabix` 指定路径。
纯文本 VCF 仍会保留。

### INV recall 或大 SV genotype F1 偏低

先确认 PAF 使用 `--mask-level 1` 生成，并在日志中确认 SynNetFilter 与
PairwiseSV 都收到相同的 `--invlen 1000`。若只有一个插入 FP 和一个 FN，且
两者序列相同、位置差等于重复单元长度，这通常是 tandem repeat 中的等价
VCF 表示。基准脚本应提供 `--ref_fasta genomes/Ref` 进行参考序列归一化；
同时保留严格坐标结果，便于区分表示差异与真正的基因型错误。

### 旧命令中的脚本名失效

version1 的 `1_Convert2Axt.pl`、`2_ChainNetSyn.pl` 等名称已由 version2 的
`01_Convert2Axt.pl`、`02_ChainNetSyn.pl` 等替换。建议直接使用
`SVGAP_pipeline.pl`，旧辅助脚本保存在 `trash/legacy_utilities/` 供追溯。
