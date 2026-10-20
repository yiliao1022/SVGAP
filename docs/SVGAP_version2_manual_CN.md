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
10. Minimap2、LAST、LASTZ、MUMmer4、AnchorWave、通用 MAF 和 AXT 入口统一
    转换为 forward AXT，并共用同一套下游算法。

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

### 4.3 Aligner 兼容性与输入命名

SVGAP 不负责运行外部全基因组 aligner。它读取 aligner 已产生的 PAF、MAF、
delta 或 AXT，将其统一转换为 forward AXT，再进入完全相同的
chain/net、SV discovery 和 population genotyping 流程。

| Aligner/入口 | SVGAP 输入 | 顶层参数 | 文件名示例 | 转换方式 |
|---|---|---|---|---|
| Minimap2 | PAF，必须包含 CIGAR（`-c`） | `--start_from align --aligner minimap2` | `RefvsSample01.chr.paf` | `paftools.js view` → MAF → AXT |
| LAST | pairwise MAF | `--start_from align --aligner last` | `RefvsSample01.chr.maf` | 补 MAF header、排序 → AXT |
| LASTZ | pairwise MAF | `--start_from align --aligner lastz` | `RefvsSample01.chr.maf` | 补 MAF header、排序 → AXT |
| MUMmer4/nucmer | delta | `--start_from align --aligner mummer4` | `RefvsSample01.chr.delta` | bundled `delta2maf` → AXT |
| AnchorWave | pairwise MAF | `--start_from align --aligner anchorwave` | `RefvsSample01.chr.maf` | 恢复 genome prefix、内部排序 → AXT |
| 其他标准 MAF | pairwise MAF | `--start_from maf --maf_dir DIR` | `RefvsSample01.maf` | 通用 MAF → AXT |
| 已有 AXT | forward AXT | `--start_from axt --axt_forward_dir DIR` | `RefvsSample01.axt` | 不再转换 |

`mummer` 和 `mummer4` 是等价参数名；aligner 名称不区分大小写。`.chr` 可省略。
若 query 使用 SVGAP 的 split-contig 坐标，文件名使用 `.ctg.paf`、`.ctg.maf`
或 `.ctg.delta`，并提供 `--align_mode ctg --align_sizes FILE`。

文件名中的 `Ref`、`Sample01` 必须与 `g.list` 和 FASTA genome ID 一致。
标准 MAF 的 sequence source 推荐为 `Ref.chr01` 和 `Sample01.chr01`。
AnchorWave 的 MAF 可以只写 `chr01`，专用转换器会依据文件名补回 genome prefix；
每个 block 的第一条 `s` 行必须是 reference，第二条必须是 query。

version2 可以利用单次 forward alignment 产生的 tnet 和 qnet，不要求制造
reciprocal alignment。简单交换 forward 文件坐标不能产生独立比对。真正独立
运行的 query-vs-reference MAF/AXT 可以通过 reciprocal 参数提供。

### 4.4 已有 aligner 输出如何交给 SVGAP

每个非参考样本准备一个 pairwise alignment 文件，并放在同一个目录。例如：

```text
alignments/
  RefvsSample01.chr.paf
  RefvsSample02.chr.paf
```

将 `.paf` 换成对应 aligner 的 `.maf` 或 `.delta` 即可。文件名左侧必须是
`--refname`，右侧必须是 `g.list` 中的 query ID。一次 pipeline run 应使用同一种
aligner/格式；若不同样本来自不同 aligner，先分别转换为 AXT，再统一从
`--start_from axt` 开始。

第 5.1–5.5 节给出了各 aligner 的完整 pipeline 命令。若只想先检查格式转换，
也可以单独运行：

```bash
perl 01_Convert2Axt.pl \
  --ali anchorwave \
  --input alignments \
  --wk axt_check \
  --tname Ref \
  --qname Query \
  --t 8
```

把 `anchorwave` 换成 `minimap2`、`last`、`lastz`、`mummer4` 或 `maf`。
转换成功后，每个 query 都应有非空文件：

```text
axt_check/Target_Ref/RefvsSample01.axt
```

正常完整运行不需要先执行此检查；`SVGAP_pipeline.pl` 会自动完成相同转换。

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

### 5.2 从 LAST MAF 开始

以下是常见模板，LAST 参数可按物种距离调整：

```bash
lastdb -P16 lastdb/Ref genomes/Ref
lastal -P16 lastdb/Ref genomes/Sample01 | last-split \
  > alignments/RefvsSample01.chr.maf

perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from align \
  --aligner last \
  --align_dir alignments \
  --outdir SVGAP_version2_run \
  --threads 16
```

FASTA header 应包含 genome prefix，例如 `Ref.chr01` 和 `Sample01.chr01`。

### 5.3 从 LASTZ MAF 开始

```bash
lastz 'genomes/Ref[multiple]' genomes/Sample01 --format=maf \
  > alignments/RefvsSample01.chr.maf

perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from align \
  --aligner lastz \
  --align_dir alignments \
  --outdir SVGAP_version2_run \
  --threads 16
```

### 5.4 从 MUMmer4 delta 开始

version2 直接接受 MUMmer4 `nucmer` 产生的 delta。为了保留倒位和重排证据，
推荐使用 many-to-many `delta-filter -m`；`delta-filter -1` 可能删除相互重叠的
真实 rearrangement alignment。

```bash
nucmer --threads 16 --maxmatch \
  -p work/RefvsSample01 genomes/Ref genomes/Sample01

delta-filter -m -i 90 -l 100 work/RefvsSample01.delta \
  > alignments/RefvsSample01.chr.delta

perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from align \
  --aligner mummer4 \
  --align_dir alignments \
  --outdir SVGAP_version2_run \
  --threads 16
```

`nucmer` 和 `delta-filter` 来自外部 MUMmer4；delta-to-MAF 转换程序已包含在
SVGAP 的 `pub/MUMMER/` 中。

### 5.5 从 AnchorWave MAF 开始

先按当前 AnchorWave 版本的 `genoAli` 或 `proali` 流程生成每个 query 的
pairwise MAF，并保存为：

```text
alignments/RefvsSample01.chr.maf
alignments/RefvsSample02.chr.maf
```

然后运行：

```bash
perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from align \
  --aligner anchorwave \
  --align_dir alignments \
  --outdir SVGAP_version2_run \
  --threads 16
```

AnchorWave MAF 中 reference 必须是每个 block 的第一条 `s` 行。允许 sequence
ID 只写 `chr01`；专用转换器会补为 `Ref.chr01` 和 `Sample01.chr01`。
AnchorWave 的 anchor/CDS 配置会影响 intergenic SV 和 rearrangement recall，
因此应检查输出 MAF 是否保留反向链和长 gap。

### 5.6 从其他标准 MAF 开始

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

通用 MAF 应为 pairwise MAF，reference 在每个 block 中位于 query 之前，且
source ID 已包含 genome prefix。AnchorWave chromosome-only MAF 应使用上一节
的 `--aligner anchorwave`，不能直接走通用 MAF 入口。

### 5.7 从 AXT 开始

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

### 5.8 Split-contig alignment 完整流程

split-contig 模式适合单条 chromosome 很长、aligner 内存压力较大，或者需要按块
并行生成 alignment 的情况。只切分 query，reference 仍使用完整 FASTA。

第一步，切分每个 query：

```bash
perl 01_Convert2Axt.pl \
  --splitfa genomes/Sample01 \
  --split_out split/Sample01 \
  --bin 10000000 \
  --slide 10000
```

该命令生成三个文件：

```text
split/Sample01/Sample01.subcontigs.fa     # 用于比对的 query 子序列
split/Sample01/Sample01.sizes             # 原始 chromosome 长度
split/Sample01/Sample01.subcontigs.map.tsv # 子序列到原 chromosome 的映射
```

不要修改 `subcontigs.fa` 的 sequence ID。ID 末尾的 `_start_end` 是恢复原始坐标
所需的信息。

第二步，用所选 aligner 将完整 reference 与 split query 比对。Minimap2 示例：

```bash
minimap2 -t 16 -x asm5 -c --cs=long --mask-level 1 \
  genomes/Ref split/Sample01/Sample01.subcontigs.fa \
  > alignments/RefvsSample01.ctg.paf
```

其他 aligner 使用相同命名规则：LAST、LASTZ、AnchorWave 或通用 MAF 输出命名为
`RefvsSample01.ctg.maf`，MUMmer4 输出命名为
`RefvsSample01.ctg.delta`。alignment 中的 query source ID 必须保留
`Sample01.chr01_1_10010000` 这样的 split ID。

第三步，运行主流程并提供对应的 `.sizes` 文件。主流程在转换 AXT 前把所有 query
坐标恢复到原始 chromosome：

```bash
perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from align \
  --aligner minimap2 \
  --align_dir alignments \
  --align_mode ctg \
  --align_sizes split/Sample01/Sample01.sizes \
  --outdir SVGAP_version2_run \
  --threads 16
```

多个 query 时，对每个样本重复切分和比对，并重复提供参数：

```bash
--align_sizes split/Sample01/Sample01.sizes \
--align_sizes split/Sample02/Sample02.sizes
```

`.ctg` 文件名会让默认的 `--align_mode auto` 自动识别 split-contig 输入，但建议
显式写 `--align_mode ctg`。若 `.sizes` 文件放在 alignment 目录，转换器也会自动
发现；显式提供更容易复现。不要给 Sample01 使用 Sample02 的 sizes 文件。

转换后检查 `01_axt/Target_Ref/RefvsSample01.axt`。其中 query sequence 名应恢复为
`Sample01.chr01`，不能再含 `_start_end` 后缀。后续 chain/net、SV calling 和
genotyping 与普通 chromosome alignment 完全相同。

### 5.9 Forward 与 reciprocal alignment

默认 `--reciprocal_strategy auto` 在只有 forward input 时选择 `qnet`，这是
所有上述 aligner 的推荐入口。若已经独立运行了 query-vs-reference MAF：

```bash
perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from maf \
  --maf_dir maf_forward \
  --maf_reciprocal_dir maf_query_vs_ref \
  --reciprocal_strategy true \
  --outdir SVGAP_version2_run \
  --threads 16
```

reciprocal 文件必须来自真正独立运行的 Query-vs-Ref alignment。若已经是
AXT，则使用 `--axt_forward_dir`、`--axt_reciprocal_dir` 和相同 strategy。

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
| `--start_from align|maf|axt` | 选择 aligner 输出、通用 MAF 或已有 AXT 入口 |
| `--aligner NAME` | `minimap2`、`last`、`lastz`、`mummer/mummer4`、`anchorwave` 或 `maf` |
| `--align_mode auto|chr|ctg` | 普通染色体或 split-contig 坐标模式 |
| `--align_sizes FILE` | `ctg` 模式的原始染色体 sizes；可重复 |
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

各入口只扫描对应后缀：Minimap2 为 `.paf`，MUMmer4 为 `.delta`，
LAST/LASTZ/AnchorWave 为 `.maf`。查看日志中的 `Found N alignment file(s)`；
若 `N` 小于 query 数量，应先修正文件名。转换成功后还应确认
`01_axt/Target_Ref/RefvsQuery.axt` 非空。

### Aligner 已运行，但倒位或复杂 SV recall 偏低

SVGAP 只能使用输入 alignment 中保留的证据。Minimap2 使用
`--mask-level 1`；MUMmer4 优先使用 `delta-filter -m`；LAST、LASTZ 和
AnchorWave 应保留 reverse-strand block 与长 gap。过早做 one-to-one 或严格
collinear 过滤会在 SVGAP 之前删除真实 rearrangement。

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

## 12. 致谢

SVGAP version2 的开发得到了 **OpenAI Codex** 的协助，包括代码审查、流程调试、
性能优化、测试验证，以及中英文手册整理。项目的科学设计、结果解释和发布决定
由项目作者负责。
