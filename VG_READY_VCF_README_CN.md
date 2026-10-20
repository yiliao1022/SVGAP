# SVGAP：生成可用于 `vg construct` 的VCF

## 1. 本次修复的核心问题

原SVGAP DEL记录由0-based half-open区间转换为symbolic `<DEL>`时使用了：

```text
VCF POS = BED_START + 1
END     = BED_END
SVLEN   = -(BED_END - BED_START)
```

这会导致：

```text
END - POS = |SVLEN| - 1
```

例如：

```text
POS=23056192; END=23062329; SVLEN=-6138; BED_START=23056191
```

`vg construct -S`因此报告：

```text
Warning: deletion END and SVLEN do not agree [canonicalize]
```

新版不再依赖不稳定的symbolic DEL坐标解释。最终vg-ready VCF会从参考基因组重新提取删除区间，并生成显式、左锚定的REF/ALT：

```text
POS   = BED_START                    # 删除前锚定碱基的1-based位置
REF   = 左锚定碱基 + 被删除序列
ALT   = 左锚定碱基
END   = POS + length(REF) - 1
SVLEN = length(ALT) - length(REF)
```

因此恒有：

```text
END - POS = |SVLEN|
```

## 2. 新增脚本

### `13_prepare_vg_vcf.pl`

将SVGAP各组件VCF合并并转换为vg-ready VCF：

- symbolic DEL重建为显式左锚定DEL；
- 旧的右锚定INS转换为左锚定显式INS；
- 显式REF必须与参考FASTA完全一致；
- `END`由REF跨度重新计算；
- `SVLEN`由`length(ALT)-length(REF)`重新计算；
- symbolic INV保留为`<INV>`，但REF改为真实参考碱基；
- 无完整替代序列的symbolic INS/CXSV不会静默写入，而是记录在report中；
- 强制VCF按参考染色体顺序和位置排序；
- 强制ID唯一。

### `14_validate_vg_vcf.pl`

验证最终VCF：

- REF是否与参考FASTA一致；
- 显式等位基因的`END`是否等于`POS + length(REF) - 1`；
- `SVLEN`是否等于`length(ALT)-length(REF)`；
- symbolic INV的POS、END和SVLEN是否一致；
- VCF是否排序；
- ID是否唯一；
- 是否仍包含未解析的symbolic INS/CXSV。

验证发现错误时退出码为2，主流程会停止。

## 3. Complex SV序列化

`12_combine_genotype_CXSV.pl`现在新增：

```text
--ref_fasta FILE
--max_explicit_allele 5000000
```

并输出：

```text
All.CXSV.vcf
```

保留旧的symbolic CXSV格式，便于兼容旧分析。

```text
All.CXSV.vg.vcf
```

为新的序列解析VCF：

- 从同一cluster中选择坐标和长度最接近cluster中心、且具有完整替代序列的代表记录；
- 正确处理forward和reverse比较；
- reverse记录会交换reference/alternate序列，并在必要时反向互补回原参考方向；
- 从参考FASTA重新提取reference-side序列；
- REF/ALT均加入共享锚定碱基；
- 纯reference-side事件标记为`DEL`；
- 纯alternate-side事件标记为`INS`；
- 两侧均有序列的replacement-like事件标记为`CXSV`。

未能恢复完整替代序列的cluster输出到：

```text
All.CXSV.vg.unresolved.tsv
```

这些事件不会进入最终构图VCF，以避免使用截断序列或伪序列构图。


## 3A. DEL、INS和small InDel模块的直接修复

除了最终整理脚本，本版还从源头修复了三个genotyping模块：

- `06_genotype_DEL.pl --ref_fasta reference.fa`：从参考基因组提取真实删除序列并生成左锚定显式DEL；如果无法恢复，symbolic fallback也严格满足`END-POS=abs(SVLEN)`。
- `07_genotype_INS.pl --ref_fasta reference.fa`：当`REFseq`为空、`ALTseq`为纯插入序列时，自动添加真实左锚定碱基，不再退回并丢失为`<INS>`。
- `08_genotype_small_indels.pl --ref_fasta reference.fa`：small DEL和small INS统一生成参考一致的显式REF/ALT。

主流程已经自动向这三个脚本传递reference FASTA。

## 4. 主流程新增输出

完整运行后，`08_vcf/`新增：

```text
All.CXSV.vg.vcf
All.CXSV.vg.unresolved.tsv

All.SV.vg.vcf
All.SV.vg.report.tsv
All.SV.vg.validation.tsv

All.SVGAP.vg.vcf
All.SVGAP.vg.report.tsv
All.SVGAP.vg.validation.tsv
All.SV.vg.vcf.gz + All.SV.vg.vcf.gz.tbi
All.SVGAP.vg.vcf.gz + All.SVGAP.vg.vcf.gz.tbi
```

其中：

- 如果系统能够找到`bgzip`和`tabix`，主流程会自动生成压缩VCF和tabix索引；若工具缺失，则保留plain VCF并在日志中给出命令。

- `All.SV.vg.vcf`：large DEL/INS/INV + small InDel + sequence-resolved CXSV；
- `All.SVGAP.vg.vcf`：在上述基础上加入SNP；
- `*.report.tsv`：每条记录的坐标转换、重建或跳过原因；
- `*.validation.tsv`：最终一致性检查。

## 5. 推荐构图命令

```bash
bgzip -f -c SVGAP_run/08_vcf/All.SVGAP.vg.vcf \
  > SVGAP_run/08_vcf/All.SVGAP.vg.vcf.gz

tabix -f -p vcf SVGAP_run/08_vcf/All.SVGAP.vg.vcf.gz

vg construct \
  -a \
  -f \
  -S \
  -r genomes/ISO1.fa \
  -v SVGAP_run/08_vcf/All.SVGAP.vg.vcf.gz \
  > ISO1.SVGAP.vg
```

说明：

- DEL、INS和CXSV在最终VCF中尽量使用显式REF/ALT序列；
- `-S`主要用于保留的symbolic `<INV>`；
- 不需要额外的insertion FASTA，因为可解析INS已经直接写入ALT；
- `-a`保存variant paths，便于后续基于原VCF genotyping；
- `-f`避免vg再次切分输入ALT序列。

## 6. 超长等位基因

默认最大显式REF或ALT长度为5 Mb：

```text
--max_vg_allele_len 5000000
```

超过该长度的记录不会直接写入vg-ready VCF，而会进入report。可根据计算资源调整，例如：

```bash
--max_vg_allele_len 10000000
```

注意：非常长的显式等位基因会显著增加VCF体积、graph构建内存和索引复杂度。

## 7. 单独整理已有SVGAP结果

不重新运行前面步骤，也可以直接整理已有VCF：

```bash
perl 13_prepare_vg_vcf.pl \
  --ref genomes/ISO1.fa \
  --vcfs All.largeSV.raw.vcf,All.smallIndel.vcf,All.CXSV.vg.vcf,All.sorted.SNP.vcf \
  --query_list query.list \
  --refname ISO1 \
  --refsize genomes/ISO1.sizes \
  --out All.SVGAP.vg.vcf \
  --report All.SVGAP.vg.report.tsv

perl 14_validate_vg_vcf.pl \
  --vcf All.SVGAP.vg.vcf \
  --ref genomes/ISO1.fa \
  --out All.SVGAP.vg.validation.tsv
All.SV.vg.vcf.gz + All.SV.vg.vcf.gz.tbi
All.SVGAP.vg.vcf.gz + All.SVGAP.vg.vcf.gz.tbi
```
