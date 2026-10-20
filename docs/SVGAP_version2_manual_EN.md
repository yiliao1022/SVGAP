# SVGAP version2 User Manual

## 1. Version naming

- **version1** is the previous workflow in the Git repository.
- **version2** is the current replacement. It includes CXSV output, vg-ready
  VCFs, adaptive filtering for short contigs, correct empty-call handling,
  a repaired `--skip_filter` path, and additional performance improvements.

Check the installed version with:

```bash
perl SVGAP_pipeline.pl --version
```

Expected output:

```text
SVGAP version2
```

## 2. Main changes in version2

1. `SVGAP_pipeline.pl` runs the complete workflow from alignment output to
   population VCFs.
2. The workflow calls large DEL, INS, INV, complex SVs, small indels, and SNPs.
3. It writes reference-validated `All.SV.vg.vcf` and `All.SVGAP.vg.vcf` files.
4. Raw qnet files retain the correct `Query.Ref` orientation and insertion evidence.
5. `--skip_filter` reads raw SynNet files and `singRaw` MAF files directly.
6. Empty variant tables skip incompatible genotypers instead of aborting the run.
7. SynNet thresholds scale automatically for short contigs.
8. Independent jobs use core Perl `fork()/wait()` parallelism and need no CPAN
   runtime module.
9. Minimap2 PAF conversion skips the unused `paftools.js call` pass by default.

## 3. Environment check

From the repository root, run:

```bash
perl util/Test_dependency.pl
```

A successful check includes:

```text
SVGAP dependency check: OK
Parallel workers: core Perl fork/wait (no CPAN module required)
```

Install the whole-genome aligner needed to create your input. For example,
Minimap2 must be available when producing PAF input. Starting from existing MAF
or AXT files does not require rerunning Minimap2.

## 4. Input files

### 4.1 Genome directory

Store one FASTA per assembly. A file may use the genome ID without an extension,
or the `.fa`, `.fasta`, or `.fna` extension:

```text
genomes/
  Ref
  Sample01
  Sample02
```

Use `GenomeID.ChromosomeID` sequence names when possible:

```text
>Ref.chr01
>Sample01.chr01
```

Use consistent chromosome IDs across assemblies.

### 4.2 Genome list

`g.list` contains the reference and all query genome IDs, one per line:

```text
Ref
Sample01
Sample02
```

The final VCF sample order follows this list after removing `--refname`.

### 4.3 Alignment file names

Minimap2 PAF:

```text
RefvsSample01.chr.paf
RefvsSample02.chr.paf
```

MAF:

```text
RefvsSample01.maf
RefvsSample02.maf
```

AXT:

```text
RefvsSample01.axt
RefvsSample02.axt
```

version2 uses tnet and qnet from the forward alignment. Do not create a false
reciprocal alignment by swapping coordinates in the forward file. If truly
independent query-vs-reference alignments are available, provide them through
the reciprocal alignment options.

## 5. Recommended commands

### 5.1 Start from Minimap2 PAF

Create one alignment per query:

```bash
minimap2 -t 16 -x asm5 -c --cs=long --mask-level 1 \
  genomes/Ref genomes/Sample01 \
  > alignments/RefvsSample01.chr.paf
```

`--mask-level 1` retains local reverse-strand chains that overlap the primary
chain. Minimap2's default value of `0.5` can suppress true inversions when
nearby variants compete for the same interval, causing both event and carrier
false negatives. This option can modestly increase PAF size and chaining time;
omit it when maximum speed matters and inversion sensitivity does not.

Run SVGAP:

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

### 5.2 Start from MAF

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

### 5.3 Start from AXT

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

## 6. Output layout

```text
SVGAP_version2_run/
  00_config/          sample lists and run configuration
  01_axt/             converted AXT files
  02_chain_forward/   chain, net, SynNet, and singRaw MAF files
  03_filter_forward/  cleaned SynNet and filtered_sing MAF files
  04_pairwise_sv/     pairwise raw SV calls
  05_combined_sv/     merged population candidate sets
  06_genotype/        DEL, INS, small-indel, INV, and CXSV genotypes
  07_snp/             SNP VCF
  08_vcf/             final VCFs
  09_annotation/      optional annotation work files
```

Main result files:

| File | Contents |
|---|---|
| `All.largeSV.raw.vcf` | Large DEL, INS, and INV records |
| `All.smallIndel.vcf` | Small insertions and deletions |
| `All.CXSV.vcf` | Complex structural variants |
| `All.sorted.SNP.vcf` | SNPs |
| `All.SV.vg.vcf` | vg-ready SV, small-indel, and CXSV records |
| `All.SVGAP.vg.vcf` | Combined vg-ready VCF including SNPs |

`All.SV.vg.validation.tsv` and `All.SVGAP.vg.validation.tsv` report checks for
reference alleles, coordinates, `END`, `SVLEN`, ordering, and unique IDs.

## 7. Common options

| Option | Meaning |
|---|---|
| `--threads N` | Maximum concurrent sample or file jobs |
| `--speciesname NAME` | Species name written to VCF metadata |
| `--skip_annotation` | Skip optional large-SV annotation |
| `--skip_filter` | Use raw SynNet and `singRaw` MAF files |
| `--skip_snp` | Skip SNP genotyping |
| `--skip_cxsv` | Skip complex-SV merging and genotyping |
| `--skip_vg_vcf` | Skip final vg-ready VCF generation |
| `--no-index-vg-vcf` | Do not run bgzip and tabix |
| `--reciprocal_strategy auto` | Select qnet or true reciprocal input automatically |

Adaptive SynNet filtering is enabled by default. To reproduce fixed thresholds,
use:

```text
--no-filter-adaptive
--filter_toplen N
--filter_nonsynlen N
--filter_synlen N
--filter_invlen N
```

`--filter_invlen` now controls both SynNet filtering and the PairwiseSV
inversion call floor. The version2 default is `1000` bp. The former 10 kb floor
discarded inversions between 1 and 10 kb.

Adaptive controls are:

```text
--filter_auto_fraction 0.10
--filter_auto_min 1000
```

The effective threshold is the smaller of the configured ceiling and the value
calculated from the target contig length.

## 8. Restarting from existing results

The following options skip completed stages:

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

Before restarting, verify that the required files exist under `--outdir`.
`--skip_filter` explicitly selects raw SynNet and `singRaw` inputs; it does not
mean “reuse existing filtered results.”

## 9. Validating a final VCF

Run the validator independently with:

```bash
perl 14_validate_vg_vcf.pl \
  --vcf SVGAP_version2_run/08_vcf/All.SVGAP.vg.vcf \
  --ref genomes/Ref \
  --out SVGAP_version2_run/08_vcf/All.SVGAP.vg.validation.tsv
```

The final line should resemble:

```text
#SUMMARY records=125;errors=0;warnings=0
```

Resolve validation errors before using the VCF to build a graph. The report
identifies records with reference-base, coordinate, or allele-length problems.

## 10. Performance guidance

1. Set `--threads` near the smaller of the available CPU count and query count.
2. Store FASTA, PAF/MAF, and output files on local SSD storage when possible.
3. Start from `--start_from axt` when validated AXT files already exist.
4. The production path omits the unused paftools diagnostic variant table.
   Run `01_Convert2Axt.pl --emit_variants` when that table is needed.
5. PairwiseSV staging uses hard links on the same filesystem and falls back to copying.
6. INS genotyping is commonly the slowest remaining version2 stage. Parallelism
   across queries usually helps more than extra threads within one sample.
7. Inversion-sensitive mode uses Minimap2 `--mask-level 1` and the default
   `--filter_invlen 1000`. For maximum speed when only inversions above 10 kb
   matter, omit `--mask-level 1` upstream and set `--filter_invlen 10000`.

## 11. Troubleshooting

### Genome FASTA is not found

Confirm that IDs in `g.list` exactly match FASTA file names and that every FASTA
contains at least one non-empty sequence.

### PAF, MAF, or AXT files are not detected

Check the `RefvsQuery` naming convention and ensure that `Ref` matches `--refname`.

### One VCF category is missing

When a category has no data records, version2 skips its genotyper normally.
Check `SVGAP_pipeline.log` for `Skip because file has no data records`.

### Final `.vcf.gz` and `.tbi` files are missing

Install the htslib `bgzip` and `tabix` programs, or provide their paths with
`--bgzip` and `--tabix`. The plain-text VCF remains available.

### INV recall or large-SV genotype F1 is low

Confirm that PAF was generated with `--mask-level 1`, then check the pipeline
log to ensure that SynNetFilter and PairwiseSV both received `--invlen 1000`.
If the only discrepancy is one insertion FP plus one FN with identical inserted
sequence and a position shift equal to the repeat unit, the records can be
equivalent tandem-repeat VCF representations. Pass `--ref_fasta genomes/Ref`
to the benchmark tool for reference-aware normalization, and retain the strict
coordinate result to distinguish representation differences from genotype errors.

### Commands using version1 script names fail

version1 names such as `1_Convert2Axt.pl` and `2_ChainNetSyn.pl` were replaced by
version2 names such as `01_Convert2Axt.pl` and `02_ChainNetSyn.pl`. Prefer the
top-level `SVGAP_pipeline.pl`. Retained legacy utilities are under
`trash/legacy_utilities/` for reference.
