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
10. Minimap2, LAST, LASTZ, MUMmer4, AnchorWave, generic MAF, and AXT entry
    points normalize to forward AXT and share the same downstream algorithms.

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

### 4.3 Aligner compatibility and file names

SVGAP does not run an external whole-genome aligner. It consumes PAF, MAF,
delta, or AXT produced by an aligner, converts the input to a common forward
AXT representation, and then uses the same chain/net, SV discovery, and
population-genotyping stages.

| Aligner or entry point | SVGAP input | Top-level options | Example file | Conversion |
|---|---|---|---|---|
| Minimap2 | PAF with CIGAR (`-c`) | `--start_from align --aligner minimap2` | `RefvsSample01.chr.paf` | `paftools.js view` → MAF → AXT |
| LAST | Pairwise MAF | `--start_from align --aligner last` | `RefvsSample01.chr.maf` | Add MAF header, sort → AXT |
| LASTZ | Pairwise MAF | `--start_from align --aligner lastz` | `RefvsSample01.chr.maf` | Add MAF header, sort → AXT |
| MUMmer4/nucmer | delta | `--start_from align --aligner mummer4` | `RefvsSample01.chr.delta` | Bundled `delta2maf` → AXT |
| AnchorWave | Pairwise MAF | `--start_from align --aligner anchorwave` | `RefvsSample01.chr.maf` | Restore genome prefixes, internal sort → AXT |
| Other standard MAF | Pairwise MAF | `--start_from maf --maf_dir DIR` | `RefvsSample01.maf` | Generic MAF → AXT |
| Existing AXT | Forward AXT | `--start_from axt --axt_forward_dir DIR` | `RefvsSample01.axt` | No conversion |

`mummer` and `mummer4` are equivalent option values, and aligner names are
case-insensitive. The `.chr` component is optional. For alignments against
SVGAP split query contigs, use `.ctg.paf`, `.ctg.maf`, or `.ctg.delta` together
with `--align_mode ctg --align_sizes FILE`.

The `Ref` and `Sample01` portions of a file name must match the genome IDs in
`g.list` and the FASTA files. Standard MAF source names should be
`Ref.chr01` and `Sample01.chr01`. AnchorWave MAF may use chromosome-only names
such as `chr01`; its dedicated converter restores genome prefixes from the file
name. The first `s` line in each AnchorWave block must be reference and the
second must be query.

version2 derives tnet and qnet from one forward alignment, so every aligner can
use the forward-only path. Swapping coordinates in a forward file does not
create an independent reciprocal alignment. Truly independent query-vs-reference
MAF or AXT can be supplied through the reciprocal options.

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

### 5.2 Start from LAST MAF

This is a common template; tune LAST parameters for the assembly divergence:

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

FASTA headers should contain genome prefixes such as `Ref.chr01` and
`Sample01.chr01`.

### 5.3 Start from LASTZ MAF

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

### 5.4 Start from MUMmer4 delta

version2 accepts delta written by MUMmer4 `nucmer`. Use many-to-many
`delta-filter -m` when rearrangement sensitivity matters. `delta-filter -1`
can discard genuine overlapping rearrangement alignments.

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

`nucmer` and `delta-filter` come from the external MUMmer4 installation. SVGAP
bundles the delta-to-MAF converter under `pub/MUMMER/`.

### 5.5 Start from AnchorWave MAF

Use the `genoAli` or `proali` workflow for your AnchorWave version to generate
one pairwise MAF per query, saved as:

```text
alignments/RefvsSample01.chr.maf
alignments/RefvsSample02.chr.maf
```

Then run:

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

Reference must be the first `s` line in each AnchorWave MAF block.
Chromosome-only source names such as `chr01` are accepted; the dedicated
converter changes them to `Ref.chr01` and `Sample01.chr01`. Anchor/CDS settings
affect intergenic-SV and rearrangement recall, so confirm that the output MAF
retains reverse-strand blocks and long gaps.

### 5.6 Start from another standard MAF

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

Generic input must be pairwise MAF with reference before query in each block,
and source names must already include genome prefixes. Use the dedicated
`--aligner anchorwave` path for chromosome-only AnchorWave MAF.

### 5.7 Start from AXT

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

### 5.8 Split-contig alignment

Split a query first:

```bash
perl 01_Convert2Axt.pl \
  --splitfa genomes/Sample01 \
  --split_out split/Sample01 \
  --bin 10000000 \
  --slide 10000
```

Align `split/Sample01/Sample01.subcontigs.fa` and use `.ctg` in the result name,
for example `RefvsSample01.ctg.paf`. The controller restores original chromosome
coordinates with:

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

`--align_sizes` is repeatable. The converter can discover `.sizes` files in the
alignment directory, but explicit paths make the run easier to reproduce.

### 5.9 Forward and reciprocal alignment

With forward input only, the default `--reciprocal_strategy auto` selects
`qnet`; this is the recommended entry point for every aligner above. If true,
independently generated query-vs-reference MAF is already available, use:

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

Reciprocal files must come from an independent Query-vs-Ref alignment. For
existing AXT, use `--axt_forward_dir`, `--axt_reciprocal_dir`, and the same
strategy.

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
| `--start_from align|maf|axt` | Select aligner output, generic MAF, or existing AXT input |
| `--aligner NAME` | `minimap2`, `last`, `lastz`, `mummer/mummer4`, `anchorwave`, or `maf` |
| `--align_mode auto|chr|ctg` | Chromosome or split-contig coordinate mode |
| `--align_sizes FILE` | Original chromosome sizes for `ctg` mode; repeatable |
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

Each entry point scans only its expected suffix: `.paf` for Minimap2, `.delta`
for MUMmer4, and `.maf` for LAST, LASTZ, or AnchorWave. Check `Found N alignment
file(s)` in the log. If `N` is smaller than the query count, fix the file names
before continuing. After conversion, confirm that
`01_axt/Target_Ref/RefvsQuery.axt` exists and is non-empty.

### The aligner ran, but inversion or complex-SV recall is low

SVGAP can only use evidence retained in its alignment input. Use
`--mask-level 1` with Minimap2 and prefer `delta-filter -m` with MUMmer4. LAST,
LASTZ, and AnchorWave output should retain reverse-strand blocks and long gaps.
An early one-to-one or strict collinearity filter can remove genuine
rearrangements before SVGAP sees them.

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
