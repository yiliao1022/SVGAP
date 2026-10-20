# SVGAP version2

SVGAP detects, merges, genotypes, and annotates structural variants, small
indels, and SNPs across assembly populations. This revision also writes
reference-consistent, vg-ready VCF files.

The previous repository workflow is designated **version1**. The current
repaired and optimized workflow is **version2**. Check an installation with:

```bash
perl SVGAP_pipeline.pl --version
```

## Quick start

Prepare one FASTA per assembly in `genomes/`. Sequence names should use the
`GenomeID.ChromosomeID` convention. Put the reference and query genome IDs in
`g.list`, one per line, and place pairwise reference-vs-query alignments in one
directory.

```bash
perl util/Test_dependency.pl

perl SVGAP_pipeline.pl \
  --refname Ref \
  --glist genomes/g.list \
  --genome_dir genomes \
  --start_from align \
  --aligner minimap2 \
  --align_dir alignments \
  --outdir SVGAP_run \
  --threads 16
```

The pipeline accepts Minimap2 PAF, MUMmer delta, MAF, or existing AXT input.
Run `perl SVGAP_pipeline.pl --help` for all modes and options.

For inversion-sensitive Minimap2 input, create PAF with
`-x asm5 -c --cs=long --mask-level 1`. SVGAP uses a 1 kb inversion floor by
default and applies it consistently in SynNet filtering and PairwiseSV calling.

Main outputs are written under `SVGAP_run/08_vcf/`:

- `All.largeSV.raw.vcf`
- `All.smallIndel.vcf`
- `All.sorted.SNP.vcf`
- `All.SV.vg.vcf`
- `All.SVGAP.vg.vcf`

`All.SV.vg.vcf` contains explicit sequence-resolved SV/small-indel/CXSV
alleles. `All.SVGAP.vg.vcf` additionally contains SNPs. Both are checked by the
included reference-aware validator.

## Documentation

- [中文手册](docs/SVGAP11_manual_CN.md)
- [version2 中文运行手册](docs/SVGAP_version2_manual_CN.md)
- [English manual](docs/SVGAP11_manual_EN.md)
- [version2 English manual](docs/SVGAP_version2_manual_EN.md)
- [VG-ready VCF 中文说明](VG_READY_VCF_README_CN.md)
- [VG-ready VCF English notes](VG_READY_VCF_README_EN.md)

## Current fixes

The 2026-10-20 revision repairs `--skip_filter`, empty-call datasets, raw qnet
orientation, and short-contig filtering. It also removes an unused paftools
pass, speeds up SNP MAF scanning, and uses core Perl process parallelism without
a CPAN runtime dependency. See [CHANGELOG_clean_package.md](CHANGELOG_clean_package.md).
