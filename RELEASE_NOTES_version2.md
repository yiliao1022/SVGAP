# SVGAP version2

SVGAP version2 is the repaired and optimized successor to the original SVGAP workflow.

## Highlights

- Unified `SVGAP_pipeline.pl` workflow from whole-genome alignment output to population VCF.
- Supports Minimap2 PAF, LAST/LASTZ MAF, MUMmer4 delta, AnchorWave MAF, generic pairwise MAF, and existing AXT input.
- Complete split-contig workflow with automatic restoration to original chromosome coordinates.
- Improved large-SV discovery and genotype accuracy, including inversion-sensitive defaults.
- Correct handling of short contigs, `--skip_filter`, and datasets with no calls.
- Core Perl process parallelism and removal of an unused conversion pass.
- Reference-validated, explicit-allele `All.SV.vg.vcf` and `All.SVGAP.vg.vcf` outputs.
- Complete Chinese and English version2 manuals.

## Version history

- The `version1` tag remains fixed at commit `4c541f4` and preserves the original workflow and README.
- The `version2` tag contains the current workflow and release documentation.
- Obsolete draft manuals were removed from version2 to avoid presenting them as version1 documentation.

## Validation

- All production Perl scripts pass syntax checks and the bundled dependency checker.
- The two-sample accuracy run recovered 44/44 repeat-aware large-SV events and 70/70 carrier genotypes.
- A real MUMmer4 delta completed the full pipeline; both vg-ready VCF outputs passed with zero validation errors and warnings.
- LAST, LASTZ, generic MAF, AnchorWave MAF, Minimap2 PAF, and split-contig conversion paths were exercised.

Download `SVGAP_version2_20261020.tar.gz` and verify it with the attached SHA256 file.
