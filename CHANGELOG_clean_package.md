# SVGAP11 clean package changelog

Release naming: the previous repository workflow is **version1**; this repaired,
optimized workflow is **version2**.

Date: 2026-07-07
Author line added to scripts: `Author: Yi Liao @ scau.edu.cn`

## Fixes

1. Fixed the recurring warning:

   `WARNING: Parallel::ForkManager is not installed; falling back to single-thread mode`

   The package already includes `lib/Parallel/ForkManager.pm` and its Perl dependencies. The problem was that scripts did not consistently add the package-local `lib/` directory to Perl `@INC`. All Perl scripts now add package-local library paths using `FindBin`, so the bundled module is found automatically when running scripts from the package directory or from subdirectories.

2. Added a consistent script header to every `.pl` file:

   `# Author: Yi Liao @ scau.edu.cn`

3. Moved non-production files to `trash/` instead of deleting them.

## Files moved to trash

See `trash/MANIFEST_trash.txt` for details.

## Validation performed

- `perl -c` was run on all main workflow scripts.
- `perl -c` was also run on optional scripts moved under `trash/`.
- A dry-run test of `02_ChainNetSyn.pl` with `--threads 2` confirmed that the Parallel::ForkManager warning no longer appears.

## 2026-07-20: vg-ready VCF revision

1. Fixed symbolic DEL anchoring so `END`, `POS` and `SVLEN` no longer differ by one base.
2. Added `13_prepare_vg_vcf.pl` to rebuild explicit DEL/INS alleles, normalize coordinates, merge components and generate final vg-ready VCFs.
3. Added `14_validate_vg_vcf.pl` to verify reference alleles, END/SVLEN consistency, sorting and unique IDs.
4. Added explicit sequence-resolved complex-SV output `All.CXSV.vg.vcf` and unresolved-sequence report.
5. Fixed sequence projection for reverse-direction CXSV records.
6. The pipeline now invokes `04_PairwiseSV.pl --exact_seq_len` so raw REF/ALT event sequences do not contain the legacy extra base.
7. Added automatic final outputs `All.SV.vg.vcf` and `All.SVGAP.vg.vcf`.
8. `06_genotype_DEL.pl`, `07_genotype_INS.pl` and `08_genotype_small_indels.pl` now use the reference FASTA to write explicit, left-anchored sequence alleles.
9. Pure insertion sequences are retained when PairwiseSV emits an empty reference-side interval.
10. The controller can automatically bgzip and tabix-index final vg-ready VCFs when htslib tools are available.

## 2026-10-20: reliability and performance revision

1. Fixed `SVGAP_pipeline.pl --skip_filter` to use raw SynNet and `singRaw` inputs.
2. Skip genotype/merge steps when an input contains only a header and no data records.
3. Added adaptive per-contig SynNet thresholds and exposed them through the controller.
4. Fixed raw qnet filename recognition so reciprocal calls keep `Query.Ref` orientation.
5. Replaced the broken dependency probe with package-root-aware executable checks.
6. Disabled the unused `paftools.js call` pass by default; `--emit_variants` restores it when wanted.
7. Reduced SNP MAF scanning overhead and switched independent job queues to core Perl `fork()/wait()`.
8. Hard-link filtered SynNet inputs into PairwiseSV staging when possible, with copy fallback.
9. Added short-contig, skip-filter, no-variant, truth-set accuracy, and vg-ready validation regressions.
