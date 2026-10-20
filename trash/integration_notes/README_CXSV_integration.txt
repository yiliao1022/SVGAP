SVGAP with CXSV integration
============================

This package keeps the original SVGAP module numbering and inserts complex SV
combination/genotyping as script 12:

  01_Convert2Axt.pl
  02_ChainNetSyn.pl
  03_SynNetFilter.pl
  04_PairwiseSV.pl        # produces *.cxsv.sort files
  05_CombineSV.pl
  06_genotype_DEL.pl
  07_genotype_INS.pl
  08_genotype_small_indels.pl
  09_SNPgenotyping.pl
  10_SVannotation.pl
  11_INVgenotyping.pl
  12_combine_genotype_CXSV.pl
  13_compare_truth_SVGAP.pl
  14_diagnose_small_indel_benchmark.pl

The updated SVGAP_pipeline.pl calls 12_combine_genotype_CXSV.pl inside the
genotyping block and writes:

  OUTDIR/06_genotype/CXSV/All.CXSV.records.projected.tsv
  OUTDIR/06_genotype/CXSV/All.CXSV.combined.tsv
  OUTDIR/06_genotype/CXSV/All.CXSV.genotype.matrix.tsv
  OUTDIR/06_genotype/CXSV/All.CXSV.vcf
  OUTDIR/08_vcf/All.CXSV.vcf

CXSV parameters exposed in SVGAP_pipeline.pl:

  --skip_cxsv
  --cxsv_break_tolerance 2000
  --cxsv_min_ro 0.50
  --cxsv_len_ratio 1.00
  --cxsv_unknown_missing

By default, non-carrier samples are genotyped as 0|0. With
--cxsv_unknown_missing, non-carriers are written as ./.
