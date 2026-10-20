# SVGAP vg-ready VCF output

This revision adds a final normalization stage for `vg construct`.

## Main changes

- Symbolic SVGAP deletions are rebuilt as explicit left-anchored alleles from the reference FASTA.
- Explicit insertions are canonicalized from SVGAP 0-based half-open coordinates.
- `END` is recalculated from the explicit REF span.
- `SVLEN` is recalculated as `length(ALT)-length(REF)`.
- Symbolic inversions are retained with the true reference base.
- Complex-SV clusters with complete sequences are exported as explicit REF/ALT alleles in `All.CXSV.vg.vcf`.
- Unresolved symbolic INS/CXSV records are excluded from the graph VCF and listed in a report.
- Final VCFs are reference-checked and validated automatically.


## Direct fixes in source genotypers

- `06_genotype_DEL.pl --ref_fasta reference.fa` writes real left-anchored explicit deletion alleles; its symbolic fallback enforces `END-POS=abs(SVLEN)`.
- `07_genotype_INS.pl --ref_fasta reference.fa` retains a pure inserted `ALTseq` when the reference-side interval is empty and adds the true reference anchor.
- `08_genotype_small_indels.pl --ref_fasta reference.fa` writes explicit reference-consistent alleles for small DEL and INS records.
- `SVGAP_pipeline.pl` passes the reference FASTA automatically and invokes `04_PairwiseSV.pl --exact_seq_len`.

## New outputs

```text
08_vcf/All.CXSV.vg.vcf
08_vcf/All.CXSV.vg.unresolved.tsv
08_vcf/All.SV.vg.vcf
08_vcf/All.SV.vg.report.tsv
08_vcf/All.SV.vg.validation.tsv
08_vcf/All.SVGAP.vg.vcf
08_vcf/All.SVGAP.vg.report.tsv
08_vcf/All.SVGAP.vg.validation.tsv
All.SV.vg.vcf.gz + All.SV.vg.vcf.gz.tbi
All.SVGAP.vg.vcf.gz + All.SVGAP.vg.vcf.gz.tbi
```

When `bgzip` and `tabix` are available, the controller creates the `.vcf.gz` and `.tbi` files automatically. Otherwise it retains the plain VCF and logs the indexing commands.

## Recommended graph construction

```bash
bgzip -f -c SVGAP_run/08_vcf/All.SVGAP.vg.vcf > SVGAP_run/08_vcf/All.SVGAP.vg.vcf.gz
tabix -f -p vcf SVGAP_run/08_vcf/All.SVGAP.vg.vcf.gz
vg construct -a -f -S -r genomes/ISO1.fa -v SVGAP_run/08_vcf/All.SVGAP.vg.vcf.gz > ISO1.SVGAP.vg
```

`-S` is retained for symbolic `<INV>` records. Resolved DEL, INS and CXSV alleles are explicit and do not require an insertion FASTA.
