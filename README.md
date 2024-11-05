# SVGAP (Structural Variants Genotyping of Assemblies on Population scales)

SVGAP is a flexible pipeline to detect, genotype, and annotate SVs in large samples of de novo genome asse
mblies. It compares each sample to a reference genome in whole-genome alignments to call SVs. The SVs iden
tified are subsequently combined across samples to generate a nonredundant call set. Each SV call in this
call set can be further re-genotyped by examining local alignment information specific to each sample to p
roduce fully genotyped VCF (Variant Call Format) files. The SVGAP pipeline consists of six main steps:

(1) processing whole genome alignments (WGA);

(2) constructing syntenic alignments;

(3) calling SVs;

(4) merging SVs from multiple samples;

(5) re-genotyping SVs for each sample; and

(6) annotating SVs.

These steps are executed sequentially using a set of perl programs.

## Step-by-step Instructional Guides
After downloading the SVGAP package from the GitHub repository (https://github.com/yiliao1022/SVGAP), in m
ost cases, you can immediately start using it in your projects without requiring any additional dependenci
es other than perl. Here, we provided a step-by-step practical instruction using rice data, which includes
 high-quality genome assemblies of the reference MingHui 63 and 48 other representative accessions.

`git clone https://github.com/yiliao1022/SVGAP`

### Step 0: Constructing a working folder
To start the project, we first constructed a folder named “rice”:
```
`mkdir rice`\
`mkdir rice/genome`
```
### Step 1: Preparing the genomes
#### Assembly quality
The SVGAP pipeline was originally designed for genome assemblies of chromosome-level quality. However, it
can also be used effectively for contig-level assemblies. Theoretically, there are no specific quality req
uirements for the genome assembly, such as N50 or contig numbers. Nevertheless, we highly recommend that u
sers only include chromosomes and long contigs when preparing the input genome assemblies. Currently, SVGA
P is optimized for haploid genomes. If you have phased genomes with multiple haplotypes, the recommended a
pproach is to treat each haplotype as an individual assembly.

#### Naming
To effectively utilize the SV-GAPS pipeline, users must adhere to specific naming guidelines for genome as
semblies and their chromosomes/contigs. Here are some key points to keep in mind:

1. Assembly names do not need to include file extensions such as ".fasta" or ".fa".

2. Avoid naming assemblies solely using numbers. For example, naming an assembly "9311" could potentially
cause issues in certain situations. Instead, consider using "N9311".

3. Do not name the query assembly in a way that completely matches or overlaps with the names of reference
 assemblies. For instance, if the reference name is "MH63", it would be advisable to avoid naming any quer
y assembly "H63" or "MH6".

4. Chromosome names should follow the format GenomeID.ChromosomeID, such as "MH63.chr01".

We named our test samples as follows:
##### Genome ID
1. Reference Name:

MH63

2. Names for 48 representative accessions:

Accession Name List1    Accession Name List2    Accession Name List3    Accession Name List4    Accession
Name List5      Accession Name List6

Basm    DHX2H   Gui630  LJ      Oooo    TUMBA

Cccc    Eeee    II32    Llll    R3551   WSSM

CJ14    FH838   J4115S  N22     R527    YueGuang

CN1B    FS32    Jjjj    N9311   R548    YX

D62B    G46     Kkkk    NAMROO  SJ18    Z02428

Dddd    Gggg    KY131   Nipp    Ssss    ZH11

DG      Guang8B Lemont  Nnnn    TSIPALA ZS97

##### Chromosome ID
Chromosome ID List1     Chromosome ID List2     Chromosome ID List3     Chromosome ID List4     Chromosome
 ID List5       Chromosome ID List6

MH63.chr01      MH63.chr02      MH63.chr03      MH63.chr04      MH63.chr05      MH63.chr06

MH63.chr07      MH63.chr08      MH63.chr09      MH63.chr10      MH63.chr11      MH63.chr12

#### Splitting (Optional)
If you are working with a plant genome with very large genome size, such as maize and pepper, you may need
 to split the genome into small pieces, and then conduct the pairwise genome alignment. Here we provide a
script called 0_SplitFa.pl to facilitate this splitting work. You just need to split the query one in the
genome assemblies folder, here is an example:

`perl ~/bin/SV-GAPS/util/0_SplitFa.pl ZH97 2000000 400000`

This script will split a FASTA file into fixed-sized windows with a specified step size. The first paramet
er "ZH97" refers to the input FASTA file, the second parameter "2000000" indicates the size of the fixed-s
ized windows, and the third parameter "400000" denotes the step size for moving along the sequence.

Moving working genomes to subfolder

1. Once the working genomes have been named and filtered for short contigs, they were moved to the "rice/g
enome" subfolder.

2. Please be note that even if you split the query genomes, the "rice/genome" should only contain unsplit
files.

### Step 2: Whole genome alignment
Currently, the SV-GAPS pipeline supports alignment results from a variety of aligners, namely LASTZ, LAST,
 MUMmer, Minimap2, GSAlign, and AnchorWave:

1. For small and medium-sized genomes (< 500 Mb), we recommend utilizing LAST and MUMmer.

2. For larger genomes, particularly in the case of plants, it is advisable to use Minimap2 or AnchorWave i
nstead.

3. If you opt for AnchorWave, please ensure that a dependable gene annotation GFF3 file is accessible for
the reference genome and that there are no considerable rearrangements between the reference and query gen
omes.

In the case of our rice data, we opted to employ MUMmer4 (with default parameters) to conduct alignments b
etween the reference genome (MH63) and each of the 48 query genomes.
Here's an example:

`/bin/mummer-4.0.0rc1/nucmer -t 4 --prefix=MH63vsBasm ./genome/MH63 ./genome/Basm`

If you would like to use other alternative aligners, you may refer the following commands:
(Minimap2)

minimap2 -c --cs=long ./genome/MH63 ./genome/Basm > MH63vsBasm.chr.paf

(LAST)

`lastdb -P4 -uNEAR MH63dat ./genome/MH63`\
`last-train -P4 --revsym -E0.05 -C2 MH63dat ./genome/Basm > MH63vsBasm.train`\
`lastal -p MH63vsBasm.train MH63dat Basm > MH63vsBasm.maf`

### Step 3: Converting the alignment output from a particular aligner to the axt alignment files
SVGAP incorporates a set of UCSC programs (e.g., axtChain, ChainNet, and netSyntenic developed by Jim Kent
) to preprocess alignment files before identifying genomic variants. Hence, it is necessary to convert ali
gnment files from other aligners with different formats, such as delta from MUMmer4, paf from minimap2, an
d maf from LAST/LASTZ/ AnchorWave, to the axt alignment format (The axt format is described in http://geno
me.ucsc.edu/goldenPath/help/axt.html) prior to processing them within the pipeline. To facilitate this con
version task, the SV-GAPS pipeline provides a script called 1_Convert2Axt.pl.

To see how to use the 1_Convert2Axt.pl script, you can run
`perl 1_Convert2Axt.pl --help`

Parameters:
--ali   [last|lastz|mummer|minimap2|anchorwave] [Required]

--input [path] to the folder store the alignment files  [Required]

--wk    [path] to a working folder save the output (default: current)   [Required]

--query verbose (produce query-as-reference alignment if invoke)        [Optional]

--tname [character] Name of the target genome assembly [default: Ref]   [Optional]

--qname [character] Name of the query genome assembly [default: Que]    [Optional]

--t     [int] Number of threads to run this script (default: 4) [Optional]

--help  Print this help information     [Optional]

In our case of rice data, we saved all delta alignment files produced by MUMmer4 to a single output folder
 "rice/alignment" and created a working folder "rice/chainnet/Target_Ref".
Rearrange output files

`mkdir  ./rice/alignment`\
`mkdir ./rice/chainnet`\
`mv *.delta ./rice/alignment`\

Run 1_Convert2Axt.pl

`perl ~/bin/SV-GAPS/1_Convert2Axt.pl -ali mummer -input ./rice/alignment -wk ./rice/chainnet`

Please note that if you split query genomes, you should provide genome index (tab-delimited form) named "g
enomeID.sizes" in the working folder where the axt file will be saved before running 1_Convert2Axt.pl. Aft
er running the above commands, the resulting axt files have been saved in the "rice/chainnet" folder and r
eady for further processing.


### Step 4: Constructing chains and nets from pairwise alignments
To accurately identify genetic variants between the reference and query genomes, the raw alignments (axt f
ile) must be processed for synteny, and paralogous regions need to be removed. This ensures that only the
orthologous regions are retained for analysis. This task can be accomplished by employing a set of UCSC pr
ograms (find necessary information from http://genomewiki.ucsc.edu/index.php/Whole_genome_alignmen_howto)
that are sequentially invoked in the script 2_ChainNetSyn.pl within the SV-GAPS pipeline which will genera
te a filtered net file for syntenic alignments only and a single-coverage (on the reference genome) collec
tion of pairwise alignments.

To see how to use the 2_ChainNetSyn.pl script, you can run:

`perl 1_Convert2Axt.pl --help`

Parameter:

--ad    Forder contains the axt alignment files [REQUIRED]

--gd    Forder contains all genome sequences    [REQUIRED]

--lst   a file list the genome IDs (one line per genome)        [REQUIRED]

--sing  [verbose] if to report the sing-coverage maf file       [Optional]

--linearGap     linearGap file for ChainNet     [Optional]

--minscore      -minScore = 1000 (default)      [Optional]

--syn   Output folder for raw synnet files default:synnet       [Optional]

--help  print this help information     [Optional]

In our case we run the command:

`perl ~/bin/SV-GAPS/2_ChainNetSyn.pl --gd ../genome/ --ad chainnet/Target_Ref/ --lst ../genome/g.lst --sin
g`

The filtered net files for syntenic alignments will be saved in ./rice/synnet/.

The sing-coverage alignment maf files will be saved in ./rice/singRaw/.


In most cases, especially when comparing two highly divergent genomes, it is recommended to utilize the ad
ditional script 3_SynNetFilter.pl to filter the syntenic alignment net files. This script enhances precisi
on while slightly reducing sensitivity, resulting in improved accuracy for such comparisons.

To see how to use the 3_SynNetFilter.pl script, you can run

`perl 1_Convert2Axt.pl --help`

Parameter:

Required:

--synnet        Folder kept the Chain/Net/Synnet files

--gd    Folder kept genome sequences

--chain Folder kept the .chain files

Optional:

--one2mul       List of sequences in query that were aligned to the reference

--rawsv output of the syntenic gaps information

--alrate        the minimum alignable proportion for chains to be processed

--toplen        the minimum length for the top syntenic chains to be processed

--nonsynlen     the minimum length for the top non-syntenic chains to be processed

--synlen        the minimum length for any lower chains that is to be processed

--invlen        the minimum length for any chains that is to be considered as inversion

--qFarCutoff    the maximum number of the qFar value in the net file when predicting an inversion event

--help  print this help information

In our case we run the command:

`perl ~/bin/SV-GAPS/3_SynNetFilter.pl --synnet synnet/ --gd ../genome/ --chain chainnet/Target_Ref`

The final filtered net files for syntenic alignments are saved in ./rice/CleanSynNet/ which will be used f
or SV calling. The cleaned single-coverage alignment maf files are saved in ./rice/filtered_sing/ which wi
ll be used for SNP calling and SV genotyping.

### Step 5: Identification of structural variants (SVs)
The SV-GAPS pipeline infers SVs directly from the final filtered net files (the net format is described in
 http://genome.ucsc.edu/goldenPath/help/net.html). Currently, five common SV types between pairwise genome
 comparisons can be routinely identified using the script 4_PairwiseSV.pl within the SV-GAPS pipeline, inc
luding insertions, deletions, tandem duplications, inversions, and translations. Additionally, local highl
y divergent regions (i.e. with gaps on both target and query sequences) are defined as a complex type.

Try:

`perl 4_PairwiseSV.pl --help for more information`

Parameter:

Required:

--syndir        folder kept the Chain/Net/Synnet files

--outdir        output folder [default: current folder]

--gd    folder kept the genome sequences

--t     Number of theads to run this script (default: 4)

Optional:

--invlen        the minimum size of inversions to report [10,000]

--cnvlen        the minimum size of CNVs (tandem duplications) to report [500]

--dellen        the maximum size of deletions to report [20,000]

--trllen        the minumum size of translocation to report [5,000]

--filltoplen    the minimum size of the top fills to be considered [20,000]

--fill2len      the minimum size of the secondary level fills to be considered [10,000]

--rate  the minimum alignable proportion used to call a local duplication [0.25]

--help  print this help information

Below is the command we run for the rice data:

`perl ~/bin/SV-GAPS/4_PairwiseSV.pl --syndir ./rice/CleanSynNet --gd ./rice/genome/ --outdir ./rice/CleanS
V`

The detected SV calls are stored in the directory "./rice/CleanSV/". For each pairwise comparison, two set
s of SVs are reported, with coordinates that correspond to both the reference and query genomes.


### Step 6: Merge and remove redundant SV calls
Since SVs are identified between the reference genome and each of the query genomes, it is essential to me
rge the SVs from all comparisons. This merging process aims to generate a non-redundant dataset of SV call
s that encompasses all the comparisons. By doing so, each SV can be further genotyped across all individua
ls involved in the study, enabling construction of the pan-genome. The script 5_Combined.pl within the SV-
GAPs pipeline can work for this task.

Try:

`perl 5_Combined.pl --help` for more information.

Parameter

--sv    folder kept pairwise SVs output [REQUIRED]

--refname       reference name [e.g. ISO1,B73, et al...]        [REQUIRED]

--outdir        output directory [default: CombinedSV]  [Optional]

--ins_tol       the maximum distance tolerance for two insertion events that are considered or tested if t
hey are the same event [default : 13]   [Optional]

--del_ov        the minimum overlap proportion for two deletion events that are considered or tested if th
ey are the same event [default : 0.90]  [Optional]

--cant  How many neighbor SV events used to test together if they belong to the same event, should increas
e when your sample number is large [default : 20]       [Optional]

--help  print this help information     [Optional]

We run the following command to generate the non-redundant SV calls for all 48 pairwise companions.

`perl ~/bin/SV-GAPS/5_Combined.pl --sv ./rice/CleanSV --refname MH63`

Executing the above command will generate a folder called "./rice/CombinedSV/". Within this folder, five t
ypes of non-redundant SVs will be generated. It is worth mentioning that the insertions and deletions have
 been categorized into two classes according to their lengths; specifically, events that are equal to or l
onger than 50 bp (> 50 bp) are in one class, while the rest are in another.


### Step 7: Re-genotyping the merged non-redundant SV calls across all samples based on pairwise (single c
overage) alignments
Although the above step has provided preliminary genotyping information, the information still lacks compl
eteness. For instance, the detection of SVs was limited to specific individuals, leaving uncertainty as to
 whether the absence of these SVs in other individuals is genuine or merely a result of inaccurate sequenc
ing and assembling, or sequence loss in those regions. Therefore, additional methods are necessary to achi
eve complete genotyping information for all SVs across all samples. The SV-GAPS utilizes pairwise single c
overage alignment files to enhance the genotyping of each SV by extracting and comparing the local alignme
nt in order to determine the SV status in each sample. The scripts 6_DELGenotyping.pl and 7_INSGenotyping.
pl within the SV-GAPS pipeline are utilized to perform genotyping for deletions and insertions (with a len
gth greater than or equal to 50 base pairs), respectively.

#### Genotyping Deletions (≥ 50 bp)
Try `perl 6_DELGenoptyping.pl --help` for more information.

Parameter

--del   the file of the merged non-redundant deletions (> 49 bp)        [REQUIRED]

--singmaf       folder kept single-coverage maf files (on the reference side)   [REQUIRED]

--refseq        reference sequence      [REQUIRED]

--refsize       chromosome sizes of the reference       [REQUIRED]

--refname       reference name, e.g. hg38, IRGSP1       [REQUIRED]

--speciesname   species name, e.g. human, Drosophila, rice , maize      [Default: Undefined]

--ucsc  path to Kent's utilities        /$Bin/pub/UCSC/

--tba   path to TBA program utilities   /$Bin/pub/TBA/

--extend        the length of sequence around the deletion breakpoints for genotyping   [Default: 1000]

--outdir        set an output dir       [Default: VCFtemp]

--t     how many threads will use       [Default: 4]

--help  print this help information

We run the tested rice data with the following command:

`perl ~/bin/SV-GAPS/6_DELGenoptyping.pl --del ./rice/CombinedSV/All.DELs.50bplarge.bed.combined.sorted.txt
 --singmaf ./rice/filtered_sing/ --refseq ./rice/genome/MH63 --refsize ./rice/genome/MH63.sizes --refname
MH63 -t 12`

This step will result in a vcf file named All.DELs.50bplarge.bed.combined.sorted.txt.vcf in the folder whe
re kept the All.DELs.50bplarge.bed.combined.sorted.txt file.

#### Genotyping Insertions (≥ 50 bp)
Try `perl 7_INSGenotyping.pl --help` for more information.

Parameter

--ins   the file of the merged non-redundant insertions (> 49 bp)       [REQUIRED]

--query folder kept all query genome sequences  [REQUIRED]

--singmaf       folder kept single-coverage maf files (on the reference side)   [REQUIRED]

--refsize       chromosome/contig sizes from the reference      [REQUIRED]

--refseq        reference genome sequence       [REQUIRED]

--stretcher     path to stretcher       [REQUIRED]

--refname       reference name, e.g. hg38, IRGSP1       [REQUIRED]

--t     how many threads used [default: 4]      [Optional]

--speciesname   species name, e.g. human, Drosophila, rice, maize       [Optional]

--ucsc  path for kent's programs        [Optional]

--outdir        set the output dir name [Optional]

--help  print this help information     [Optional]

We run the tested rice data with the following command:

`perl ~/bin/SV-GAPS/7_INSGenoptyping.pl --ins ./rice/CombinedSV/All.INTs.50bplarge.bed.combined.sorted.txt
 --query ../genome/ -sing filtered_sing/ --stretcher ../../SV-GAPS/pub/Stretcher/ --refseq ../genome/MH63
--refsize ../genome/MH63.sizes --refname MH63`

This step will result in a vcf file named All.INTs.50bplarge.bed.combined.sorted.txt.vcf in the folder whe
re kept the All.INTs.50bplarge.bed.combined.sorted.txt file.

### Step 8: Re-genotyping InDels
Small deletions and insertions (< 50 bp) are separately genotyped with 8_InDelGenotyping.pl.
Try `perl 8_InDelGenotyping.pl --help` for more information.

Parameter

--indel raw merged InDels file  [Required]

--indeltype     InDel types [ins|del]   [Required]

--singmaf       path to pairwise genome algnment maf files      [Required]

--refseq        reference genome sequence       [Required]

--refsize       reference chromosome sizes file [Required]

--refname       reference name, e.g. hg38, IRGSP1       [Required]

--query folder kept all query genomes   [Required]

--t     how many threads to use [default:4]     [Required]

--speciesname   species name, e.g. human, Drosophila, rice, maize       [Optional]

--outdir        output dir      [Optional]

--head  print VCF header if invoked     [Optional]

--help  print this help information     [Optional]

#### Genotyping small deletions in our case:

`perl ~/bin/SV-GAPS/8_InDelGenoptyping.pl --indel ./rice/CombinedSV/All.DELs.50bpsmall.bed.combined.sorted
.txt --indeltype del --refseq ./rice/genome/MH63 --refsize ./rice/genome/MH63.sizes --refname MH63 --query
 ./rice/genome/ --singmaf ./rice/filtered_sing/ --outdir MH63_del --t 12`

This step will result in a vcf file named All.DELs.50bpsmall.bed.combined.sorted.txt.vcf in the folder whe
re kept the All.DELs.50bpsmall.bed.combined.sorted.txt file.

#### Genotyping small insertions in our case:

`perl ~/bin/SV-GAPS/8_InDelGenoptyping.pl --indel ./rice/CombinedSV/All.INTs.50bpsmall.bed.combined.sorted
.txt --indeltype ins --refseq  ./rice/genome/MH63 --refsize ./rice/genome/MH63.sizes --refname MH63 --quer
y ./rice/genome/ --singmaf ./rice/filtered_sing/ --outdir MH63_ins --t 12`

This step will result in a vcf file named All.INTs.50bpsmall.bed.combined.sorted.txt.vcf in the folder whe
re kept the All.INTs.50bpsmall.bed.combined.sorted.txt file.

### Step 9: Genotyping SNPs
The SV-GAPs also offer a script called 9_SNPgenotyping.pl, which allows for the direct calling of SNPs fro
m the single-coverage maf files. It generates a vcf file that includes all samples and their respective SN
P information.

Try `perl 9_SNPgenotyping.pl --help` for more information.

Parameter

--refseq        reference sequence      [REQUIRED]

--refsize       reference sizes [REQUIRED]

--refname       reference name, e.g. hg38, IRGSP1       [REQUIRED]

--singmaf       path to single coverage maf file        [REQUIRED]

--t     how many processors will invoke [default: 4]

--speciesname   species name, e.g. human, Drosophila, rice, maize       [default: Unknown]

--outdir        define an output dir    [default: VCFtemp]

--help  print this help information

The command we used for the test data is:

`perl ~/bin/SV-GAPS/9_SNPgenotyping.pl --refseq ./rice/genome/MH63 --refsize ./rice/genome/MH63.sizes --re
fname MH63 --singmaf ../rice/filtered_sing/ --outdir SNPVCF --t 8`

The SNP genotyped file All.SNP.vcf will report in the SNPVCF folder.

### Step 10: Annotating Structural Variants
When we acquire precise catalogs of structural variants, it is crucial to elucidate the potential mechanis
ms driving the formation of SVs and their functional impacts. We offer a script called 10_SVannotation.pl
to infer the primary common factors in plant genomes, such as tandem duplication, inverted tandem duplicat
ion, transposon elements (TEs) insertion, nest TEs insertion, and complete/partial gene duplication, that
may be involved in SV formation.

Try `perl 10_SVannotation.pl --help` for more information.

Parameter

--ref   the reference genome                                            [REQUIRED]


--vcf   the SV vcf file                                                 [REQUIRED]


--svtype        INS or DEL                                                      [REQUIRED]


--te    a TE liabrary (e.g. output from EDTA pipeline)                  [REQUIRED]

--cds   a cds file from this species or a closely related species       [REQUIRED]

--out           annotated vcf file                                              [REQUIRED]


--help  print this help information

The command we used for the test data is:

`perl ~/bin/SV-GAPS/10_SVannotation.pl --ref ./rice/genome/MH63 --vcf All.INTs.50bplarge.bed.combined.sort
ed.txt.vcf --svtype INS --te TE.EDTA.lab.fa --cds MSU7.cds.fa --out INS.annotated.vcf`

`perl ~/bin/SV-GAPS/10_SVannotation.pl --ref ./rice/genome/MH63 --vcf All.DELs.50bplarge.bed.combined.sort
ed.txt.vcf --svtype DEL --te TE.EDTA.lab.fa --cds MSU7.cds.fa --out INS.annotated.vcf`

The annotated .vcf file(s) will report in the current folder.
