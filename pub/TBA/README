--modified on 11/12/08
Add refsrc=src option to maf2fasta.  The default is to use the
source of the first component of the first alignment block in the
maf file.

-- See README2 for more updated programs.  10/12/05

-- big revisions 1/13/05

--modified on 11/23/04
add a switch to maf_order to keep single-row alignment or not.

--modified on 11/22/04
added program descriptions.

--modified on 11/02/04
maf_order is added. 

README version 10 -- modified on 11/01/04
add a small tool "get_standard_headers", which return in format of
"1:end:+:srcSize" assuming the sequence contig starts at position
1. The start position in fasta header is inclusive while the end
position is exclusive, so the srcSize is one larger than end. 

README version 10 -- modified on 10/26/04
-- maf2lav and maf2fasta allow the species with more than one contigs. 
Parameters are the same as before.    

README version 10 -- modified on 10/25/04
-- Multiz maf_Project maf_project_simple programs are removed.
-- all_bz multiz programs no longer need mapping file. 

README version 10 -- modofied on 10/09/04
-- Var_multiz program becomes Multiz version 10.
-- var_multiz program becomes multiz version 10. 


README version3 -- modified on 01/13/05
_________________________________________________________
---------------------------------------------------------
| 	Examples commands to run all_bz, tba, multiz	|
---------------------------------------------------------
all_bz "(((((human chimp) galago)(mouse rat)) chicken)(tetra (fugu zebrafish)))" blastz_spec_file
tba "(((((human chimp) galago)(mouse rat)) chicken)(tetra (fugu zebrafish)))" *.maf destination-file
multiz human.chimp.galago.maf human.mouse.rat.maf 1 
maf_project tba.maf zebrafish 

_________________________________________________________
---------------------------------------------------------
|		Some requirements			|
---------------------------------------------------------
*1. format of sequence header:

>string1:string2:int1:char:int3

Note:
string1 is usually species name, string2 is usually chromosome 
information. "string1.string2" is used as src field in maf struct.

int1 is start position(1 based, inclusive).
int3 is src size. 
char is +/-.

It's allowed not to contain standard header as long as there is 
only one contig in the sequence file.


*2. MSA header can also be accepted.
">${COMMON_NAME}|${ENCODE_REGION}|${FREEZE_DATE}|${NCBI_TAXON_ID}
|${ASSEMBLY_PROVIDER}|${ASSEMBLY_DATE}|${ASSEMBLY_ID}|${CHROMOSOME}
|${CHROMOSOME_START}|${CHROMOSOME_END}|${CHROM_LENGTH}|${STRAND}
|${ACCESSION}.${VERSION}|${NUM_BASES}|${NUM_N}|${THIS_CONTIG_NUM}
|${TOTAL_NUM_CONTIGS}|${OTHER_COMMENTS}

empty field shall be represented by ".".

*3. mafFile component positions start at 0 instead of 1. 

*4. species-guid-tree is used as arugment in many of the 
following programs, it consists of double quotes, parenthesis,
and species names, e.g.  "((HUMAN CHIMP)(RAT MOUSE))"

_________________________________________________________
---------------------------------------------------------
| 	 	program descriptions			|
---------------------------------------------------------

---------------------< all_bz >--------------------------
all_bz [-+] species-guid-tree [blastz_spec]

No mapping file is needed. all_bz is a wrapper for blastzWrapper,
it generates blastzWrapper/blastz commands for pairs of specified
sequences. Optional blastz-spec file contains command-lind
options for the blastz runs. 

A '+' tells all_bz to echo each blastz command to stdout 
before executing it.
A '-' tells all_bz to simply echo the commands 

---------------------< blastzWrapper >--------------------------
blastzWrapper seq-file1 seq-file2 [options]

[options] are the same ones as defined in "blastz" program.
blastzWrapper processes two sequence files each containing one or 
more contigs and runs blastz. Each sequence contig must 
follow the sequence header format described above.

---------------------< lav2maf >------------------------
lav2maf blastz-file seq-file1 seq-file2

The program transform a blastz output format file blastz-file
 into maf format. seq-file1 and seq-file2 are sequence files
used to run blastz.

---------------------< maf2lav >------------------------
maf2lav align.maf seq1 seq2

The program transform a maf format file align.maf into blastz 
output format. seq1 and seq2 are source sequences.

-------------------------------< maf2fasta >---------------------------------
maf2fasta refseq-file maf-file [beg end] [fasta[2]][?] [iupac2n] [refsrc=src]

The program transform a maf format file maf-file into fasta format with
reference to refseq-file. The maf-file shall not have inversion or overlap
regions. The maf-file shall be referenced(top row) by species of refseq-file.

[beg end] option limits the resulted fasta region to be within beg and end in
respect to the reference.

[fasta[2]][?]: The default output format is that used for text alignments by
MultiPipMaker; appending the command-line argument "fasta" or "fasta2"
requests FastA-format output; with "fasta2", alignment rows are split into
rows of length <= 50 (set by COL_WIDTH). To identify gaps as either "within
-alignment" or "between-alignment", append a character to the word "fasta" or
"fasta2" that will replace '-' between two local alignments, as in "fasta@"

[iupac2n]: If a nucleotide character other than one of "ACGTNacgtn" is found
in the refseq-file it is mapped to either "N" or "n" depending on its case.

[refsrc=src]: Specify that the component source for the reference species is
src.  By default, fasta2maf will assume that the source of the first
component of the first alignment block in the maf-file is the source of the
reference species.

--------------------< maf_order >----------------------
maf_order maf-file species1 species2 ... [nohead] [all] 

maf-file is the maf file to be ordered. It is followed
by species names need to be included in the ordered
file, the ordering of components in a maf block follows
the ordering of the species names in the argments.
Species not in the arguments are excluded.

[nohead]: if nohead is specified, maf header is not
shown in the result projection file.

[all]: if all is specified, single-row blocks are also
included, otherwise excluded from resulted file.

---------------------< maf_project >--------------------
maf_project maf-file reference [from to] [filename-for-other-mafs] [species-guid-tree] [nohead]

maf_project is able to process maf file where there might
be more than one contigs for any species.

Arguments:
reference: the sequence to which the maf file is projected. 
The result maf blocks always have reference sequence in the
top row, and maf blocks are ordered by the starting position
of the top row.

[from to]: the output maf blocks are limited to [from to] area
in respect to positions on reference sequence.

[filename-for-other-mafs]: collect maf blocks not contained
in projected file. WHEN THIS ARGUMENT IS NOT SPECIFIED, THE
PROJECTED MAF BLOCKS ARE BEAUTIFIED.

[species-guid-tree]: species not specified in this argument 
are screened out.

[nohead]: if nohead is specified, maf header is not shown in 
the result projection file.

---------------------< mafFind >-------------------------------
mafFind file.maf beg end [species-prefix] [slice]

mafFind finds mafs intersecting a particular interval. The mafs 
whose first row intersects positions beg-end are printed. For
non-reference species, a command like
      mafFind file.maf beg end mm3
asks for mafs that have a row where"mm3" is a prefix of the 
"source"(e.g., "mm3.chr7") and which intersects positions beg-end. 
Finally, the last argument "slice" asks that ends of the reported 
mafs be trimmed to make them precisely match beg-end, as in:
      mafFind file.maf beg end mm3 slice

---------------------< multiz version 10 >---------------------
multiz [R=?] [M=?] maf-file1 maf-file2 v [out1] [out2] [nohead]

multiz version 10 allows each species containing more than one 
contigs without reqirement for mapping-file. 
 
maf-file1 and maf-file2 are two maf files to be aligned, each 
topped by a same reference sequence. The alignment of reference 
sequence with other components might be just for purpose of 
determing approximate alignment between two files, thus the 
alignment might be fixed or not, this is specified by v value, 
which can be only 0 or 1.

0 - neither alignment of reference in each file is fixed.
1 - the alignment of reference in the first file is fixed.

[R=?] species radius values in dynamic programming, by default 30
[M=?] species minimum output width, by default 1, which means 
output all blocks.

[out1] collects unused blocks from maf-file1
[out2] collects unused blocks from maf-file2
[nohead] specifies not to have maf header for output


---------------------< pair2tb >------------------------
pair2tb pairwise.maf seq-file1 seq-file2

The program assumes the input file pairwise.maf does not 
have overlapping blocks. It also assumes the top components
correspond to seq-file1, and second components correspond to
seq-file2. seq-file1 and seq-file2 allow for more than one 
contig.


---------------------< single_cov2 >--------------------
single_cov2 pairwise.maf [F=deleted.maf]

This program removes overlapped regions from pairwise.maf.
mapping file is not needed. Optional [F=deleted.maf] specifies filename
collecting removed regions. 

---------------------< tba >----------------------------
tba [+-] [R=?] [M=?] species-guid-tree maf_source destination-file

tba program passes R M values to multiz program if any of 
them are provided. This version requires destination filename
to be provided.

arg '-' means no execute, '+' means verbose.

---------------------< maf_checkThread >-----------------
maf_checkThread projected-maf-file

this tool check the threading condition of the reference. It assumes the sequence starts at position 0. So there is one error when the sequence starts at a position other than 0. The alignment file must be projected onto reference. And the threading condition is checked for reference only. It has to be done for all species in TBA alignments to make sure the whole alignment satisifies threading condition.
