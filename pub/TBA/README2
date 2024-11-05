blastzWrapper replaced Blastz.
roast replaced autoMZ.            -- 05/07/2008 

This README2 shall replace the prvious README. -- 10/12/2005

Notice to previous users: 
1. Compiling is the same.
2. Input sequence files format is the same including naming issues.
3. New parameters are added to existing programs, and their default values shall give the same programs as before. e.g. the behavior of the programs with the same arguments you used before shall give the same result as before, otherwise please report to mhou@cse.psu.edu.
4. A new program "multic" is added intending to be new "multiz" eventually, and currently they have the same parameter list, see details below.
5. A new program called "roast" is added to run a series of multiz/multic to get reference-only multiple alignment by using a guide tree. roast and tba have similar parameters, see below for details.


Description and usage: values in parenthesis are defaults.
1. blastzWrapper 
blastzWrapper.v10:  -- wrapper of blastz, passing all arguments to blastz.
args: seqfile1 seqfile2 [options]

2. single_cov2
single_cov2.v11: -- screening out overlapped regions.
args: pairwise.maf [R=species] [F=deleted.maf]
By default, single coverage is done for both species in the pairwise alignment; if R=species specified, single coverage is done for the specified species only. Blocks do not need to be sorted, but the first rows of all blocks must be of the same species; the second rows of all blocks must be of the same species. This can be ensured by running maf_project, or maf_order, see below.

3. maf_project 
maf_project.v10:  -- extract maf-file entries that name a given reference sequence.
args: file.maf reference [from to] [filename-for-other-mafs] [species-guid-tree] [nohead]

4. maf_order 
maf_order.v10:  -- order rows according to a given list. Rows of species not specified are excluded.
args: maf-file species1 species2 .. [nohead] [all]
When "all" is not turned on, single-row blocks are not in output.

5. lav2maf
lav2maf.v11:  -- convert blastz output to maf file.
 args: blastz.output seq-file1 seq-file2

6. maf2lav
maf2lav.v10:  -- convert maf file to blastz output
 args = align.maf seq1 seq2

7. all_bz
all_bz.v12: -- generate all blastz commands for pairs of specified sequences.
args: [-+] [b=?] [A=?] [F=reference] [h=?] [q=?] [t=?] [D=?] species-guid-tree [blastz_specfile]
        +(off) verbose
        -(off) output command only.
        b(2) 0: run post-process only 1: run blastzWrapper only, transform to maf 2: run both
        A(1) 0: toast 1: single_cov2
        F(null) null: single coverage is done for both species; reference: single coverage is done for reference only, effective in single_cov2
        h(300) minimum chaining size, effective in toast
        q(600) minimum cluster size, effective in toast
        t(0) 0: full size for toast 1: reduced size for toast.
        D(1) 0: run all_bz for roast 1: run all_bz for TBA.

8. multiz
multiz.v10.6:  -- aligning two files of alignment blocks where top rows are always the reference, reference in both files cannot have duplicats
args: [R=?] [M=?] file1 file2 v? [out1 out2] [nohead] [all]
        R(30) radius in dynamic programming.
        M(1) minimum output width.
        out1 out2(null) null: stdout; out1 out2: file names for collecting unused input.
        nohead(null) null: output maf header; nohead: not to output maf header.
        all(null) null: not to output single-row blocks; all: output all blocks.

9. multic
multic.v11.1:  -- aligning two files of alignment blocks where top rows are always the reference, reference in both files can contain duplicats (different to multiz)
args: [R=?] [M=?] file1 file2 v? [out1 out2] [nohead] [all]
        R(30) radius in dynamic programming.
        M(1) minimum output width.
        out1 out2(null) null: stdout; out1 out2: file names for collecting unused input.
        nohead(null) null: output maf header; nohead: not to output maf header.
        all(null) null: not to output single-row blocks; all: output all blocks.

10. tba 
tba.v12: TBA -- threaded block alignment.
args: [+-] [R=?] [M=?] [E=?] [P=?] [X=?] species-guid-tree maf-source destination
        R(30) dynamic programming radius.
        M(1) minimum block length of output.
        E(null) null: no reference centric alignment, single coverage is guaranteed for every species; reference: refernece centric alignment, singe coverage is guaranteed for reference species.
        P(null) null: run multiz; P=multic specifies to run multic.
        X(0) utilize maf files with different suffix from differnt post processing.
                0: .sing.maf from single coverage pairwise alignment
                1: .toast.maf from full size toast
                2: .toast2.maf from reduced size toast

11. roast 
roast.v1: roast -- reference guided multiple alignment.
args: [+-] [R=?] [M=?] [P=?] [X=?] E=reference-species species-guid-tree maf-source destination
        R(30) dynamic programming radius.
        M(1) minimum block length of output.
        P(multiz) multiz: single coverage for reference row multic: no requirement on single coverage.
        X(0) utilize maf files with different suffix from differnt post processing.
                0: .sing.maf from single coverage pairwise alignment
                1: .toast.maf from full size toast
                2: .toast2.maf from reduced size toast

12. maf_checkThread -- Check threading condition or check overlap of the reference row in alignment. maf_project needs to be run first.
args: projected-maf-flie [overlap]
By default the threading condition is checked; if "overlap" is specified, only overlapping is checked.

Examples
1. ENCODE multiple alignment with reference to human
1) using tba
	all_bz + F=human guid-tree blastz-spec-file
	tba + E=human P=multic guid-tree maf-source destination
With above commands, all_bz produces pairwise alignments using blastz, these alignments are in files with postfix of ".orig.maf"; and then all_bz reads files with ".orig.maf" and runs single_cov2 with reference of human, which means the single coverage property is only enforced on human when human is aligned to other species; while for pairwise alignment between two species other than human, reciprocal single coverage is enforced. This step produces pairwise alignment files with postfix of ".sing.maf". The purpose for keeping blastz output is to avoid running blastz repeatedly when the pairwise post processing program needs to be changed, e.g. choose single_cov2, reference-specified single_cov2, or TOAST in future. In those cases, you can run all_bz with "b=0" to specify not to run blastz again as long as the sequence files are the same and parameters for blastz are the same.

tba in this example utilizes multic. In tba, species including human are used as references to guide the alignment procedure. Since single coverage is only ensured for human, not for other species, multic instead of multiz shall be used.

2) using roast 
	all_bz + F=human guid-tree blastz-spec-file  // can be skipped
	roast + E=human guid-tree maf-source destination
The all_bz procedure can be skipped if all_bz is already done for running tba. But it does not hold for tba if roast is run first, since tba requires more pairwise alignments. 

By default, roast uses multiz instead of multic. roast only utilizes one reference, that is human in this case. Input alignments with human are ensured single coverage on human, so multiz can be used.
