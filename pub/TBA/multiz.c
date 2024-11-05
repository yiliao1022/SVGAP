/* multiz.c   version 11.2
 *
 * Change argument row2 to all, make it default not to output
 * single-row blocks. -- 10/06/05
 * Add option not to output single-row blocks. -- 09/28/05
 *
 * most recent modification: previously not to break some large
 * blocks when aligned to very small blocks, the unused small blocks
 * overlap with large blocks when they are both output to the same
 * file. Fixed by enforcing aligning for any size. Version output
 * are made consistent. -- 09/26/05
 *
*  when input block only contain few bases
*  length, it was ignored and regarded as unused. If both input
*  files contain the short blocks at the same positions, overlapping
*  happens. Fixed by enforcing aligning for any size. -- 09/15/05
*
*  Modified on Dec. 22, if [out1][out2] are specified, unused
*  blocks are output to out1 and out2, if not specified,
*  unused blocks are output to stdout.
*
*  Modified on Aug. 31 to fix the bug where blocks might
*  be lost when column width is below a certain value.
*
*  Modified on Aug. 31 to allow for arguments of R and M.
*
*  Variant to multiz program. It aligns two files of
*  alignment blocks where top row is always the reference,
*  assuming blocks are increasing ordered based on the
*  start position on the refernece seqence. Single-coverage
*  on reference is required at this stage.
*
*  Three arguments are required: char* arg1, char* arg2,
*  int v. arg1 and arg2 are two files need to be aligned
*  together. The blocks are always topped by a reference,
*  for being used to determine the approximate alignment,
*  but the referece is fixed in that block or not, depending
*  on the argument of reference value of v:
*     0:  neither is fixed
*     1:  the first block is fixed
*/

#include <stdio.h>

#include <stdlib.h>
#include "maf.h"
#include "util.h"
#include "multi_util.h"
#include "mz_preyama.h"
#include "mz_scores.h"

static const char rcsid[] = "$Id: multiz.c 148 2009-01-21 22:30:05Z rico $";

extern int row2;
extern int radius;
extern int MIN_OUTPUT_WID;
extern int LRG_BREAK_WID;
extern int SML_BREAK_WID;

#define VERSION 11.2

int multiz(struct mafAli** wk_list1, struct mafAli** wk_list2, FILE* fpw1, FILE* fpw2, int v) {
    struct mafAli *a1, *a2, *new_ali;
    int beg1, end1, beg2, end2, beg, end, col_beg, col_end;
    int test=0;

    a1 = retrieve_first(wk_list1);
    a2 = retrieve_first(wk_list2);

    while (1) {
        while ( a1 != NULL && ( a2 == NULL || a1->components->start + a1->components->size -1 < a2->components->start)) {
            if ( a1->components->size >= MIN_OUTPUT_WID && fpw1 != NULL && (row2==0 || a1->components->next!=NULL))
                mafWrite(fpw1, a1);
            mafAliFree(&a1);
            a1 = retrieve_first(wk_list1);
        }
        while ( a2 != NULL && ( a1 == NULL || a2->components->start + a2->components->size -1 < a1->components->start)) {
            if ( a2->components->size >= MIN_OUTPUT_WID && fpw2 != NULL && (row2==0 || a2->components->next!=NULL))
                mafWrite(fpw2, a2);
            mafAliFree(&a2);
            a2 = retrieve_first(wk_list2);
        }
        if ( a1 == NULL && a2 == NULL)
            break;
        if ( a1 == NULL || a2 == NULL)
            continue;
        if ( a1->components->start + a1->components->size - 1 < a2->components->start)
            continue;
        if ( a2->components->start + a2->components->size - 1 < a1->components->start)
            continue;

        if ( a1->components->start == 11305 )
            test++;


        beg1 = a1->components->start;                           // at this point, a1 a2 overlap or
        end1 = a1->components->start + a1->components->size - 1;// cover, print uncovered front part
        beg2 = a2->components->start;                           // pre_yama middle/covered part,
        end2 = a2->components->start + a2->components->size - 1;// then add another ali to process end part

        /* // the unused positions overlap with the unbroken part
        if ( beg1 > beg2 && end1 < end2 && end2-beg2+1 > LRG_BREAK_WID && end1-beg1 + 1 < SML_BREAK_WID) { // not to break a2
          if ( a1->components->size >= MIN_OUTPUT_WID && fpw1 != NULL )
        mafWrite(fpw1, a1);
          mafAliFree(&a1);
          a1 = retrieve_first(wk_list1);
          continue;
        }

        if ( beg2 > beg1 && end2 < end1 && end1-beg1+1 > LRG_BREAK_WID && end2-beg2 + 1 < SML_BREAK_WID) { // not to break a1
          if ( a2->components->size >= MIN_OUTPUT_WID && fpw2 !=NULL)
            mafWrite(fpw2, a2);
          mafAliFree(&a2);
          a2 = retrieve_first(wk_list2);
          continue;
        }
        */

        if ( beg1 < beg2 && beg2-beg1 >= MIN_OUTPUT_WID && fpw1 != NULL ) {
            col_beg = mafPos2Col(a1->components, beg1, a1->textSize);
            for (; col_beg>0 && a1->components->text[col_beg-1]=='-'; col_beg--)
                ;
            col_end = mafPos2Col(a1->components, beg2-1, a1->textSize);
            for (; col_end<a1->textSize-1 && a1->components->text[col_end+1]=='-'; col_end++)
                ;
            print_part_ali_col(a1, col_beg, col_end, fpw1);
        } else if ( beg2 < beg1 && beg1-beg2 >= MIN_OUTPUT_WID && fpw2 != NULL) {
            col_beg = mafPos2Col(a2->components, beg2, a2->textSize);
            for (; col_beg>0 && a2->components->text[col_beg-1]=='-'; col_beg--)
                ;
            col_end = mafPos2Col(a2->components, beg1-1, a2->textSize);
            for (; col_end<a2->textSize-1 && a2->components->text[col_end+1]=='-'; col_end++)
                ;
            print_part_ali_col(a2, col_beg, col_end, fpw2);
        }


        beg = beg1 > beg2 ? beg1 : beg2;
        end = end1 < end2 ? end1 : end2;

        if ( beg == beg1 ) { // for gaps in front
            col_beg = mafPos2Col(a1->components, beg1, a1->textSize);
            if ( col_beg != 0 && fpw1 != NULL)
                print_part_ali_col(a1, 0, col_beg-1, fpw1);
        }
        if ( beg == beg2) { // for gaps in front
            col_beg = mafPos2Col(a2->components, beg2, a2->textSize);
            if ( col_beg != 0 && fpw2 != NULL)
                print_part_ali_col(a2, 0, col_beg-1, fpw2);
        }
        new_ali = pre_yama(a1, a2, beg, end, radius, v, fpw2);

        if ( new_ali != NULL && new_ali->components->size >= MIN_OUTPUT_WID)
            mafWrite(stdout, new_ali);
        mafAliFree(&new_ali);

        if ( end1 < end2)
            a2 = keep_ali(a2, end1+1);

        if ( end2 < end1)
            a1 = keep_ali(a1, end2+1);

        if ( end1 <= end2 ) {
            col_end = mafPos2Col(a1->components, end1, a1->textSize);
            if ( col_end < a1->textSize-1 && fpw1 != NULL)
                print_part_ali_col(a1, col_end+1, a1->textSize-1, fpw1);
            mafAliFree(&a1);
            a1 = retrieve_first(wk_list1);
        }
        if ( end2 <= end1 ) {
            col_end = mafPos2Col(a2->components, end2, a2->textSize);
            if ( col_end < a2->textSize-1 && fpw2 != NULL)
                print_part_ali_col(a2, col_end+1, a2->textSize-1, fpw2);
            mafAliFree(&a2);
            a2 = retrieve_first(wk_list2);
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    char USAGE[10000];
    struct mafFile *maf1, *maf2;
    struct mafAli *ali, *cp_list1, *cp_list2, *wk_list1, *wk_list2;
    int x, nohead=0, v, i;
    FILE *fpw1=NULL, *fpw2=NULL;
    char ref_chr[200], cmd[200], args[2000];

    sprintf(cmd, "multiz.v%.1f", VERSION);
    argv0 = cmd;

    strcpy(USAGE, "args: [R=?] [M=?] file1 file2 v? [out1 out2] [nohead] [all]\n");
    strcat(USAGE, "\tR(30) radius in dynamic programming.\n");
    strcat(USAGE, "\tM(1) minimum output width.\n");
    strcat(USAGE, "\tout1 out2(null) null: stdout; out1 out2: file names for collecting unused input.\n");
    strcat(USAGE, "\tnohead(null) null: output maf header; nohead: not to output maf header.\n");
    strcat(USAGE, "\tall(null) null: not to output single-row blocks; all: output all blocks.\n");

    strcpy(args, cmd);
    strcat(args, " ");
    for (i=1; i<argc; i++) {
        strcat(args, argv[i]);
        strcat(args, " ");
    }

    // ---- process arguments
    while (argc > 1 && strchr("RMLS", (x=argv[1][0])) && argv[1][1] == '=') {
        if ( x== 'R') {
            radius = atoi(argv[1]+2);
            if (radius < 0)
                fatal("radius cannot be negative");
        } else if ( x=='M') {
            MIN_OUTPUT_WID = atoi(argv[1]+2);
            if (MIN_OUTPUT_WID < 0)
                fatal("MIN_OUTPUT_WID cannot be negative");
        } else if ( x=='L') {
            LRG_BREAK_WID = atoi(argv[1]+2);
            if (LRG_BREAK_WID < 0)
                fatal("LRG_BREAK_WID cannot be negative");
        } else if ( x=='S') {
            SML_BREAK_WID = atoi(argv[1]+2);
            if (SML_BREAK_WID < 0)
                fatal("SML_BREAK_WID cannot be negative");
        } else
            fatalf("illegal flag %c", x);
        ++argv;
        --argc;
    }

    if ( strcmp(argv[argc-1], "all") == 0 ) {
        row2=0;
        argc--;
    }

    if ( strcmp(argv[argc-1], "nohead")==0 ) {
        argc--;
        nohead = 1;
    }

    if ( argc != 4 && argc !=6)
        fatalf(" -- aligning two files of alignment blocks where top rows are always the reference, reference in both files cannot have duplicats\n%s", USAGE);

    if ( argc == 6) {
        fpw1 = fopen(argv[4], "w");
        fpw2 = fopen(argv[5], "w");
    } else
        fpw1 = fpw2 = stdout;

    v = atoi(argv[3]);   // v is used later to determine the way to use yama
    if ( v != 0 && v != 1 )
        fatal("v can only be value of 0, 1 ");

    // ---- end of processing arguments

    if ( nohead == 0) {
        mafWriteStart(stdout, "multiz");
        printf("# %s\n", args);
    }
    init_scores70();

    maf1 = mafReadAll(argv[1], 1);
    maf2 = mafReadAll(argv[2], 1);
    cp_list1 = maf1->alignments;
    cp_list2 = maf2->alignments;
    maf1->alignments = maf2->alignments = NULL;
    mafFileFree(&maf1);
    mafFileFree(&maf2);

    while ( cp_list1 != NULL && cp_list2 != NULL) {
        wk_list1 = wk_list2 = NULL;
        strcpy(ref_chr, cp_list1->components->src);

        seperate_cp_wk(&cp_list1, &wk_list1, ref_chr);
        seperate_cp_wk(&cp_list2, &wk_list2, ref_chr);

        multiz(&wk_list1, &wk_list2, fpw1, fpw2, v);
    }

    if ( cp_list1 != NULL && fpw1 != NULL)
        for (ali=cp_list1; ali!=NULL; ali=ali->next)
            if ( row2==0 || ali->components->next != NULL )
                mafWrite(fpw1, ali);

    if ( cp_list2 != NULL && fpw2 != NULL)
        for (ali=cp_list2; ali!=NULL; ali=ali->next)
            if ( row2==0 || ali->components->next != NULL )
                mafWrite(fpw2, ali);

    if (fpw1 != NULL)
        fclose(fpw1);
    if ( fpw2 != NULL)
        fclose(fpw2);

    mafWriteEnd(stdout);
    return 0;
}
