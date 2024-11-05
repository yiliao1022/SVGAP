/*
 *  get_covered.c version 10
 */

#include "util.h"
#include "multi_util.h"
#include "maf.h"
#include "mz_scores.h"

static const char rcsid[] = "$Id: get_covered.c 142 2008-11-12 18:55:23Z rico $";


int get_covered(struct mafAli** wk_list1, struct mafAli** wk_list2) {
    struct mafAli *a1, *a2;
    int beg1, beg2, beg, end1, end2, end;

    init_scores70();
    a1 = retrieve_first(wk_list1);
    a2 = retrieve_first(wk_list2);

    while (1) {
        while ( a1 != NULL && a2 != NULL && a1->components->start + a1->components->size-1 < a2->components->start) {
            //mafWrite(stdout, a1);
            mafAliFree(&a1);
            a1 = retrieve_first(wk_list1);
        }
        while ( a1 != NULL && a2 != NULL && a2->components->start + a2->components->size-1 < a1->components->start) {
            mafAliFree(&a2);
            a2 = retrieve_first(wk_list2);
        }
        if ( a1==NULL )
            return 0;
        if ( a2==NULL ) {
            //mafWrite(stdout, a1);
            mafAliFree(&a1);
            while ( (a1=retrieve_first(wk_list1)) != NULL ) {
                //mafWrite(stdout, a1);
                mafAliFree(&a1);
            }
            return 0;
        }
        if ( a1->components->start + a1->components->size - 1 < a2->components->start)
            continue;
        if ( a2->components->start + a2->components->size - 1 < a1->components->start)
            continue;

        beg1 = a1->components->start;                           // at this point, a1 a2 overlap or
        end1 = a1->components->start + a1->components->size - 1;// cover, print uncovered front part
        beg2 = a2->components->start;                           // pre_yama middle/covered part,
        end2 = a2->components->start + a2->components->size - 1;// then add another ali to process end part

        //if ( beg1 < beg2 )
        //print_part_ali(a1, beg1, beg2-1, stdout);

        beg = beg1 > beg2 ? beg1 : beg2;
        end = end1 < end2 ? end1 : end2;

        print_part_ali(a1, beg, end, stdout);

        if ( end1 < end2 ) {                       //  --     -    -
            mafAliFree(&a1);                          //   --   ---   --
            a1 = retrieve_first(wk_list1);
            a2 = keep_ali(a2, end1+1);
        } else if ( end2 < end1 ) {
            mafAliFree(&a2);
            a2 = retrieve_first(wk_list2);
            a1 = keep_ali(a1, end2+1);
        } else {
            mafAliFree(&a1);
            mafAliFree(&a2);
            a1 = retrieve_first(wk_list1);
            a2 = retrieve_first(wk_list2);
        }
    }
    return 0;
}

// second file contain only discontinuous block for a single seq
// or for its top seq
// first file's top seq is the same sa the second file
// seq may have different contigs
int main(int argc, char** argv) {
    char cur_chr[200];
    struct mafFile *maf1, *maf2;
    struct mafAli *cp_list1, *cp_list2, *wk_list1, *wk_list2;

    if ( argc != 3 )
        fatal("arguments: file1 file2");

    mafWriteStart(stdout, "get_covered");

    maf1 = mafReadAll(argv[1], 1);
    maf2 = mafReadAll(argv[2], 1);

    cp_list1 = maf1->alignments;
    cp_list2 = maf2->alignments;
    maf1->alignments = maf2->alignments = NULL;
    mafFileFree(&maf1);
    mafFileFree(&maf2);

    while ( cp_list1 != NULL && cp_list2 != NULL ) {
        wk_list1 = wk_list2 = NULL;
        strcpy(cur_chr, cp_list2->components->src);

        seperate_cp_wk(&cp_list2, &wk_list2, cur_chr);
        seperate_cp_wk(&cp_list1, &wk_list1, cur_chr);

        get_covered(&wk_list1, &wk_list2);
    }
    /*
    if ( cp_list1 != NULL ) {
      while ( (ali=retrieve_first(&cp_list1))!=NULL ) {
        mafWrite(stdout, ali);
        mafAliFree(&ali);
      }
    }
    */
    mafWriteEnd(stdout);
    return 0;
}
