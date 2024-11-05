// pair2tb version 11
// pair2tb.c -- convert a pairwise .maf file to a threaded blockset
// The given pairwise alignment must be single-cov in both sequences.
// the next two arguments are specie names, followed by a string or not:
// "contigs+" or "contigs-", no mapping files are needed

#define VERSION 11

#include "util.h"
#include "maf.h"
#include "mz_scores.h"
#include "multi_util.h"
#include "maftop2tb.h"
#include "seq.h"

static const char rcsid[] = "$Id: pair2tb.c 142 2008-11-12 18:55:23Z rico $";


int main(int argc, char **argv) {
    struct mafFile *mf;
    struct mafAli *a;
    char cmd[500];

    sprintf(cmd, "pair2tb.v%d", VERSION);
    argv0 = cmd;
    if (argc != 4)
        fatal("-- convert a pairwise .maf file to a threaded blockset.\nargs: pairwise.maf seqfile1 seqfile2");

    mafWriteStart(stdout, cmd);

    mf = mafReadAll(argv[1], 1);
    for (a = mf->alignments ; a != NULL; a = a->next)
        mafWrite(stdout, a);

    //mf->alignments = maf_sort_top(mf->alignments);
    mf->alignments = getMafBetween(mf, argv[2], stdout);

    flip_comps(mf->alignments);
    //mf->alignments = maf_sort_top(mf->alignments);
    mf->alignments = getMafBetween(mf, argv[3], stdout);

    mafFileFree(&mf);
    mafWriteEnd(stdout);
    return 0;
}
