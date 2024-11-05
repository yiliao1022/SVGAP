// maf_checkThread.c -- check threading condition of an alignment for a sequence
// first parameter is testing file, second par is name
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "util.h"
#include "maf.h"

int main(int argc, char **argv) {
    int lastEnd, plusStart, totalError;
    struct mafFile *mf;
    struct mafAli *b;
    struct mafComp *c;

    if ( argc < 2)
        fatal("args: maf-flie");

    totalError = 0;

    mf = mafOpen(argv[1], 0);
    lastEnd = -1;
    while ((b = mafNext(mf)) != NULL) {
        c = b->components;
        ////plusStart = getPlusStart(c->strand, c->start, c->size, c->srcSize);
        plusStart = c->start;
        if ( plusStart < lastEnd + 1) {
            printf("%s not threaded at %d ", c->src, c->start);
            totalError++;
        }
        lastEnd = plusStart + c->size - 1;
        mafAliFree(&b);
    }
    mafFileFree(&mf);
    printf("Total Errors: %d\n", totalError);

    return 0;
}
