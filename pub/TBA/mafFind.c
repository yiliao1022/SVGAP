/* mafFind -- find mafs intersecting a particular interval
*
* The basic command syntax is
*	mafFind file.maf beg end
* The mafs whose first row intersects positions beg-end are printed. For
* non-reference species, a command like
*	mafFind file.maf beg end mm3
* asks for mafs that have a row where"mm3" is a prefix of the "source"
* (e.g., "mm3.chr7") and which intersects positions beg-end. Finally, the last
* argument "slice" asks that ends of the reported mafs be trimmed to make them
* precisely match beg-end, as in:
*	mafFind file.maf beg end mm3 slice
*/

#define VERSION 1

#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "mz_scores.h"

int main(int argc, char **argv) {
    struct mafFile *mf;
    struct mafAli *a, *A;
    struct mafComp *c;
    int i, b, e, B, E, bcol, ecol, sz, slice = 0, species_len=-1;
    char cmd[100], *species = NULL;

    sprintf(cmd, "mafFind.v%d", VERSION);
    argv0 = cmd;

    if ((i = argc) > 4 && same_string(argv[i-1], "slice")) {
        slice = 1;
        --i;
    }
    if (i == 5) {
        species = argv[--i];
        species_len = strlen(species);
    }
    if (i != 4)
        fatal(" -- find mafs intersecting a particular interval.\nargs: file.maf beg end [species-prefix] [slice]");
    B = atoi(argv[2]);
    E = atoi(argv[3]);
    mf = mafOpen(argv[1], 0);
    mafWriteStart(stdout, cmd);
    printf("# %s", cmd);
    for (i = 1; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    init_scores70();
    while ((a = mafNext(mf)) != NULL) {
        c = a->components;
        if (species != NULL)
            while (c != NULL &&
                    strncmp(c->src, species, species_len))
                c = c->next;
        if (c != NULL) {
            b = c->start;
            e = b + c->size - 1;
            if (e >= B && b <= E) {
                if (slice) {
                    sz = a->textSize;
                    bcol = mafPos2Col(c, MAX(b,B), sz);
                    ecol = mafPos2Col(c, MIN(e,E), sz);
                    A = mafSlice(a, bcol, ecol+1);
                    if ( A==NULL)
                        continue;
                    A = mafRowDashRm(A);
                    if ( A==NULL)
                        continue;
                    A->score = mafScoreRange(A, 0, A->textSize);
                    mafWrite(stdout, A);
                    mafAliFree(&A);
                } else
                    mafWrite(stdout, a);
            }
        }
        mafAliFree(&a);
    }
    mafWriteEnd(stdout);
    return 0;
}
