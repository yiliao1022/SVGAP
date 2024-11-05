// single_cov2.c version 10
// -- convert a pairwise .maf file to single-cov in both sequences
// -- suppose each row contain chrs for the same species, e.g. output from lav2maf

#include <stdlib.h>
#include "util.h"
#include "maf.h"
#include "mz_scores.h"
#include "multi_util.h"
#include "seq.h"

static const char rcsid[] = "$Id: single_cov2.c 142 2008-11-12 18:55:23Z rico $";


#define IGNORE 6	/* don't make alignments smaller than this */
#define STOP_CRITERIA 0.99
#define VERSION 11

struct mafAli *mafOverlap(struct mafAli *a, FILE *fp) {
    struct mafAli **A, *b, *alast;
    struct mafComp *c, *d;
    int nali, i, j, x, *min_later_beg, L, R, c_end_pos, d_end_pos, col_beg, col_end;
    double lost1, lost2;

    if (a==NULL)
        return NULL;
    A = mafArray(a, &nali);
    min_later_beg = ckalloc(nali*sizeof(int));

    min_later_beg[nali-1] = A[nali-1]->components->start;
    for (i = nali-2; i >= 0; --i)
        min_later_beg[i] =
            MIN(A[i]->components->start, min_later_beg[i+1]);

    for (i = 0; i < nali - 1; ++i) {
        for (j = i+1; j < nali; ++j) {
            if (A[i] == NULL)
                break;
            c = A[i]->components;
            c_end_pos = c->start + c->size - 1;
            if (min_later_beg[j] > c_end_pos)
                break;
            if (A[j] == NULL)
                continue;
            d = A[j]->components;
            d_end_pos = d->start + d->size - 1;
            if (d_end_pos < c->start)
                continue;
            if (c_end_pos < d->start)
                continue;
            // there is an overlap
            L = d->start  > c->start  ? d->start  : c->start;
            R = c_end_pos < d_end_pos ? c_end_pos : d_end_pos;

            if (c_end_pos >= d_end_pos) {
                if (fp != NULL)
                    fprintf(fp, "deleted %s:%d-%d\n",
                            c->src, d->start, d_end_pos);
                mafAliFree(&A[j]);
                continue;
            }

            /*
            fprintf(stderr, "overlap: L = %d, R = %d\n", L, R);
            mafWrite(stderr, A[i]);
            mafWrite(stderr, A[j]);
            */
            col_beg = mafPos2Col(c, L, A[i]->textSize);
            lost1 = mafScoreRange(A[i], col_beg, A[i]->textSize-col_beg);

            col_end = mafPos2Col(d, R, A[j]->textSize);
            lost2 = mafScoreRange(A[j], 0, col_end+1);

//fprintf(stderr, "%s:%d-%d lost1 = %g, lost2 = %g\n", c->src, L, R, lost1, lost2);

            if (fp != NULL)
                fprintf(fp, "deleting %s:%d-%d\n", c->src, L, R);
            if (lost1 <= lost2) {
                x = A[i]->textSize;
                b = mafSlice(A[i], 0, mafPos2Col(c, L, x));
                mafAliFree(&A[i]);
                A[i] = b;
            } else {
                x = A[j]->textSize;
                b = mafSlice(A[j], mafPos2Col(d, c_end_pos+1, x), A[j]->textSize);
                mafAliFree(&A[j]);
                A[j] = b;
            }
            if (A[i]->textSize < IGNORE)
                mafAliFree(&A[i]);
            if (A[j]->textSize < IGNORE)
                mafAliFree(&A[j]);
        }
    }

    // put the modified alignments back into a linked list
    a = alast = NULL;
    for (i = 0; i < nali; ++i)
        if ((b = A[i]) != NULL) {
            if (alast == NULL)
                a = b;
            else
                alast->next = b;
            alast = b;
        }

    alast->next = NULL;
    free(A);
    free(min_later_beg);
    return a;
}

void iterative_single_cov(struct mafAli **head, FILE* fp) {
    char ref_chr[500];
    struct mafAli *a, *b, *prev, *last, *cp_list, *wk_list, *res_last;

    cp_list = *head;
    *head = res_last = NULL;

    while ( cp_list != NULL) {
        strcpy(ref_chr, cp_list->components->src);
        prev = a = cp_list;
        last = wk_list = NULL;
        for (; a != NULL; ) {
            if ( strcmp(ref_chr, a->components->src)==0) { // add to wk_ls
                if (a == cp_list) {
                    cp_list = a->next;
                    a->next = NULL;
                    b = a;
                    prev = a = cp_list;
                } else {
                    prev->next = a->next;
                    a->next = NULL;
                    b=a;
                    a=prev->next;
                }
                if ( wk_list == NULL)
                    wk_list = last = b;
                else {
                    last->next = b;
                    last = b;
                }
            } else {
                prev = a;
                a=a->next;
            }
        }
        a = wk_list;
        a = mafOverlap(a, fp);
        if ( *head == NULL)
            *head = a;
        else
            res_last->next = a;
        if ( a!=NULL) {
            for (; a->next!=NULL; a=a->next)
                ;
            res_last = a;
        }
    }
    return;
}

int main(int argc, char **argv) {
    char cmd[500];
    FILE *fp = NULL;
    struct mafFile *mf;
    struct mafAli *a, **A;
    char* reference=NULL;
    int i, nali, orig1, orig2, res1, res2;

    sprintf(cmd, "single_cov2.v%d", VERSION);
    argv0 = cmd;

    if ( argc < 2)
        fatal("-- screening out overlapped regions.\nargs: pairwise.maf [R=species] [F=deleted.maf]\nBy default, single coverage is done for both species; if S=species specified, single coverage is done for the specified species only.\nThe first rows of all blocks must be of the same species; the second rows of all blocks must be of the same species.\n");

    mafWriteStart(stdout, "single_cov2");
    printf("# %s", cmd);
    for (i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    if (strncmp(argv[argc-1], "F=", 2) == 0) {
        fp = ckopen(argv[argc-1]+2, "w");
        argc--;
    }
    if (strncmp(argv[argc-1], "R=", 2) == 0) {
        reference = copy_string(argv[argc-1]+2);
        argc--;
    }

    init_scores70();

    mf = mafReadAll(argv[1], 1);
    a = mf->alignments;
    if ( a==NULL) {
        mafWriteEnd(stdout);
        return 0;
    }
    orig1 = orig2 = res1 = res2 = 0;
    for (; a!=NULL; a=a->next) {
        orig1 += a->components->size;
        orig2 += a->components->next->size;
    }
    a = mf->alignments;
    mf->alignments = NULL;
    mafFileFree(&mf);

    if ( a == NULL ) {
        if ( fp != NULL )
            fclose(fp);
        mafWriteEnd(stdout);
        return 0;
    }
    if ( reference == NULL || strcmp(a->components->name, reference)==0 )
        iterative_single_cov(&a, fp);
    flip_comps(a);
    if ( reference == NULL || strcmp(a->components->name, reference)==0 )
        iterative_single_cov(&a, fp);

    if (fp != NULL)
        fclose(fp);
    flip_comps(a);
    A = mafArray(a, &nali);
    for (i = 0; i < nali; ++i) {
        a = A[i];
        a = mafRowDashRm(a);
        if ( a == NULL || a->components->next == NULL )
            continue;
        a->score = mafScoreRange(a, 0, a->textSize);
        mafWrite(stdout, a);
        res1 += a->components->size;
        res2 += a->components->next->size;
    }
    for (i=0; i<nali; i++)
        A[i] = NULL;
    free(A);

    if ( (double)(res1+res2)/(orig1+orig2) < STOP_CRITERIA )
        fprintf(stderr, "%d bases loss out of %d\n", orig1+orig2-res1-res2, orig1+orig2);
    mafWriteEnd(stdout);
    return 0;
}
