/* maf_project version 12
* maf_project.c -- extract maf-file entries that name a given reference
*  sequence.
*
* Put reference sequence in the top row, order alignments by starting position
* in the reference sequence, and reorder rows as suggested by the phylogenetic
* tree.
*
* With compile-time options, human-readable output (i.e., when no file is
* named for the non-reference mafs) is adjusted to fuse each narrow alignment
* (i.e., having at most FUSE_SIZE columns) to one of its immediate neighbors,
* provided that it meets one of the following three conditions:
*
* 1. Every species with a row in the narrow block has a row in the neighbor
*    and the rows are adjacent in the species; rows of dashes are added to
*    the narrow block where it lacks a species in the neighbor.
* 2. The two blocks satisfy condition 1 after a single block of at most
*    DISCARD_SIZE columns is wedged between them; a block not involving the
*    reference sequence is examined to see if completely fills in the hole
*    between the narrow block and it neighbor.
* 3. The narrow block has at most DISCARD_SIZE columns, and the two blocks
*    can be made to satisfy condition 1 by removing <= 3 rows from the narrow
*    block. If the narrow block has at most MUST_FUSE columns, this is done
*    regardless of how many rows must be discarded.
*/

// compile with -DSTATS to get statistics about fusings of narrow blocks

#define VERSION 12
#define FUSE_SIZE 30	// try to fuse blocks with at most this many columns
#define DISCARD_SIZE 20	// can discard rows of blocks at most this wide
#define MUST_FUSE 10	// must fuse blocks with at most this many columns

static const char rcsid[] = "$Id: maf_project.c 142 2008-11-12 18:55:23Z rico $";


#define CAUTIOUS	// do extra work to check for internal errors

static const char USAGE[] = "args: file.maf reference [from to]\
                            [filename-for-other-mafs] [species-guid-tree] [nohead]";

#ifdef FUSE_SIZE
static int nfuse, nwedge, ndistroy, ncompress;
#endif

#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "maf_order.h"
#include "mz_scores.h"

struct name_list {
    char *name;
    struct name_list *next;
};

/* information at each node of a phylogenetic tree */
struct tree_node {
    int target, type;
    struct name_list *names;	/* leaf-species names */
};


static int abut(struct mafAli *a, struct mafAli *b) {
    struct mafComp *c, *d;

    for (c = a->components; c != NULL; c = c->next) {
        for (d = b->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src))
                break;
        if (d == NULL || d->paralog != c->paralog || c->strand != d->strand ||
                c->start + c->size != d->start)
            return 0;
    }
    for (c = b->components; c != NULL; c = c->next) {
        for (d = a->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src))
                break;
        if (d == NULL || d->paralog != c->paralog || c->strand != d->strand ||
                d->start + d->size != c->start)
            return 0;
    }
    return 1;
}

#ifdef FUSE_SIZE

// maybe fusing left a seam that can be closed
static void accordion(struct mafAli *a, int n1) {
    int i, sp, min_space, n = a->textSize;
    struct mafComp *c;

    for (min_space = n, c = a->components; c != NULL; c = c->next) {
        sp = 0;
        for (i = n1-1; i >= 0 && c->text[i] == '-'; --i)
            ++sp;
        for (i = n1; i < n && c->text[i] == '-'; ++i)
            ++sp;
        min_space = MIN(sp, min_space);
    }
    if (min_space > 0) {
//printf("------------- min_space = %d, n1 = %d, n = %d\n", min_space, n1, n);
//mafWrite(stdout, a);
        for (c = a->components; c != NULL; c = c->next) {
            for (i = n1; i > 0 && c->text[i-1] == '-'; --i)
                ;
            for ( ; i+min_space <= n; ++i)
                c->text[i] = c->text[i+min_space];
        }
        a->textSize -= min_space;
        ++ncompress;
//mafWrite(stdout, a);
    }
}

// fuse alignment b to the end of a
static void fuse(struct mafAli *a, struct mafAli *b) {
    struct mafComp *c, *d, *unmatched, *x;
    int  n1, n2, n, i;
    char *new;

    ++nfuse;
    n1 = a->textSize;
    n2 = b->textSize;
    a->textSize = n = n1 + n2;
    // extend each component of a
    for (c = a->components; c != NULL; c = c->next) {
        new = ckalloc((n+1)*sizeof(char));
        for (i = 0; i < n1; ++i)
            new[i] = c->text[i];
        new[n] = '\0';
        for (d = b->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src))
                break;
        if (d != NULL) { // same species occurs in b
#ifdef CAUTIOUS
            if (d->strand != c->strand || d->start != c->start + c->size)
                fatal("possible use of unprojected alignment");
#endif
            for (i = 0; i < n2; ++i)
                new[n1+i] = d->text[i];
            c->size += d->size;
        } else	// does not occur in b
            for (i = 0; i < n2; ++i)
                new[n1+i] = '-';
        //free(c->text);
        c->text = new;
    }
    // find components of b that are not in a, and add dashes on the left
    unmatched = NULL;
    for (d = b->components; d != NULL; d = d->next) {
        for (c = a->components; c != NULL; c = c->next)
            if (same_string(c->src, d->src))
                break;
        if (c == NULL) {
            x = mafCpyComp(d);
            x->text = new = ckalloc((n+1)*sizeof(char));
            for (i = 0; i < n1; ++i)
                new[i] = '-';
            for (i = 0; i < n2; ++i)
                new[i+n1] = d->text[i];
            new[n] = '\0';
            x->next = unmatched;
            unmatched = x;
        }
    }
    // add these new components to the end of a
    for (c = a->components; c->next != NULL; c = c->next)
        ;
    c->next = unmatched;

    accordion(a, n1);

    a->score = mafScoreRange(a, 0, a->textSize);
}

// fuse, can wedge or add all-dash rows to a, else return 0
static int fuseLeft(struct mafAli *a, struct mafAli *b, struct mafAli *orph) {
    struct mafComp *c, *d, *w, *xw;
    struct mafAli *wedge, *xwedge;
    int sep, compEnd, abut;

    // check each row of a for a row in b from the same species
    abut = 1;
    for (c = a->components; c != NULL; c = c->next) {
        for (d = b->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src))
                break;
        if (d == NULL || c->strand != d->strand)
            return 0;
        sep = d->start - c->start - c->size;
        if (sep < 0 || sep > DISCARD_SIZE)
            return 0;
        if (sep > 0)
            abut = 0;
    }
    // if wedging is not necessary ..
    if (abut) {
        fuse(a,b);
        return 1;
    }
    // find a row of a that does not abut the corresponding row of b
    for (c = a->components; c != NULL; c = c->next) {
        for (d = b->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src))
                break;
        if (d == NULL || c->start + c->size < d->start)
            break;
    }
    if (c == NULL || d == NULL)
        fatal("bad left wedge");
    for (wedge = orph; wedge != NULL; wedge = wedge->next) {
        for (w = wedge->components; w != NULL; w = w->next)
            if (same_string(w->src, c->src))
                break;
        if (w != NULL && w->strand == c->strand &&
                w->start == c->start + c->size)
            break;
    }
    if (wedge == NULL)
        return 0;
    // does wedge completely fill in the hole?
    for (c = a->components; c != NULL; c = c->next) {
        compEnd = c->start + c->size;
        for (w = wedge->components; w != NULL; w = w->next)
            if (same_string(c->src, w->src)) {
                if (c->strand != w->strand ||
                        compEnd != w->start)
                    return 0;
                compEnd += w->size;
                break;
            }
        for (d = b->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src)) {
                if (d->start != compEnd)
                    return 0;
                break;
            }
    }
    /*
    printf( "-----------------------wedging left \n");
    mafWrite(stdout, a);
    mafWrite(stdout, wedge);
    mafWrite(stdout, b);
    */

    // make a copy of wedge with only the rows from b
    xwedge = ckalloc(sizeof(struct mafAli));
    xwedge->components = NULL;
    xwedge->score = 0.0;
    xwedge->textSize = wedge->textSize;
    for (d = b->components; d != NULL; d = d->next) {
        for (w = wedge->components; w != NULL; w = w->next)
            if (same_string(d->src, w->src))
                break;
        if (w != NULL && w->strand == d->strand &&
                w->start + w->size == d->start) {
            xw = mafCpyComp(w);
            xw->text = (char*)malloc( (strlen(w->text)+1)*sizeof(char));
            strcpy(xw->text, w->text);
            xw->next = xwedge->components;
            xwedge->components = xw;
        }
    }
    fuse(xwedge, b);
    fuse(a, xwedge);
    /*
    mafWrite(stdout, a);
    */
    ++nwedge;
    return 1;
}

// fuse, can wedge or add all-dash rows to b, else return 0
static int fuseRight(struct mafAli *a, struct mafAli *b, struct mafAli *orph) {
    struct mafComp *c, *d, *w, *xw;
    struct mafAli *wedge, *xwedge;
    int sep, compStart, abut;

    // check each row of b for a row in a from the same species
    abut = 1;
    for (c = b->components; c != NULL; c = c->next) {
        for (d = a->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src))
                break;
        if (d == NULL || c->strand != d->strand)
            return 0;
        sep = c->start - d->start - d->size;
        if (sep < 0 || sep > DISCARD_SIZE)
            return 0;
        if (sep > 0)
            abut = 0;
    }
    // if wedging is not necessary ..
    if (abut) {
        fuse(a,b);
        return 1;
    }
    // find a row of b that does not abut the corresponding row of a
    for (c = b->components; c != NULL; c = c->next) {
        for (d = a->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src))
                break;
        if (d == NULL || c->start > d->start + d->size)
            break;
    }
    if (c == NULL || d == NULL)
        fatal("bad right wedge");
    for (wedge = orph; wedge != NULL; wedge = wedge->next) {
        for (w = wedge->components; w != NULL; w = w->next)
            if (same_string(w->src, c->src))
                break;
        if (w != NULL && w->strand == d->strand &&
                w->start == d->start + d->size)
            break;
    }
    if (wedge == NULL)
        return 0;
    // does wedge completely fill in the hole?
    for (c = b->components; c != NULL; c = c->next) {
        compStart = c->start;
        for (w = wedge->components; w != NULL; w = w->next)
            if (same_string(c->src, w->src)) {
                if (w->start + w->size != compStart)
                    return 0;
                compStart = w->start;
                break;
            }
        for (d = a->components; d != NULL; d = d->next)
            if (same_string(c->src, d->src)) {
                if (d->start + d->size != compStart)
                    return 0;
                break;
            }
    }
    /*
    printf( "-----------------------wedging right \n");
    mafWrite(stdout, a);
    mafWrite(stdout, wedge);
    mafWrite(stdout, b);
    */

    // make a copy of wedge with only the rows from a
    xwedge = ckalloc(sizeof(struct mafAli));
    xwedge->components = NULL;
    xwedge->score = 0.0;
    xwedge->textSize = wedge->textSize;
    for (c = a->components; c != NULL; c = c->next) {
        for (w = wedge->components; w != NULL; w = w->next)
            if (same_string(c->src, w->src))
                break;
        if (w != NULL && w->strand == c->strand &&
                w->start == c->start + c->size) {
            xw = mafCpyComp(w);
            xw->text = (char*)malloc( (strlen(w->text)+1)*sizeof(char));
            strcpy(xw->text, w->text);
            xw->next = xwedge->components;
            xwedge->components = xw;
        }
    }
    fuse(a, xwedge);
    fuse(a, b);
    /*
    mafWrite(stdout, a);
    */
    ++nwedge;
    return 1;
}

struct mafAli *beautify(struct mafAli *projection, struct mafAli *orphans) {
    struct mafAli *A, *B, *C;
    int fuseResult;

    A = NULL;
    B = projection;
    while (B != NULL) {
        C = B->next;
        if (B->textSize > FUSE_SIZE) {
            A = B;
            B = C;
        } else if (A != NULL && fuseRight(A, B, orphans)) {
            A->next = B = C;
        } else if (C != NULL && fuseLeft(B, C, orphans)) {
            B->next = C->next;
#ifdef DISCARD_SIZE
        } else if (B->textSize <= DISCARD_SIZE && A!=NULL) {
            struct mafComp *c, *d, *n;
            int i, j;
            // B can be fused to A if we delete i of B's rows

            for (i = 0, c = B->components; c != NULL; c = c->next) {
                for (d = A->components; d != NULL; d = d->next)
                    if (same_string(c->src, d->src))
                        break;
                if (d == NULL || d->strand != c->strand ||
                        d->start + d->size != c->start)
                    ++i;
            }
            // B can be fused to C if we delete j of B's rows
            if (C == NULL)
                break;
            for (j = 0, c = B->components; c != NULL;
                    c = c->next) {
                for (d = C->components; d != NULL; d = d->next)
                    if (same_string(c->src, d->src))
                        break;
                if (d == NULL || d->strand != c->strand ||
                        c->start + c->size != d->start)
                    ++j;
            }
            if (i <= j && (i <= 3 || B->textSize <= MUST_FUSE)) {
#ifdef DEBUG_DESTROY
                printf("fuse after tossing %d rows\n", i);
                mafWrite(stdout, A);
                mafWrite(stdout, B);
#endif
                // first row is always the target species
                for (c = B->components; (n = c->next) != NULL; ) {
                    for (d = A->components; d != NULL; d = d->next)
                        if (same_string(n->src, d->src))
                            break;
                    if (d == NULL || d->strand != n->strand ||
                            d->start + d->size != n->start)
                        c->next = n->next;
                    else
                        c = n;
                }
                mafColDashRm(B);
                fuseResult = fuseRight(A, B, orphans);
#ifdef DEBUG_DESTROY
                mafWrite(stdout, A);
#endif
                if ( fuseResult > 0 ) {
                    B = A->next = C; // deleted B
                    ++ndistroy;
                } else {
                    A=B;
                    B=C;
                }
            } else if (j <= 3 || B->textSize <= MUST_FUSE) {
#ifdef DEBUG_DESTROY
                printf("fuse before tossing %d rows\n", j);
                mafWrite(stdout, B);
                mafWrite(stdout, C);
#endif
                for (c = B->components; (n = c->next) != NULL; ) {
                    for (d = C->components; d != NULL; d = d->next)
                        if (same_string(n->src, d->src))
                            break;
                    if (d == NULL || d->strand != n->strand ||
                            n->start + n->size != d->start)
                        c->next = n->next;
                    else
                        c = n;
                }
                mafColDashRm(B);
                fuseResult = fuseLeft(B, C, orphans);
                if ( fuseResult > 0 ) {
                    B->next = C->next; // deleted C
                    ++ndistroy;
                } else {
                    A=B;
                    B=C;
                }
#ifdef DEBUG_DESTROY
                mafWrite(stdout, B);
#endif
            } else {
                A = B;
                B = C;
            }
#endif
        } else {
            A = B;
            B = C;
        }
    }

#ifdef DEBUG_DESTROY
    exit(0);
#endif

    return projection;
}
#endif

static struct name_list *get_names(char *target, char *tree_spec) {
    char *q, buf[500], *p;
    int top = -1;
    struct tree_node tree[1000];
    struct name_list *n;

    for (q = tree_spec; *q != '\0'; ++q) {
        if (*q == '(') {
            if (++top >= 1000)
                fatal("parse_tree: stack overflow");
            tree[top].type = '(';
        } else if (*q == ')') {
            if (top < 1 || tree[top].type != 0 ||
                    tree[top-1].type != '(') {
                q[1] = '\0';
                fatalf("parse error: %s", tree_spec);
            }
            tree[top-1] = tree[top];
            --top;
        } else if (isalpha(*q)) {
            /* leaf species */
            if (++top >= 1000)
                fatal("parse_tree: stack overflow");
            p = buf;
            while (isalpha(*q) || isdigit(*q) || *q == '_' ||
                    *q == '.')
                *p++ = *q++;
            --q;
            *p = '\0';
            n = ckalloc(sizeof(struct name_list));
            n->name = copy_string(buf);
            n->next = NULL;
            tree[top].target = same_string(buf, target);
            tree[top].type = 0;
            tree[top].names = n;
        } else if (*q != ' ')
            fatalf("improper character in tree specification: %c",
                   *q);
        if (top > 0 && tree[top-1].type == 0 && tree[top].type == 0) {
            /* merge the lists of names */
            if (tree[top-1].target && tree[top].target)
                fatal("both children have the target species");
            if (tree[top].target) {
                // names for subtree with target go first
                n = tree[top].names;
                if (n == NULL)
                    fatal("empty list of names");
                while (n->next != NULL)
                    n = n->next;
                n->next = tree[top-1].names;
                tree[top-1].names = tree[top].names;
            } else {
                n = tree[top-1].names;
                if (n == NULL)
                    fatal("empty list of names");
                while (n->next != NULL)
                    n = n->next;
                n->next = tree[top].names;
            }
            tree[top-1].target |= tree[top].target;
            --top;
        }
    }
    return tree[top].names;
}

int main(int argc, char **argv) {
    struct mafFile *mf;
    struct mafAli *a, *prev, *next, *A = NULL, *B, **align, *projection;
    int nohead=0;
#ifdef FUSE_SIZE
    struct mafAli *orphans = NULL;
#ifdef STATS
    int total_width = 0, nremain = 0, nofuse = 0;
#ifdef DISCARD_SIZE
    int nodiscard = 0;
#endif
#ifdef MUST_FUSE
    int nomust = 0;
#endif
#endif
#endif

    struct mafComp *c, *b, **location;
    int i, nali, beg, end, nspecies, orig_argc = argc;
    FILE *fp;
    char cmd[500], *target, **species, *x;
    struct name_list *n=NULL, *m=NULL;
    char  ref_chr[200];

    if ( strcmp(argv[argc-1], "nohead")==0) {
        nohead = 1;
        argc--;
    }

    sprintf(cmd, "maf_project.v%d", VERSION);
    argv0 = cmd;
    x = argv[argc-1];
    if (argc == 6 || (argc == 4 && strchr(argv[3], '('))) {
        n = get_names(argv[2], argv[argc-1]);
        for (m = n, nspecies = 0; m != NULL; m = m->next)
            ++nspecies;
        species = ckalloc(nspecies*sizeof(char *));
        for (i = 0, m = n; m != NULL; ++i, m = m->next) {
            //	species[i] = seq_ident(m->name);
            species[i] = (char*)malloc(500*sizeof(char));
            strcpy(species[i], m->name);
        }
        --argc;
    } else {
        n = NULL;
        nspecies = 0;
        location = NULL;
        species = NULL;
    }
    if (argc == 5 && (beg = atoi(argv[3])) >= 0 &&
            (end = atoi(argv[4])) > beg)
        argc = 3;
    else
        beg = end = -1;
    if (argc != 3 && argc != 4)
        fatalf(" -- extract maf-file entries that name a given reference sequence.\n%s", USAGE);
    //	target = seq_ident(argv[2]);
    target=copy_string(argv[2]);
    fp = (argc == 4 ? ckopen(argv[3], "w") : NULL);

    if ( nohead==0) {
        mafWriteStart(stdout, cmd);
        printf("# %s", argv0);
        for (i = 1; i < orig_argc; ++i)
            printf(" %s", argv[i]);
        putchar('\n');
    }

    mf = mafReadAll(argv[1], 1);
    /*
    for(a=mf->alignments; a!=NULL; a=a->next) 
      for(c=a->components; c!=NULL; c=c->next) { 
        for(chr=chrs, spe=spes; chr!=NULL; chr=chr->next, spe=spe->next)
          if ( strcmp(c->src, chr->str)==0)
    	break;
        if ( chr != NULL )
          c->name = copy_string(spe->str);
        else
          c->name = copy_string(c->src);
      }
    */

    a=mf->alignments;
    mf->alignments=NULL;
    mafFileFree(&mf);
    for (; a!=NULL; a=next) {
        next = a->next;
        a->next = NULL;
        for (c = a->components; c != NULL; c = c->next) {
            //parseSrcName(c->src, name, src);
            if (same_string(c->name, target) || same_string(c->src, target))
                break;
        }
        if (c != NULL) {
            if (c != a->components) {
                for (b = a->components;
                        b != NULL && b->next != c; b = b->next)
                    ;
                if (b == NULL)
                    fatal("maf_project: cannot happen");
                b->next = c->next;
                c->next = a->components;
                a->components = c;
            }
            if (c->strand == '-')
                rc(a);
            a->next = A;
            A = a;
        } else if (fp != NULL) {
            mafWrite(fp, a);
            mafAliFree(&a);
        }
#ifdef FUSE_SIZE
        else {
            a->next = orphans;
            orphans = a;
        }
#endif
    }

    init_scores70();

    while (A!=NULL) {
        strcpy(ref_chr, A->components->src);
        B=NULL;
        for (prev=a=A; a!=NULL;a=next ) {
            next=a->next;
            if ( strcmp(ref_chr, a->components->src)!=0) {
                prev->next = next;
                a->next = B;
                B=a;

            } else
                prev = a;
        }

        align = mafArray(A, &nali);
        if (nali == 0)
            fatal("no alignments in the projection");

        projection = align[0];
        for (i = 1; i < nali; ++i)
            align[i-1]->next = align[i];
        align[nali-1]->next = NULL;
        free(align);

        for (a = projection; (A = a->next) != NULL; a = A)
            if (abut(a,A)) {
                fuse(a,A);
                a->next = A->next;
                A = a;
            }
#ifdef FUSE_SIZE
        if (fp == NULL)
            projection = beautify(projection, orphans);

        for (a = projection; (A = a->next) != NULL; a = A)
            if (abut(a,A)) {
                fuse(a,A);
                a->next = A->next;
                A = a;
            }

#ifdef STATS

        fprintf(stderr,
                "%d blocks fused; %d involved wedging, %d involved removing rows, %d compressions\n",
                nfuse, nwedge, ndistroy, ncompress);
        for (a = projection; a != NULL; a = a->next) {
            ++nremain;
            i = a->textSize;
            total_width += i;
            if (i <= FUSE_SIZE)
                ++nofuse;
#ifdef DISCARD_SIZE
            if (i <= DISCARD_SIZE)
                ++nodiscard;
#endif
#ifdef MUST_FUSE
            if (i <= MUST_FUSE)
                ++nomust;
#endif
        }
        fprintf(stderr, "%d blocks remain (average width %4.1f):\n",
                nremain, (float)total_width/(float)nremain);
        fprintf(stderr, "  %d of width <= %d", nofuse, FUSE_SIZE);
#ifdef DISCARD_SIZE
        fprintf(stderr, ", %d of width <= %d", nodiscard, DISCARD_SIZE);
#endif
#ifdef MUST_FUSE
        fprintf(stderr, ", %d of width <= %d", nomust, MUST_FUSE);
#endif
        fprintf(stderr, "\n");
#endif
#endif

#ifdef NEVER
#ifdef CAUTIOUS
        for (a = projection; (A = a->next) != NULL; a = A) {
            c = a->components;
            if (c->start + c->size != A->components->start) {
                mafWrite(stderr, a);
                mafWrite(stderr, A);
                fatal("not adjacent");
            }
        }
#endif
#endif
        if ( n!=NULL)
            init_maf_order(nspecies, species);
        while ( projection != NULL ) {
            a = projection;
            projection = projection->next;
            a->next = NULL;
            c = a->components;
            if (beg < 0 || (c->start <= end && c->start + c->size > beg)) {
                if (n != NULL)
                    a = maf_order_ali(a);
                if ( a != NULL) {
                    mafWrite(stdout, a);
                }
            }
            if ( a != NULL )
                mafAliFree(&a);
        }

        A=B;
    }

    mafWriteEnd(stdout);
    return 0;
}
