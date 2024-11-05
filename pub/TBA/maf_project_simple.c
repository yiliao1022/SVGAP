/* maf_project version 10
* maf_project.c -- extract maf-file entries that name a given reference
*  sequence.
*
* Put reference sequence in the top row, order alignments by starting position
* in the reference sequence, and reorder rows as suggested by the phylogenetic
* tree.
*
*/
// compile with -DSTATS to get statistics about fusings of narrow blocks

#define VERSION 20
#define FUSE_SIZE 30	// try to fuse blocks with at most this many columns
#define DISCARD_SIZE 20	// can discard rows of blocks at most this wide
#define MUST_FUSE 10	// must fuse blocks with at most this many columns

static const char rcsid[] = "$Id: maf_project.c,v 1.12 2006-02-03 21:02:16 mhou Exp $";


#define CAUTIOUS	// do extra work to check for internal errors

static const char USAGE[] = "args: file.maf reference [from to]\
                            [filename-for-other-mafs] [species-guid-tree] [nohead]";

#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "maf_order.h"

struct name_list {
    char *name;
    struct name_list *next;
};

/* information at each node of a phylogenetic tree */
struct tree_node {
    int target, type;
    struct name_list *names;	/* leaf-species names */
};

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

    a=mf->alignments;
    mf->alignments=NULL;
    mafFileFree(&mf);
    for (; a!=NULL; a=next) {
        next = a->next;
        a->next = NULL;
        for (c = a->components; c != NULL; c = c->next) {
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
