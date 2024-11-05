/* maf2fasta.c.c -- convert an ordered reference-sequence maf file to a few
*    long rows of nucleotides and dashes.  The default output format is that
*    used for text alignments by MultiPipMaker; appending the command-line
*    argument "fasta" or "fasta2" requests FastA-format output; with "fasta2",
*    alignment rows are split into rows of length <= 50 (set by COL_WIDTH).
*    To identify gaps as either "within-alignment" or "between-alignment",
*    append a character to the word "fasta" or "fasta2" that will replace '-'
*    between two local alignments, as in "fasta@".
*
*    Some maf files contain rows from a number of sources that should all be
*    placed in the same row. For example, all sources "mm3.chr1", "mm3.chr2",
*    etc. might be from mouse, and "rn3.chr1", "rn3.chr2" etc. from rat.
*    Unless otherwise specified, maf2fasta.c will put sequence from each of
*    these sources in its own row. If the last command-line argument (after
*    "fasta" or "fasta2") is "alias=filename", then the file named "filename"
*    will be read. It should contains lines like "mm3 mouse" or "rn3 rat".
*    The result will be that source names starting with the letters "mm3"
*    will be replaced by "mouse" (and hence be put in a single row named
*    "mouse"), and names starting with "rn3" will be replaced by "rat".
*
*    It is also possible to specify begin and end positions (starts at 0) in
*    the reference sequence, to be spanned by the generated alignment. If you
*    want alignments to start at a particular point and extend to the end of
*    the sequence, you can use a number known to exceed the sequence length
*    (e.g., for a human chromosome use 300000000 -- three hundred million) for
*    the second position. Thus, typical uses include:
*	maf2fasta.c human human.maf
*	maf2fasta.c human human.maf fasta#
*	maf2fasta.c human human.maf 0 1000 fasta# alias=row_names
*    Here, "human" names a file containing the reference sequence (used to
*    fill in the first row of the alignment between entries of the .maf file)
*    and "human.maf" names a maf-format file of local alignments where
*    (1) each alignment starts with a row from "human" and (2) the alignments
*    are sorted by start position in "human". The file named "row_names"
*    should contain lines like "mm3 -> mouse".
*/

#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "seq.h"

#include <limits.h>

#define VERSION 3

#define COL_WIDTH 50

struct edge {
    struct name_rec *name;
    struct edge *next;
};

#define WHITE 0
#define GRAY 1
#define BLACK 2
struct name_rec {
    char *name;
    int color;
    struct edge *follows;
    struct name_rec *next;
};

char **row_name;
int next_row_name;

// depth-first search
static void dfs(struct name_rec *n) {
    struct edge *e;

    if (n==NULL)
        return;
    if (n->color != WHITE)
        fatal("please apply the maf_order program");
    n->color = GRAY;
    for (e = n->follows; e != NULL; e = e->next)
        if (e->name->color != BLACK)
            dfs(e->name);
    if (next_row_name <= 0)
        fatal("underflow in row names");
    row_name[--next_row_name] = copy_string(n->name);
    n->color = BLACK;
}

// compare reference row of a block with the reference sequence
static void check_ref(struct mafComp *c, uchar *seq, int seq_len, int ncol, int start, int iupac2n) {
    char *s, x;
    int pos, col;

    pos = c->start - (start-1);
    s = c->text;

    for (col = 0; col < ncol; ++col) {
        if ((x = s[col]) != '-') {
            if (pos >= seq_len)
                fatalf("maf position %d >= fa size %d", pos, seq_len);

            if (iupac2n && (strchr("ACGTNacgtn", seq[pos]) == NULL))
                seq[pos] = (seq[pos] == toupper(seq[pos])) ? 'N' : 'n';

            if (toupper(x) != toupper(seq[pos]))
                fatalf("ref-seq mismatch at position %d", pos);
            ++pos;
        }
    }
}

int main(int argc, char **argv) {
    struct mafFile *mf;
    struct mafAli *a, *A, *last_a;
    struct mafComp *c, *d;
    SEQ *sf;
    uchar *s;
    char **row, star = '-', cmd[100], species[200], chrName[200], strand, species_name[200], chr_name[200];
    char *first_comp_src = NULL, *ref_src = NULL;
    struct name_rec *names, *m, *n;
    struct edge *e;
    int i, j, seq_len, nrow, ncol, next_pos, col, fasta = 0, BEG = 0,
            END = INT_MAX, start, tail, srcSize, beg, end, copyback = 0,
                    iupac2n = 0;

    sprintf(cmd, "maf2fasta.v%d", VERSION);
    argv0 = cmd;

    for (i = 1; i < argc; i++) {
        if (strncmp(argv[i], "fasta2", 6) == 0) {
            if (strlen(argv[i]) > 6)
                star = argv[i][6];
            fasta = 2;
            copyback++;
            continue;
        }

        if (strncmp(argv[i], "fasta", 5) == 0) {
            if (strlen(argv[i]) > 5)
                star = argv[i][5];
            fasta = 1;
            copyback++;
            continue;
        }

        if (strncmp(argv[i], "iupac2n", 7) == 0) {
            iupac2n = 1;
            copyback++;
            continue;
        }

        if (strncmp(argv[i], "refsrc=", 7) == 0) {
            ref_src = copy_string(argv[i] + 7);
            if (strlen(argv[i]) == 7)
                fatal("refsrc argument missing");
            copyback++;
            continue;
        }

        argv[i - copyback] = argv[i];
    }
    argc -= copyback;

    if (argc != 3 && argc != 5)
        fatal("args = refseq-file maf-file [beg end] [fasta[2]][?] [iupac2n] [refsrc=src]");

    if (argc == 5) {
        BEG = MAX(0,atoi(argv[3]));
        END = atoi(argv[4]);
        if (BEG > END)
            fatalf("BEG = %d > END = %d", BEG, END);
    }

    sf = seq_open(argv[1]);

    while ( seq_read(sf)) {
        beg = BEG;
        end = END;
        s = SEQ_CHARS(sf);
        seq_len = SEQ_LEN(sf);
        if (parseHeader(argv[1], sf, species, chrName, &start, &tail, &strand, &srcSize) != 0) {
            start = 1;
            tail = srcSize = seq_len;
            strand = '+';
        }
        if ( beg < start - 1)
            beg = start - 1;
        if ( end > tail - 1)
            end = tail - 1;


        // extract the relevant maf entries, chopping at beg and end;
        // also chop overlaps
        mf = mafOpen(argv[2], 0);
        A = last_a = NULL;
        next_pos = beg;
        if (ref_src != NULL)
            first_comp_src = copy_string(ref_src);
        while ((a = mafNext(mf)) != NULL) {
            if ((c = a->components) == NULL)
                fatal("empty maf entry");
            if (ref_src == NULL && first_comp_src == NULL)
                first_comp_src = copy_string(c->src);
            if (!same_string(c->src, first_comp_src))
                continue;
            if ( (c->strand=='+' && (c->start+1 > end || c->start+c->size < start )) || (c->strand=='-' && (c->srcSize-(c->start+c->size-1) > end || c->srcSize - c->start < start))) {
                mafAliFree(&a);
                continue;
            }
            if (c->start + c->size <= next_pos) {
                mafAliFree(&a);
            } else if (c->start > end)
                break;
            else {	// intersects interval next_pos..end
                if (c->start < next_pos) {
                    j = mafPos2Col(c, next_pos, a->textSize);
                    a = mafSlice(a, j, a->textSize);
                    c = a->components;
                }
                next_pos = c->start + c->size;
                if (c->start + c->size > end+1) {
                    j = mafPos2Col(c, end, a->textSize);
                    a = mafSlice(a, 0, j+1);
                }
                //do_alias(a, alias);
                a->next = NULL;
                if (last_a == NULL)
                    A = a;
                else
                    last_a->next = a;
                last_a = a;
            }
        }
        ckfree(first_comp_src);
        first_comp_src = NULL;
        if ( A==NULL)
            continue;

        mafFileFree(&mf);
        // find and properly order the species names
        nrow = ncol = 0;
        next_pos = beg;
        names = NULL;
        for (a = A; a != NULL; a = a->next) {
            c = a->components;
            if (c->start < next_pos)
                fatalf("alignments out of order at pos %d", c->start);
            check_ref(c, s, seq_len, a->textSize, start, iupac2n);
            ncol += (c->start - next_pos + a->textSize);
            next_pos = c->start + c->size;
            if (nrow == 0) {
                names = ckalloc(sizeof(struct name_rec));
                parseSrcName(c->src, species_name, chr_name);
                names->name = copy_string(species_name);
                names->follows = NULL;
                names->next = NULL;
                nrow = 1;
            } else {
                parseSrcName(c->src, species_name, chr_name);
                if (!same_string(species_name, names->name))
                    fatalf("conflicting ref-seq names: %s and %s",
                           names->name, species_name);
            }
            for (m = names ; (d = c->next) != NULL; c = d, m = n) {
                parseSrcName(d->src, species_name, chr_name);
                for (n = names; n != NULL; n = n->next)
                    if (same_string(n->name, species_name))
                        break;
                if (n == NULL) {
                    n = ckalloc(sizeof(struct name_rec));
                    parseSrcName(d->src, species_name, chr_name);
                    n->name = copy_string(species_name);
                    n->next = m->next;
                    n->follows = NULL;
                    m->next = n;
                    ++nrow;
                }
                for (e = m->follows; e != NULL; e = e->next)
                    if (same_string(e->name->name, n->name))
                        break;
                if (e == NULL) {
                    e = ckalloc(sizeof(struct edge));
                    e->name = n;
                    e->next = m->follows;
                    m->follows = e;
                }
            }
        }
        ncol += (end - next_pos + 1);

        for (n = names; n != NULL; n = n->next)
            n->color = WHITE;
        next_row_name = nrow;
        row_name = ckalloc(nrow*sizeof(char *));
        dfs(names);
        if (next_row_name != 0)
            fatal("not enough row names");

        // allocate the rectangular matrix for the alignment
        row = ckalloc(nrow*sizeof(char *));
        for (i = 0; i < nrow; ++i)
            row[i] = ckalloc((ncol+1)*sizeof(char));
        col = 0;
        next_pos = beg;
        // fill in the matrix
        for (a = A; a != NULL; a = a->next) {
            c = a->components;
            // put in the gap before the start of the local alignment
            for (j = next_pos; j < c->start; ++j, ++col) {
                row[0][col] = s[j-start+1];
                for (i = 1; i < nrow; ++i)
                    row[i][col] = star;
            }
            // put in the local alignment
            for (i = 0; i < nrow; ++i) {
                for (c=a->components; c!=NULL; c=c->next) {
                    parseSrcName(c->src, species_name, chr_name);
                    if ( same_string(species_name, row_name[i])) {
                        for (j = 0; j < a->textSize; ++j)
                            row[i][col+j] = c->text[j];
                        break;
                    }
                }
                if ( c==NULL)
                    for (j = 0; j < a->textSize; ++j)
                        row[i][col+j] = star;
            }

            c = a->components;
            next_pos = c->start + c->size;
            col += a->textSize;
        }
        // put in the gap after the last local alignment
        for (j = 0; j < ncol - col; ++j)
            row[0][col+j] = s[next_pos+j-start+1];
        for (i = 1; i < nrow; ++i)
            for (j = col; j < ncol; ++j)
                row[i][j] = star;

        // print the alignment
        if (beg != 0 || end != seq_len - 1) {
            char new_name[1000];
            sprintf(new_name, "%s:%d-%d", row_name[0], beg, end);
            row_name[0] = copy_string(new_name);
        }
        if (fasta == 1) {
            for (i = 0; i < nrow; ++i) {
                row[i][ncol] = '\0';
                printf(">%s\n%s\n", row_name[i], row[i]);
            }
        } else if (fasta == 2) {
            for (i = 0; i < nrow; ++i) {
                printf(">%s\n", row_name[i]);
                for (col = j = 0; j < ncol; ++j) {
                    putchar(row[i][j]);
                    if (++col == COL_WIDTH) {
                        putchar('\n');
                        col = 0;
                    }
                }
                if (col != 0)
                    putchar('\n');
            }
        } else {
            printf("%d %d\n", nrow, ncol);
            for (i = 0; i < nrow; ++i)
                printf("%s\n", row_name[i]);
            for (i = 0; i < nrow; ++i) {
                row[i][ncol] = '\0';
                printf("%s\n", row[i]);
            }
        }
    }

    return 0;
}
