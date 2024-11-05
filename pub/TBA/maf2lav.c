// maf2lav version 11
// maf2lav -- convert two rows of a .maf file to a lav file
#include "util.h"
#include "seq.h"
#include "maf.h"
#include "multi_util.h"
#include "mz_scores.h"

static const char rcsid[] = "$Id: maf2lav.c 142 2008-11-12 18:55:23Z rico $";


#define END_RUN(x) (x == '-' || x == '\0')
#define VERSION 11

struct pair {
    int textSize, score;
    struct mafComp *c1;
    struct mafComp *c2;
    struct pair *next;
}
*forward, *backward, *endfor, *endback;

void print_pair(struct pair *p) {
    int i, b1, b2, e1, e2, matches, gap;
    char *t1, *t2;

    for ( ; p != NULL; p = p->next) {
        b1 = p->c1->start + 1;
        b2 = p->c2->start + 1;
        printf("a {\n  s %d\n  b %d %d\n  e %d %d\n",
               p->score, b1, b2, b1 + p->c1->size-1, b2 + p->c2->size-1);
        t1 = p->c1->text;
        t2 = p->c2->text;
        gap = 1;
        e1 = b1 - 1;
        e2 = b2 - 1;
        matches = 0;
        for (i = 0; i <= p->textSize; ++i) {
            if (gap == 0 && (i == p->textSize ||
                             t1[i] == '-' || t2[i] == '-')) {
                printf("  l %d %d %d %d %d\n", b1, b2, e1, e2,
                       (100*matches)/(e1-b1+1));
                gap = 1;
            } else if (gap && t1[i] != '-' && t2[i] != '-') {
                b1 = e1+1;
                b2 = e2+1;
                matches = gap = 0;
            }
            if (i == p->textSize)
                break;
            if (t1[i] != '-')
                ++e1;
            if (t2[i] != '-')
                ++e2;
            if (gap == 0 && toupper(t1[i]) == toupper(t2[i]))
                ++matches;
        }
        printf("}\n");
    }
}

void make_lav(char *name1, int len1, char *head1,
              char *name2, int len2, char *head2, int n) {

    printf("#:lav\ns {\n");
    printf("  \"%s\" 1 %d 0 1\n  \"%s\" 1 %d 0 %d\n}\n",
           name1, len1, name2, len2, n);
    printf("h {\n  \"%s\"\n  \"%s\"\n}\n", head1, head2);
    print_pair(forward);
    printf( "x {\n  n 0\n}\n#:lav\ns {\n");
    printf("  \"%s\" 1 %d 0 1\n  \"%s-\" 1 %d 1 %d\n}\n",
           name1, len1, name2, len2, n);
    printf("h {\n  \"%s\"\n  \"%s (reverse complement)\"\n}\n",
           head1, head2);
    print_pair(backward);
    printf("x {\n  n 0\n}\n");
}

void disconnect(struct mafAli *a, struct mafComp *c) {
    struct mafComp *b;

    if (a->components == c)
        a->components = c->next;
    else {
        for (b = a->components; b != NULL && b->next != c; b = b->next)
            ;
        if (b == NULL)
            fatal("failure to disconnect");
        b->next = c->next;
    }
}

void record (struct mafAli *a, struct mafComp *c1, struct mafComp *c2) {
    struct pair *new = ckalloc(sizeof(struct pair));
    struct mafAli A;
    char *s = c1->text, *t = c2->text;
    int i, j;

    disconnect(a, c1);
    disconnect(a, c2);
    for (i = j = 0; i < a->textSize; ++i)
        if (s[i] != '-' || t[i] != '-') {
            if (j < i) {
                s[j] = s[i];
                t[j] = t[i];
            }
            ++j;
        }
    s[j] = t[j] = '\0';
    new->textSize = j;
    new->c1 = c1;
    new->c2 = c2;
    new->next = NULL;

    A.components = c1;
    c1->next = c2;
    c2->next = NULL;
    A.textSize = j;
    A.score = 0.0;
    A.next = NULL;
    new->score = mafScoreRange(&A, 0, A.textSize);
    if (c2->strand == '+') {
        if (endfor == NULL)
            forward = new;
        else
            endfor->next = new;
        endfor = new;
    } else {
        if (endback == NULL)
            backward = new;
        else
            endback->next = new;
        endback = new;
    }
}

int main(int argc, char **argv) {
    SEQ *sf1, *sf2;
    char name1[500], name2[500], *head1, *head2, species1[200], species2[200], strand1, strand2, cmd[200], chr1[200], chr2[200];
    struct mafFile *mf;
    struct mafAli *a;
    struct mafComp *c, *c1, *c2;
    int offset1, len1, offset2, len2, end1, end2, srcSize1, srcSize2, n, c2s, c_len;


    sprintf(cmd, "maf2lav.v%d", VERSION);
    argv0 = cmd;

    if (argc != 4)
        fatal("args = align.maf seq1 seq2");

    init_scores70();

    mf = mafReadAll(argv[1], 0);

    printf("#:lav\nd {\n  \"mav2lav %s %s %s\"\n}\n",
           argv[1], argv[2], argv[3]);

    sf1 = seq_open(argv[2]);

    while (seq_read(sf1)) {
        if (parseHeader(argv[2], sf1, species1, chr1, &offset1, &end1, &strand1, &srcSize1)==0) {
            strcpy(name1, species1);
            if ( strcmp(species1, chr1) != 0 ) {
                strcat(name1, ".");
                strcat(name1, chr1);
            }
        } else {
            strcpy(name1, argv[2]);
            offset1 = 0;
        }
        len1 = SEQ_LEN(sf1);
        head1 = copy_string(SEQ_HEAD(sf1));

        n = 0;
        sf2 = seq_open(argv[3]);
        while (seq_read(sf2)) {
            if (parseHeader(argv[3], sf2, species2, chr2, &offset2, &end2, &strand2, &srcSize2)==0) {
                strcpy(name2, species2);
                if ( strcmp(species2, chr2) != 0) {
                    strcat(name2, ".");
                    strcat(name2, chr2);
                }
            } else {
                strcpy(name2, argv[3]);
                offset2 = 0;
            }
            ++n;
            len2 = SEQ_LEN(sf2);
            head2 = copy_string(SEQ_HEAD(sf2));
            forward = backward = endfor = endback = NULL;
            for (a = mf->alignments; a != NULL; a = a->next) {
                c1 = c2 = NULL;
                for (c = a->components; c != NULL; c = c->next) {
                    c_len = strlen(c->src);
                    if (c1==NULL && strncmp(c->src, name1, c_len)==0 && ((c->strand=='+' && c->start+1 >= offset1 && c->start+c->size-1 < end1) || (c->start=='-' && c->srcSize-c->start-c->size+1 >= offset1 && c->srcSize - c->start < end1)))
                        c1 = c;
                    else if (strncmp(c->src, name2, c_len)==0 && ((c->strand=='+' && c->start+1 >= offset2 && c->start+c->size-1 < end2) || (c->strand=='-' && c->srcSize - c->start - c->size+1 >= offset2 && c->srcSize - c->start < end2)))
                        c2 = c;
                }
                if (c1 != NULL && c2 != NULL) {
                    if (c1->strand == '-') {
                        c1->start = c1->srcSize - (c1->start + c1->size);
                        c1->strand = '+';
                        do_revcompl(c1->text, a->textSize);
                        c2->start = c2->srcSize - (c2->start + c2->size);
                        c2->strand = (c2->strand == '-' ? '+' : '-');
                        do_revcompl(c2->text, a->textSize);
                    }
                    if (c2->strand == '+')
                        c2s = c2->start - offset2 +1;
                    else {
                        c2s = c2->srcSize -
                              (c2->start + c2->size);
                        c2s -= (offset2-1);
                        c2s = len2 - (c2s + c2->size);
                    }
                    //fprintf(stderr, "testing start = %d\n", c2->start);
                    if (c2s >= 0 && c2s < len2) {
                        c1->start -= (offset1-1);
                        c2->start = c2s;
                        record(a, c1, c2);
                    }
                }
            }
            make_lav(argv[2], len1, head1, argv[3], len2, head2, n);
        }
        seq_close(sf2);
    }
    printf("m {\n  n 0\n}\n#:eof\n");

    return 0;
}
