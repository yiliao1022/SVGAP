/* seq.c -- procedures for reading a file of DNA sequences */

static const char rcsid[] =
    "$Id: seq.c 142 2008-11-12 18:55:23Z rico $";

#include "util.h"
#include "seq.h"
#include "nib.h"
#include "charvec.h"
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>

#define _ (-1)

int KEEP_SEQ = 0;

extern const signed char fasta_encoding[256];
const signed char fasta_encoding[] = {
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, 0, _, 1, _, _, _, 2, _, _, _, _, _, _, _, _,
                                         _, _, _, _, 3, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                     };

extern const signed char nfasta_encoding[256];
const signed char nfasta_encoding[] = {
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, 0, _, 1, _, _, _, 2, _, _, _, _, _, _, _, _,
                                          _, _, _, _, 3, _, _, _, _, _, _, _, _, _, _, _,
                                          _, 0, _, 1, _, _, _, 2, _, _, _, _, _, _, _, _,
                                          _, _, _, _, 3, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                          _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                                      };

static const unsigned char nfasta_ctype[256];
static const unsigned char nfasta_ctype[] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 2, 2, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
            0, 1, 3, 1, 3, 0, 0, 1, 3, 0, 0, 3, 0, 3, 1, 0,
            0, 0, 3, 3, 1, 0, 3, 3, 1, 3, 0, 0, 0, 0, 0, 0,
            0, 1, 3, 1, 3, 0, 0, 1, 3, 0, 0, 3, 0, 3, 1, 0,
            0, 0, 3, 3, 1, 0, 3, 3, 1, 3, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };

static const unsigned char dna_complement[] =
    "                                                                "
    " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
    "                                                                "
    "                                                                ";
/* ................................................................ */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */
/* ................................................................ */

static const unsigned char fasta_decoding[] = "ACGT";

static void Fatalfr(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fflush(stdout);
    print_argv0();
    (void)vfprintf(stderr, fmt, ap);
    (void)fprintf(stderr, ": %s\n", strerror(errno));
    va_end(ap);
    exit(1);
}

void init_KEEP_SEQ() {
    KEEP_SEQ = 1;
}

void close_KEEP_SEQ() {
    KEEP_SEQ = 0;
}


static SEQ *seq_mask_inplace(SEQ *seq);
static int getpair(FILE *fp, int *a, int *b);
static char *byte_fill_range(uchar *p, int l, int c, int a, int b);
static int ws(int c);
static int getnwc(FILE *fp);
static void un_getc(int c, FILE *fp);
static void char_append(charvec_t *s, int c);

enum {
    Nfasta_bad=0, Nfasta_nt=1, Nfasta_ws=2, Nfasta_amb=3
};

static int ws(int c) {
    return (c == ' ' || c == '\t' || c == '\n' || c == '\f' || c == '\r');
}

static int getnwc(FILE *fp) {
    int c = EOF;
    if (!feof(fp)) do c = getc(fp);
        while (c != EOF && ws(c));
    return c;
}

static void un_getc(int c, FILE *fp) {
    if (c != EOF)
        if (ungetc(c, fp) == EOF)
            fatalf("cannot ungetc '%c'", c);
}

static void char_append(charvec_t *s, int c) {
    if (!charvec_append(s, (char)c))
        fatal("cannot append");
}

int seq_read_fasta(SEQ *seq) {
    int b, c;
    charvec_t shdr = charvec_INIT(ckrealloc,ckfree);
    charvec_t sseq = charvec_INIT(ckrealloc,ckfree);

    if (feof(seq->fp) || ferror(seq->fp))
        return 0;

    if (seq->count > 0) {
        if (seq->flags & SEQ_IS_SUBRANGE) {
            return 0;
        } else {
            seq->from = 1;
            seq->slen = -1; /* will be computed below */
        }
    }

    if (KEEP_SEQ==0 && seq->header) ZFREE(seq->header);
    if (KEEP_SEQ==0 && seq->seq) ZFREE(seq->seq);

    seq->offset = ftell(seq->fp);

    /* --- header --- */
    c = getnwc(seq->fp);
    if (c == '>') {
        //while (c != '\n' && c != '\r' && c != EOF) {
        while (c != '\n' && c != EOF) {
            char_append(&shdr, c);
            c = getc(seq->fp);
        }
    } else {
        un_getc(c, seq->fp);
    }
    if (ferror(seq->fp))
        Fatalfr("seq_read(%s)", seq->fname);
    char_append(&shdr, 0);
    seq->header = shdr.a;
    seq->hlen = shdr.len;

    /* --- seq --- */
    b = '\n';
    c = getnwc(seq->fp);
    //while ((c != EOF) && !((b == '\n' || b == '\r') && c == '>')) {
    while ((c != EOF) && !(b == '\n' && c == '>')) {
        switch (nfasta_ctype[c]) {
        case Nfasta_nt:
            char_append(&sseq, c);
            break;
        case Nfasta_ws:
            /* skip space */
            break;
        case Nfasta_amb:
            if (seq->flags & SEQ_ALLOW_AMB) {
                char_append(&sseq, c);
                break;
            }
            /* FALLTHRU */
        default:
            fatalf("non-DNA character '%c' in sequence '%s'",
                   c, seq->fname);
            break;
        }
        b = c;
        c = getc(seq->fp);
    }
    un_getc(c, seq->fp);
    if (ferror(seq->fp))
        Fatalfr("seq_read(%s)", seq->fname);

    /* check conformance */
    if (SEQ_LEN(seq) == -1) {
        char_append(&sseq, 0);
        charvec_fit(&sseq);
        seq->seq = (uchar*)sseq.a;
        seq->slen = sseq.len;
        if (seq->slen > 0) --seq->slen;  /* don't include '\0' */
    } else {
        charvec_t ssub = charvec_INIT(ckrealloc,ckfree);
        int i;

        if (SEQ_FROM(seq) < 1 ||
                (int)sseq.len < SEQ_FROM(seq) ||
                SEQ_TO(seq) < 1 ||
                (int)sseq.len < SEQ_TO(seq) ||
                SEQ_TO(seq) < SEQ_FROM(seq))
            fatalf("range [%d,%d] incommensurate with sequence [%d,%d]",
                   SEQ_FROM(seq), SEQ_TO(seq), 1, sseq.len);

        for (i = SEQ_FROM(seq); i <= SEQ_TO(seq); ++i)
            char_append(&ssub, sseq.a[i-1]);
        char_append(&ssub, 0);
        charvec_fini(&sseq);
        seq->seq = (uchar*)ssub.a;
    }

    seq->flags = seq->flags &~ SEQ_IS_REVCOMP;
    if (seq->flags & SEQ_DO_REVCOMP) {
        (void)seq_revcomp_inplace(seq);
    }
    if (seq->flags & SEQ_HAS_MASK) {
        (void)seq_mask_inplace(seq);
    }

    seq->count++;
    return 1;
}

static SEQ *seq_mask_inplace(SEQ *seq) {
    int a, b;
    FILE *fp = fopen(seq->maskname, "r");

    if (fp == 0) {
        Fatalfr("cannot open '%s'", seq->maskname);
        return 0;
    } else {
        while (getpair(fp, &a, &b))
            byte_fill_range(SEQ_CHARS(seq),SEQ_LEN(seq),'X',a,b);
        fclose(fp);
        seq->flags |= SEQ_HAS_MASK;
        return seq;
    }
}

#define BUF 128
static int getpair(FILE *fp, int *a, int *b) {
    char buf[BUF];

    /* XXX - should handle comments, etc */
    if (fgets(buf, (int)sizeof(buf), fp) == 0)
        return 0;
    if (sscanf(buf, "%d%d", a, b) != 2)
        return 0;
    return 1;
}

static char *byte_fill_range(uchar *p, int l, int c, int a, int b) {
    /* fill [a,b] (1-indexed) in p with c */

    a--;
    b = b-a; /* make it into a 0-indexed, open interval */
    return (b < 0 || l < b) ? 0 : memset(p+a, c, (size_t)b);
}

static int parse_fname(const char* arg,
                       char **fname, int *from, int *len, char **maskfile) {
    char *p = 0;
    int flags = 0;

    /* "seqfile{maskfile}[from,to]-" */
    *fname = copy_string(arg);

    p = (*fname)+strlen(*fname)-1;
    if (*p == '-') {
        *p = 0;
        flags |= SEQ_DO_REVCOMP;
    }

    if ((p = strchr(*fname, '['))) {
        int to;

        if (sscanf(p+1, "%d,%d", from, &to) != 2)
            return -1;
        if (*from <= 0 || *from > to)
            return -1;
        *p = '\0';
        *len = to - *from + 1;
        flags |= (SEQ_DO_SUBRANGE|SEQ_IS_SUBRANGE);
    } else {
        *from = 1;
        *len = -1;
    }

    if ((p = strchr(*fname, '{'))) {
        char *q = strchr(p+1, '}');
        if (q) {
            *p = *q = 0;
            if (maskfile) {
                *maskfile = copy_string(p+1);
                flags |= SEQ_DO_MASK;
            }
        }
    } else {
        *maskfile = copy_string(""); /* XXX ugh */
    }
    return flags;
}

static int check_flags(int flags) {
    switch (flags & (SEQ_DISALLOW_AMB|SEQ_ALLOW_AMB)) {
    case 0:
        /* default is to allow ambiguious */
        flags |= SEQ_ALLOW_AMB;
        break;
    case SEQ_ALLOW_AMB:
    case SEQ_DISALLOW_AMB:
        break;
    case SEQ_DISALLOW_AMB|SEQ_ALLOW_AMB:
        fatalf("seq_open: contradictory flags: SEQ_DISALLOW_AMB|SEQ_ALLOW_AMB");
    }
    return flags;
}


SEQ* seq_open_type(const char *fname, int type) {
    SEQ *s = ckallocz(sizeof(SEQ));
    int r, flags = 0;

    r = parse_fname(fname,
                    &(s->fname), &(s->from), &(s->slen), &(s->maskname));
    if (r == -1)
        fatalf("improper positions specification: %s", fname);

    s->type = type;
    s->flags = check_flags(r|flags);
    s->fp = ckopen(s->fname, "rb");
    s->count = 0;
    s->offset = 0;
    s->contigs = NULL;
    s->numContigs = -1;
    s->flipped = 'y';
    return s;
}

SEQ* seq_open(const char *fname) {
    return seq_open_type(fname, SEQ_TYPE_GUESS);
}

SEQ* seq_get_all(const char* fname) {
    SEQ* sf;
    int count=0, i=0;

    sf = seq_open(fname);
    while ( seq_read(sf) != 0 )
        count++;
    seq_close(sf);

    sf = seq_open(fname);
    sf->contigs = (SEQ*)malloc(count*sizeof(SEQ));
    KEEP_SEQ = 1;
    sf->numContigs = count;
    while ( seq_read(sf) != 0 ) {
        sf->contigs[i].header = copy_string(sf->header);
        sf->contigs[i].seq = (uchar*)(copy_string((char*)sf->seq));
        sf->contigs[i].fp = sf->fp;
        sf->contigs[i].slen = sf->slen;
        sf->contigs[i].flipped = 'n';
        i++;
    }
    return sf;
}

SEQ* seq_close_all(SEQ* sf) {
    int i;

    if (sf != NULL )
        for (i=0; i<sf->numContigs; i++) {
            free(sf->contigs[i].header);
            free(sf->contigs[i].seq);
        }
    if ( sf->header != NULL )
        free(sf->header);
    if ( sf->seq != NULL )
        free(sf->seq);
    if ( sf->fp )
        fclose(sf->fp);
    return NULL;
}

SEQ *seq_copy(const SEQ *s) {
    SEQ *ss = ckallocz(sizeof(SEQ));
    *ss = *s;
    ss->seq = (uchar*)copy_string((const char*)s->seq);
    ss->header = copy_string(s->header);
    ss->fname = copy_string(s->fname);
    ss->maskname = copy_string(s->fname);
    ss->fp = 0; /* XXX - no subsequent seq_read operations allowed */
    ss->offset = 0; /* XXX - no subsequent seq_read operations allowed */
    return ss;
}

SEQ* seq_close(SEQ *s) {
    if (s) {
        if (!(s->flags & SEQ_IS_SUBSEQ)) {
            if (s->fp)
                fclose(s->fp);
            if (s->fname) ckfree(s->fname);
            if (s->header) ckfree(s->header);
            if (s->seq) ckfree(s->seq);
            if (s->maskname) ckfree(s->maskname);
        }
        memset(s, 0, sizeof(SEQ));
        ckfree(s);
    }
    return 0;
}

uchar dna_cmpl(uchar ch) {
    /* XXX - assumes ascii, returns space on error. */
    return dna_complement[ch];
}

void do_revcomp(uchar *s, int len) {
    uchar *p = s + len - 1;

    while (s<=p) {
        uchar c;

        c = dna_cmpl(*s);
        *s = dna_cmpl(*p);
        *p = c;
        ++s, --p;
    }
}

SEQ *seq_revcomp_inplace(SEQ *seq) {
    do_revcomp(SEQ_CHARS(seq), SEQ_LEN(seq));
    seq->flags ^= SEQ_IS_REVCOMP;
    return seq;
}

int seq_read_nib(SEQ *seq) {
    int n;
    if (seq->header) ZFREE(seq->header);
    if (seq->seq) ZFREE(seq->seq);

    seq->offset = ftell(seq->fp);
    seq->seq = seq_freadnib(seq->fp, seq->from-1, seq->slen, &(seq->slen));

    if (seq->seq == 0) return 0;

    n = strlen(seq->fname)+1+22+1+22+1;
    seq->header = ckalloc(n);
    seq->hlen = snprintf(seq->header, n, "%s:%d-%d",
                         seq->fname, seq->from, seq->from+seq->slen-1);

    seq->count++;
    return 1;
}

int seq_read(SEQ *seq) {
    if (seq->type == SEQ_TYPE_GUESS) {
        if (strstr(seq->fname, ".nib")) seq->type = SEQ_TYPE_NIB;
        else if (strstr(seq->fname, ".fa")) seq->type = SEQ_TYPE_FASTA;
    }

    switch (seq->type) {
    case SEQ_TYPE_GUESS:
    case SEQ_TYPE_FASTA:
        return seq_read_fasta(seq);
    case SEQ_TYPE_NIB:
        return seq_read_nib(seq);
    default:
        fatalf("unknown file type (%d)", seq->type);
    }
    /*NOTREACHED*/
    return 0;
}

SEQ *seq_get(const char *fname) {
    SEQ *s = seq_open(fname);
    int r = seq_read(s);
    if (r < 0)
        Fatalfr("could not read from %s", fname);
    else if (r == 0)
        return 0;
    else
        return s;
    /*NOTREACHED*/
    return 0;
}

