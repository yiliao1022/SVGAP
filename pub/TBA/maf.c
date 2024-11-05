/*
 *   maf.c version 12
 */

// mini-maf.c - Read/write maf format.  Stolen from Jim Kent & seriously abused
#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "mz_scores.h"

static const char rcsid[] = "$Id: maf.c 142 2008-11-12 18:55:23Z rico $";


struct mafFile *mafOpen(char *fileName, int verbose) {
    struct mafFile *mf;
    FILE *fp;
    char buf[500], *s;

    mf = ckalloc(sizeof(struct mafFile));
    mf->next = NULL;
    fp = mf->fp = ckopen(fileName, "r");
    if (fgets(buf, 500, fp) == NULL) {
        if ( ferror(fp) )
            fprintf(stderr, "error file reading\n");
        fatalf("empty file %s", fileName);
    }
    if (sscanf(buf, "##maf version=%d", &(mf->version)) != 1)
        fatalf("improper maf header line: %s", buf);
    if ((s = strstr(buf, "scoring=")) != NULL)
        mf->scoring = copy_string(s+8);
    else
        mf->scoring = NULL;
    mf->alignments = NULL;
    mf->fileName = copy_string(fileName);
    mf->line_nbr = 0;
    mf->verbose = verbose;
    return mf;
}

/* ------- start of code to read arbitrarily long lines and echo comments --- */

//static int need(unsigned long n, char **linep, unsigned long *np) {
int need(unsigned long n, char **linep, unsigned long *np) {
    char *p;

    if (*np <= n) {
        *np += (*np >> 5) + 16;
        if ((p = realloc(*linep, *np)) == 0)
            return -1;
        *linep = p;
    }
    return 0;
}

//static unsigned long get_line(char **linep, unsigned long *np, FILE *fp) {
unsigned long get_line(char **linep, unsigned long *np, FILE *fp) {
    int ch;
    unsigned long n = 0;

    while ((ch = fgetc(fp)) != -1) {
        if (need(n+1, linep, np))
            return -1;
        (*linep)[n++] = ch;
        if (ch == '\n')
            break;
    }
    if (need(n+1, linep, np))
        return -1;
    (*linep)[n] = 0;
    if (n == 0 && ch == -1)
        return -1;
    return n;
}

//static unsigned long get_maf_line(char **linep, unsigned long *np, FILE *fp,
//  struct mafFile *mf) {
unsigned long get_maf_line(char **linep, unsigned long *np, FILE *fp,
                           struct mafFile *mf) {
    long nn;

    while ((nn = get_line(linep, np, fp)) > 1) {
        mf->line_nbr++;
        if  ((*linep)[0] == '#') {
            if (mf->verbose && strstr(*linep, "eof")==NULL)
                printf("%s", *linep);
        } else
            break;
    }
    return nn;
}

// parse line "a ..." to get score
int parseScoreLine(char* line, struct mafAli* ali) {
    char *ptr1=line+1, *ptr2, *buffer=NULL;
    int len, row, curr=0;
    struct mafComp* comp=ali->components;

    ali->score = (double)MIN_INT;
    while ( *ptr1 != '\0') {
        while ( *ptr1 == ' '  || *ptr1 == '\t')
            ptr1++;
        if ( *ptr1 == '\n' || *ptr1 == '\0')
            break;
        ptr2 = ptr1;
        while ( *ptr2 != ' ' && *ptr2 != '\t' && *ptr2 != '\n' && *ptr2 != '\0')
            ptr2++;
        len = ptr2 - ptr1;
        buffer = (char*)malloc((len+1)*sizeof(char));
        strncpy(buffer, ptr1, len);
        buffer[len] = '\0';
        if ( strncmp(buffer, "score=", 6) == 0 )
            ali->score = atof(buffer+6);
        else if ( strncmp(buffer, "amplifier=", 10) == 0 ) {
            row = atoi(buffer+10);
            while (curr < row ) {
                comp = comp->next;
                curr++;
            }
            comp->paralog = 'a';
        } else if ( strncmp(buffer, "copy=", 5) == 0) {
            row = atoi(buffer+5);
            while (curr < row) {
                comp = comp->next;
                curr++;
            }
            comp->paralog = 'c';
        }
        free(buffer);
        buffer = NULL;
        ptr1 = ptr2+1;
    }
    return 0;
}

/* ------- end of code to read arbitrarily long lines and echo comments --- */
struct mafAli *mafNext(struct mafFile *mf) {
    FILE *fp = mf->fp;
    struct mafAli *a = ckalloc(sizeof(struct mafAli));
    struct mafComp *c, *last;
    char buf[500], blockHeaderLine[1000];
    static char *line = NULL;
    static unsigned long nb = 0;
    int i, len;

    while ((len=get_maf_line(&line, &nb, fp, mf)) != -1)
        if ( line[0] != '#' && line[0] != '\n' && line[0] != ' ')
            break;
    if (len == -1) {
        fclose(fp);
        mf->fp = NULL;
        return NULL;
    }
    /*
    if (strncmp(line, "a score=", 8) == 0)
           	a->score = atof(line+8);
    else if (line[0] == 'a')
    	a->score = MIN_INT;
    */
    if ( strncmp(line, "a", 1) == 0)
        strcpy(blockHeaderLine, line);
    else
        fatalf("Expecting 'a (score=xxx)' in file %s, line %d:\n%s",
               mf->fileName, mf->line_nbr, line);
    a->textSize = 0;
    last = a->components = NULL;
    a->next = NULL;
    while ((len = get_maf_line(&line, &nb, fp, mf)) != -1 &&
            line[0] != '\n' && line[0] != ' ' && line[0] != '#') {
        c = ckalloc(sizeof(struct mafComp));

        c->text = ckalloc(len*sizeof(char));
        if ( line[0] != 's' )
            continue;
        if (sscanf(line, "s %s %d %d %c %d %s",
                   buf, &(c->start), &(c->size), &(c->strand),
                   &(c->srcSize), c->text) != 6)
            fatalf("bad component in file %s, line %d:\n%s",
                   mf->fileName, mf->line_nbr, buf);
        c->src = copy_string(buf);
        parseSrcName2(c);
        c->paralog = 's';        // singleton by default
        c->mafPosMap = NULL;
        //c->name = NULL;
        c->next = NULL;
        if (a->components == NULL) {
            a->textSize = strlen(c->text);
            a->components = c;
        } else {
            if (a->textSize != (signed)strlen(c->text))
                fatalf("line %d of %s: inconsistent row size",
                       mf->line_nbr, mf->fileName);
            last->next = c;
        }
        last = c;
        /* Do some sanity checking. */
        if (c->srcSize <= 0 || c->size <= 0)
            fatalf("Size <= 0 at line %d of file %s:\n%s",
                   mf->line_nbr, mf->fileName, line);
        if (c->start < 0 || c->start + c->size > c->srcSize) {
            if (c != a->components) {
                c = a->components;
                fprintf(stderr,
                        "in maf entry with top row %s:%d len = %d,\n",
                        c->src, c->start, c->size);
            }
            fatalf("Bad coordinates at line %d of file %s:\n%s",
                   mf->line_nbr, mf->fileName, line);
        }
        for (i = len = 0; i < a->textSize; ++i)
            if (c->text[i] != '-')
                ++len;
        if (len != c->size)
            fatalf("Actual size %d, claimed size %d at line %d of file %s:\n%s", len, c->size, mf->line_nbr, mf->fileName, line);

    }
    parseScoreLine(blockHeaderLine, a);
    mf->line_nbr++;
    return a;
}


struct mafFile *mafReadAll(char *fileName, int verbose) {
    struct mafFile *mf = mafOpen(fileName, verbose);
    struct mafAli *a, *last;

    for (last = NULL; (a = mafNext(mf)) != NULL; last = a)
        if (last == NULL)
            mf->alignments = a;
        else
            last->next = a;
    return mf;
}

void mafWriteStart(FILE *f, char *scoring) {
    fprintf(f, "##maf version=1 scoring=%s\n", scoring);
}

void mafWriteEnd(FILE *f) {
    fprintf(f, "##eof maf\n");
}

//static int digitsBaseTen(int x) {
int digitsBaseTen(int x) {
    int digCount;

    if (x < 0)
        fatalf("digitsBaseTen: negative argument %d", x);
    for (digCount = 1; x >= 10; x /= 10, ++digCount)
        ;
    return digCount;
}


void mafWrite(FILE *f, struct mafAli *a) {
    struct mafComp *c;
    int srcChars = 0, startChars = 0, sizeChars = 0, srcSizeChars = 0, row;
    char src[500], name[200], chr[200];

    fprintf(f, "a");
    if ( a->score != MIN_INT )
        fprintf(f, " score=%3.1f", a->score);
    for (row = 0, c = a->components; c != NULL; c = c->next, row++) {
        switch ( c->paralog ) {
        case 's':
            break;
        case 'a':
            fprintf(f, " amplifier=%d", row);
            break;
        case 'c':
            fprintf(f, " copy=%d", row);
            break;
        default:
            fatalf("Wrong character: \'%c\'", c->paralog);
        }
    }
    fprintf(f, "\n");

    /* Figure out length of each field. */
    for (c = a->components; c != NULL; c = c->next) {
        srcChars = MAX(srcChars, (signed int)strlen(c->src));
        startChars = MAX(startChars, digitsBaseTen(c->start));
        sizeChars = MAX(sizeChars,digitsBaseTen(c->size));
        srcSizeChars = MAX(srcSizeChars, digitsBaseTen(c->srcSize));
    }
    for (c = a->components; c != NULL; c = c->next) {
        parseSrcName(c->src, name, chr);
        strcpy(src, name);
        if ( strcmp( name, chr) != 0) {
            strcat(src, ".");
            strcat(src, chr);
        }
        fprintf(f, "s %-*s %*d %*d %c %*d %s\n",
                srcChars, src, startChars, c->start, sizeChars, c->size,
                c->strand, srcSizeChars, c->srcSize, c->text);
    }
    fprintf(f, "\n");	// blank separator line
}

void mafCompFree(struct mafComp **pComp) {
    struct mafComp *c = *pComp;

    if (c != NULL) {
        //if ( c->name != NULL)
        //free(c->name);
        free(c->src);
        free(c->text);
        free(c->contig);
        free(c->name);
        if ( c->mafPosMap != NULL )
            free(c->mafPosMap);
        free(c);
        *pComp = NULL;
    }
}

void mafCompsFree(struct mafComp* pComp) {
    struct mafComp* c, *next;
    for (c = pComp; c != NULL; c = next) {
        next = c->next;
        mafCompFree(&c);
    }
}

void mafAliFree(struct mafAli **pAli) {
    struct mafAli *a;

    if ( pAli == NULL || *pAli == NULL)
        return;
    a = *pAli;
    mafCompsFree(a->components);
    a->components = NULL;
    /*
    for (c = a->components; c != NULL; c = next) {
      next = c->next;
      mafCompFree(&c);
    }
    */
    free(a);
    *pAli = NULL;
}

void mafFileFree(struct mafFile **pmf) {
    struct mafFile *mf = *pmf;
    struct mafAli *a, *next;

    if (mf->fp != NULL && strcmp(mf->fileName, "/dev/stdin") != 0 )
        fclose(mf->fp);
    if (mf->scoring != NULL)
        free(mf->scoring);
    free(mf->fileName);
    for (a = mf->alignments; a != NULL; a = next) {
        next = a->next;
        mafAliFree(&a);
    }
    free(mf);
    *pmf = NULL;
}

// rm columns which contain all dashes.
struct mafAli* mafColDashRm(struct mafAli *a) {
    int i, col;
    struct mafComp *c;

    if (a==NULL)
        return NULL;

    for (i = col = 0; col < a->textSize; ++col) {
        for (c = a->components; c != NULL; c = c->next)
            if (c->text[col] != '-')
                break;
        if (c != NULL) {
            if (i < col)
                for (c = a->components; c != NULL; c = c->next)
                    c->text[i] = c->text[col];
            ++i;
        }
    }
    if (i < a->textSize) {
        a->textSize = i;
        for (c = a->components; c != NULL; c = c->next)
            c->text[i] = '\0';
    }
    return a;
}

// rm components which contain all dashes
struct mafAli* mafRowDashRm(struct mafAli *ali) {
    struct mafComp* comp, *prev;
    char* s;

    if ( ali==NULL)
        return NULL;
    for (comp=prev=ali->components; comp != NULL;) {
        for (s=comp->text; *s != '\0'; s++)
            if ( *s != '-' )
                break;
        if ( *s == '\0' ) { // delete this comp
            if ( comp == ali->components ) {
                ali->components = comp->next;
                comp->next = NULL;
                mafCompFree(&comp);
                comp = prev = ali->components;
            } else {
                prev->next = comp->next;
                comp->next = NULL;
                mafCompFree(&comp);
                comp = prev->next;
            }
            continue;
        }
        prev = comp;
        comp = comp->next;
    }

    if ( ali->components == NULL) {
        mafAliFree(&ali);
        return NULL;
    }
    return ali;
}

struct mafComp *mafNewComp(char *src, int start, int size, char strand, int srcSize, int len, char paralog, char *name, char* contig) {
    struct mafComp *c = ckalloc(sizeof(struct mafComp));

    c->src = copy_string(src);
    c->name = copy_string(name);
    c->contig = copy_string(contig);
    c->start = start;
    c->size = size;
    c->strand = strand;
    c->srcSize = srcSize;
    c->text = ckalloc((len+1)*sizeof(char));
    c->paralog = paralog;
    c->mafPosMap = NULL;
    c->next = NULL;
    return c;
}


struct mafComp* mafCpyComp(struct mafComp* t) {
    struct mafComp* c = ckalloc(sizeof(struct mafComp));

    c->src = copy_string(t->src);
    c->name = copy_string(t->name);
    c->contig = copy_string(t->contig);
    c->start = t->start;
    c->size = t->size;
    c->strand = t->strand;
    c->srcSize = t->srcSize;
    c->paralog = t->paralog;
    c->next = NULL;
    c->mafPosMap = NULL;
    return c;
}

struct mafAli* mafNewAli(double score, int textSize) {
    struct mafAli *a = ckalloc(sizeof(struct mafAli));

    a->components = NULL;
    a->next = NULL;
    a->score = score;
    a->textSize = textSize;
    return a;
}

struct mafAli* duplicate_ali(struct mafAli* template) {
    struct mafAli* newAli;
    struct mafComp *ncomp, *tcomp, *pcomp;

    newAli = (struct mafAli*)malloc(sizeof(struct mafAli));
    newAli->components = NULL;

    for ( tcomp = template->components; tcomp != NULL; tcomp = tcomp->next) {
        ncomp = mafCpyComp(tcomp);
        ncomp->text = copy_string(tcomp->text);

        if ( newAli->components == NULL )
            newAli->components = ncomp;
        else {
            for (pcomp = newAli->components; pcomp->next != NULL; pcomp=pcomp->next)
                ;
            pcomp->next = ncomp;
        }
    }
    newAli->next = NULL;
    newAli->textSize = template->textSize;
    newAli->score = template->score;
    return newAli;
}

struct mafAli* make_part_ali(struct mafAli* template, int cbeg, int cend) {
    struct mafAli* retAli;
    struct mafComp *ncomp, *tcomp, *pcomp;
    int i;

    retAli = (struct mafAli*)malloc(sizeof(struct mafAli));
    retAli->components = NULL;

    for ( tcomp=template->components; tcomp!=NULL; tcomp=tcomp->next) {
        ncomp = mafCpyComp(tcomp);
        ncomp->text = (char*)malloc( (cend-cbeg+2)*sizeof(char));
        for (i=cbeg; i<=cend; i++)
            ncomp->text[i-cbeg] = tcomp->text[i];
        ncomp->text[cend-cbeg+1] = '\0';
        for (i=0, ncomp->start=tcomp->start-1; i<cbeg; i++)
            if (tcomp->text[i] != '-')
                ncomp->start++;
        ncomp->start++;
        for (i=cbeg, ncomp->size=0; i<=cend; i++)
            if (tcomp->text[i] != '-')
                ncomp->size++;
        if ( retAli->components == NULL)
            retAli->components = ncomp;
        else {
            for (pcomp=retAli->components; pcomp->next!=NULL; pcomp=pcomp->next)
                ;
            pcomp->next = ncomp;
        }
    }
    retAli->next = NULL;
    retAli->textSize = cend - cbeg+1;
    retAli = mafRowDashRm(retAli);
    if ( retAli != NULL )
        retAli->score = mafScoreRange(retAli, 0, cend-cbeg+1);
    return retAli;
}
