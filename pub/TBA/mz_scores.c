// mz_scores.c -- handles some details of alignment scores for multiz

#include "util.h"
#include "mz_scores.h"

static const char rcsid[] = "$Id: mz_scores.c 142 2008-11-12 18:55:23Z rico $";


int **ss70, **ss85, *gop70, *gop85;

static const uchar nchars[] = "ACGT";
static const int HOX70[4][4] = {	// for human-rodent
                                   {  91, -114,  -31, -123 },
                                   {-114,  100, -125,  -31 },
                                   { -31, -125,  100, -114 },
                                   {-123,  -31, -114,   91 },
                               };

static const int HOX85[4][4] = {	// for mouse-rat alignments
                                   {  86, -135,  -68, -157 },
                                   {-135,  100, -148,  -68 },
                                   { -68, -148,  100, -135 },
                                   {-157,  -68, -135,   86 },
                               };

#define GAP_OPEN70	400
#define GAP_EXTEND70	30

#define GAP_OPEN85	600
#define GAP_EXTEND85	50

#define UNSPECIFIED	-100	// score for an unspecified aligned pair

#define CLEN(s) (sizeof((s))-1)
#define NACHARS 128

static void init_scores(const int s[4][4], int **ss, int *gop,  int filler,
                        int gap_op, int gap_ex) {
    int i, j, a, b, A, B, X, D;

    // first, fill in ss, the substitution-score matrix
    for (i = 0; i < NACHARS; ++i)
        for (j = 0; j < NACHARS; ++j)
            ss[i][j] = filler;
    for (i = 0; i < (signed)CLEN(nchars); ++i) {
        A = nchars[i];
        a = tolower(A);
        for (j = 0; j < (signed)CLEN(nchars); ++j) {
            B = nchars[j];
            b = tolower(B);
            ss[A][B] = ss[a][B] = ss[A][b] = ss[a][b] =
                                                 s[i][j];
        }
    }
    for (i = 0; i < NACHARS; ++i)
        ss['-'][i] = ss[i]['-'] = -gap_ex;
    ss['-']['-'] = 0;

    // initialize the "gap-open penalty" array, for quasi-natural gap costs
    for (i = 0; i < 16; ++i)
        gop[i] = 0;
    D = 1;	// dash
    X = 0;	// anything other than '-'
    /* The six gap-open configurations are:
    	xx	x-	x-	-x	-x	--
    	x-	xx	-x	x-	--	-x
    	 1	 2	 3	 4	 5	 6
       The last two may result in a gap being counted twice, as with:
    	xx-xx
    	x---x
    	 1 5
       (A digit indicates the second column of a column-pair that conforms
       to the indicated rule, 1-6.)
       However, it permits a gap to be detected by inspection of adjacent
       columns in a case like:
    	x-----xx
    	x------x
    	      5
       GAP's arguments give column 1, followed by column 2.
    */
    GAP(X,X,X,D) = GAP(X,X,D,X) = GAP(X,D,D,X) =  GAP(D,X,X,D) =
                                      GAP(D,D,X,D) = GAP(D,D,D,X) = gap_op;
    gap_extend = gap_ex;
}

static int **alloc_scores(int **gop_ptr) {
    int **matrix, i;

    *gop_ptr = ckalloc(16*sizeof(int));
    matrix = ckalloc(NACHARS*sizeof(int *));
    matrix[0] = ckalloc(NACHARS*NACHARS*sizeof(int));
    for (i = 1; i < NACHARS; ++i)
        matrix[i] = matrix[i-1] + NACHARS;
    return matrix;
}

void init_scores70() {
    static int init = 1;

    if (init) {
        init = 0;
        ss70 = alloc_scores(&gop70);
        init_scores(HOX70, ss70, gop70, UNSPECIFIED, GAP_OPEN70,
                    GAP_EXTEND70);
    }
    ss = ss70;
    gop = gop70;
    gap_open = GAP_OPEN70;
    gap_extend = GAP_EXTEND70;
}

void init_scores85() {
    static int init = 1;

    if (init) {
        init = 0;
        ss85 = alloc_scores(&gop85);
        init_scores(HOX85, ss85, gop85, UNSPECIFIED, GAP_OPEN85,
                    GAP_EXTEND85);
    }
    ss = ss85;
    gop = gop85;
    gap_open = GAP_OPEN85;
    gap_extend = GAP_EXTEND85;
}

double mafScoreRange(struct mafAli *maf, int start, int size) {
    uchar ai, ar, bi, br;
    int i;
    double score;
    struct mafComp *c1, *c2;

    if (start < 0 || size <= 0 || start+size > maf->textSize)
        fatalf("mafScoreRange: start = %d, size = %d, textSize = %d\n",
               start, size, maf->textSize);
    if (ss == NULL)
        fatal("mafScoreRange: scores not initialized");

    score = 0.0;
    for (i = start; i < start+size; ++i) {
        for (c1 = maf->components; c1 != NULL; c1 = c1->next) {
            br = c1->text[i];
            for (c2 = c1->next; c2 != NULL; c2 = c2->next) {
                bi = c2->text[i];
                score += SS(br,bi);
                if (i > 0) {
                    ar = c1->text[i-1];
                    ai = c2->text[i-1];
                    score -= GAP2(ar,ai,br,bi);
                }
            }
        }
    }
    return score;
}

