/*
 *  mz_preyama.c version 13
 *
 *  revision on Feb. 2, 2006
 */


#include "util.h"
#include "multi_util.h"
#include "mz_preyama.h"
#include "mz_yama.h"
#include "maf.h"
#include "mz_scores.h"
#include "align_util.h"

static const char rcsid[] = "$Id: mz_preyama.c 142 2008-11-12 18:55:23Z rico $";

// smooth out the left and right bounds for the dynamic-programming region
void smooth(int *LB, int *RB, int M, int N, int radius) {
    int i, j, radi = MIN(M, radius);

    // make LB and RB monotone
    for (i = j = 0; i <= M; ++i)
        LB[i] = j = MAX(j, LB[i]);
    for (i = M, j = N; i >= 0; --i)
        RB[i] = j = MIN(j, RB[i]);

    // expand line into a sausage
    for (i = M; i > radi; --i)
        LB[i] = MIN(MAX(LB[i]-radi,0),LB[i-radi]);
    for ( ; i >= 0; --i)
        LB[i] = 0;
    for (i = 0; i < M - radi; ++i)
        RB[i] = MAX(MIN(RB[i]+radi,N),RB[i+radi]);
    for ( ; i <= M; ++i)
        RB[i] = N;
}

// construct a new mafAli structure with the block returned from yama procedure
struct mafAli *mafBuild(uchar **A_new, int nrow, int ncol, struct mafAli *a2, int cbeg2, struct mafAli *a3, int cbeg3, int top) {
    struct mafAli *A = ckalloc(sizeof(struct mafAli));
    struct mafComp *nc, *d, *prev_nc = NULL;
    int i, j, beg, strt;

    A->textSize = ncol;
    A->next = NULL;
    A->components = NULL;
    d = a2->components;
    beg = cbeg2;
    for (i = 0; i < nrow; ++i, d = d->next) {
        if (d == NULL) {
            if ( top == 0 )
                d = a3->components->next;
            else
                d = a3->components;
            beg = cbeg3;   // if [cbeg3] happens to be non'-', it's fine
        }
        for (strt = d->start-1, j = 0; j < beg; ++j)
            if (d->text[j] != '-')
                ++strt;
        nc = mafCpyComp(d);
        nc->start = strt+1;
        nc->size = 0;
        nc->text = (char*)malloc( (ncol+1)*sizeof(char));
        for (j = 0; j < ncol; ++j)
            if ((nc->text[j] = A_new[j+1][i]) != '-')
                ++(nc->size);
        nc->text[ncol] = '\0';
        if ( nc->size == 0) {
            mafCompFree(&nc);
            continue;
        }
        if (prev_nc == NULL)
            A->components = nc;
        else
            prev_nc->next = nc;
        prev_nc = nc;
    }
    if ( A->components==NULL)
        return NULL;
    A->score = mafScoreRange(A, 0, ncol);
    return A;
}

// rm columns with all dashes, update N value also
// X[m][n] where m is col, n is row, m starts 1, n starts 0
// after this procedure, X and N are reassigned values so that no all-gaps col
// below ref row
int* rmColDash(uchar** X, int* N, int row) { // x starts from 1
    int i, j, k;
    int* mapArray;
    mapArray = (int*)malloc((*N + 1)*sizeof(int));
    for (i=1; i<=*N; i++)   // starts from 1
        mapArray[i] = -1;
    for (i=k=1; i<=*N; i++) { // X starts from 1
        for (j=0; j<row; j++)
            if ( X[i][j]!='-') //********!!!**********
                //if ( X[i][j+1]  != '-' ) // the non-ref seq
                break;
        if (j!=row) { // non-all-dash
            if ( i!=k)
                for (j=0; j<row; j++)
                    X[k][j] = X[i][j];
            mapArray[i] = k;
            k++;
        }
    }
    *N = k-1;
    return mapArray;  // maps from old columns of X to new columns
}

// maps from A to B
int* mapping(uchar** A, int a_row1, int a_row2, int a_col1, int a_col2, uchar** B, int b_row1, int b_row2, int b_col1, int b_col2) {
    int *mapArray, i, j, k, l;

    if ( a_row2 - a_row1 != b_row2 - b_row1 )
        fatalf("not equal rows:!\n");
    mapArray = (int*)malloc((a_col2-a_col1+2)*sizeof(int));
    for (i=a_col1; i<=a_col2; i++)
        mapArray[i-a_col1+1] = -1;

    i=a_col1;
    k=b_col1;
    while (i<=a_col2 && k<=b_col2) {
        j=l=1000000;   // assign j and l large values
        while (i<=a_col2) {
            for (j=a_row1; j<=a_row2; j++)
                if ( A[i][j] != '-')
                    break;
            if ( j > a_row2 )  // non-all-dashes
                i++;
            else
                break;
        }
        while (k<=b_col2) {
            for (l=b_row1; l<=b_row2; l++)
                if ( B[k][l] != '-')
                    break;
            if ( l > b_row2 )
                k++;
            else
                break;
        }
        if ( j<=a_row2 && l<=b_row2)
            mapArray[i] = k;
        i++;
        k++;
    }
    return mapArray;
}

// determine upper and lower bounds for yama procedure
// v can only be 0, or 1
struct mafAli* pre_yama(struct mafAli *a1, struct mafAli* a2, int beg, int end,int radius, int v, FILE* fpw2) {
    int cbeg1, cend1, cbeg2, cend2, M, N, K, L, i, j, M_new, M3,N3;
    int *LB, *RB, *LB2, *RB2;
    uchar **A, **B, **AL_new, **AL_new2, **A2, **B2;
    struct mafComp* comp;
    struct mafAli* val;
    int tmp1, tmp2, M_cp, N_cp;
    int *map1, *map2, *map3, *map4;
    int curr1, curr2;

    for (K=0, comp=a1->components; comp!=NULL; comp=comp->next, K++)
        ;
    for (L=0, comp=a2->components->next; comp!=NULL; comp=comp->next, L++)
        ;

    cbeg1 = mafPos2Col(a1->components, beg, a1->textSize); // starts from 0
    cend1 = mafPos2Col(a1->components, end, a1->textSize);
    cbeg2 = mafPos2Col(a2->components, beg, a2->textSize);
    cend2 = mafPos2Col(a2->components, end, a2->textSize);

    M = cend1 - cbeg1 + 1;  // length
    N = cend2 - cbeg2 + 1;
    B = (uchar**)malloc(N*sizeof(uchar*)) - 1;
    B[1] = (uchar*)malloc(L*N*sizeof(uchar));
    for (i=2; i<=N; i++)
        B[i] = B[i-1] + L;
    for (i=1; i<=N; i++)
        for (j=0, comp=a2->components->next; j<L; j++, comp=comp->next)
            B[i][j] = comp->text[cbeg2 + i - 1];    // starts from 1
    N_cp = N;
    map2 = rmColDash(B, &N, L);   // map2 starts from 1
    if ( N < 1) {
        free(B[1]);
        free(B+1);
        free(map2);
        return NULL;
    }

    A = (uchar**)malloc(M*sizeof(uchar*)) - 1;
    if ( v == 0 )                              // if top of A is not fixed
        K--;                                     // K needs to be subtracted
    if ( K==0 ) {
        if (fpw2 != NULL)
            print_part_ali_col(a2, cbeg2, cend2, fpw2);
        free(B[1]);
        free(B+1);
        free(map2);
        free(A+1);
        return NULL;
    }

    A[1] = (uchar*)malloc(K*M*sizeof(uchar));
    for (i=2; i<=M; i++)
        A[i] = A[i-1] + K;

    for (i=1; i<=M; i++) {
        if ( v==0)                               // to exclude the top ref
            comp = a1->components->next;
        else                                      // to include the top ref
            comp = a1->components;                  // since it's going to be fixed
        for (j= 0; j<K; j++, comp=comp->next)
            A[i][j] = comp->text[cbeg1 + i - 1];
    }
    M_cp = M;

    if (v==0) {
        map1 = rmColDash(A, &M, K);   // starts from 1
        if (M < 1) {
            free(B[1]);
            free(B+1);
            free(map2);
            free(A[1]);
            free(A+1);
            free(map1);
            return NULL;
        }
    } else {
        map1 = (int*)malloc((M+1)*sizeof(int));
        for (i=1; i<=M; i++)
            map1[i] = i;
    }
    LB = (int*)malloc((M+1)*sizeof(int));       // Initialize LB, RB
    RB = (int*)malloc((M+1)*sizeof(int));       // they're tuned by two top refs
    for (i=0; i<=M; i++) {
        LB[i] = 0;
        RB[i] = N;
    }

    for (i=cbeg1, j=cbeg2; i<=cend1; ) {       // revision on feb.1 2006
        while ( a1->components->text[i] == '-' )
            i++;
        while ( a2->components->text[j] == '-' )
            j++;
        curr1 = map1[i-cbeg1+1];
        curr2 = map2[j-cbeg2+1];
        if ( curr1 == -1 || curr2 == -1 ) {
            i++;
            j++;
            continue;
        }
        if ( LB[curr1] == 0 || LB[curr1] > curr2 )
            LB[curr1] = curr2;
        if ( RB[curr1] == N || RB[curr1] < curr2 )
            RB[curr1] = curr2;
        i++;
        j++;
    }
    smooth(LB, RB, M, N, radius);
    yama(A, K, M, B, L, N, LB, RB, &AL_new, &M_new);
    if ( v == 1) {
        val = mafBuild(AL_new, K+L, M_new, a1, cbeg1, a2, cbeg2, 0);
        free(LB);
        free(RB);
    } else {  // need to do one more yama to align ref and AL_new
        free(LB);
        free(RB);

        A2 = (uchar**)malloc(M_cp*sizeof(uchar*)) - 1;   // Creat single ref block A
        A2[1] = (uchar*)malloc(M_cp*sizeof(uchar));
        for (i=2; i<=M_cp; i++)
            A2[i] = A2[i-1] + 1;
        for (i=1; i<=M_cp; i++)
            A2[i][0] = a1->components->text[cbeg1 + i - 1];

        M3 = M_cp;
        map3 = rmColDash(A2, &M3, 1);   // mapping from top ref to single seq

        map4 = mapping(A, 1, K, 1, M, AL_new, 0, K-1, 1, M_new); // mapping from A to AL_new

        LB2=(int*)malloc((M3+1)*sizeof(int));  // M3 is the length for seq
        RB2=(int*)malloc((M3+1)*sizeof(int));
        for (i=0; i<=M3; i++) {
            LB2[i] = 0;
            RB2[i] = M_new;
        }
        for (i=1; i<=M_cp; i++) { // walking along the original top row in a1
            tmp1 = map3[i];
            if ( map1[i] == -1)
                continue;
            tmp2 = map4[map1[i]];  // mapping from (mapping from original align) to AL_new
            if ( tmp1 == -1 || tmp2 == -1 )
                continue;
            if ( LB2[tmp1] == 0 || LB2[tmp1] > tmp2)
                LB2[tmp1] = tmp2;
            if ( RB2[tmp1] == M_new || RB2[tmp1] < tmp2)
                RB2[tmp1] = tmp2;
        }
        smooth(LB2, RB2, M3, M_new, radius);

        free(map3);
        free(map4);
        B2 = (uchar**)malloc(N_cp*sizeof(uchar*)) - 1;   // Creat single ref block A
        B2[1] = (uchar*)malloc(N_cp*sizeof(uchar));
        for (i=2; i<=N_cp; i++)
            B2[i] = B2[i-1] + 1;
        for (i=1; i<=N_cp; i++)
            B2[i][0] = a2->components->text[cbeg2 + i - 1];
        N3=N_cp;
        map3 = rmColDash(B2, &N3, 1);
        map4 = mapping(B, 0, L-1, 1, N, AL_new, K, K+L-1, 1, M_new);
        LB=(int*)malloc((N3+1)*sizeof(int));
        RB=(int*)malloc((N3+1)*sizeof(int));
        for (i=0; i<=N3; i++) {
            LB[i]=0;
            RB[i]=M_new;
        }
        for (i=1; i<=N_cp; i++) { // walking along the original top row in a2
            tmp1 = map3[i];
            tmp2 = map4[map2[i]];
            if ( tmp1 == -1 || tmp2 == -1 )
                continue;
            if ( LB[tmp1] == 0 || LB[tmp1] > tmp2)
                LB[tmp1] = tmp2;
            if ( RB[tmp1] == M_new || RB[tmp1] < tmp2)
                RB[tmp1] = tmp2;
        }
        smooth(LB, RB, N3, M_new, radius);
        if ( M3 != N3)
            fatal("M3 not equals N3!!\n");
        for ( i=0; i<=M3; i++) {
            LB[i] = MIN(LB[i], LB2[i]);
            RB[i] = MAX(RB[i], RB2[i]);
        }
        yama(A2, 1, M3, AL_new, K+L, M_new, LB, RB, &AL_new2, &M_new);
        val = mafBuild(AL_new2, K+L+1, M_new, a1, cbeg1, a2, cbeg2, 0);
        free(map4);
        free(map3);
        free(A2[1]);
        free(A2+1);
        free(B2[1]);
        free(B2+1);
        free(LB);
        free(RB);
        free(LB2);
        free(RB2);
        free(AL_new2[1]);
        free(AL_new2+1);
    }
    free(A[1]);
    free(A+1);
    free(B[1]);
    free(B+1);
    free(AL_new[1]);
    free(AL_new+1);
    free(map1);
    free(map2);
    return val;
}

// verify equality of characters in what are presumed to be identical sequences.
static void must_equal(char x, char y, int pos, // starting position
                       struct mafComp *c1, int cbeg1, int cend1,	// begin and end columns
                       struct mafComp *c2, int cbeg2, int cend2,
                       int file2) {					// 2 or 3
    int i;
    char z;

    if (toupper(x) != toupper(y)) {
        fprintf(stderr, "%c != %c\n", x, y);
        fprintf(stderr, "in file 1, positions %d... of %s are:\n  ",
                pos, c1->src);
        for (i = cbeg1; i <= cend1; ++i)
            if ((z = c1->text[i]) != '-')
                fputc(z, stderr);
        fputc('\n', stderr);
        fprintf(stderr, "while in file %d they are:\n  ", file2);
        for (i = cbeg2; i <= cend2; ++i)
            if ((z = c2->text[i]) != '-')
                fputc(z, stderr);
        fputc('\n', stderr);
        exit(1);
    }
}


struct mafAli* pre_yama2(struct mafAli* a1, struct mafAli* a2, struct mafAli* a3, int beg1, int end1, int begN, int endN, int radius, struct pwuAliFiles* pws) {

    struct mafAli *val;
    struct mafComp *c, *c1, *c2, *c3, *d;
    int cbeg1, cend1, cbeg2, cend2, cbeg3, cend3, i, j, k,
    M, N, K, L, M_new, *LB, *RB;
    uchar **A, **B, **AL_new, x, y, z;

    if ((c = a1->components) == NULL || (c1 = c->next) == NULL)
        fatal("pre_yama: cannot find c and c1");
    if (c1->next != NULL)
        fatal("pre_yama: a1 is not a pairwise alignment");
    // count K = number of rows in a2
    for (d = c2 = a2->components, K = 0; d != NULL; ++K, d = d->next)
        ;
    // count L = number of rows in a3
    for (d = c3 = a3->components, L = 0; d != NULL; ++L, d = d->next)
        ;
    if (K == 0 || L == 0)
        fatal("pre_yama: an alignment has 0 rows");
    if (!same_string(c->src, c2->src))
        fatalf("pre_yama: first rows for sequences %s != %s",
               c->src, c2->src);
    if (c->srcSize != c2->srcSize)
        fatalf("pre_yama: first row srcSizes %d != %d",
               c->srcSize, c1->srcSize);
    if (c->strand != c2->strand)
        fatalf("pre_yama: first rows on opposite strands");
    if (!same_string(c1->src, c3->src))
        fatalf("pre_yama: first rows (2) for sequences %s != %s",
               c1->src, c3->src);
    if (c1->srcSize != c3->srcSize)
        fatalf("pre_yama: first row (2) srcSizes %d != %d",
               c1->srcSize, c3->srcSize);
    if (c1->strand != c3->strand)
        fatalf("pre_yama: first rows (2) on opposite strands");

    cbeg1 = mafPos2Col(c, beg1, a1->textSize);
    cend1 = mafPos2Col(c, end1, a1->textSize);
    // check that a1 aligns beg1 to end1, and begN to endN
    if (cbeg1 != mafPos2Col(c1, begN, a1->textSize))
        fatalf("pre_yama: mismatch of beg1 and begN");
    if (cend1 != mafPos2Col(c1, endN, a1->textSize))
        fatalf("pre_yama: mismatch of end1 and endN");
    cbeg2 = mafPos2Col(c2, beg1, a2->textSize);
    cend2 = mafPos2Col(c2, end1, a2->textSize);
    cbeg3 = mafPos2Col(c3, begN, a3->textSize);
    cend3 = mafPos2Col(c3, endN, a3->textSize);

    if ( connectionAgreement2(a2, a3, cbeg2, cend2, cbeg3, cend3, pws) == 0 )
        return 0;

    // determine and check sizes of slices of a2 and a3
    M = cend2 - cbeg2 + 1;
    N = cend3 - cbeg3 + 1;
    if (MAX(M, N) < 2)
        return NULL;

    // prepare first alignment for yama
    A = (uchar **)ckalloc(M*sizeof(uchar *)) - 1;
    A[1] = ckalloc(K*M*sizeof(uchar));
    for (i = 2; i <= M; ++i)
        A[i] = A[i-1] + K;
    for (i = 1; i <= M; ++i)
        for (j = 0, d = c2; j < K; ++j, d = d->next)
            A[i][j] = d->text[cbeg2 + i - 1];

    // prepare second alignment for yama
    B = (uchar **)ckalloc(N*sizeof(uchar *)) - 1;
    B[1] = ckalloc(L*N*sizeof(uchar));
    for (i = 2; i <= N; ++i)
        B[i] = B[i-1] + L;
    for (i = 1; i <= N; ++i)
        for (j = 0, d = c3; j < L; ++j, d = d->next)
            B[i][j] = d->text[cbeg3 + i - 1];

    // Determine sausage
    LB = ckalloc((M+1)*sizeof(int));
    RB = ckalloc((M+1)*sizeof(int));
    for (i = 0; i <= M; ++i) {
        LB[i] = 0;
        RB[i] = N;
    }

    // Column k of c and c1 corresponds to column i of A and column j of B.
    // Several "reality checks" guard against data/internal errors.
    for (i = j = 0, k = cbeg1; k <= cend1; ++k) {
        if ((x = c->text[k]) != '-') {
            do {
                if (++i > M)
                    fatal("pre_yama: bad scene");
            } while ((z = A[i][0]) == '-');
#ifdef RM_HACK
            // RepeatMasker .align files sometimes give incorrect
            // sequence positions
            if (toupper(x) != toupper(z))
                return NULL;
#else
            must_equal(x, z, beg1, c, cbeg1, cend1,
                       c2, cbeg2, cend2, 2);
#endif
        }
        if ((y = c1->text[k]) != '-') {
            do {
                if (++j > N)
                    fatal("pre_yama: ouch");
            } while ((z = B[j][0]) == '-');
            must_equal(y, z, begN, c1, cbeg1, cend1,
                       c3, cbeg3, cend3, 3);
        }
        // For fixed i, LB[i] and RB[i] are the minimum and maximum
        // corresponding j, resp.
        if (LB[i] == 0)
            LB[i] = j;
        RB[i] = j;
    }
    if (i != M || j != N)
        fatalf("pre_yama: i = %d, M = %d, j = %d, N = %d", i, M, j, N);

    smooth(LB, RB, M, N, radius);

    yama(A, K, M, B, L, N, LB, RB, &AL_new, &M_new);

    val = mafBuild(AL_new, K+L, M_new, a2, cbeg2, a3, cbeg3, 1);

    free(AL_new[1]);
    free(AL_new+1);
    free(A[1]);
    free(A+1);
    free(B[1]);
    free(B+1);
    free(LB);
    free(RB);
    return val;
}
