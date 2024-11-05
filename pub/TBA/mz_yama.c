/* mz_yama.c -- dynamic programming procedure to align two alignments
*
* The algorithm uses sum-of-pairs substitution scores and the "quasi-natural
* gap costs" developed by:
*	Altschul, S. (1989) Gap costs for multiple sequence alignment.
*	J. theor. Biol. 138, 297-309.
* The graph model for aligning using affine gap scores, which employs an
* (M+1)-by-(N+1) grid (when the sequence lengths are M and N) with nodes called
* C, D and I at each grid point, is described by:
*	E. Myers and W. Miller (1989) Approximate matching of regular
*	expressions. Bull. Math. Biol. 51, 5-37.
* A slight modification to make the model appropriate for quasi-natural gap
* costs is explained by:
*	K.-M. Chao, R. Hardison and W. Miller (1994) Recent developments in
*	linear-space alignment methods: a survey. J. Comp. Biol. 1, 271-291.
* End-gaps (overhangs) are not charged a gap-open penalty.
*/

#include "util.h"
#include "maf.h"
#include "mz_scores.h"
#include "mz_yama.h"

static const char rcsid[] = "$Id: mz_yama.c 142 2008-11-12 18:55:23Z rico $";


#define FLAG_C (0)
#define FLAG_I (1<<0)
#define FLAG_D (1<<1)
#define SELECT_CID (FLAG_I | FLAG_D | FLAG_C)
// MININT is a hugely negative number away from the underflow threshold
#define MININT INT_MIN/2

/* Dynamic Programming structure. */
typedef struct DP {
    int D, C, I;
    int unused;       /* padding for efficiency with linux/GCC */
}
dp_node;

// populate a column from a column of height K and a column of height L
static void new_col(uchar *from1, int K, uchar *from2, int L, uchar *to_col) {
    //void new_col(uchar *from1, int K, uchar *from2, int L, uchar *to_col) {
    int i;

    for (i = 0; i < K; ++i)
        to_col[i] = from1[i];
    for (i = 0; i < L; ++i)
        to_col[i+K] = from2[i];
}

// Entry point: (documented in mz_yama.h)
void yama(uchar **A, int K, int M, uchar **B, int L, int N, int *LB, int *RB,
          uchar ***OAL, int *OM) {
    int C, D, I, i, j, k, m, m_new, n, row, col, s, t, u, v, x, y, z,
    diag_c, diag_d, diag_i, tback_size, nedit;
    uchar **tback_row, *tback, *tbp, st, flag_d, flag_i, flag_c, node,
    *dashes, *script, **al_new;
    dp_node *dp;

    if (LB[0] != 0 || RB[M] != N)
        fatalf("LB and RB not terminated properly: %d %d %d", LB[0], RB[M], N);
    for (tback_size = row = 0; row <= M; row++) {
        j = RB[row] - LB[row];
        // temporary sanity check:
        if (j < MIN(N,10)) {
            fatalf("RB[%d] - LB[%d] < %d, %d %d %d", row, row, MIN(N,10), RB[row], LB[row], N);
        }
        tback_size += (j+1);
        if (row > 0 && LB[row] < LB[row-1])
            fatal("LB not monotonic");
        if (row > 0 && RB[row] < RB[row-1])
            fatal("RB not monotonic");
    }

    dashes = ckalloc(MAX(K,L)*sizeof(uchar));
    for (i = 0; i < MAX(K,L); ++i)
        dashes[i] = '-';

    tback = ckalloc(tback_size * sizeof(uchar));
    tback_row = ckalloc((M+1)*sizeof(uchar *));
    tbp = tback_row[0] = tback;
    *(tbp++) = 0;	/* unused */

    dp = ckalloc((N+1)*sizeof(dp_node));
    dp[0].C = dp[0].D = dp[0].I = 0;
    flag_i = FLAG_I;
    for (col = 1; col <= RB[0]; ++col) {
        dp[col].C = dp[col].D = MININT;
        for (j = n = 0; j < L; ++j)
            if (B[col][j] != '-')
                ++n;
        dp[col].I = dp[col-1].I - n*K*gap_extend;
        *(tbp++) = (flag_i << 4);
    }
    for ( ; col <= N; ++col)
        dp[col].C = dp[col].D = dp[col].I = MININT;

    C = D = I = MININT;
    for (row = 1; row <= M; ++row) {
        tback_row[row] = tbp - LB[row];

        col = LB[row] - 1;
        if (LB[row-1] <= col) {
            diag_c = dp[col].C;
            diag_d = dp[col].D;
            diag_i = dp[col].I;
        } else
            diag_c = diag_d = diag_i = MININT;

        C = D = I = MININT;
        while (++col <= RB[row]) {
            // C, D, and I values are for the previous column
            // diag_c, diag_d and diag_i are for (row-1,col-1)

            // compute I and flag_i
            if (col > LB[row]) {
                x = C;
                y = D;
                z = I;
                /* The main complexity is for the quasi-natural
                *  gap costs, which depend on the types (C,D,I)
                *  of the last two edges in a path, i.e., on
                *  the types of the last two nodes.
                */
                if (row < M)	// not for end-gaps
                    for (i = 0; i < K; ++i) {
                        s = (A[row][i] == '-');
                        u = 1;
                        for (j = 0; j < L; ++j) {
                            t = (col > 1 &&
                                 B[col-1][j] == '-');
                            v = (B[col][j] == '-');
                            if (col > LB[row-1]+1)
                                x -= GAP(s,t,u,v);
                            y -= GAP(s,1,u,v);
                            if (col > LB[row]+1)
                                z -= GAP(1,t,u,v);
                        }
                    }
                if (x >= y && x >= z) {
                    /* It is best to precede the insertion
                    *  edge with a substitution edge, i.e.,
                    *  to reach the I node at grid point
                    *  (row,col) from the C node at
                    *  (row-1,col).
                    */
                    I = x;
                    flag_i = FLAG_C;
                } else if (y > z) {
                    // It is best to come from the D node.
                    I = y;
                    flag_i = FLAG_D;
                } else {
                    I = z;
                    flag_i = FLAG_I;
                }
                /* Each of n non-dashes in the inserted column
                * is aligned to K dashes.
                */
                for (j = n = 0; j < L; ++j)
                    if (B[col][j] != '-')
                        ++n;
                I -= n*K*gap_extend;

            } else {
                I = MININT;
                flag_i = 0;	// unused
            }

            // compute C and flag_c
            if (col > LB[row-1]) {
                x = diag_c;
                y = diag_d;
                z = diag_i;
                if (col > 1)	// no gap-open penalty at start
                    for (i = 0; i < K; ++i) {
                        s = (row > 1 && A[row-1][i] == '-');
                        u = (A[row][i] == '-');
                        for (j = 0; j < L; ++j) {
                            t = (B[col-1][j] == '-');
                            v = (B[col][j] == '-');
                            if (row > 1 &&
                                    col > LB[row-2] + 1)
                                x -= GAP(s,t,u,v);
                            if (row > 1)
                                y -= GAP(s,1,u,v);
                            if (col > LB[row-1] + 1)
                                z -= GAP(1,t,u,v);
                        }
                    }
                if (x >= y && x >= z) {
                    C = x;
                    flag_c = FLAG_C;
                } else if (y > z) {
                    C = y;
                    flag_c = FLAG_D;
                } else {
                    C = z;
                    flag_c = FLAG_I;
                }
                for (i = 0; i < K; ++i)
                    for (j = 0; j < L; ++j)
                        C += SS(A[row][i],B[col][j]);
            } else {
                C = MININT;
                flag_c = 0;
            }

            // compute D and flag_d
            x = dp[col].C;
            y = dp[col].D;
            z = dp[col].I;
            if (0 < col && col < N)	// not for end-gaps
                for (i = 0; i < K; ++i) {
                    s = (row > 1 && A[row-1][i] == '-');
                    u = (A[row][i] == '-');
                    v = 1;
                    for (j = 0; j < L; ++j) {
                        t = B[col][j] == '-';
                        if (row > 1 && col > LB[row-2])
                            x -= GAP(s,t,u,v);
                        if (row > 1)
                            y -= GAP(s,1,u,v);
                        if (col > LB[row-1])
                            z -= GAP(1,t,u,v);
                    }
                }
            if (x >= y && x >= z) {
                D = x;
                flag_d = FLAG_C;
            } else if (y > z) {
                D = y;
                flag_d = FLAG_D;
            } else {
                D = z;
                flag_d = FLAG_I;
            }
            // Each of n non-dashes in the deleted column
            // is aligned to L dashes.

            for (j = n = 0; j < K; ++j)
                if (A[row][j] != '-')
                    ++n;
            D -= n*L*gap_extend;

            // complete processing of the grid point (row,col)
            diag_c = dp[col].C;
            diag_d = dp[col].D;
            diag_i = dp[col].I;
            dp[col].C = C;
            dp[col].D = D;
            dp[col].I = I;

            // pack three flags into one byte of traceback memory
            *(tbp++) = flag_c | (flag_d<<2) | (flag_i << 4);
        }
    }

    script = ckalloc((M+N)*sizeof(uchar));
    nedit = 0;
    // start with best node at grid point (M,N)
    row = M;
    col = N;
    if (C >= D && C >= I)
        node = FLAG_C;
    else if (D >= I)
        node = FLAG_D;
    else
        node = FLAG_I;
    while (row > 0 || col > 0) {
        /* Trace back from a C, D or I node at grid point (row, col).
        *  The type (C, D or I) tells the previous grid point on an
        *  optimal path.  The path's node type at that previous grid
        *  point is packed into the traceback byte.
        */
        if (row < 0 || col < 0) {
            fatal("Error generating edit script.");
        }
        st = tback_row[row][col];	// three flags in one byte
        script[nedit++] = node;
        if (node == FLAG_I) {
            col--;
            node = (st>>4);
        } else if (node == FLAG_D) {
            row--;
            node = (st>>2) & SELECT_CID;
        } else if (node == FLAG_C) {
            row--;
            col--;
            node = st & SELECT_CID;
        } else
            fatal("illegal node type in traceback");
    }

    *OM = m_new = nedit;
    al_new = (uchar **)ckalloc(m_new*sizeof(uchar *)) - 1;
    al_new[1] = (uchar *)ckalloc(m_new*(K+L)*sizeof(uchar));
    for (i = 2; i <= m_new; ++i)
        al_new[i] = al_new[i-1] + (K+L);

    i = j = m = 0;
    while (--nedit >= 0) {
        if ((k = script[nedit]) == FLAG_C)
            new_col(A[++i], K, B[++j], L, al_new[++m]);
        else if (k == FLAG_I)
            new_col(dashes, K, B[++j], L, al_new[++m]);
        else if (k == FLAG_D)
            new_col(A[++i], K, dashes, L, al_new[++m]);
        else
            fatalf("Illegal edit op: %d", k);
    }
    if (i != M || j != N || m  != m_new)
        fatalf("new_align: i=%d, j=%d, m=%d, M=%d, N=%d, M_new=%d\n",
               i, j, j, M, N, m_new);
    *OAL = al_new;

    free(tback_row);
    free(tback);
    free(dp);
    free(dashes);
    free(script);
}
