#ifndef MZ_YAMA_H
#define MZ_YAMA_H

/* yama performs dynamic-programming global alignment on two alignments.
*  Gaps at the extreme ends (overhangs) are not charged a gap-open penalty.
*
* A is a K-by-M matrix representing a multiple alignment, with rows indexed
* by 0 <= k < K, and columns indexed by 1 <= i <= M.  Entries are in the list
* "acgtACGT-". A is stored by columns, so the entry in row k and column i is
* A[i][k].
*
* B is a L-by-N matrix, with structure analogous to A.
*
* LB[i] and RB[i] are left (lower) and right (upper) bounds in row i
* (0 <= i <= M) in the dynamic programming grid.
*
* AL_new is the (K+L)-by-M_new matrix resulting from a global alignment of
* A and B.  It can be freed by: "free(AL_new[1]); free(AL_new+1);".
*
*/

void yama(uchar **A, int K, int M, uchar **B, int L, int N, int *LB, int *RB, uchar ***OAL, int *OM);

#endif
