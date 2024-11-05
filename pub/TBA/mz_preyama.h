#ifndef MZ_PREYAMA
#define MZ_PREYAMA
#include "align_util.h"

/*
 * mz_preyama.h version 10
 *
 * The top sequence of two blocks a1 and a2 must share a same region
 * (same specie, overlapped regions). a1 and a2 are two blocks to be
 * aligned from beg to end, which are positions(not in respect to the
 * block) on the top reference sequence for each block.

 * pre_yama aligns columns of a1 and a2 from beg to end corresponding
 * to positions beg-end to the columns of a1 corresponding to positions
 * beg-end to the columns of a2, using the top reference sequence to
 * determine an approximation to the result.

 * The block is always topped by a reference, for being used to determine
 * the approximate alignment, but the referece is fixed in that block
 * or not, depending on the argument of reference:
 *     0:  neither is fixed
 *     1:  the first block is fixed
 */

struct mafAli *pre_yama(struct mafAli *a1, struct mafAli *a2, int beg, int end, int radius, int reference, FILE* fpw2);

struct mafAli* pre_yama2(struct mafAli* a1, struct mafAli* a2, struct mafAli* a3, int beg1, int end1, int begN, int endN, int radius, struct pwuAliFiles* pws);


#endif
