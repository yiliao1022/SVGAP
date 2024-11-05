#include "maf.h"

void init_maf_order(int num, char** strVec);
void free_maf_order();

// input ali is deallocated if return is NULL
struct mafAli* maf_order_ali(struct mafAli* ali);

// the ordering of alignment blocks from root is reversed in return list
struct mafAli* maf_order_list(struct mafAli* root);
