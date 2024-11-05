#include "maf.h"

struct mafAli* get_unused();
int compar_start(const void* a, const void* b);
struct mafAli* maf_sort_list(struct mafAli* root, char* ref, int unVal);
struct mafAli* maf_sort_top(struct mafAli* root);
