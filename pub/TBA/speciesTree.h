#ifndef SPECIESTREE_H
#define SPECIESTREE_H
#include "align_util.h"

void do_cmd(const char *fmt, ...);

/* parse treeStr, tree is allocated, assuming enough size */
int  parseSpeciesTree(char* treeStr, TreeNode* tree, int nbz, char** bz_file, int* id, char* buf, int (*operation)(TreeNodePtr, TreeNodePtr, int, int, char**));

//NameListPtr formNameList(char* str);

#endif
