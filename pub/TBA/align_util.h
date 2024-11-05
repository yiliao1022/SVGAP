#ifndef ALIGN_UTIL_H
#define ALIGN_UTIL_H

/*
 *  align_util.h version 1
 *    -- structures and precedures shared by aligners
 */

#include "maf.h"

struct position {
    int x, y;
};

// transform a single ali to ordered uAli with a given list of species
struct ordereduAli {
    char** contigs;
    int* begins, *ends;  // positively oriented
};

//----< wrapper of ali structure with mark array >------
struct uAli {
    struct mafAli* ali; // ali contains textSize field
    char* usedArr;      // 'u' unused; 'o' occupied
    struct uAli* next;
    char* sortContig;
    struct ordereduAli *oali;
    int index, start, end, cbeg, cend, rows, *indexes, incoming;
    char flipped;
};

//----< sorted array of uAlis >-----
struct sortuAlis {
    struct uAli** sortuAliArr;
    int *fronts, *ends;
    char* sortedSpeciesName;
    int sortuAliArrSize;
};

//----< collection of uAlis from an alignment file >-----
struct uAliFile {
    struct uAli** uAliArr;
    struct sortuAlis** sorted;
    char** speciesNames;
    int*   speciesAliCount;
    char* filename;
    int uAliCount, speciesCount;
};

//----< collection of pariwise alignment files >------
struct pwuAliFiles {
    struct uAliFile** pwuAliFileArrs;
    int pairK;  // number of pairwise alignment files
};


typedef struct name_list_ {
    char *name;
    struct name_list_ *next;
}
NameList, *NameListPtr;

typedef struct tree_node_ {
    int id, type;
    NameListPtr names;      /* leaf-species names */
}
TreeNode, *TreeNodePtr;

struct uAli* ali2uAli(struct mafAli* ali);
struct uAli** aliList2uAliArr(struct mafAli* subtreeAli, int treeCount);
void rc_uAli(struct uAli* uali);
void mark_uAli(struct uAli* A, int beg, int end, struct mafAli* nalilist);
void print_unused_ali(struct uAli* iAli, FILE* fpw);

void initialize_uAliFile(struct mafAli* head, struct uAliFile* tAli);
struct uAliFile* alilist2uAliFile(struct mafAli* head);
struct uAliFile* create_uAliFile(char* filename);
struct pwuAliFiles* create_pws(NameListPtr leftNames, NameListPtr rightNames, char* postfix);
NameListPtr formNameList(char* str);

char** compose_namelist(struct uAli* list, int* num);
//int connectionAgreement(struct mafAli* nlist, struct pwuAliFiles* pws, char** leftNames, int leftK, char** rightNames, int rightK);
int connectionAgreement2(struct mafAli* a2, struct mafAli* a3, int cbeg2, int cend2, int cbeg3, int cend3, struct pwuAliFiles* pws);
int mark_infered_pws(struct mafAli* nlist, struct pwuAliFiles* pws);
int sort_uAli_contigs(struct uAli** uAliArr, int arrSize);

struct sortuAlis* do_sortuAlis(struct uAli** uAliArr, int totalSize, char* name, int subSize);
struct uAli* Find_Exemplar(struct sortuAlis* sAli, char* src, int beg, int end, int* startIndex);

#endif
