#ifndef MULTI_UTIL_H
#define MULTI_UTIL_H

/*
 *   multi_util.h  version 12
 *    -- procedures shared among components
 */
#include "maf.h"
#include "seq.h"

struct str_node {
    char* str;
    struct str_node* next;
};

/*
struct uAli {
  struct mafAli* ali;
  char* used;
  struct uAli* next;
  int start, end, size, index;
  char flipped, *topcontig;
};

struct pwuAli {
  struct uAli* pw;
  struct pwuAli* next;
  int cbeg, cend;
};

struct pwAlis {
  struct pwuAli ***bArrs;
  int *bCountArr, *countArr;
  char** topSpeciesArr, **botSpeciesArr, **fns;
  int pairK;
};

typedef struct name_list_ {
  char *name;
  struct name_list_ *next;
} NameList, *NameListPtr;

typedef struct tree_node_ {
  int id, type;
  NameListPtr names;      // leaf-species names
} TreeNode, *TreeNodePtr;

*/

void do_revcompl(char *s, int len);
void rev_comp(struct mafComp* c, int textSize);
void rc(struct mafAli *a);
//void rc_uAli(struct uAli *a);
void flip_comp(struct mafAli *a);
void flip_comps(struct mafAli *a);


int parse_fasta(char *line, char **ident, int *start, int *end, int *srcSize);
int parseHeader(char* fn, SEQ* sf, char* name, char* chr, int* start, int* end, char* strnd, int* size);
// inspect fasta header line to determine a sequence's name in maf entries
char *seq_ident(char *filename);


struct mafAli **mafArray(struct mafAli *a, int *nali);

void order_rows_species(struct mafAli **a, struct mafComp **location, char **species, int nspecies);

// keep alignment starting from position beg on reference
struct mafAli* keep_ali(struct mafAli* ali, int beg);

// print part of an alignment ali from beg to end which are
// on the first component, according to col number
int print_part_ali_col(struct mafAli* ali, int beg, int end, FILE* fp);
struct mafAli* make_part_ali_col(struct mafAli* ali, int cbeg, int cend);

//  print part of an alignment ali from beg to end which are
//  on the first component, not col number, starts at 0.
int print_part_ali(struct mafAli* ali, int beg, int end, FILE* fp);

// locate column of struct mafAli with a given position in sequence 1
int mafPos2Col(struct mafComp *c, int pos, int textSize);

// form a maf from a run of consecutive columns of a given maf
struct mafAli *mafSlice(struct mafAli *a, int start_col, int beyond_col);

/* ------------ start of code to manage the output list of alignments ----------
*  The output alignments are not guarateed to be found in increasing order of start
*  position in the reference sequence.  As they are generated, we buffer them into a
*  properly sorted list, and frequently flush the list.
*/

// flush output-list entries that start before a specified position
//struct mafAli *mafFlush(struct mafAli *output_list, int pos, FILE* fp);
struct mafAli *mafFlush(struct mafAli *output_list, int pos);

// insert an alignment into the output list, sorting by components->start
struct mafAli *mafInsert(struct mafAli *output_list, struct mafAli *a);

void insert_ali(struct mafAli** ali_head, struct mafAli* ali);
int compar_ali_top_start(const void* a, const void* b);
int change_neg_pos(struct mafComp* comp, int pos);

// currently break_ali cut ali with two rows
// bp is on top component maf position, ali2 is already allocated
void break_ali(struct mafAli* ali, int bp, struct mafAli* ali2);

struct mafAli* retrieve_first(struct mafAli** head);
void seperate_cp_wk(struct mafAli** cp_list, struct mafAli** wk_list, char* chr);

int colPos2Maf_after(struct mafComp* comp, int col);
int colPos2Maf_before(struct mafComp* comp, int col);

int test_ali(struct mafAli* ali);
void parseSrcName(char* srcName, char* name, char* src);
void parseSrcName2(struct mafComp*);

int overlap(int beg1, int end1, int beg2, int end2);
void print_ali_list(struct mafAli* root, FILE* fpw);

int mafPos2Col_v2(struct mafComp *c, int pos);

//void mark_uAli(struct uAli* A, int beg, int end);

//struct uAli** aliList2uAliArr(struct mafAli* subtreeAli, int treeCount);
//struct uAli* ali2uAli(struct mafAli* head);

//int sort_uAli_contigs(struct uAli** uAliArr, int arrSize, int structSize);
//struct pwAlis* setup_pwAlis(NameListPtr leftNames, NameListPtr rightNames, char* postfix);

int compute_ss(struct mafAli* ali, int cbeg, int cend);

int y_intercept(struct mafAli* pw, int x0);

#endif
