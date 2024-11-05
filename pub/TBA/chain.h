#ifndef CHAIN_H
#define CHAIN_H
#include "maf.h"

typedef double big_score_t;

struct hsp_ali {
    int b1, b2, e1, e2, pct_id, score, used, chain, succ, contribute;
    big_score_t total_score;
    struct mafAli* ali;
    struct hsp_ali* next;
};

typedef struct {
    int *perm;
    int *rev_perm;
    struct hsp_ali *start, *query;
    int X, Y;		/* target point */
}
gp_t;

struct kdnode {
    int	bucket;
    int	cutval;
    int rightmost;
    big_score_t max_score;
    struct	kdnode *loson, *hison;
    int	lopt, hipt;
};

struct best_predecessor {
    int num;
    big_score_t contrib;
};

struct chain {
    char* contigA, *contigB;
    struct hsp_ali* segments;
    int b1, b2, e1, e2, score, keep_len, len1, len2;
    char orient;
    struct chain* next;
    char* used_array_1;
    char* used_array_2;
    char strand, marker;
};

struct chain* all_chains(int n);

struct chain* process_chain(struct mafAli* head);

struct chain* cutoff_chain(struct chain* chn, int beg, int end, struct chain* nchain, int row);

void print_chain(struct chain*);
void print_chains(struct chain*);

int compar_chain_score(const void* a, const void* b);

void insert_hali_head_chain(struct chain* chn, struct hsp_ali* hali);

void insert_head_chain(struct chain** head, struct chain* target);

int part_chain_score(struct chain* chn, int beg, int end, int* spb);
int top2botChainPos_after(struct chain* chn, int pos);
int top2botChainPos_before(struct chain* chn, int pos);
int bot2topChainPos_after(struct chain* chn, int pos);
int bot2topChainPos_before(struct chain* chn, int pos);

struct chain* mafs2chain(int num, char** strvec);
struct chain* duplicate_chain(struct chain* template);
int query_chain_length(struct chain* chn, int beg, int end);
struct chain** chainList2Arr(struct chain* chainhead, int* count);
void flip_chains(struct chain* head);
void print_unmarked_region(struct chain* chain_ptr, char* markArr, int refSeqStart);

#endif
