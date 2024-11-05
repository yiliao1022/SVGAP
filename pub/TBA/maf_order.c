// maf_order.c -- order rows in a ali list or ali according to a given list
// version 1

#include "util.h"
#include "maf.h"
#include "multi_util.h"

#define MAX_ENTRIES 100
#define VERSION 11

static struct mafComp** location;
static char** species;
static int nspecies;

// num shall be larger than 0
void init_maf_order(int num, char** strVec) {
    nspecies = num;
    species = strVec;
    if ( num < 1)
        location = NULL;
    else
        location = ckalloc(nspecies*sizeof(struct mafComp*));
}

void free_maf_order() {
    free(location);
}

// input ali is deallocated if NULL returned
struct mafAli* maf_order_ali(struct mafAli* a) {
    struct mafComp *prev, *curr, *last;
    int i;

    for (i = 0; i < nspecies; ++i)
        location[i] = NULL;
    prev = curr = a->components;
    while (curr != NULL) {
        for (i = 0; i < nspecies; i ++) {
            if (same_string(curr->name, species[i]))
                break;
        }
        if (i == nspecies) {
            if (prev == curr) {
                a->components = curr->next;
                mafCompFree(&curr);
                prev = curr = a->components;
            } else {
                prev->next = curr->next;
                mafCompFree(&curr);
                curr = prev->next;
            }
            continue;
        }
        if (location[i] != NULL) {
            mafWrite(stderr, a);
            fatalf("species `%s' named more than once", curr->src);
        }
        location[i] = curr;
        prev = curr;
        curr = curr->next;
    }

    for (last = NULL, i = 0; i < nspecies; ++i)
        if (location[i] != NULL) {
            if (last == NULL)
                a->components = location[i];
            else
                last->next = location[i];
            last = location[i];
        }
    if (last == NULL) {
        mafAliFree(&a);
        return NULL;
    }
    last->next = NULL;

    a = mafColDashRm(a);
    if ( a != NULL && a->components->strand == '-')
        rc(a);
    return a;
}

// the ordering of alignment blocks in root is reversed in output
struct mafAli* maf_order_list(struct mafAli* root) {
    struct mafAli* ret, *ali;

    ret = NULL;
    while ( root != NULL ) {
        ali = root;
        root = root->next;
        ali->next = NULL;
        ali = maf_order_ali(ali);
        if (ali != NULL ) {
            ali->next = ret;
            ret = ali;
        }
    }
    return ret;
}
