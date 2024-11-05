#include "util.h"
#include "maf.h"
#include "multi_util.h"

static char* ref;
static int unused;
static struct mafAli* unList;

static void init_maf_sort(char* refStr, int unVal) {
    ref = refStr;
    unused = unVal;
    unList = NULL;
}

struct mafAli* get_unused() {
    return unList;
}

int compar_start(const void* a, const void* b) {
    return (*((struct mafAli **)a))->components->start - (*((struct mafAli **)b))->components->start;
}

struct mafAli* ref_mvto_top(struct mafAli* root) {
    struct mafAli* ali, *ret=NULL;
    struct mafComp* comp, *prev;

    while ( root != NULL ) {
        ali = root;
        root = root->next;
        ali->next = NULL;
        for (comp=ali->components; comp!=NULL; comp=comp->next)
            if ( strcmp(comp->name, ref)==0 || strcmp(comp->src, ref)==0)
                break;
        if ( comp!=NULL) {
            if ( comp!=ali->components) {
                for (prev=ali->components; prev->next != comp; prev=prev->next)
                    ;
                prev->next = comp->next;
                comp->next = ali->components;
                ali->components = comp;
            }
            if ( ali->components->strand != '+')
                rc(ali);
            ali->next = ret;
            ret = ali;
        } else if ( unused == 1) {
            ali->next = unList;
            unList = ali;
        } else
            mafAliFree(&ali);
    }
    return ret;
}

// sort list according to first component starting position
// input list value changed
struct mafAli* maf_sort_top(struct mafAli* root) {
    struct mafAli** aliArr, *ali, *ret=NULL;
    int i, count;

    for (count=0, ali=root; ali != NULL; ali=ali->next)
        count++;

    aliArr = (struct mafAli**)malloc(count*sizeof(struct mafAli*));
    for (i=0, ali=root; i<count; i++, ali=ali->next)
        aliArr[i] = ali;
    for (i=0; i<count; i++)
        aliArr[i]->next = NULL;

    qsort((void*)aliArr, count, sizeof(struct mafAli*), compar_start);

    for (i=count-1; i>=0; i--) {
        aliArr[i]->next = ret;
        ret = aliArr[i];
        aliArr[i] = NULL;
    }
    return ret;
}

struct mafAli* maf_sort_list(struct mafAli* root, char* ref, int unusedVal) {

    init_maf_sort(ref, unusedVal);
    root = ref_mvto_top(root);
    root = maf_sort_top(root);
    return root;
}
