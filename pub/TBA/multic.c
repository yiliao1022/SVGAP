/*****************************************\
 *  multic - aligning two blocksets which
 *           are topped by reference, no
 *           single coverage requirement.
 *  ver 12.1
 *
 *  Language: C
 *  Platform: Linux Debian
 *  Application: independent executable, TBA
 *  Author:   Minmei Hou, Ph.D Candidate
 *            Penn State University
 *            mhou@cse.psu.edu
 *            April, 2005
 *
\*****************************************/
/*
  Neither of the blocksets need to be single-
  coveraged. But both blocksets need to be
  sorted according to top-row start position.
*/


#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "mz_scores.h"
#include "mz_preyama.h"
#include "mz_yama.h"

static const char rcsid[] = "$Id: multic.c 142 2008-11-12 18:55:23Z rico $";

extern int radius;
extern int MIN_OUTPUT_WID;
extern int row2;
extern int CONNECTION_THRESHOLD;

#define VERSION 12.1
//#define SHORTEST_BLOCKSIZE 30

int ALIGN_CATE = 0;
char* COLOR_ROW_NAME = NULL;

struct aliNode {
    struct mafAli* ali;
    char* unused;
    struct aliNode* next;
};

int any_identical_species(struct mafComp* A, struct mafComp* B) {
    struct mafComp* compA, *compB;

    for (compA=A; compA != NULL; compA=compA->next)
        for (compB=B; compB != NULL; compB=compB->next)
            if ( strcmp(compA->name, compB->name)==0)
                return 1;
    return 0;
}

void overlap_wrapper(struct aliNode* A, struct aliNode* B, int v) {
    struct mafComp* compA, *compB, *comp;
    struct mafAli* nali;
    int over_beg, over_end, a_end, b_end, A_cbeg, A_cend, B_cbeg, B_cend, i;

    compA = A->ali->components;
    compB = B->ali->components; // top reference

    a_end = compA->start + compA->size - 1;
    b_end = compB->start + compB->size - 1;
    over_beg = compA->start > compB->start ? compA->start : compB->start;
    over_end = a_end < b_end ? a_end : b_end;

    if ( over_beg > over_end )
        fatalf("there is no overlapping! %d %d\n", over_beg, over_end);
    nali = pre_yama(A->ali, B->ali, over_beg, over_end, radius, v, NULL);
    if ( nali==NULL)
        return;

    if ( A->ali->components->paralog == B->ali->components->paralog ) //both black/red
        nali->components->paralog = A->ali->components->paralog;

    else if ( A->ali->components->paralog == 'a' && B->ali->components->paralog == 'c') {
        for (comp=A->ali->components->next; comp!=NULL; comp=comp->next)
            if ( comp->paralog == 'c')
                break;
        if ( comp != NULL ) // there is non-reference red
            nali->components->paralog = 'a';
        else                // all-black
            nali->components->paralog = 'c';
    }
    else if ( A->ali->components->paralog == 'c' && B->ali->components->paralog == 'a') {
        for (comp=B->ali->components->next; comp!=NULL; comp=comp->next)
            if ( comp->paralog == 'c')
                break;
        if ( comp != NULL ) // there is non-reference red
            nali->components->paralog = 'a';
        else                // all-black
            nali->components->paralog = 'c';
    }

    over_beg = nali->components->start;
    over_end = nali->components->start + nali->components->size - 1;
    if ( nali->textSize >= MIN_OUTPUT_WID )
        mafWrite(stdout, nali);
    mafAliFree(&nali);

    // ---- mark the used region -------//
    if ( over_beg < compA->start || over_beg > compA->start+compA->size-1
            || over_end < compA->start || over_end > compA->start+compA->size-1 )
        fatalf("index out of boundary: %d-%d, %d-%d", over_beg, over_end, compA->start, compA->start+compA->size-1);
    A_cbeg = mafPos2Col(A->ali->components, over_beg, A->ali->textSize);
    A_cend = mafPos2Col(A->ali->components, over_end, A->ali->textSize);

    for (i=A_cbeg; i<=A_cend; i++)
        A->unused[i] = 'o';

    if ( over_beg < compB->start || over_beg > compB->start+compB->size-1
            || over_end < compB->start || over_end > compB->start+compB->size-1 )
        fatalf("index out of boundary: %d-%d, %d-%d", over_beg, over_end, compB->start, compB->start+compB->size-1);
    B_cbeg = mafPos2Col(B->ali->components, over_beg, B->ali->textSize);
    B_cend = mafPos2Col(B->ali->components, over_end, B->ali->textSize);

    for (i=B_cbeg; i<=B_cend; i++)
        B->unused[i] = 'o';
}

int multih(struct aliNode* A, struct aliNode* B, int v) {
    struct aliNode *a=A, *b=B, *bk;
    struct mafComp* comp, *compA, *compB;
    int a_beg, a_end, b_end, copyA, copyB;

    for (a=A, bk=B; a!=NULL; a=a->next) {

        if ( a->ali == NULL )
            continue;
        if ( ALIGN_CATE==2 && a->ali->components->paralog == 'a' )
            continue;
        for (comp = a->ali->components->next, copyA=0; comp!=NULL; comp=comp->next)
            if ( comp->paralog == 'c' )
                copyA++;
        if ( ALIGN_CATE != 0 && copyA > 1)
            fatalf("A: each block shall contain at most one copy paralog: %d", copyA);
        a_beg = a->ali->components->start;
        a_end = a_beg + a->ali->components->size - 1;
        for (; bk != NULL; bk=bk->next ) {
            if ( bk->ali == NULL )
                continue;
            comp=bk->ali->components;
            if ( comp->start + comp->size - 1 >= a_beg )
                break;
        }
        if ( bk==NULL)
            return 0;
        for (b=bk; b != NULL; b=b->next) {
            if ( b->ali == NULL )
                continue;
            if ( ALIGN_CATE==2 && b->ali->components->paralog == 'a')
                continue;
            for (comp=b->ali->components->next, copyB=0; comp!=NULL; comp=comp->next)
                if ( comp->paralog == 'c')
                    copyB++;
            if ( ALIGN_CATE != 0 && copyB > 1)
                fatalf("B: each block shall contain at most one copy paralog: %d", copyB);
            if ( ALIGN_CATE != 0 && copyA > 0 && copyB > 0) //two non-ref red rows
                continue;
            if ( b->ali->components->start > a_end )
                break;

            compA = a->ali->components;
            compB = b->ali->components->next;
            if ( v== 0)
                compA = compA->next;
            if ( any_identical_species(compA, compB) ) { // do not align bc species conflict
                if ( ALIGN_CATE != 0 && (copyA==0 && copyB==0)) { // all non-reference rows are black
                    b->ali->components->paralog = 'a';
                    if ( COLOR_ROW_NAME == NULL )
                        fatal("No COLOR_ROW_NAME specified!");
                    for (comp=b->ali->components->next; comp!=NULL; comp=comp->next)
                        if ( strcmp(comp->name, COLOR_ROW_NAME) == 0 )
                            break;
                    if ( comp == NULL )
                        fatal("COLOR_ROW_NAME specified wrong!");
                    comp->paralog = 'c';
                }
                continue;
            }

            b_end = b->ali->components->start + b->ali->components->size - 1;
            //if ( a->ali->components->size >= SHORTEST_BLOCKSIZE
            //   && b->ali->components->size >= SHORTEST_BLOCKSIZE
            //   && b_end - a->ali->components->start >= SHORTEST_BLOCKSIZE
            //   && a_end - b->ali->components->start >= SHORTEST_BLOCKSIZE )
            if (a->ali->components->start > b_end || b->ali->components->start > a_end )
                continue;
            overlap_wrapper(a, b, v);  // force alignment even for 1 base to avoid overlapping
        }
    }
    return 0;
}

struct aliNode* create_aliNode_list(struct mafAli* A) {
    struct aliNode **aliN_array, *aliN_root, *aliN;
    struct mafAli* ali;
    int i, count, j;

    for (count=0, ali=A; ali!=NULL; ali=ali->next)
        count++;
    if ( count == 0)
        return NULL;
    aliN_array = (struct aliNode**)malloc(count*sizeof(struct aliNode*));
    for (i=0, ali=A; i<count; i++, ali=ali->next) {
        aliN = (struct aliNode*)malloc(sizeof(struct aliNode));
        aliN_array[i] = aliN;
        aliN_array[i]->ali = ali;
        aliN_array[i]->unused = (char*)malloc(ali->textSize*sizeof(char));
        for (j=0; j<ali->textSize; j++)
            aliN_array[i]->unused[j] = 'u';
    }
    for (i=0; i<count; i++)
        aliN_array[i]->ali->next = NULL;
    for (i=0; i<count-1; i++)
        aliN_array[i]->next = aliN_array[i+1];
    aliN_array[count-1]->next = NULL;
    aliN_root = aliN_array[0];
    for (i=0; i<count; i++)
        aliN_array[i] = NULL;
    free(aliN_array);
    return aliN_root;
}

static void print_unused_ali_multic(struct aliNode* aliN, FILE* fpw) {
    struct mafAli* ali=aliN->ali, *nali;
    char* unused=aliN->unused;
    int i, j, size;

    if ( fpw==NULL || aliN->ali==NULL)
        return;

    size = ali->textSize;
    for (i=j=0; i<size && j<size;) {
        while ( i<size && unused[i]=='o') //used
            i++;
        if ( i>=size)
            break;
        j=i; // i locates on '0'
        while ( j<size && unused[j]=='u') //notused
            j++;
        j--; // next is >0 or j is at last pos
        nali = make_part_ali(ali, i, j);
        //if ( rm_repeated_block(&root, nali) != 1 )
        if ( nali!=NULL) {
            mafWrite(fpw, nali);
            mafAliFree(&nali);
        }
        i=j+1;
    }
}



// parses arguments
// args: [s=?] [R=?] [M=?] file1 file2 v [out1 out2] [nohead]
int main(int argc, char** argv) {
    char USAGE[10000], args[2000], cmd[200], ref_chr[200], x;
    struct mafFile* maf;
    struct mafAli *ali, *cp_list1, *cp_list2, *wk_list1, *wk_list2;
    struct aliNode *A, *B, *aliN, *nextN;
    int v, nohead=0, i;
    FILE* fpw[2];
    struct aliNode* arrN[2];

    sprintf(cmd, "multic.v%.1f", VERSION);
    argv0 = cmd;

    strcpy(USAGE, "args: [R=?] [M=?] [C=?] file1 file2 v? [out1 out2] [nohead] [all]\n");
    strcat(USAGE, "\tR(30) radius in dynamic programming.\n");
    strcat(USAGE, "\tM(1) minimum output width.\n");
    strcat(USAGE, "\tout1 out2(null) null: stdout; out1 out2: file names for collecting unused input.\n");
    strcat(USAGE, "\tnohead(null) null: output maf header; nohead: not to output maf header.\n");
    strcat(USAGE, "\tall(null) null: not to output single-row blocks; all: output all blocks.\n");

    if ( argc < 2 )
        fatalf("%s\n", USAGE);

    strcpy(args, cmd);
    strcat(args, " ");
    for (i=1; i<argc; i++) {
        strcat(args, argv[i]);
        strcat(args, " ");
    }

    // ---- process arguments
    while ( argc > 1 && strchr("sRMC", (x=argv[1][0])) && argv[1][1] == '=') {
        switch (x) {
        case 's':
            ALIGN_CATE = atoi(argv[1]+2);
            break;
        case 'R':
            radius = atoi(argv[1]+2);
            if ( radius < 0)
                fatal("radius cannot be negative");
            break;
        case 'M':
            MIN_OUTPUT_WID = atoi(argv[1]+2);
            if ( MIN_OUTPUT_WID < 0)
                fatal("MIN_OUTPUT_WID cannot be negative");
            break;
            //case 'c': COLOR_ROW_NAME = copy_string(argv[1]+2);
        case 'C':
            CONNECTION_THRESHOLD = atoi(argv[1]+2);
            if (CONNECTION_THRESHOLD < 0 || CONNECTION_THRESHOLD > 100)
                fatalf("%s\n", USAGE);
            break;
        default:
            fatalf("illegal flag %c", x);
        }
        ++argv;
        --argc;
    }

    if ( strcmp(argv[argc-1], "all") == 0 ) {
        row2=0;
        argc--;
    }

    if ( strcmp(argv[argc-1], "nohead")==0 ) {
        argc--;
        nohead = 1;
    }

    if ( argc != 4 && argc != 6)
        fatalf(" -- aligning two files of alignment blocks where top rows are always the reference, reference in both files can contain duplicats\n%s", USAGE);

    if ( argc == 6 ) {
        fpw[0] = fopen(argv[4], "w");
        fpw[1] = fopen(argv[5], "w");
    } else
        fpw[0] = fpw[1] = stdout;

    v = atoi(argv[3]);
    if ( v != 0 && v != 1)
        fatal("v can only be value of 0 or 1");

    // ---- end of processing args ------

    if ( nohead == 0) {
        mafWriteStart(stdout, "multih.c");
        printf("# %s\n", args);
    }

    init_scores70();
    maf = mafReadAll(argv[1], 1);
    cp_list1 = maf->alignments;
    maf->alignments = NULL;
    mafFileFree(&maf);
    maf = mafReadAll(argv[2], 1);
    cp_list2 = maf->alignments;
    maf->alignments = NULL;
    mafFileFree(&maf);

    while ( cp_list1 != NULL && cp_list2 != NULL ) {
        wk_list1 = wk_list2 = NULL;
        strcpy(ref_chr, cp_list1->components->src);

        seperate_cp_wk(&cp_list1, &wk_list1, ref_chr);
        seperate_cp_wk(&cp_list2, &wk_list2, ref_chr);

        A = create_aliNode_list(wk_list1);
        B = create_aliNode_list(wk_list2);

        if ( A!=NULL && B!=NULL)
            multih(A, B, v);

        arrN[0] = A;
        arrN[1] = B;

        for (i=0; i<2; i++)
            for (aliN=arrN[i]; aliN!=NULL; ) {
                nextN = aliN->next;
                aliN->next = NULL;
                if ( fpw[i] != NULL && aliN->ali->textSize >= MIN_OUTPUT_WID )
                    print_unused_ali_multic(aliN, fpw[i]);
                mafAliFree(&(aliN->ali));
                free(aliN->unused);
                free(aliN);
                aliN = nextN;
            }
    }

    if ( cp_list1 != NULL && fpw[0] != NULL)
        for (ali=cp_list1; ali!=NULL; ali=ali->next)
            if ( row2==0 || ali->components->next != NULL )
                mafWrite(fpw[0], ali);

    if ( cp_list2 != NULL && fpw[1] != NULL)
        for (ali=cp_list2; ali!=NULL; ali=ali->next)
            if ( row2==0 || ali->components->next != NULL )
                mafWrite(fpw[1], ali);

    for (i=0; i<2; i++)
        if (fpw[i] != NULL )
            fclose(fpw[i]);

    mafWriteEnd(stdout);
    return 0;
}
