// tba version 12
// tba.c -- execute commands to run multiz and satellite programs, creating a
// threaded-block alignment

#define VERSION 12

// arg '-' means no execute, '+' means verbose

#include <unistd.h>
#include <sys/types.h>
#include "util.h"
#include "multi_util.h"
#include "maf.h"
#include "mz_scores.h"
#include "speciesTree.h"

static const char rcsid[] = "$Id: tba.c 142 2008-11-12 18:55:23Z rico $";

#define MZ "multiz"
#define MC "multic"
#define MP "maf_project"
#define P2 "pair2tb"
#define GCD "get_covered"
#define DEFAULT_MIN_WIDTH "1"

#define ORIG_SUFFIX ".orig.maf"
#define SING_SUFFIX ".sing.maf"
#define TOAST_SUFFIX ".toast.maf"
#define REDUCE_SUFFIX ".toast2.maf"
#define NAMESIZE 500
#define BIG_BUF_SIZE 50000

static char* REF = NULL;
static char SUFFIX[200];

char _tba_A[NAMESIZE], _tba_B[NAMESIZE], _tba_C[NAMESIZE], _tba_D[NAMESIZE], _tba_E[NAMESIZE], _tba_F[NAMESIZE], _tba_H[NAMESIZE], _tba_L[NAMESIZE], _tba_T[NAMESIZE], _tba_V[NAMESIZE], _tba_W[NAMESIZE], _tba_U[NAMESIZE], _tba_X[NAMESIZE], _tba_Y[NAMESIZE], _tba_REF[NAMESIZE];

char DESTINATION[500];
char project_cmd[500];

static char mz[500];	// multiz command perhaps with R= and M= options
extern int execute;
extern int verbose;
extern int force;
extern char* PREFIX;
extern char* OPERAT;

int get_rid_of_top(char* input_maf, char* output_maf) {
    struct mafFile* mf;
    struct mafAli* ali;
    struct mafComp* comp;
    FILE *fpw;

    if ( execute == 0)
        return 0;
    fpw = fopen(output_maf, "w");
    mf = mafOpen(input_maf, 1);
    while ( (ali=mafNext(mf))!=NULL) {
        comp=ali->components;
        if  (comp->next != NULL) {
            ali->components = comp->next;
            comp->next = NULL;
            mafCompFree(&comp); // need to recalculate the score
            ali->score = mafScoreRange(ali, 0, ali->textSize);
            mafWrite(fpw, ali);
        }
        mafAliFree(&ali);
    }
    mafFileFree(&mf);
    fclose(fpw);
    return 0;
}

#ifdef DEAD_CODE
// added 9/5/03
void iterate_porder(int id, char *name_list) {
    char cmd1[10000], cmd2[10000], cmd3[1000];
    int iterate = 1;

    sprintf(cmd1, "%s %s%s%d %s > %stemp",
            MD, PREFIX, OPERAT, id, name_list, PREFIX);
    sprintf(cmd2, "%s %stemp > %s%s%d", MD, PREFIX, PREFIX, OPERAT, id);
    sprintf(cmd3, "mv %stemp %s%s%d", PREFIX, PREFIX, OPERAT, id);

    while (iterate) {
        iterate = 0;
        if (verbose)
            printf("%s\n", cmd1);
        if (execute && system(cmd1) != 0)
            fatalf("command '%s' failed", cmd1);
        if (verbose)
            printf("%s\n", cmd2);
        if (execute && (iterate = system(cmd2)) != 0) {
            if (verbose)
                printf("%s\n", cmd3);
            system(cmd3);
        }
    }
}
#endif

static char *bz_cmd(char *x, char *y, int nbz, char **bz_file) {
    int i;
    char buf[500];

    sprintf(buf, "%s.%s%s", x, y, SUFFIX);
    for (i = 0; i < nbz; ++i)
        if (same_string(buf, bz_file[i]))
            return copy_string(buf);
    return NULL;
}

/* merge_trees - generate commands to merge the threaded-block alignments for
*  the two offspring of a node in the phylogenetic tree.
*/
static int tba_merge(TreeNodePtr x, TreeNodePtr y, int id, int nbz, char **bz_file) {
    NameListPtr n=NULL, n1, n2, prev;
    TreeNodePtr z;
    int found_cmd, single_left, single_right, single1, single2, swap=0;
    char  tmp_cmd[200], _right_maf[200], _left_maf[200], _middle_maf[200], *cmd, *mp = MP, *p2 = P2, *gcd = GCD;
    char _tba_head_F[NAMESIZE], _tba_head_E[NAMESIZE];

    if ( (n1=x->names) == NULL || (n2=y->names) == NULL )
        fatal("merge_tree: empty sub tree");

    if ( n1->next == NULL && n2->next == NULL) {
        if ((cmd = bz_cmd(n1->name, n2->name, nbz, bz_file)) == NULL)
            fatalf("no alignment found for %s and %s", n1->name,
                   n2->name);
        do_cmd("%s %s %s > %s", mp, cmd, n1->name, _tba_X);
        do_cmd("%s %s %s %s > %s%s%d", p2, _tba_X, n1->name, n2->name,
               PREFIX, OPERAT, id);
        return 0;
    }

    if ( REF != NULL )
        for (prev=n=x->names; n != NULL; prev=n,n=n->next )   // left-side
            if ( strcmp(n->name, REF)==0 ) {
                if ( n!=x->names ) {
                    prev->next = n->next;
                    n->next = x->names;
                    x->names = n;
                }
                break;
            }

    sprintf(_middle_maf, "%smiddle.maf", PREFIX);
    sprintf(_left_maf, "%sleft.maf%d", PREFIX, id);
    sprintf(_right_maf, "%sright.maf%d", PREFIX, id);


    if ( REF != NULL && n == NULL ) { // right-side
        for (prev=n=y->names; n != NULL; prev=n,n=n->next)
            if ( strcmp(n->name, REF)==0 ) {
                if ( n!=y->names ) {
                    prev->next = n->next;
                    n->next = y->names;
                    y->names = n;
                }
                break;
            }
        if ( n != NULL ) {
            z = x;
            x = y;
            y = z;
            do_cmd("mv %s %s", _right_maf, _middle_maf);
            do_cmd("mv %s %s", _left_maf, _right_maf);
            do_cmd("mv %s %s", _middle_maf, _left_maf);
            swap = 1;
        }
    }

    n1 = x->names;
    n2 = y->names;
    single_left = (n1->next == NULL);
    single_right = (n2->next == NULL);
    found_cmd = 0;
    single1=single2=0;
    for (n1 = x->names; n1 != NULL; n1 = n1->next)
        for (n2 = y->names; n2 != NULL; n2 = n2->next) {
            if ((cmd = bz_cmd(n1->name, n2->name, nbz, bz_file))!= NULL) {
                do_cmd("rm -f %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", _tba_Y, _tba_X, _tba_U, _tba_W, _tba_T, _tba_E, _tba_F, _tba_H, _tba_L, _tba_V, _tba_A, _tba_B, _tba_C, _tba_D, _tba_REF);
                do_cmd("touch %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", _tba_Y, _tba_X, _tba_U, _tba_W, _tba_T, _tba_E, _tba_F, _tba_H, _tba_L, _tba_V, _tba_A, _tba_B, _tba_C, _tba_D, _tba_REF);
                if ( !single_left) {
                    do_cmd(project_cmd, mp, _left_maf, n1->name, _tba_A, _tba_B);
                    if ( !single_right) {
                        do_cmd("%s %s %s 1 %s %s > %s", mz, _tba_B, cmd, _tba_Y, _tba_X, _tba_H);
                        do_cmd(project_cmd, mp, _tba_H, n2->name, _tba_U, _tba_B);
                        do_cmd(project_cmd, mp, _right_maf, n2->name, _tba_C, _tba_D);
                        if ( REF!=NULL && ((strcmp(REF, n1->name)==0 && n2->next!=NULL) || (strcmp(REF, n2->name)==0 && n1->next!=NULL)) )
                            do_cmd("%s %s %s 1 %s %s nohead > %s", mz, _tba_D, _tba_B, _tba_E, _tba_F, _tba_REF);
                        else
                            do_cmd("%s %s %s 1 %s %s nohead >> %s%s%d", mz, _tba_D, _tba_B, _tba_E, _tba_F, PREFIX, OPERAT, id);
                        if ( REF==NULL || strcmp(REF, n1->name) != 0 ) {
                            do_cmd("cat %shead %s > %shead_F", PREFIX, _tba_F, PREFIX);
                            sprintf(_tba_head_F, "%shead_F", PREFIX);
                            do_cmd(project_cmd, mp, _tba_head_F, n2->name, _tba_W, _tba_H);
                            get_rid_of_top(_tba_H, _tba_F);
                        }
                    } else { // single_right
                        do_cmd("%s %s %s > %s", mp, cmd, n1->name, _tba_X);
                        do_cmd("%s %s %s %s > %s", p2, _tba_X, n1->name, n2->name, _tba_D); // _D is pairwise tba including single-row
                        do_cmd(project_cmd, mp, _tba_D, n2->name, _tba_V, _tba_H);       // _H is _D projected onto seq2 includig single-row
                        if ( single2 == 0) {
                            do_cmd(project_cmd, mp, _tba_H, n1->name, _tba_C, _tba_D);
                            single2++;
                        } else {
                            do_cmd(project_cmd, mp, _right_maf, n2->name, _tba_V, _tba_D); // _D is projected from alignment file onto seq2
                            do_cmd("%s %s %s > %s", gcd, _tba_H, _tba_D, _tba_L);     // _H first file, _D second file, _L shall have parts from _H contained in _D, _D contains positions not aligned so far, If positions are not in _D, it means they're aligned already.
                            do_cmd(project_cmd, mp, _tba_L, n1->name, _tba_C, _tba_D);
                        }
                        do_cmd("%s %s %s 1 %s %s nohead >> %s%s%d", mz, _tba_B, _tba_D, _tba_F, _tba_E, PREFIX, OPERAT, id);
                        if ( REF==NULL || strcmp(REF, n1->name) != 0) {
                            do_cmd("cat %shead %s > %shead_E", PREFIX, _tba_E, PREFIX);
                            sprintf(_tba_head_E, "%shead_E", PREFIX);
                            do_cmd(project_cmd, mp, _tba_head_E, n1->name, _tba_T, _tba_H);
                            get_rid_of_top(_tba_H, _tba_E);
                        }
                    }
                }
                else if ( single_left) {
                    do_cmd("%s %s %s > %s", mp, cmd, n1->name, _tba_X);
                    do_cmd("%s %s %s %s > %s", p2, _tba_X, n1->name, n2->name, _tba_B);
                    do_cmd(project_cmd, mp, _tba_B, n1->name, _tba_V, _tba_H);
                    if ( single1 == 0) {
                        do_cmd(project_cmd, mp, _tba_H, n2->name, _tba_A, _tba_B);
                        single1++;
                    } else {
                        do_cmd(project_cmd, mp, _left_maf, n1->name, _tba_V, _tba_B);
                        do_cmd("%s %s %s > %s", gcd, _tba_H, _tba_B, _tba_L);
                        do_cmd(project_cmd, mp, _tba_L, n2->name, _tba_A, _tba_B);
                    }
                    do_cmd(project_cmd, mp, _right_maf, n2->name, _tba_C, _tba_D);
                    do_cmd("%s %s %s 1 %s %s nohead >> %s%s%d", mz, _tba_D, _tba_B, _tba_E, _tba_F, PREFIX, OPERAT, id);
                    if ( REF==NULL || strcmp(REF, n1->name) != 0) {
                        do_cmd("cat %shead %s > %shead_F", PREFIX, _tba_F, PREFIX);
                        sprintf(_tba_head_F, "%shead_F", PREFIX);
                        sprintf(tmp_cmd, "%stmp_F", PREFIX);
                        do_cmd(project_cmd, mp, _tba_head_F, n2->name, _tba_W, tmp_cmd);
                        get_rid_of_top(tmp_cmd, _tba_F);
                    }
                }
                force = 1;
                do_cmd("grep -v -h eof %shead %s %s %s %s %s > %s", PREFIX, _tba_A, _tba_Y, _tba_U, _tba_F, _tba_W, _left_maf);
                do_cmd("grep -v -h eof %shead %s %s %s > %s", PREFIX, _tba_C, _tba_E, _tba_T, _right_maf);
                force = 0;
                if ( REF!=NULL && !single_left && !single_right ) {
                    force = 0;
                    if ( strcmp(REF, n1->name)==0 && n2->next!=NULL )
                        do_cmd("grep -v eof %s >> %s", _tba_REF, _left_maf);
                    else if ( strcmp(REF, n2->name)==0 && n1->next!=NULL)
                        do_cmd("grep -v eof %s >> %s", _tba_REF, _right_maf);
                    force = 1;
                }
                found_cmd = 1;
            }
        }
    if (!found_cmd) {
        fprintf(stderr, "Warning! No alignments connect tree with leaves:\n");
        for (n = x->names; n; n = n->next)
            fprintf(stderr, "  %s", n->name);
        putc('\n', stderr);
        fprintf(stderr, "and tree with leaves:\n");
        for (n = y->names ; n; n = n->next)
            fprintf(stderr, "  %s", n->name);
        putc('\n', stderr);
        //exit(1); // in case the sequence of a species in guide tree does not exist
    }
    if ( swap == 1 ) {
        //z = x;  //swap by values
        //x = y;
        //y = z;
        do_cmd("mv %s %s", _right_maf, _middle_maf);
        do_cmd("mv %s %s", _left_maf, _right_maf);
        do_cmd("mv %s %s", _middle_maf, _left_maf);
    }
    return 0;
}

int main(int argc, char **argv) {
    char name_list[10000], USAGE[10000], cmd[500], big_buf[BIG_BUF_SIZE], *bz_file[10000], buf[500], mzPar[500], *mzOpt, *p;
    int i, top = -1, id = 0, nbz, X;
    TreeNode tree[1000];
    pid_t pid;

    sprintf(cmd, "tba.v%d", VERSION);
    argv0 = cmd;

    strcpy(USAGE, "args: [+-] [R=?] [M=?] [E=?] [P=?] [X=?] species-guid-tree maf-source destination\n");
    strcat(USAGE, "\tR(30) dynamic programming radius.\n");
    strcat(USAGE, "\tM(1) minimum block length of output.\n");
    strcat(USAGE, "\tE(null) null: no reference centric alignment, single coverage is guaranteed for every species; reference: refernece centric alignment, singe coverage is guaranteed for reference species.\n");
    strcat(USAGE, "\tP(null) null: run multiz; P=multic specifies to run multic.\n");
    strcat(USAGE, "\tX(0) utilize maf files with different suffix from differnt post processing.\n\t\t0: .sing.maf from single coverage pairwise alignment\n\t\t1: .toast.maf from full size toast\n\t\t2: .toast2.maf from reduced size toast\n");

    if (argc < 4)
        fatalf("TBA -- threaded block alignment.\n%s", USAGE);

    strcpy(DESTINATION, argv[argc-1]);

    pid = getpid();
    OPERAT = copy_string("tba");
    PREFIX = (char*)malloc(NAMESIZE*sizeof(char));
    sprintf(PREFIX, "/tmp/_%s_%d_", OPERAT, pid);
    sprintf(_tba_A, "%sA", PREFIX);
    sprintf(_tba_B, "%sB", PREFIX);
    sprintf(_tba_C, "%sC", PREFIX);
    sprintf(_tba_D, "%sD", PREFIX);
    sprintf(_tba_E, "%sE", PREFIX);
    sprintf(_tba_F, "%sF", PREFIX);
    sprintf(_tba_H, "%sH", PREFIX);
    sprintf(_tba_L, "%sL", PREFIX);
    sprintf(_tba_T, "%sT", PREFIX);
    sprintf(_tba_W, "%sW", PREFIX);
    sprintf(_tba_U, "%sU", PREFIX);
    sprintf(_tba_V, "%sV", PREFIX);
    sprintf(_tba_X, "%sX", PREFIX);
    sprintf(_tba_Y, "%sY", PREFIX);
    sprintf(_tba_REF, "%sREF", PREFIX);

    sprintf(big_buf, "# %s", cmd);
    for (i = 1; i < argc; ++i) {
        if (strlen(big_buf) + strlen(argv[1]) + 2 >= BIG_BUF_SIZE)
            fatal("overflow = big_buf");
        strcat(big_buf, " ");
        strcat(big_buf, argv[i]);
    }

    strcpy(project_cmd, "%s %s %s %s > %s");
    mzOpt = MZ;
    strcpy(mzPar, " ");
    verbose = 0;
    execute = 1;
    buf[0] = '\0';
    strcpy(SUFFIX, SING_SUFFIX);

    //----------------< process arguments >---------
    if (argc > 1 && same_string(argv[1], "-")) {
        strcat(buf, "-");
        execute = 0;
        verbose = 1;
        --argc;
        ++argv;
    } else if (argc > 1 && same_string(argv[1], "+")) {
        strcat(buf, "+");
        verbose = 1;
        --argc;
        ++argv;
    }

    while (argc > 1 && strchr("RMEPX",argv[1][0]) && argv[1][1] == '=') {
        if ( argv[1][0]=='E' )
            REF = copy_string(argv[1]+2);
        else if ( argv[1][0]=='P' ) {
            if (strstr(MC, argv[1]+2))
                mzOpt = MC;
            else if (strstr(MZ, argv[1]+2)==NULL)
                fatalf("the optional multiple aligner can be multiz or multic only.\n%s", USAGE);
        } else if ( argv[1][0]=='X') {
            X = atoi(argv[1]+2);
            switch (X) {
            case 1:
                strcpy(SUFFIX, TOAST_SUFFIX);
                break;
            case 2:
                strcpy(SUFFIX, REDUCE_SUFFIX);
                break;
            default:
                if ( X!=0 )
                    fatalf("Parameter X can only be 0, 1, 2, 3.\n%s", USAGE);
            }
        } else {
            strcat(mzPar, argv[1]);
            strcat(mzPar, " ");
        }
        --argc;
        ++argv;
    }

    if (strstr(mzPar, "M=") == NULL) {
        strcat(mzPar, "M=");
        strcat(mzPar, DEFAULT_MIN_WIDTH);
        strcat(mzPar, " ");
    }

    strcpy(mz, mzOpt);
    strcat(mz, mzPar);

    if (argc == 5 && same_string(argv[2], "-f")) {
        char buf2[500], *s;
        FILE *fp = ckopen(argv[3], "r");

        for (nbz = 0; fgets(buf2, 500, fp) != NULL; ++nbz) {
            if ((s = strchr(buf2, '\n')) != NULL)
                *s = '\0';
            bz_file[nbz] = copy_string(buf2);
        }
        fclose(fp);
    } else {
        nbz = argc - 3;
        for (i = 0; i < nbz; ++i)
            bz_file[i] = argv[i+2];
    }
    //----------------< end of processing arguments >-----------

    init_scores70();

    do_cmd("rm -f %s", DESTINATION);
    do_cmd("echo \"##maf version=%d scoring=multiz\" > %shead", VERSION, PREFIX);

    do_cmd("echo \"##maf version=%d scoring=%s\" > %s", VERSION, cmd, DESTINATION);
    do_cmd("echo \"%s\" >> %s", big_buf, DESTINATION);

    top = parseSpeciesTree(argv[1], tree, nbz, bz_file, &id, buf, tba_merge);

    if ((p = strchr(name_list, '+')) != NULL)
        *p = ' ';
    //	do_cmd("%s %stba%d %s nohead >> %s", mo, PREFIX, id-1, name_list,
    //	  DESTINATION);
    force = 1;
    do_cmd("grep -v eof %s%s%d >> %s", PREFIX, OPERAT, id-1, DESTINATION);
    force = 0;
    do_cmd("rm %s*", PREFIX);
    if (top > 0)
        fatal("tree specification contains too many '('");
    if (top != 0 || tree[0].type != 0)
        fatal("tree specification is improper");
    do_cmd("echo \"##eof maf\" >> %s", DESTINATION);
    return 0;
}
