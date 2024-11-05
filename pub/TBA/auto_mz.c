#include <unistd.h>
#include "util.h"
#include "multi_util.h"
#include "speciesTree.h"
#include "mz_scores.h"

static const char rcsid[] = "$Id: auto_mz.c,v 1.5 2006-02-06 20:58:53 mhou Exp $";

#define VERSION 3

#define MC "multic"
#define MZ "multiz"
#define MP "maf_project"
#define DEFAULT_MIN_WIDTH "1"

#define SING_SUFFIX ".sing.maf"
#define TOAST_SUFFIX ".toast.maf"
#define REDUCE_SUFFIX ".toast2.maf"

#define NAMESIZE 500
#define BIG_BUF_SIZE 50000

static char* REF = NULL;
static char SUFFIX[200];

char _MZ_U1[NAMESIZE], _MZ_U2[NAMESIZE], _MZ_O1[NAMESIZE], _MZ_O2[NAMESIZE], DESTINATION[500], project_cmd[500];
static char mz[500];
extern int execute;
extern int verbose;
extern int force;
extern char* PREFIX;
extern char* OPERAT;


int contain_species(NameListPtr nlist, char* species) {
    NameListPtr n;
    for (n=nlist; n!=NULL; n=n->next)
        if ( strcmp(n->name, species)==0 )
            return 1;
    return 0;
}

//int mz_merge(TreeNodePtr x, TreeNodePtr y, int id, int nbz, char** bz_file) {
int mz_merge(TreeNodePtr x, TreeNodePtr y, int id, int nbz, char** bz_file) {
    char _right_maf[200], _left_maf[200], _middle_maf[200];
    NameListPtr n1, n2;
    int left, right;


    nbz = nbz;
    bz_file = bz_file;

    if ( x->names == NULL || y->names == NULL )
        fatal("mz_merge:; emtpy sub-tree");

    sprintf(_middle_maf, "%smiddle.maf", PREFIX);
    sprintf(_left_maf, "%sleft.maf%d", PREFIX, id);
    sprintf(_right_maf, "%sright.maf%d", PREFIX, id);

    n1 = x->names;
    n2 = y->names;

    if ( n1->next == NULL && strcmp(n1->name, REF)==0 ) {
        force = 1;
        if ( n2->next == NULL )
            do_cmd("grep -v eof %s.%s%s >> %s%s%d", n1->name, n2->name, SUFFIX, PREFIX, OPERAT, id);
        else
            do_cmd("grep -v eof %s >> %s%s%d", _right_maf, PREFIX, OPERAT, id);
        do_cmd("rm -f %s %s", _right_maf, _left_maf);
        return 0;
    }
    if ( n2->next == NULL && strcmp(n2->name, REF)==0 ) {
        force = 1;
        if ( n1->next == NULL )
            do_cmd("grep -v eof %s.%s%s >> %s%s%d", n2->name, n1->name, SUFFIX, PREFIX, OPERAT, id);
        else
            do_cmd("grep -v eof %s >> %s%s%d", _left_maf, PREFIX, OPERAT, id);
        do_cmd("rm -f %s %s", _left_maf, _right_maf);
        return 0;
    }


    if ( n1->next == NULL )
        do_cmd("cp %s.%s%s %s", REF, n1->name, SUFFIX, _left_maf);
    if ( n2->next == NULL )
        do_cmd("cp %s.%s%s %s", REF, n2->name, SUFFIX, _right_maf);

    do_cmd(project_cmd, MP, _left_maf, REF, _MZ_O1, _MZ_U1);
    do_cmd(project_cmd, MP, _right_maf, REF, _MZ_O2, _MZ_U2);
    do_cmd("mv %s %s", _MZ_U1, _left_maf);
    do_cmd("mv %s %s", _MZ_U2, _right_maf);

    left = contain_species(n1, REF);
    right = contain_species(n2, REF);

    if ( left==0 && right==0 ) {
        do_cmd("%s %s %s 0 %s %s >> %s%s%d", mz, _left_maf, _right_maf, _MZ_U1, _MZ_U2, PREFIX, OPERAT, id);
        if ( n1->next == NULL && n2->next == NULL ) {
            force = 1;
            do_cmd("grep -v -h eof %s %s >> %s%s%d", _MZ_U1, _MZ_U2, PREFIX, OPERAT, id);
            return 0;
        }
    } else {
        if ( right == 1 ) {
            do_cmd("mv %s %s", _right_maf, _middle_maf);
            do_cmd("mv %s %s", _left_maf, _right_maf);
            do_cmd("mv %s %s", _middle_maf, _left_maf);
        }
        do_cmd("%s %s %s 1 %s %s >> %s%s%d", mz, _left_maf, _right_maf, _MZ_U1, _MZ_U2, PREFIX, OPERAT, id);
    }
    do_cmd("mv %s %s", _MZ_U1, _left_maf);
    do_cmd("mv %s %s", _MZ_U2, _right_maf);
    return 0;
}

int main(int argc, char **argv) {
    char name_list[10000], USAGE[10000], cmd[500], big_buf[BIG_BUF_SIZE], *bz_file[10000], buf[500], mzPar[500], *mzOpt, *p;
    int i, top = -1, id = 0, nbz, X, Cvalue=-1;
    TreeNode tree[1000];
    pid_t pid;
    char *tmpDir = "/tmp";        /* default for T= */

    sprintf(cmd, "roast.v%d", VERSION);
    argv0 = cmd;

    strcpy(USAGE, "args: [+-] [R=?] [M=?] [P=?] [T=?] [X=?] [C=?] E=reference-species species-guid-tree maf-source destination\n");
    strcat(USAGE, "\tR(30) dynamic programming radius.\n");
    strcat(USAGE, "\tM(1) minimum block length of output.\n");
    strcat(USAGE, "\tP(multiz) multiz: single coverage for reference row multic: no requirement on single coverage.\n");
    strcat(USAGE, "\tT(/tmp) specify alternate temp directory\n");
    strcat(USAGE, "\tX(0) utilize maf files with different suffix from differnt post processing.\n\t\t0: .sing.maf from single coverage pairwise alignment\n\t\t1: .toast.maf from full size toast\n\t\t2: .toast2.maf from reduced size toast\n");

    if (argc < 5)
        fatalf("roast -- reference guided multiple alignment.\n%s", USAGE);

    strcpy(DESTINATION, argv[argc-1]);



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

    while (argc > 1 && strchr("RMEPXCT",argv[1][0]) && argv[1][1] == '=') {
        switch ( argv[1][0] ) {
            //if ( strncmp(argv[1], "E=", 2)== 0)
        case 'E' :
            REF = copy_string(argv[1]+2);
            break;
            //else if ( argv[1][0]=='P' ) {
        case 'P':
            if (strstr(MC, argv[1]+2))
                mzOpt = MC;
            else if (strstr(MZ, argv[1]+2)==NULL)
                fatalf("the optional multiple aligner can be multiz or multic only.\n%s", USAGE);
            break;
            //}
            //else if ( strncmp(argv[1], "T=", 2)== 0)
        case 'T':
            tmpDir = argv[1]+2;
            break;
            //else if ( argv[1][0]=='X') {
        case 'X':
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
            break;
        case 'C':   // contintue to next case handling
            Cvalue = atoi(argv[1]+2);
            if ( Cvalue < 0 || Cvalue > 100 )
                fatalf("%s\n", USAGE);
        case 'R':
        case 'M':
            //else {
            strcat(mzPar, argv[1]);
            strcat(mzPar, " ");
            //}
            break;
        default:
            break;
        }
        --argc;
        ++argv;
    }


    if (REF == NULL)
        fatalf("fatal -- reference is not specified.\n%s", USAGE);

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
        for (i = 0; i < nbz; ++i) {
            bz_file[i] = argv[i+2];
        }
    }
    //----------------< end of processing arguments >-----------
    OPERAT = copy_string("MZ");
    pid = getpid();
    PREFIX = (char*)malloc(NAMESIZE*sizeof(char));
    sprintf(PREFIX, "%s/_%s_%d_", tmpDir, OPERAT, pid);
    sprintf(_MZ_U1, "%sU1", PREFIX);
    sprintf(_MZ_U2, "%sU2", PREFIX);
    sprintf(_MZ_O1, "%sO1", PREFIX);
    sprintf(_MZ_O2, "%sO2", PREFIX);

    init_scores70();

    do_cmd("rm -f %s", DESTINATION);
    do_cmd("echo \"##maf version=1 scoring=multiz.%d\" > %shead", VERSION, PREFIX);
    do_cmd("echo \"##maf version=1 scoring=%s.%d\" > %s", cmd, VERSION, DESTINATION);
    do_cmd("echo \"%s\" >> %s", big_buf, DESTINATION);

    top = parseSpeciesTree(argv[1], tree, nbz, bz_file, &id, buf, mz_merge);

    if ((p = strchr(name_list, '+')) != NULL)
        *p = ' ';
    do_cmd("%s %s%s%d %s %s > %s", MP, PREFIX, OPERAT, id-1, REF, _MZ_O1, _MZ_U1);
    force = 1;
    do_cmd("grep -v eof %s >> %s", _MZ_U1, DESTINATION);
    force = 0;
    do_cmd("rm %s*", PREFIX);
    if (top > 0)
        fatal("tree specification contains too many '('");
    if (top != 0 || tree[0].type != 0)
        fatal("tree specification is improper");
    do_cmd("echo \"##eof maf\" >> %s", DESTINATION);
    return 0;
}
