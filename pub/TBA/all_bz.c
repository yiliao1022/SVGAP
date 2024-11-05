/* all_bz version 12

*  -- generate all blastz commands for pairs of specified sequences with
*     post-processing.
*
*  The command syntax is:
*	all_bz [-+] string-with-names [blastz-specfile]
*  where "string-with-names" contains names of sequence files (and perhaps
*  non-name characters), such as:
*	all_bz "(((human chimp) baboon) (mouse rat))"
*  A '+' tells all_bz to echo each blastz command to stdout before executing it.
*  A '-' tells all_bz to simply echo the commands (which can be captured in a
*  file, edited, then executed).  The optional last argument contains
*  command-line options for the blastz runs.  For instance, the "specfile":

# The is a sample file specifying command options for blastz
#define MAMMAL human baboon

* : *
	Y=3400
frog : *
	G=11
MAMMAL : zfish
	Q=HoxD50
MAMMAL : MAMMAL
	C=2 B=0
human : zfish
	O=300

*  says e.g. that (1) all blastz runs will specify Y=3400, (2) an alignment
*  involving frog will also specify G=11, (3) an alignment of human and zfish
*  will specify all of: Y=3400 Q=HoxD55 O=300, and (4) a human-vs-baboon
*  alignment will use: Y=3400 C=2 B=0
*
*  CURRENTLY WE ALWAYS SPECIFY: Y=3400 H=2000
*/

#include <stdio.h>
#include "util.h"

#define VERSION 15
#define TOAST_SUFFIX "toast.maf"
#define TOAST2_SUFFIX "toast2.maf"
#define NON_NAME " ()"

static const char rcsid[] = "$Id: all_bz.c 142 2008-11-12 18:55:23Z rico $";

#define BZ_CMD "blastzWrapper %s %s Y=9000 H=0 %s | lav2maf /dev/stdin %s %s | maf_sort /dev/stdin %s > %s.%s.orig.maf"
#define BZ_T2_CMD "blastzWrapper %s %s Y=9000 H=0 T=2 %s | lav2maf /dev/stdin %s %s | maf_sort /dev/stdin %s > %s.%s.orig.maf"
#define SIN_CMD_PRE "single_cov2 %s.%s.orig.maf "
#define SIN_CMD_POST " > %s.%s.sing.maf"
#define CLEAN_CMD "blastz_clean %s %s.%s.orig.maf | maf_sort /dev/stdin %s > %s.%s.clean.maf"
#define TOAST_CMD "toast %s %s %s %s %s %s.%s.clean.maf %s.%s.clean.maf %s.%s.clean.maf | maf_sort /dev/stdin %s > %s.%s.%s"
#define TOAST2_CMD "chain R=%s %s %s.%s.toast.maf | maf_sort /dev/stdin %s > %s.%s.%s"

static char buf[500], cmd[10000], min_cluster_chain_str[100], singleton_cmd[100], min_chain[100], inflation_cmd[100], *reference=NULL, close_alignment[500];
static char annotation_file[500];
static int execute=1;
static int verbose=0;
static int RUN_BZ=2;
static int RUN_TBA=1;
static int POST_PROC=1;
static int reduce=0;
static char *x[100], *y[100], *z[100];
static int nrules;


struct name_record {
    char *name;
    struct name_record *next;
};

int is_comment(char *s) {
    while (*s != '\0' &&isspace(*s))
        ++s;
    return (*s == '\0' || (*s == '#' && strncmp(s, "#define ", 8)));
}

char *spec_line(char *buf, int buf_size, FILE *fp) {
    int n = 0;

    while (fgets(buf+n, buf_size-n, fp) != NULL) {
        n = strlen(buf);
        if (n > 1 && buf[n-2] != '\\')
            return buf;
        --n;
        buf[n-1] = ' ';
        buf[n] = '\0';
    }
    return NULL;
}

// parse a blastz specfile; results go into globals x, y, z;
static void get_specs(char *filename) {
    FILE *fp = ckopen(filename, "r");
    char buf[5000], *name[100], *val[100], *s, *t, *u,
    *v;
    int nmacro = 0, i;

    while (spec_line(buf, 5000, fp) != NULL) {
        if (is_comment(buf))	// blank line
            continue;
        if (strncmp(buf, "#define ", 8) == 0) { // macro definition
            if (nmacro >= 100)
                fatal("too many macros");
            for (s = buf+8; isspace(*s); ++s)
                ;
            if ((t = strchr(s, ' ')) == NULL &&
                    (t = strchr(s, '\t')) == NULL)
                fatalf("end of name: %s", buf);
            for (u = t+1; *u != '\0' && isspace(*u); ++u)
                ;
            if (*u == '\0')
                fatalf("start of definition: %s", buf);
            if ((v = strchr(u, '\n')) != NULL)
                *v = '\0';
            *t = '\0';
            name[nmacro] = copy_string(s);
            val[nmacro] = copy_string(u);
            ++nmacro;
        } else { // must look like " name1 : name2 "
            for (s = buf; *s != '\0'; ++s)
                if (*s == ':' || isspace(*s))
                    break;
            if ((u = strchr(s, ':')) == NULL)
                fatalf("needs ':' in %s", buf);
            for (++u; *u != '\0' && isspace(*u); ++u)
                ;
            for (v = u; *v != '\0' && !isspace(*v); ++v)
                ;
            if (u == v)
                fatalf("confused by %s", buf);
            *s = *v = '\0';
            // name1 (buf)  might be a macro
            for (i = 0; i < nmacro; ++i)
                if (same_string(buf, name[i]))
                    break;
            x[nrules] = copy_string(i < nmacro ? val[i] : buf);
            // name2 (u) might be a macro
            for (i = 0; i < nmacro; ++i)
                if (same_string(u, name[i]))
                    break;
            y[nrules] = copy_string(i < nmacro ? val[i] : u);
            // next line gives command options
            while (fgets(buf, 500, fp) != NULL && is_comment(buf))
                ;
            // it must start with white-space
            if (!isspace(buf[0]))
                fatalf("missing space at start of %s", buf);
            // trim white-space from the front
            for (s = buf+1; isspace(*s); ++s)
                ;
            if ((t = strchr(s, '\n')) == NULL)
                fatalf("missing newline in %s", buf);
            // trim white-space from the rear
            while (t > s && isspace(t[-1]))
                --t;
            *t = '\0';
            // record the options in z
            z[nrules] = copy_string(s);
            ++nrules;
        }
    }
    fclose(fp);
}

// does name s match pattern t from the blastz-specifile?
static int match(char *s, char *t) {
    return (same_string(t, "*") || strstr(t, s));
}

// put the blastz command options for two species into buf
static void options(char *buf, char *name1, char *name2) {
    int i;

    buf[0] = '\0';
    for (i = 0; i < nrules; ++i)
        if ((match(name1, x[i]) && match(name2, y[i])) ||
                (match(name1, y[i]) && match(name2, x[i]))) {
            if (buf[0] != '\0')
                strcat(buf, " ");
            strcat(buf, z[i]);
        }
}

//------< run blastz for a pair of species >------
void command_bz(char* mname, char* nname, const int t2) {

    if ( RUN_BZ != 0) {
        options(buf, mname, nname);
        if (strstr(buf, "NOALIGN")) {
            if (verbose)
                printf("do not align %s and %s\n", mname, nname);
            return;
        }
        if ( t2 == 0)
            sprintf(cmd, BZ_CMD, mname, nname, buf, mname, nname, mname, mname, nname);
        else
            sprintf(cmd, BZ_T2_CMD, mname, nname, buf, mname, nname, mname, mname, nname);
        if ( verbose)
            printf("%s\n", cmd);
        if ( execute && system(cmd) != 0)
            fatal("all_bz BZ quits");
    }

    //sprintf(cmd, CLEAN_CMD, close_alignment, mname, nname, mname, mname, nname);
    //if ( verbose)
    //fprintf(stderr, "%s\n", cmd);
    //system(cmd);


#ifdef MAF2LAV
    sprintf(ifcmd, "maf2lav %s.%s.orig.maf %s %s > %s.%s.orig.bz", mname, nname, mname, nname, mname, nname);
    system(ifcmd);
#endif
}

//------< run post-process for a pair of species >---
void command_pp(char* mname, char* nname) {
    static char tmpStr[500];

    options(buf, mname, nname);
    if (strstr(buf, "NOALIGN")) {
        if (verbose)
            printf("do not align %s and %s\n", mname, nname);
        return;
    }

    if ( POST_PROC == 1 ) { // single_cov2
        sprintf(cmd, SIN_CMD_PRE, mname, nname);
        sprintf(tmpStr, SIN_CMD_POST, mname, nname);
        if ( reference != NULL && (strcmp(mname, reference)==0 || strcmp(nname, reference)==0) ) {
            strcat(cmd, "R=");
            strcat(cmd, reference);
        }
        strcat(cmd, tmpStr);
        if ( verbose)
            printf("%s\n", cmd);
        if ( execute && system(cmd) != 0)
            fatal("all_bz post-process quits");
    } else {
        sprintf(cmd, CLEAN_CMD, close_alignment, mname, nname, mname, mname, nname);
        if ( verbose)
            fprintf(stderr, "%s\n", cmd);
        system(cmd);

        sprintf(cmd, TOAST_CMD, singleton_cmd, min_cluster_chain_str, min_chain, inflation_cmd, annotation_file, mname, nname, mname, mname, nname, nname, mname, mname, nname, TOAST_SUFFIX);
        if ( verbose)
            fprintf(stderr, "%s\n", cmd);
        if ( execute && system(cmd) != 0 )
            fatal("all_bz post-process quits");

        if ( POST_PROC == 2 ) { // POST_RPOC == 2  // toast2
            if ( reference == NULL )
                fatal("To use A=2, reference must be specified\n");
            sprintf(cmd, TOAST2_CMD, reference, inflation_cmd, mname, nname, mname, mname, nname, TOAST2_SUFFIX);

            if ( verbose)
                printf("%s\n", cmd);
            if ( execute && system(cmd) != 0)
                fatal("all_bz post-process quits");
        }
    }

#ifdef MAF2LAV
    sprintf(ifcmd, "maf2lav %s.%s.toast.maf %s %s > %s.%s.toast.bz", mname, nname, mname, nname, mname, nname);
    system(ifcmd);
#endif
}


int main(int argc, char **argv) {
    char *s, *t, x, mainCmd[500], USAGE[10000];
    struct name_record *m, *n, *names;
    int TWO=2;

    sprintf(mainCmd, "all_bz.v%d", VERSION);
    argv0 = copy_string(mainCmd);

    strcpy(USAGE, "-- generate all blastz commands for pairs of specified sequences.\n");
    strcat(USAGE, "args: [-+] [b=?] [A=?] [F=reference] [T=annotation-file] [h=?] [q=?] [D=?] [f=?] species-guid-tree [blastz_specfile]\n");
    strcat(USAGE, "\t+(off) verbose\n\t-(off) output command only.\n");
    strcat(USAGE, "\tb(2) 0: run post-process only 1: run blastzWrapper only, transform to maf 2: run both\n");
    strcat(USAGE, "\tA(1) 0: toast 1: single_cov2 2: toast, following by chain and single cov on reference\n");
    strcat(USAGE, "\tF(null) null: single coverage is done for both species; reference: single coverage is done for reference only, effective in single_cov2\n");
    strcat(USAGE, "\tT(null): annotation file path and name, used for running toast and chaining procedure\n");
    strcat(USAGE, "\th(300) minimum chaining size, effective in toast\n");
    strcat(USAGE, "\tq(600) minimum cluster size, effective in toast\n");
    strcat(USAGE, "\tD(1) 0: run all_bz for roast 1: run all_bz for TBA.\n");
    strcat(USAGE, "\tc(500): parameter transfered to blastz_clean, alignments closer than c are subjected to be cleaned.\n");
    strcat(USAGE, "\tf(2) x% is used for determine in-paralogs, effective in toast.\n");

    //-----< process arguments >--------------
    if (argc > 1 && same_string(argv[1], "-")) {
        execute = 0;
        verbose = 1;
        --argc;
        ++argv;
    } else if (argc > 1 && same_string(argv[1], "+")) {
        verbose = 1;
        --argc;
        ++argv;
    }
    min_cluster_chain_str[0]=singleton_cmd[0]=min_chain[0]=inflation_cmd[0]=annotation_file[0]=close_alignment[0]=' ';
    min_cluster_chain_str[1]=singleton_cmd[1]=min_chain[1]=inflation_cmd[0]=annotation_file[1]=close_alignment[1]='\0';

    while ( argc > 1 && argv[1][1]=='=' && strchr("bAFThqscDf", (x=argv[1][0])) ) {
        switch (x) {
        case 'b':
            RUN_BZ = atoi(argv[1]+2);
            if ( RUN_BZ != 0 && RUN_BZ != 1 && RUN_BZ != 2 )
                fatalf("argument b can only be 0, 1, 2.\n%s", USAGE);
            break;
        case 'A':
            POST_PROC = atoi(argv[1]+2);
            if ( POST_PROC != 0 && POST_PROC != 1 && POST_PROC != 2)
                fatalf("argument A can only be 0, 1 or 2.\n%s", USAGE);
            break;
        case 'F':
            reference = copy_string(argv[1]+2);
            break;
        case 'T':
            sprintf(annotation_file, "A=%s", argv[1]+2);
            break;
        case 'f':
            strcpy(inflation_cmd, argv[1]);
            break;
        case 'h':
            strcpy(min_chain, argv[1]);
            break;
        case 'q':
            strcpy(min_cluster_chain_str, argv[1]);
            break;
        case 's':
            strcpy(singleton_cmd, argv[1]);
            reduce = atoi(argv[1]+2);
            if ( reduce != 0 && reduce != 1 )
                fatalf("argument t can only be 0, 1.\n%s", USAGE);
            break;
        case 'D':
            RUN_TBA=atoi(argv[1]+2);
            if ( RUN_TBA != 0 && RUN_TBA != 1 )
                fatalf("argument D can only be 0, 1.\n%s", USAGE);
            break;
        case 'c':
            strcpy(close_alignment, argv[1]);
            break;
        default:
            break;
        }
        --argc;
        ++argv;
    }

    if (argc == 3) {
        get_specs(argv[argc-1]);
        argc--;
    } else if (argc != 2)
        fatalf("%s", USAGE);

    m = names = NULL;
    for (s = argv[1]; *s != '\0'; s = t) {
        while (*s != '\0' && strchr(NON_NAME,*s))
            ++s;
        if (*s == '\0')
            break;
        for (t = s+1; *t != '\0' && !strchr(NON_NAME,*t); ++t)
            ;
        x = *t;	// save last character
        *t = '\0';
        n = ckalloc(sizeof(*n));
        n->name = copy_string(s);
        *t = x;	// replace last character
        n->next = NULL;
        if (m == NULL)
            names = n;
        else
            m->next = n;
        m = n;
    }
    //--------< end of processing arguments >-----

    //--------< run blastzWrapper >----------------------
    //if ( RUN_BZ != 0) {
    if ( RUN_TBA == 0 )
        command_bz(reference, reference, TWO);
    for (m = names; m != NULL; m = m->next) {
        if ( RUN_TBA == 0 ) {
            if ( reference == NULL )
                fatalf("reference must be specified for running roast  and its all_bz.\n%s", USAGE);
            if ( strcmp(m->name, reference)==0 )
                continue;
            command_bz(reference, m->name, 0);
            command_bz(m->name, m->name, 2);
        } else {
            n = ( POST_PROC == 0 ? m : m->next );
            for (; n != NULL; n = n->next) {
                if ( strcmp(m->name, n->name)==0 )
                    command_bz(m->name, n->name, 2);
                else
                    command_bz(m->name, n->name, 0);
            }
        }
    }
    //}

    if ( RUN_BZ == 1)  // run bz only
        return 0;

    //--------< post process >--------------------
    for (m = names; m != NULL; m = m->next) {
        if ( RUN_TBA == 0 ) {
            if ( reference == NULL )
                fatalf("reference must be specified for running roast and its all_bz.\n%s", USAGE);
            if ( strcmp(m->name, reference)==0 )
                continue;
            command_pp(reference, m->name);
        } else
            for (n = m->next; n != NULL; n = n->next) {
                command_pp(m->name, n->name);
#ifdef MAF2LAV
                sprintf(ifcmd, "laj %s.%s.orig.bz %s.%s.toast.bz &", m->name, n->name, m->name, n->name);
                system(ifcmd);
#endif
            }
    }


    return 0;
}
