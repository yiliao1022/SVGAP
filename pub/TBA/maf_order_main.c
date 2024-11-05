// maf_order_main.c -- order rows in a maf-file entries according to a given list
// args: maf-file species1 species2 .. [nohead] [all]
// By default, only blocks with more than a single row are output; if "all"
// is specified, then even single-row blocks are output.

#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "maf_order.h"

#define MAX_ENTRIES 100
#define VERSION 10

int main(int argc, char **argv) {
    char buf[500];
    struct mafFile* maf;
    struct mafAli *a;
    int head = 1, all=0; // default

    sprintf(buf, "maf_order.v%d", VERSION);
    argv0 = buf;

    if (argc < 3)
        fatal(" -- order rows according to a give list.\nargs: maf-file species1 species2 .. [nohead] [all]\n\t[nohead] if nohead is turned on, there is no maf header\n\t[all] if all is turned on, single-row blocks are also in ouput\n");

    if ( same_string(argv[argc-1], "all")) {
        all = 1;
        --argc;
    }

    if (same_string(argv[argc-1], "nohead")) {
        head = 0;
        --argc;
    }

    if ( head )
        mafWriteStart(stdout, "maf_order");

    init_maf_order(argc-2, argv+2);
    maf = mafOpen(argv[1], 1);
    while ( (a=mafNext(maf))!=NULL) {
        a = maf_order_ali(a);
        if ( a==NULL )
            continue;
        if ( all==1 || a->components->next != NULL )
            mafWrite(stdout, a);
        mafAliFree(&a);
    }
    free_maf_order();
    mafFileFree(&maf);
    mafWriteEnd(stdout);
    return 0;
}
