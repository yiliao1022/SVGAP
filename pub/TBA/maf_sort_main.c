#include "util.h"
#include "maf.h"
#include "multi_util.h"
#include "maf_sort.h"

int main(int argc, char** argv) {
    struct mafFile* maf;
    struct mafAli *root;
    int unused=0;
    FILE* fpw=NULL;

    if ( argc < 3 )
        fatal("args: maf-file species-name [unused-ali-file]");

    if (argc==4) {
        unused = 1;
        fpw = fopen(argv[3], "w");
    }

    maf = mafReadAll(argv[1], 0);
    root = maf->alignments;
    maf->alignments = NULL;
    mafFileFree(&maf);
    root = maf_sort_list(root, argv[2], unused);
    mafWriteStart(stdout, "maf_project_simple");
    print_ali_list(root, stdout);

    if ( fpw!=NULL) {
        root = get_unused();
        print_ali_list(root, fpw);
        fclose(fpw);
    }
    mafWriteEnd(stdout);
    return 0;
}
