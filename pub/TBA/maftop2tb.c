// maftop2tb version 10
// maftop2tb.c -- convert a .maf file to a threaded blockset for top ref.
// the maf file must be projected to the reference sequence.
// the next two arguments are specie names, followed by a string or not:
// "contigs+" or "contigs-", no mapping files are needed

#define VERSION 10

#include "util.h"
#include "maf.h"
#include "mz_scores.h"
#include "multi_util.h"
#include "seq.h"
#include "maftop2tb.h"

static const char rcsid[] = "$Id: maftop2tb.c 142 2008-11-12 18:55:23Z rico $";


struct mafAli* getMafBetween(struct mafFile* mf, char* seqfile, FILE* fpw) {
    int i, j, start, end, size, nali, flag1, flag2;
    SEQ * sf;
    char src[500], name[200], strand, comp_name[500];
    unsigned char* s;
    struct mafAli **A, X, *a, *cp_list, *wk_list, *prev, *last, *b, *bkup=NULL, *it;
    struct mafComp C;

    X.next = NULL;
    X.score = 0.0;
    X.components = &C;
    C.next = NULL;
    C.strand = '+';
    C.paralog = 's';
    C.mafPosMap = NULL;

    sf = seq_open(seqfile);
    cp_list = mf->alignments;
    //mf->alignments = NULL;
    while (seq_read(sf)) {
        s = SEQ_CHARS(sf);
        if ( parseHeader(seqfile, sf, name, src, &start, &end, &strand, &size) != 0)
            fatalf("Wrong format: %s", sf->header);
        //end--;   // since the end itself is not included in the seq  <-- inclusive now
        if ( cp_list == NULL) { // no alignment
            strcpy(comp_name, name);
            strcat(comp_name, ".");
            strcat(comp_name, src);
            C.src = comp_name;
            C.name = name;
            C.contig = src;
            C.srcSize = size;
            C.start = start-1;
            X.textSize = C.size = end-start+1;
            C.text = (char*)malloc((C.size+1)*sizeof(char));
            for (j=0; j<C.size; j++)
                C.text[j] = s[j];
            C.text[C.size] = '\0';
            mafWrite(fpw, &X);
            free(C.text);
            break;
        }

        prev = a = cp_list;
        wk_list = last = NULL;
        while ( a != NULL) {
            if ( strcmp(name, a->components->name)==0 && strcmp(src, a->components->contig)==0 && a->components->start >= start-1 && a->components->start + a->components->size <= end) { // mv to wk_list
                if ( a == cp_list) {
                    cp_list = cp_list->next;
                    a->next = NULL;
                    b = a;
                    a = prev = cp_list;
                } else {
                    prev->next = a->next;
                    a->next = NULL;
                    b = a;
                    a = prev->next;
                }
                if ( wk_list == NULL)
                    wk_list = last = b;
                else {
                    last->next = b;
                    last = b;
                }
            } else {
                prev = a;
                a = a->next;
            }
        }

        if ( wk_list == NULL ) {// no align with this seq
            strcpy(comp_name, name);
            strcat(comp_name, ".");
            strcat(comp_name, src);
            C.src = comp_name;
            C.name = name;
            C.contig = src;
            C.srcSize = size;
            C.start = start-1;
            X.textSize = C.size = end-start+1;
            C.text = (char*)malloc((C.size+1)*sizeof(char));
            for (j=0; j<C.size; j++)
                C.text[j] = s[j];
            C.text[C.size] = '\0';
            mafWrite(fpw, &X);
            free(C.text);
            continue;
        }

        A = mafArray(wk_list, &nali);

        C.name = copy_string(name);
        C.contig = copy_string(src);
        strcpy(comp_name, name);
        strcat(comp_name, ".");
        strcat(comp_name, src);
        C.src = comp_name;

        C.srcSize = size;

        flag1 = start-1;
        for (i=0; i<nali; i++) {
            flag2 = A[i]->components->start;

            if ( flag2 > flag1 ) {
                C.start = flag1;
                C.paralog = 's';
                X.textSize = C.size = flag2-flag1;
                C.text = (char*)malloc((C.size+1)*sizeof(char));
                for (j=0; j<C.size; j++)
                    C.text[j] = s[flag1-start+1+j];
                C.text[C.size] = '\0';
                mafWrite(fpw, &X);
                free(C.text);
            }
            if ( A[i]->components->start + A[i]->components->size > flag1)
                flag1 = A[i]->components->start + A[i]->components->size;
        }
        if ( flag1 <= end-1 ) { // the last piece
            C.start = flag1;
            C.paralog = 's';
            X.textSize = C.size = end-flag1;
            C.text = (char*)malloc((C.size+1)*sizeof(char));
            for (j=0; j<C.size; j++)
                C.text[j] = s[flag1-start+1+j];
            C.text[C.size] = '\0';
            mafWrite(fpw, &X);
            free(C.text);
        }
        for (i=0; i<nali; i++)
            A[i] = NULL;
        free(A);
        if (bkup == NULL)
            bkup = wk_list;
        else {
            for (it=wk_list; it->next!=NULL; it=it->next)
                ;
            it->next = bkup;
            bkup = wk_list;
        }
    }
    seq_close(sf);
    if ( cp_list != NULL ) {
        if ( bkup == NULL )
            bkup = cp_list;
        else {
            for (it=cp_list; it->next!=NULL; it=it->next)
                ;
            it->next = bkup;
            bkup = cp_list;
        }
    }
    return bkup;
}

void getMafBetweenFileInput(char* inFile, char* seqfile, char* outFile) {
    struct mafFile* maf;
    FILE* fpw;

    maf = mafReadAll(inFile, 1);
    fpw = fopen(outFile, "w");

    getMafBetween(maf, seqfile, fpw);
    mafFileFree(&maf);
    fclose(fpw);
}
