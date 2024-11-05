// lav2maf.c -- convert a blastz output file to a .maf file

#define VERSION 13

#include "util.h"
#include "maf.h"
#include "mz_scores.h"
#include "seq.h"
#include "multi_util.h"
//#include "contig.h"

static const char rcsid[] = "$Id: lav2maf.c 146 2009-01-21 22:02:30Z rico $";


#define MAX_ALIGN_LEN 10000000
#define MAX_CONTIG    1000000

char t1[MAX_ALIGN_LEN], t2[MAX_ALIGN_LEN];

int main(int argc, char **argv) {
    char buf[500], cmd[500], chr1[500], chr2[500],  name1[200], name2[200], src1[500], src2[500], strand1, strand2, cur_file1[500], cur_file2[500];
    int i, j, k, b1, b2, e1, e2, old_e1, old_e2, len1=-1, len2=-1, dir1,dir2, contig1, contig2;
    int srcSize1, srcSize2, tmp, start1, start2, end1, end2;
    struct mafAli B, *a;
    struct mafComp C1, C2, *c1, *c2;
    FILE *fp;
    SEQ *sf1=NULL, *sf2=NULL;
    unsigned char *s1=NULL, *s2=NULL;

    sprintf(cmd, "lav2maf.v%d", VERSION);
    argv0 = cmd;

    if (argc != 4)
        fatal(" -- convert blastz output to maf file.\n args: blastz.output seq-file1 seq-file2");

    a = &B;
    c1 = &C1;
    c2 = &C2;
    a->components = c1;
    a->score = 0.0;
    c1->next = c2;
    c2->next = NULL;
    c1->strand = c2->strand = '+';
    c1->text = t1;
    c2->text = t2;
    c1->paralog = c2->paralog = 's';

    init_scores70();
    init_KEEP_SEQ();
    //cur_contig1 = cur_contig2 = MAX_CONTIG;
    mafWriteStart(stdout, cmd);
    fp = ckopen(argv[1], "r");
    if (fgets(buf, 500, fp) == NULL || !same_string(buf, "#:lav\n"))
        fatalf("%s is not a blastz output file", argv[1]);
    while (fgets(buf, 500, fp) && !same_string(buf, "#:lav\n"))
        if (same_string(buf, "d {\n")) {
            printf("#\n");
            while (fgets(buf, 500, fp) && buf[0] != '}')
                printf("#%s", buf+1);
        }

    //contigsA = seq2contigs(argv[2], &numA);
    //contigsB = seq2contigs(argv[3], &numB);
    sf1 = seq_get_all(argv[2]);
    sf2 = seq_get_all(argv[3]);

    while (fgets(buf, 500, fp)) {
        if (same_string(buf, "s {\n")) {
            fgets(buf, 500, fp);
            if (sscanf(buf, "  \"%s %*d %*d %d %d", cur_file1, &dir1, &contig1) != 3)
                fatalf("Wrong format, cannot find seq file or orient or contig in %s", buf);
            contig1--;
            tmp = strlen(cur_file1);
            if (dir1 == 0)
                cur_file1[tmp-1] = '\0';
            else
                cur_file1[tmp-2]='\0';

            fgets(buf, 500, fp);
            if (sscanf(buf, "  \"%s %*d %*d %d %d", cur_file2, &dir2, &contig2) != 3)
                fatalf("Wrong format, Cannot find seq file or orient or contig in %s", buf);

            contig2--;
            tmp = strlen(cur_file2);
            if ( dir2 == 0)
                cur_file2[tmp-1]='\0';
            else
                cur_file2[tmp-2]='\0';

            if ( (dir1 == 0 && sf1->contigs[contig1].flipped=='y') || (dir1 == 1 && sf1->contigs[contig1].flipped=='n') ) {
                do_revcompl( (char*)(sf1->contigs[contig1].seq), sf1->contigs[contig1].slen);
                sf1->contigs[contig1].flipped = ( sf1->contigs[contig1].flipped == 'n' ? 'y' : 'n' );
            }
            if ( (dir2 == 0 && sf2->contigs[contig2].flipped=='y') || (dir2 == 1 && sf2->contigs[contig2].flipped=='n') ) {
                do_revcompl( (char*)(sf2->contigs[contig2].seq), sf2->contigs[contig2].slen);
                sf2->contigs[contig2].flipped = ( sf2->contigs[contig2].flipped == 'n' ? 'y' : 'n' );
            }
            s1 = sf1->contigs[contig1].seq - 1;
            len1 = sf1->contigs[contig1].slen;
            s2 = sf2->contigs[contig2].seq - 1;
            len2 = sf2->contigs[contig2].slen;

            /*
            if ( cur_contig1 > contig1) {
              if ( sf1 != NULL)
                seq_close(sf1);
              //sf1 = seq_open(cur_file1);
              sf1 = seq_open(argv[2]);
              seq_read(sf1);
              reversed1 = 0;
              cur_contig1 = 1;
            }
            while (cur_contig1 < contig1) {
              if ( seq_read(sf1) == 0)
                fatal("EOF in seq1");
              reversed1 = 0;
              ++cur_contig1;
            }
            if ( cur_contig2 > contig2) {
              if ( sf2 != NULL)
                seq_close(sf2);
              //sf2 = seq_open(cur_file2);
              sf2 = seq_open(argv[3]);
              seq_read(sf2);
              reversed2 = 0;
              cur_contig2 = 1;
            }
            while (cur_contig2 < contig2) {
              if ( seq_read(sf2) == 0)
                fatal("EOF in seq2");
              reversed2 = 0;
              ++cur_contig2;
            }			 
            s1 = SEQ_CHARS(sf1) - 1;  // to make seq start at pos 1
            len1 = SEQ_LEN(sf1);
            s2 = SEQ_CHARS(sf2) - 1;
            len2 = SEQ_LEN(sf2);

            if ( dir1==1 && reversed1==0) {
              sf1 = seq_revcomp_inplace(sf1);  // only the text is reversed
              reversed1 = 1;
            }
            if ( dir1==0 && reversed1==1) {
              sf1 = seq_revcomp_inplace(sf1);
              reversed1 = 0;
            }
              
            if ( dir2==1 && reversed2==0) {
              sf2 = seq_revcomp_inplace(sf2);  // the header does not know
              reversed2 = 1;
            }
            if ( dir2==0 && reversed2==1) {
              sf2 = seq_revcomp_inplace(sf2);
              reversed2 = 0;
            }
            */

        } else if (same_string(buf, "h {\n")) {
            fgets(buf, 500, fp);
            if ( parseHeader(argv[2], &(sf1->contigs[contig1]), name1, chr1, &start1, &end1, &strand1, &srcSize1) != 0)
                fatalf("Wrong format: %s", buf);
            fgets(buf, 500, fp);
            if ( parseHeader(argv[3], &(sf2->contigs[contig2]), name2, chr2, &start2, &end2, &strand2, &srcSize2) != 0)
                fatalf("Wrong format: %s", buf);

            //need to change strand according to dir1,2, since header didn't change

            c1->srcSize = srcSize1;
            c2->srcSize = srcSize2;
            if ( strcmp(name1, chr1) == 0)
                strcpy(src1, name1);
            else
                sprintf(src1, "%s.%s", name1, chr1);
            if ( strcmp(name2, chr2) == 0)
                strcpy(src2, name2);
            else
                sprintf(src2, "%s.%s", name2, chr2);
            c1->src = src1;
            c2->src = src2;
            start1--;
            start2--;
            end1--;
            end2--;
            if (strand1=='+' && dir1==0)
                c1->strand = '+';
            else if (strand1=='-' && dir1==1) {
                c1->strand = '+';
                start1 = srcSize1 - 1 - end1;
            } else if (strand1=='+' && dir1==1) {
                c1->strand = '-';
                start1 = srcSize1 - 1 - end1;
            } else // strand1=='-' && dir1==0
                c1->strand = '-';
            if (strand2=='+' && dir2==0)
                c2->strand = '+';
            else if (strand2=='-' && dir2==1) {
                c2->strand = '+';
                start2 = srcSize2 - 1 - end2;
            } else if (strand2=='+' && dir2==1) {
                c2->strand = '-';
                start2 = srcSize2 - 1 - end2;
            } else // strand2=='-' && dir2==0
                c2->strand = '-';
        } else if (same_string(buf, "a {\n")) {
            fgets(buf, 500, fp);	// score
            fgets(buf, 500, fp); // can't trust this b1, b2;
            fgets(buf, 500, fp); // ditto for e1, e2
            if (sscanf(buf, "  e %d %d", &e1, &e2) != 2)
                fatalf("cannot parse: %s", buf);
            if (e1 > len1)
                fatal("first sequence length is incorrect");
            if (e2 > len2)
                fatal("second sequence length is incorrect");
            old_e1 = old_e2 = -1;
            k = 0;
            while (fgets(buf, 500, fp) && buf[0] != '}') {
                if (sscanf(buf, "  l %d %d %d %d",
                           &b1, &b2, &e1, &e2) != 4)
                    fatalf("cannot read end-points: %s",
                           buf);
                if (old_e1 == -1) {
                    c1->start = start1 + b1 - 1;
                    c2->start = start2 + b2 - 1;
                } else {
                    if (k + MAX(e1+start1-1-c1->start,e2+start2-1-c2->start)
                            >= MAX_ALIGN_LEN)
                        fatal("laf2maf: alignment too long");
                    for (j = old_e1 + 1; j < b1; ++j, ++k) {
                        t1[k] = s1[j];
                        t2[k] = '-';
                    }
                    for (j = old_e2 + 1; j < b2; ++j, ++k) {
                        t1[k] = '-';
                        t2[k] = s2[j];
                    }
                }
                for (i = b1, j = b2; i <= e1; ++i, ++j, ++k) {
                    t1[k] = s1[i];
                    t2[k] = s2[j];
                }
                old_e1 = e1;
                old_e2 = e2;
            }
            c1->size = start1 + e1 - c1->start;
            c2->size = start2 + e2 - c2->start;
            t1[k] = t2[k] = '\0';
            a->textSize = k;
            a->score = mafScoreRange(a, 0, k);
            a->next = NULL;
            if ( a->components->start == a->components->next->start
                    && a->components->size  == a->components->next->size
                    && a->components->srcSize == a->components->next->srcSize
                    && strcmp(a->components->src, a->components->next->src)==0
                    && a->components->size == sf1->contigs[contig1].slen
                    && a->components->next->size == sf2->contigs[contig2].slen )
                continue;
            mafWrite(stdout, a);
        }
    }
    //if ( sf1 != NULL)
    //seq_close(sf1);
    //if ( sf2 != NULL)
    //seq_close(sf2);
    seq_close_all(sf1);
    seq_close_all(sf2);

    mafWriteEnd(stdout);
    return 0;
}
