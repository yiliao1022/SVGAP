/* extract a range of positions from a nib file
 * reverse "to" and "from" to extract from reverse strand
 * For end of sequence, use negative number or huge number.
*/

#define VERSION 2

#include "seq.h"

char dna_compl[256] =
    "                                             -                  "
    " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
    "                                                                "
    "                                                                ";
/* .............................................-.................. */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */

int main(int argc, char **argv) {
    SEQ *sf;
    uchar *s;
    char cmd[100];
    int len, i, col, from, to;
    int start, end;

    sprintf(cmd, "dna_nib.v%d", VERSION);
    argv0 = cmd;
    if (argc != 5)
        fatal("args = nib-file from to fast-header");

    sf = seq_get(argv[1]);
    s = SEQ_CHARS(sf);
    len = SEQ_LEN(sf);

    from = atoi(argv[2]);
    if (from < 0 || from >= len)
        from = len-1;
    to = atoi(argv[3]);
    if (to < 0 || to >= len)
        to = len-1;
    printf(">%s:%d-%d:+:%d\n", argv[4], from, to, len);

    col = 0;
    if (from < to) {
        start = from;
        end = to;
    } else {
        /* reverse strand */
        start = to;
        end = from;
    }
    for (i = start; i <= end; ++i) {
        putchar(from < to ? s[i] : dna_compl[(uchar)s[start-i+end]]);
        if (++col == 50) {
            putchar('\n');
            col = 0;
        }
    }
    if (col != 0)
        putchar('\n');
    seq_close(sf);

    return 0;
}
