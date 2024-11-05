/****************************************************\
 * A small program to get sequence length to make   *
 * part of a standard header, useful for a sequence *
 * file not following                               *
 * ">name:chr:start-end:strand:srcSize" where start *
 * and end are 1 based.                             *
 *                                                  *
 * Result:                                          *
 * start position is assigned to be 1;              *
 * strand is assigned to be '+';                    *
 * end and srcSize are assigned to be sequence len. *
 *                                                  *
 * Any fasta sequence file shall start with '>'.    *
\****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "seq.h"

int main(int argc, char** argv) {
    SEQ* sf;

    if ( argc < 2)
        fatal("args: seq-file");
    sf = seq_open(argv[1]);
    while (seq_read(sf)) {
        printf("%s ==>\n", sf->header);
        printf("1-%d:+:%d\n", sf->slen, sf->slen);
    }
    seq_close(sf);
    return 0;
}


