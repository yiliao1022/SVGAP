#include <unistd.h>
#include "speciesTree.h"
#include "util.h"
#include "multi_util.h"

static const char rcsid[] = "$Id: speciesTree.c 142 2008-11-12 18:55:23Z rico $";


extern int force;
extern int execute;
extern int verbose;

extern char* PREFIX;
extern char* OPERAT;

#define BIG_BUF_SIZE 50000

// do_cmd -- print and possible execute a command
void do_cmd(const char *fmt, ...) {
    char buf[BIG_BUF_SIZE];
    va_list ap;

    // I restructured this to support the x86_64 architecture
    // without using va_copy from C99. -rico

    va_start(ap, fmt);
    (void)vsprintf(buf, fmt, ap);
    va_end(ap);

    if (verbose) {
        (void)printf("%s\n", buf);
        (void)fflush(stdout);
    }
    if (execute)
        if (system(buf) != 0)
            if (force == 0)
                fatalf("command '%s' failed", buf);
}

int  parseSpeciesTree(char* treeStr, TreeNode* tree, int nbz, char** bz_file, int* id, char* buf, int (*operation)(TreeNodePtr, TreeNodePtr, int, int, char **)) {
    char name_list[10000], *q, *p;
    NameListPtr n;
    int top=-1, i, j;

    for (q = treeStr; *q != '\0'; ++q) {
        if (*q == '(') {
            if (++top >= 1000)
                fatal("parse_tree: stack overflow");
            tree[top].type = '(';
        } else if (*q == ')') {
            if (top < 1 || tree[top].type != 0 || tree[top-1].type != '(') {
                q[1] = '\0';
                fatalf("parse error: %s", treeStr);
            }
            tree[top-1] = tree[top];
            --top;
        } else if (isalpha(*q)) {
            /* leaf species */
            if (++top >= 1000)
                fatal("parse_tree: stack overflow");
            p = buf;
            while (isalpha(*q) || isdigit(*q) || *q == '_' || *q == '.')
                *p++ = *q++;
            --q;
            *p = '\0';
            n = ckalloc(sizeof(NameList));
            n->name = copy_string(buf);
            n->next = NULL;
            tree[top].id = -1;
            tree[top].type = 0;
            tree[top].names = n;
        } else if (*q != ' ')
            fatalf("improper character in tree specification: %c", *q);
        if (top > 0 && tree[top-1].type == 0 && tree[top].type == 0) {
            // Process an internal node by merging its children.
            // One optimization is to identify nodes where both
            // children are leaves, to avoid a use of multiz.
            i = tree[top-1].id;
            j = tree[top].id;
            if (i >= 0)
                do_cmd("mv %s%s%d %sleft.maf%d", PREFIX, OPERAT, i, PREFIX, *id);
            if (j >= 0)
                do_cmd("mv %s%s%d %sright.maf%d", PREFIX, OPERAT, j, PREFIX, *id);
            do_cmd("cp %shead %s%s%d", PREFIX, PREFIX, OPERAT, *id);
            operation(&(tree[top-1]), &(tree[top]), *id, nbz, bz_file);
            if (i >= 0 || j >= 0) {
                if ( execute == 1) {
                    force = 1;
                    do_cmd("grep -v maf %sleft.maf%d >> %s%s%d", PREFIX, *id, PREFIX, OPERAT, *id);
                    do_cmd("grep -v maf %sright.maf%d >> %s%s%d", PREFIX, *id, PREFIX, OPERAT, *id);
                    force = 0;
                }
                name_list[0] = '\0';
                for (n = tree[top-1].names; n != NULL; n = n->next) {
                    strcat(name_list, " ");
                    strcat(name_list, n->name);
                }
                strcat(name_list, " +");
                for (n = tree[top].names; n != NULL; n = n->next) {
                    strcat(name_list, " ");
                    strcat(name_list, n->name);
                }
            }
            /* merge the lists of names */
            n = tree[top-1].names;
            if (n == NULL)
                fatal("empty list of names");
            while (n->next != NULL)
                n = n->next;
            n->next = tree[top].names;
            --top;
            tree[top].id = (*id)++;
        }
    }
    return top;
}

/*
NameListPtr formNameList(char* str) {
  char *ptr, tail='a';
  NameListPtr nname, ret=NULL;

  while ( *str == ' ')
    str++;

  while ( (ptr=strchr(str, ' ')) != NULL || (ptr=strchr(str, '\0')) != NULL ) {
    if ( *ptr == '\0' )
      tail = '\0';
    *ptr = '\0';
    nname = (NameListPtr)malloc(sizeof(NameListPtr));
    nname->name = copy_string(str);
    nname->next = ret;
    ret = nname;
    if ( tail == '\0' )
      break;
    *ptr = ' ';
    str = ptr+1;
    while (*str == ' ')
      str++;
  }
  return ret;
}
*/
