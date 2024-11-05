#include "maf.h"
#include "util.h"
#include "multi_util.h"
#include "align_util.h"

static int MAX_SPECIES = 1000;
int CONNECTION_THRESHOLD = 50;
int SAME_CONNECTION = 30;

#ifdef MYDEBUG
//#define strcmp(x1, x2) mydebug_strcmp(__FILE__, __LINE__, x1, x2)
#define free(x) mydebug_free(__FILE__, __LINE__, x)
#endif

struct uAli* ali2uAli(struct mafAli* ali) {
    struct uAli* wrapper;
    int i, size;

    if ( ali == NULL )
        return NULL;

    wrapper = (struct uAli*)malloc(sizeof(struct uAli));
    wrapper->ali = ali;
    size = ali->textSize;
    wrapper->usedArr = (char*)malloc(size);
    for (i=0; i<size; i++)
        wrapper->usedArr[i] = 'u';

    wrapper->sortContig = NULL;
    wrapper->next = NULL;
    wrapper->indexes = NULL;
    wrapper->oali = NULL;
    wrapper->flipped = 'n';
    wrapper->index = wrapper->start = wrapper->end = wrapper->cbeg = wrapper->cend = -1;

    return wrapper;
}

struct uAli** aliList2uAliArr(struct mafAli* subtreeAli, int treeCount) {
    struct mafAli *ali, *next;
    char* unused;
    struct uAli** subtree;
    int i, textSize, j;

    if ( treeCount == 0 )
        return NULL;

    subtree = (struct uAli**)malloc(treeCount*sizeof(struct uAli*));
    subtree[0] = (struct uAli*)malloc(treeCount*sizeof(struct uAli));
    for (i=0, ali=subtreeAli; i<treeCount; i++, ali=next) {
        next = ali->next;
        ali->next = NULL;
        subtree[i] = subtree[0] + i;
        subtree[i]->ali = ali;
        subtree[i]->index = i;
        subtree[i]->flipped = 'n';
        subtree[i]->sortContig = NULL;
        subtree[i]->next = NULL;
        subtree[i]->oali = NULL;
        subtree[i]->indexes = NULL;
        textSize = ali->textSize;
        unused = subtree[i]->usedArr = (char*)malloc(textSize*sizeof(char));
        for (j=0; j<textSize; j++)
            unused[j] = 'u';
    }
    return subtree;
}

void rc_uAli(struct uAli *a) {
    int len, i, j;
    char used;

    rc(a->ali);
    len = a->ali->textSize;
    for (i=0, j=len-1; i<=j; i++, j--) {
        used = a->usedArr[i];
        a->usedArr[i] = a->usedArr[j];
        a->usedArr[j] = used;
    }
    a->flipped = ( a->flipped == 'n' ? 'y' : 'n' );
}

void mark_uAli(struct uAli* A, int Beg, int End, struct mafAli* nalilist) {
    struct mafComp* comp=A->ali->components, *ncomp;
    char* used=A->usedArr;
    struct mafAli* nali;
    int cbeg, cend, i, nbeg, nend, beg, end, nstart, ustart, uend, cond;

    for ( nali=nalilist; nali!=NULL; nali=nali->next) {
        cond = 0;
        for (ncomp=nali->components; ncomp!=NULL; ncomp=ncomp->next) {
            for (comp=A->ali->components; comp!=NULL; comp=comp->next)
                if ( strcmp(ncomp->name, comp->name) == 0 )
                    break;
            if ( comp != NULL ) {
                nstart = ncomp->start;
                nend = ncomp->start + ncomp->size - 1;
                ustart = comp->start;
                uend = comp->start + comp->size - 1;
                if ( ncomp->strand != comp->strand ) {
                    uend = comp->srcSize - comp->start - 1;
                    ustart = uend - (comp->size-1);
                }
                if ( strcmp(ncomp->src, comp->src)!=0 || nstart > uend || ustart > nend ) {
                    cond = 1;
                    break;
                }
            }
        }
        if ( cond == 1 )
            continue;

        comp = A->ali->components;
        for (ncomp=nali->components; ncomp!=NULL; ncomp=ncomp->next)
            if ( strcmp(ncomp->name, comp->name) == 0 )
                break;
        if ( ncomp == NULL )
            continue;

        if ( ncomp->strand == '+' )
            nbeg = ncomp->start;
        else
            nbeg = ncomp->srcSize - ncomp->start - ncomp->size - 1;

        nend   = nbeg + ncomp->size - 1;

        beg = ( Beg >= nbeg ? Beg : nbeg );
        end = ( End <= nend ? End : nend );

        nbeg = (beg > comp->start ? beg : comp->start);
        nend = (end < comp->start+comp->size-1 ? end : comp->start+comp->size-1);
        if ( nbeg > nend)
            continue;

        cbeg = mafPos2Col(comp, nbeg, A->ali->textSize);
        cend = mafPos2Col(comp, nend, A->ali->textSize);

        for (i=cbeg; i<=cend; i++)
            used[i] = 'o';
    }

}

void print_unused_ali(struct uAli* iAli, FILE* fpw) {
    struct mafAli* ali=iAli->ali, *nali;
    char* unused=iAli->usedArr;
    int i, j, size;

    if ( fpw==NULL || ali==NULL)
        return;

    size = ali->textSize;
    for (i=j=0; i<size && j<size;) {
        while ( i<size && unused[i]=='o') //used
            i++;
        if ( i>=size)
            break;
        j=i; // i locates on '0'
        while ( j<size && unused[j]=='u') //notused
            j++;
        j--; // next is >0 or j is at last pos
        nali = make_part_ali(ali, i, j);
        if ( nali!=NULL) {
            mafWrite(fpw, nali);
            mafAliFree(&nali);
        }
        i=j+1;
    }
}


int compar_uAli_start(const void* a, const void* b) {
    return (*((struct uAli**)a))->start - (*((struct uAli**)b))->start;
}

int compar_uAli_src(const void* a, const void* b) {
    return strcmp((*((struct uAli**)a))->sortContig, (*((struct uAli**)b))->sortContig);
}

int sort_uAli_contigs(struct uAli** uAliArr, int arrSize) {
    int i, front;
    char* prevStr;

    if ( arrSize == 0 )
        return 0;

    qsort(uAliArr, arrSize, sizeof(struct uAli*), compar_uAli_src);
    for (i=1, front=0, prevStr=uAliArr[0]->sortContig; i<arrSize; i++) {
        if ( strcmp(uAliArr[i]->sortContig, prevStr) != 0 ) { // end of a contig
            prevStr = uAliArr[i]->sortContig;
            qsort(uAliArr+front, i-front, sizeof(struct uAli*), compar_uAli_start);
            front = i;
        }
    }
    qsort(uAliArr+front, i-front, sizeof(struct uAli*), compar_uAli_start);  // the last segment
    return 0;
}


struct sortuAlis* do_sortuAlis(struct uAli** uAliArr, int totalSize, char* name, int subSize) {
    struct sortuAlis* suAlis;
    struct mafComp* comp;
    int i, k;

    suAlis = (struct sortuAlis*)malloc(sizeof(struct sortuAlis));
    if ( subSize > 0 )
        suAlis->sortuAliArr = (struct uAli**)malloc(subSize*sizeof(struct uAli*));
    else
        suAlis->sortuAliArr = NULL;
    suAlis->sortuAliArrSize = subSize;
    suAlis->sortedSpeciesName = copy_string(name);

    for (i=k=0; i<totalSize; i++) {
        for (comp = uAliArr[i]->ali->components; comp!=NULL; comp=comp->next)
            if ( strcmp(comp->name, name)==0 )
                break;
        if ( comp != NULL ) {
            suAlis->sortuAliArr[k++] = uAliArr[i];
            if ( uAliArr[i]->sortContig != NULL ) {
                free( uAliArr[i]->sortContig );
                uAliArr[i]->sortContig = NULL;
            }
            uAliArr[i]->sortContig = copy_string(comp->contig);
            if ( comp->strand == '+' )
                uAliArr[i]->start = comp->start;
            else
                uAliArr[i]->start = comp->srcSize - comp->start - comp->size;
            uAliArr[i]->end = uAliArr[i]->start + comp->size - 1;
        }
    }
    if ( k != subSize )
        fatalf("Inconsistent sizes k: %d and subSize: %d\n", k, subSize);
    sort_uAli_contigs(suAlis->sortuAliArr, subSize);

    suAlis->fronts = (int*)malloc(subSize*sizeof(int));
    suAlis->ends   = (int*)malloc(subSize*sizeof(int));
    for (k=0; k<subSize; k++) {
        suAlis->fronts[k] = suAlis->sortuAliArr[k]->start;
        suAlis->ends[k] = suAlis->sortuAliArr[k]->end;
    }
    return suAlis;
}

void initialize_uAliFile(struct mafAli* head, struct uAliFile* tAli) {
    char* names[MAX_SPECIES];
    struct mafComp* comp;
    struct mafAli* ali;
    int speciesCount=0, uAliCount=0, i;

    for (ali=head; ali != NULL; ali = ali->next) {
        for (comp=ali->components; comp != NULL; comp=comp->next) {
            for (i=0; i<speciesCount && i<1000; i++)
                if ( strcmp(comp->name, names[i]) == 0 )
                    break;
            if ( i == speciesCount )
                names[speciesCount++] = copy_string(comp->name);
        }
        uAliCount++;
    }

    tAli->speciesCount = speciesCount;
    if ( speciesCount > 0 ) {
        tAli->speciesAliCount = (int*)malloc(speciesCount*sizeof(int));
        tAli->speciesNames = (char**)malloc(speciesCount*sizeof(char*));
        tAli->sorted = (struct sortuAlis**)malloc(speciesCount*sizeof(struct sortuAlis*));
    } else {
        tAli->speciesAliCount = NULL;
        tAli->speciesNames = NULL;
        tAli->sorted = NULL;
    }
    tAli->uAliCount = uAliCount;

    for (i=0; i<speciesCount; i++) {
        (tAli->speciesNames)[i] = names[i];
        (tAli->speciesAliCount)[i] = 0;
    }

    for (ali=head; ali!=NULL; ali=ali->next)
        for (comp=ali->components; comp!=NULL; comp=comp->next) {
            for (i=0; i<speciesCount; i++)
                if ( strcmp(comp->name, names[i])==0 )
                    break;
            if ( i==speciesCount )
                fatalf("non-included species: %d\n", i);
            (tAli->speciesAliCount)[i]++;
        }

    tAli->uAliArr = aliList2uAliArr(head, uAliCount);

    for (i=0; i<speciesCount; i++)
        (tAli->sorted)[i] = do_sortuAlis(tAli->uAliArr, uAliCount, (tAli->speciesNames)[i], (tAli->speciesAliCount)[i]);
}

struct uAliFile* alilist2uAliFile(struct mafAli* head) {
    struct uAliFile* tAli;

    tAli = (struct uAliFile*)malloc(sizeof(struct uAliFile));
    initialize_uAliFile(head, tAli);
    tAli->filename = NULL;

    return tAli;
}

struct uAliFile* create_uAliFile(char* filename) {
    struct uAliFile* tAli;
    struct mafFile* maf;
    struct mafAli* head;

    maf = mafReadAll(filename, 0);
    head = maf->alignments;
    maf->alignments = NULL;
    mafFileFree(&maf);
    tAli = alilist2uAliFile(head);
    tAli->filename = copy_string(filename);
    return tAli;
}

struct pwuAliFiles* create_pws(NameListPtr leftNames, NameListPtr rightNames, char* postfix) {
    char filename[2000];
    struct pwuAliFiles* pws;
    NameListPtr leftName, rightName;
    int i, pairK, leftCount, rightCount;

    leftCount = rightCount = 0;
    for (leftName = leftNames; leftName != NULL; leftName = leftName->next)
        leftCount++;
    for (rightName = rightNames; rightName != NULL; rightName = rightName->next)
        rightCount++;
    pws = (struct pwuAliFiles*)malloc(sizeof(struct pwuAliFiles));
    pairK = pws->pairK = leftCount*rightCount;
    pws->pwuAliFileArrs = (struct uAliFile**)malloc(pairK*sizeof(struct uAliFile*));
    for (i=0, leftName = leftNames; leftName != NULL; leftName = leftName->next)
        for (rightName = rightNames; rightName != NULL; rightName = rightName->next) {
            sprintf(filename, "%s.%s.%s", leftName->name, rightName->name, postfix);
            pws->pwuAliFileArrs[i++] = create_uAliFile(filename);
        }
    if ( i != pairK )
        fprintf(stderr, "i and pairK not equal: %d %d\n", i, pairK);

    return pws;
}

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

int inSameList(char* nameA, char* nameB, char** nameList, int K) {
    int resultA=0, resultB=0, k;

    for (k=0; k<K; k++) {
        if ( strcmp(nameList[k], nameA)==0 ) {
            resultA = 1;
            continue;
        }
        if ( strcmp(nameList[k], nameB)==0 ) {
            resultB = 1;
            if (resultA == 1)
                return 1;
        }
    }
    if ( resultA == resultB )
        return 1;
    return 0;
}

char** compose_namelist(struct uAli* list, int* num) {
    struct uAli* uali;
    struct mafComp* comp;
    char** namelist;
    int index=0, i;

    namelist = (char**)malloc(MAX_SPECIES*sizeof(char*));
    for (uali=list; uali!=NULL; uali=uali->next)
        for (comp=uali->ali->components; comp!=NULL; comp=comp->next) {
            for (i=0; i<index; i++)
                if ( strcmp(comp->name, namelist[i])==0 )
                    break;
            if (i==index)
                namelist[index++] = copy_string(comp->name);
        }
    *num = index;
    return namelist;
}

/*
int connectionAgreement(struct mafAli* nlist, struct pwuAliFiles* pws, char** leftNames, int leftK, char** rightNames, int rightK) {
  struct mafAli *nali, *pw;
  struct mafComp *compA, *compB, *compa=NULL, *compb=NULL;
  char *topspecies, *botspecies;
  int i, j, k, pairK=pws->pairK, line1, line2, overbeg, overend, cbeg, cend, beg1, end1, beg2, end2;
  int *existConnections, expectConnection, existConnection=0;

  existConnections = (int*)malloc(pairK*sizeof(int));
  for (i=0; i<pairK; i++)
    existConnections[i] = 0;
  expectConnection = leftK*rightK;

  for (nali=nlist; nali!=NULL; nali=nali->next) {
    for (compA=nali->components; compA!=NULL; compA=compA->next) {
      for (compB=compA->next; compB!=NULL; compB=compB->next) {
	if ( inSameList(compA->name, compB->name, leftNames, leftK) )
	  continue;
        for (i=0; i<pairK; i++) {
	  if ( pws->pwuAliFileArrs[i]->uAliCount == 0 )
	    continue;
	  if ( pws->pwuAliFileArrs[i]->speciesCount < 2 )
	    fatal("pairwise alignment species number less than 2\n");
	  topspecies = pws->pwuAliFileArrs[i]->speciesNames[0];
	  botspecies = pws->pwuAliFileArrs[i]->speciesNames[1];
	  if ( ( strcmp(compA->name, topspecies) == 0 && strcmp(compB->name, botspecies) == 0 )
	       || ( strcmp(compA->name, botspecies) == 0 && strcmp(compB->name, topspecies) == 0 ) )
            break;
	}
        if ( i == pairK)
          continue;

        line1 = 0;
        if ( compA->strand != '+' ) {
          rev_comp(compA, nali->textSize);
          rev_comp(compB, nali->textSize);
          line1 = 1;
        }

	for (k=0; k<pws->pwuAliFileArrs[i]->speciesCount; k++)
	  if ( strcmp(compA->name, pws->pwuAliFileArrs[i]->sorted[k]->sortedSpeciesName) == 0 )
	    break;
	if ( k == pws->pwuAliFileArrs[i]->speciesCount )
	  fatalf("no sorted species: %s\n", compA->name);

        for (j=0; j<pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArrSize; j++) {
	  if ( pws->pwuAliFileArrs[i]->sorted[k]->fronts[j] > (compA->start + compA->size - 1 ) )
	    continue;
	  if ( pws->pwuAliFileArrs[i]->sorted[k]->ends[j] < compA->start )
	    continue;
	  pw = pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArr[j]->ali;
	  if ( strcmp(compA->name, pw->components->name)==0 ) {
	    compa = pw->components;
	    compb = compa->next;
	  }
	  else {
	    compb = pw->components;
	    compa = compb->next;
	  }
	  if ( strcmp(compa->contig, compA->contig) != 0 || strcmp(compb->contig, compB->contig)!=0 )
	    continue;
	  if ( compa->strand == '+')
	    if ( compb->strand != compB->strand )
	      continue;
	  line2 = 0;
	  if ( compa->strand == '-') {
	    if ( compb->strand == compB->strand )
	      continue;
	    rev_comp(compa, pw->textSize);
	    rev_comp(compb, pw->textSize);
	    line2 = 1;
	  }

          overbeg = (compA->start > compa->start ? compA->start : compa->start )
	    ;
          overend = (compA->start+compA->size < compa->start+compa->size ? compA->start+compA->size-1 : compa->start+compa->size-1 );
          cbeg = mafPos2Col(compA, overbeg, nali->textSize);
          cend = mafPos2Col(compA, overend, nali->textSize);
          beg1 = colPos2Maf_after(compB, cbeg);
          end1 = colPos2Maf_before(compB, cend);
          cbeg = mafPos2Col(compa, overbeg, pw->textSize);
          cend = mafPos2Col(compa, overend, pw->textSize);
          beg2 = colPos2Maf_after(compb, cbeg);
          end2 = colPos2Maf_before(compb, cend);
          if ( overlap(beg1, end1, beg2, end2) == 1 )
            existConnections[i] = 1;

	  if ( line2 == 1 ) {
	    rev_comp(compa, pw->textSize);
	    rev_comp(compb, pw->textSize);
	  }
        }

        if ( line1 == 1 ) {
          rev_comp(compA, nali->textSize);
          rev_comp(compB, nali->textSize);
        }
      }
    }
  }
  for (i=0; i<pairK; i++)
    if ( existConnections[i] == 1 )
      existConnection++;

  if ( existConnection*100/expectConnection >= CONNECTION_THRESHOLD )
    return 1;
  return 0;
}
*/

/* Given: intervals cbeg1-cend1 and cbegN-cendN, where
*    (1) positions cbeg1-cend1 are contained in the top row of alignment leftali, must be positive orient
*    (3) positions cbegN-cendN are contained in the top row of alignment rightali. might be negative orient
*/
int connectionAgreement2(struct mafAli* leftali, struct mafAli* rightali, int cbeg1, int cend1, int cbegN, int cendN, struct pwuAliFiles* pws) {
    struct mafAli* pw;
    struct mafComp *compA, *compB, *compa=NULL, *compb=NULL;
    char *topspecies, *botspecies;
    struct position a, b, c, d;
    int i, j, k, pairK=pws->pairK, marker1, marker2, overbeg, overend, cbeg, cend, beg1, end1, beg2, end2, leftK, rightK, tmp;
    int *existConnections, expectConnection, existConnection=0;
    int ab_mid_y, cd_mid_y, overmid;

    if ( leftali->components->strand == '-' )
        fatalf("left top component is not positive orientation: %s\n", leftali->components->name);

    for (leftK=0, compA=leftali->components; compA != NULL; compA=compA->next)
        leftK++;
    for (rightK=0, compB=rightali->components; compB != NULL; compB=compB->next)
        rightK++;

    existConnections = (int*)malloc(pairK*sizeof(int));
    for (i=0; i<pairK; i++)
        existConnections[i] = 0;
    expectConnection = leftK*rightK;

    for (compA=leftali->components; compA!=NULL; compA=compA->next) {
        marker1 = 0;
        if ( compA->strand == '-' ) {
            rev_comp(compA, leftali->textSize);
            for (compB=rightali->components; compB!=NULL; compB=compB->next)
                rev_comp(compB, rightali->textSize);
            tmp = cendN;
            cendN = rightali->textSize - cbegN - 1;
            cbegN = rightali->textSize - tmp - 1;
            tmp = cend1;
            cend1 = leftali->textSize - cbeg1 -1;
            cbeg1 = leftali->textSize - tmp - 1;
            marker1 = 1;
        }
        for (compB=rightali->components; compB!=NULL; compB=compB->next) {
            for (i=0; i<pairK; i++) {
                if ( pws->pwuAliFileArrs[i]->uAliCount == 0 )
                    continue;
                if ( pws->pwuAliFileArrs[i]->speciesCount < 2 )
                    fatal("pairwise alignment species number less than 2\n");
                topspecies = pws->pwuAliFileArrs[i]->speciesNames[0];
                botspecies = pws->pwuAliFileArrs[i]->speciesNames[1];
                if ( ( strcmp(compA->name, topspecies) == 0 && strcmp(compB->name, botspecies) == 0 )
                        || ( strcmp(compA->name, botspecies) == 0 && strcmp(compB->name, topspecies) == 0 ) )
                    break;
            }
            if ( i == pairK)
                continue;

            for (k=0; k<pws->pwuAliFileArrs[i]->speciesCount; k++)
                if ( strcmp(compA->name, pws->pwuAliFileArrs[i]->sorted[k]->sortedSpeciesName) == 0 )
                    break;
            if ( k == pws->pwuAliFileArrs[i]->speciesCount )
                fatalf("no sorted species: %s\n", compA->name);

            for (j=0; j<pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArrSize; j++) {
                if ( pws->pwuAliFileArrs[i]->sorted[k]->fronts[j] > (compA->start + compA->size - 1 ) )
                    continue;
                if ( pws->pwuAliFileArrs[i]->sorted[k]->ends[j] < compA->start )
                    continue;
                pw = pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArr[j]->ali;
                if ( strcmp(compA->name, pw->components->name)==0 ) {
                    compa = pw->components;
                    compb = compa->next;
                } else {
                    compb = pw->components;
                    compa = compb->next;
                }
                if ( strcmp(compa->contig, compA->contig) != 0 || strcmp(compb->contig, compB->contig)!=0 )
                    continue;
                if ( compa->strand == '+')
                    if ( compb->strand != compB->strand )
                        continue;
                marker2 = 0;
                if ( compa->strand == '-') {
                    if ( compb->strand == compB->strand )
                        continue;
                    rev_comp(compa, pw->textSize);
                    rev_comp(compb, pw->textSize);
                    marker2 = 1;
                }

                a.x = beg2 = colPos2Maf_after(compA, cbeg1);
                b.x = end2 = colPos2Maf_before(compA, cend1);

                overbeg = (beg2 > compa->start ? beg2 : compa->start );
                overend = (end2 < compa->start+compa->size-1 ? end2 : compa->start+compa->size-1);
                if ( overbeg > overend )
                    continue;

                a.y = beg1 = colPos2Maf_after(compB, cbegN);
                b.y = end1 = colPos2Maf_before(compB, cendN);

                cbeg = mafPos2Col(compa, overbeg, pw->textSize);
                cend = mafPos2Col(compa, overend, pw->textSize);
                beg2 = colPos2Maf_after(compb, cbeg);
                end2 = colPos2Maf_before(compb, cend);

                if ( overlap(beg1, end1, beg2, end2) == 1 ) {
                    c.x = compa->start;
                    c.y = compb->start;
                    d.x = compa->start + compa->size - 1;
                    d.y = compb->start + compb->size - 1;
                    overbeg = ( a.x > c.x ? a.x : c.x );
                    overend = ( b.x < d.x ? b.x : d.x );
                    overmid = (overbeg + overend)/2;
                    ab_mid_y = b.y - (b.x - overmid)*(b.y - a.y)/(double)(b.x - a.x);
                    cd_mid_y = d.y - (d.x - overmid)*(d.y - c.y)/(double)(d.x - c.x);
                    if ( ab_mid_y - cd_mid_y >= -SAME_CONNECTION && ab_mid_y - cd_mid_y <= SAME_CONNECTION )
                        existConnections[i] = 1;
                }
                if ( marker2 == 1 ) {
                    rev_comp(compa, pw->textSize);
                    rev_comp(compb, pw->textSize);
                }
            }
        } // for compB
        if ( marker1 == 1 ) {
            rev_comp(compA, leftali->textSize);
            for (compB=rightali->components; compB!=NULL; compB=compB->next)
                rev_comp(compB, rightali->textSize);
            tmp = cendN;
            cendN = rightali->textSize - cbegN - 1;
            cbegN = rightali->textSize - tmp - 1;
            tmp = cend1;
            cend1 = leftali->textSize - cbeg1 -1;
            cbeg1 = leftali->textSize - tmp - 1;
        } // if (marker1==1)
    } // for compA

    for (i=0; i<pairK; i++)
        if ( existConnections[i] == 1 )
            existConnection++;

    if ( existConnection*100/expectConnection >= CONNECTION_THRESHOLD )
        return 1;
    return 0;
}


int mark_infered_pws(struct mafAli* nlist, struct pwuAliFiles* pws) {
    struct uAli* upw;
    struct mafAli *nali, *pw;
    struct mafComp *compA, *compB, *compa=NULL, *compb=NULL;
    char *topspecies, *botspecies;
    int i, j, k, m, line1, line2, overbeg, overend, cbeg, cend, beg1, end1, beg2, end2;

    for (nali=nlist; nali!=NULL; nali=nali->next) {
        for (compA=nali->components; compA!=NULL; compA=compA->next) {
            for (compB=compA->next; compB!=NULL; compB=compB->next) {
                for (i=0; i<pws->pairK; i++) {
                    if ( pws->pwuAliFileArrs[i]->uAliCount == 0 )
                        continue;
                    if ( pws->pwuAliFileArrs[i]->speciesCount < 2 )
                        fatal("pairwise alignment species number less than 2\n");
                    topspecies = pws->pwuAliFileArrs[i]->speciesNames[0];
                    botspecies = pws->pwuAliFileArrs[i]->speciesNames[1];
                    if ( ( strcmp(compA->name, topspecies) == 0 && strcmp(compB->name, botspecies) == 0 )
                            || ( strcmp(compA->name, botspecies) == 0 && strcmp(compB->name, topspecies) == 0 ) )
                        break;
                }
                if ( i == pws->pairK)
                    continue;


                line1 = 0;
                if ( compA->strand != '+' ) {
                    rev_comp(compA, nali->textSize);
                    rev_comp(compB, nali->textSize);
                    line1 = 1;
                }

                for (k=0; k<pws->pwuAliFileArrs[i]->speciesCount; k++)
                    if ( strcmp(compA->name, pws->pwuAliFileArrs[i]->sorted[k]->sortedSpeciesName) == 0 )
                        break;
                if ( k == pws->pwuAliFileArrs[i]->speciesCount )
                    fatalf("no sorted species: %s\n", compA->name);

                for (j=0; j<pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArrSize; j++) {
                    if ( pws->pwuAliFileArrs[i]->sorted[k]->fronts[j] > (compA->start + compA->size - 1 ) )
                        continue;
                    if ( pws->pwuAliFileArrs[i]->sorted[k]->ends[j] < compA->start )
                        continue;
                    pw = pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArr[j]->ali;
                    if ( strcmp(compA->name, pw->components->name)==0 ) {
                        compa = pw->components;
                        compb = compa->next;
                    } else {
                        compb = pw->components;
                        compa = compb->next;
                    }
                    if ( strcmp(compa->contig, compA->contig) != 0 || strcmp(compb->contig, compB->contig)!=0 )
                        continue;
                    if ( compa->strand == '+')
                        if ( compb->strand != compB->strand )
                            continue;
                    line2 = 0;
                    if ( compa->strand == '-') {
                        if ( compb->strand == compB->strand )
                            continue;
                        rev_comp(compa, pw->textSize);
                        rev_comp(compb, pw->textSize);
                        line2 = 1;
                    }

                    overbeg = (compA->start > compa->start ? compA->start : compa->start )
                              ;
                    overend = (compA->start+compA->size < compa->start+compa->size ? compA->start+compA->size-1 : compa->start+compa->size-1 );
                    cbeg = mafPos2Col(compA, overbeg, nali->textSize);
                    cend = mafPos2Col(compA, overend, nali->textSize);
                    beg1 = colPos2Maf_after(compB, cbeg);
                    end1 = colPos2Maf_before(compB, cend);
                    cbeg = mafPos2Col(compa, overbeg, pw->textSize);
                    cend = mafPos2Col(compa, overend, pw->textSize);
                    beg2 = colPos2Maf_after(compb, cbeg);
                    end2 = colPos2Maf_before(compb, cend);
                    if ( overlap(beg1, end1, beg2, end2) == 1 ) {
                        upw = pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArr[j];
                        for (m=cbeg; m<=cend; m++)
                            upw->usedArr[m] = 'o';
                    }

                    if ( line2 == 1 ) {
                        rev_comp(compa, pw->textSize);
                        rev_comp(compb, pw->textSize);
                    }
                }

                if ( line1 == 1 ) {
                    rev_comp(compA, nali->textSize);
                    rev_comp(compB, nali->textSize);
                }
            }
        }
    }
    return 0;
}

int retrieve_exemplar(struct sortuAlis* sAli, int length, char* src, int index, int pos, int limit) {
    int* starts, *ends;
    struct mafComp* comp;
    int remain, max_remain, max_index, exist;

    starts = sAli->fronts;
    ends = sAli->ends;
    max_remain = remain = 0;
    max_index = -1;
    for (comp=sAli->sortuAliArr[index]->ali->components; comp!=NULL; comp=comp->next)
        if ( strcmp(comp->src, src) == 0 )
            break;
    if ( comp != NULL && starts[index] > limit )
        return -1;

    for (exist=0; index < length; index++) {
        for ( comp=sAli->sortuAliArr[index]->ali->components; comp!=NULL; comp=comp->next)
            if ( strcmp(comp->src, src) == 0 )
                break;
        if ( comp == NULL ) {
            if ( exist == 0)
                continue;
            else
                break;
        }
        exist = 1;
        if ( starts[index] > pos )
            break;
        if ( ends[index] < pos )
            continue;
        remain = ends[index] - pos + 1;
        if ( remain > max_remain ) {
            max_remain = remain;
            max_index = index;
        }
    }
    if ( max_index == -1 && index != length )
        return index;

    return max_index;
}



// aliArr is sorted according to the species row's starting point
// assuming sAli has all positions
struct uAli* Find_Exemplar(struct sortuAlis* sAli, char* src, int beg, int end, int* startIndex) {
    struct uAli *puAli, *head, *tail;
    int index=*startIndex, pos=beg, currEnd, max_index=-1, length;

    length = sAli->sortuAliArrSize;
    head = tail = NULL;
    while ( index < length ) {
        max_index = retrieve_exemplar(sAli, length, src, index, pos, end);
        if ( max_index == -1 )
            break;
        puAli = sAli->sortuAliArr[max_index];
        if ( head == NULL )
            head = tail = puAli;
        else {
            tail->next = puAli;
            tail = puAli;
        }
        currEnd = sAli->ends[max_index];
        if ( currEnd >= end )
            break;
        pos = currEnd + 1;
        index = max_index + 1;
    }
    if ( max_index > *startIndex )
        *startIndex = max_index;
    return head;
}

//#define TEST_ALIGN
#ifdef TEST_ALIGN
int main(int argc, char** argv) {
    struct uAliFile* treeAli;
    struct pwuAliFiles* pws;
    NameListPtr leftnames, rightnames;

    if (argc < 4)
        printf("args: maf-file left-names right-names postfix\n");

    treeAli = create_uAliFile(argv[1]);

    leftnames = formNameList(argv[2]);
    rightnames = formNameList(argv[3]);
    pws = create_pws(leftnames, rightnames, argv[4]);


    return 0;
}

#endif
