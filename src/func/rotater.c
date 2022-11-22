#include <malloc.h>
#include "rotater.h"

void makeOandC(const int* bwt, int n, int** O, int* C, int alphabetSize) {

    int* firstRow = calloc(--alphabetSize,sizeof *firstRow);
    if(bwt[0]) {
        firstRow[bwt[0]-1] = 1;
    }
    C[bwt[0]]++;
    O[0] = firstRow;
    int* prevRow;

    for(int i=1; i<n; i++) {
        int* row = malloc(alphabetSize*sizeof *row);
        prevRow = O[i-1];
        int curChar = bwt[i];
        for(int j=0; j<alphabetSize; j++) {
            row[j] = prevRow[j] + (curChar-1==j);
        }
        C[curChar]++;
        O[i] = row;
    }

    int accum = 0;
    for(int i=0; i<alphabetSize+1; i++) {
        int c = C[i];
        C[i] = accum;
        accum += c;
    }
}

int oLookUp(int** o, int searchChar, int i) {
    if(!searchChar || !i) return 0;
    return o[i-1][searchChar-1];
}

int jump(int bwtIndex, int jumpChar, int* c, int** o) {
    return c[jumpChar] + oLookUp(o, jumpChar, bwtIndex);
}

void rotateString(const int* string, int stringLen, int* c, int** o, int btwLen, struct Range* r) {
    int jumpChar;
    r->start = 0;
    r->end = btwLen;
    for(int i=stringLen-1; i>=0; i--) {
        if(r->start>=r->end) return;
        jumpChar = string[i];
        r->start = jump(r->start, jumpChar, c, o);
        r->end = jump(r->end, jumpChar, c, o);
    }

}
