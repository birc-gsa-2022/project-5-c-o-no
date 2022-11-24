#include <malloc.h>
#include <stdio.h>
#include "approx.h"
#include "rotater.h"
#define forAlphabet(code) for(int sym=1; sym<5; sym++)code(sym, rec, com); //This is weird, I love it!
#define resetRange rec->r->start = rStart; rec->r->end = rEnd;

void makeD(int* D, int* C, int** RO, int* pattern, int n, int m, struct Range* r) {
    int jumpChar;
    r->start = 0;
    r->end = n+1;
    for(int i=0; i<m; i++) {
        jumpChar = pattern[i];
        r->start = jump(r->start, jumpChar, C, RO);
        r->end = jump(r->end, jumpChar, C, RO);
        D[i] = i ? D[i-1] : 0;
        if(r->start>=r->end) {
            D[i]++;
            r->start = 0;
            r->end = n+1;
        }
    }
}

void recurseM(int sym, struct Recur* rec, struct CommenRec* com) {
    //TODO
    //editString[editIndex] = 'M';
}
void recurseI(int sym, struct Recur* rec, struct CommenRec* com) {
    //TODO
}
void recurseD(int sym, struct Recur* rec, struct CommenRec* com) {
    //TODO
}

//TODO make not recursive
void recurseApprox(struct Recur* rec, struct CommenRec* com) {
    if(com->D[rec->patIndex] < rec->k) return;
    if(rec->patIndex == com->m) {
        //TODO report
        //TODO will need to copy editString if it is not printed directly
    }
    int patChar = com->pattern[rec->patIndex];
    int rStart = rec->r->start;
    int rEnd = rec->r->end;

    forAlphabet(recurseI);
    resetRange;
    if(rec->patIndex) { //Removes initial deletions
        forAlphabet(recurseD);
    }
    resetRange;
    forAlphabet(recurseM);
    resetRange;
}

void runApprox(int* pattern, int patIndex, int n, int m, int* D, int* C, int** O, int k, char* editString, int editIndex, struct Range* r) {
    struct CommenRec* com = malloc(sizeof *com);
    com->pattern = pattern;
    com->n = n;
    com->m = m;
    com->D = D;
    com->C = C;
    com->O = O;
    struct Recur* rec = malloc(sizeof *rec);
    rec->patIndex = patIndex;
    rec->k = k;
    rec->editString = editString;
    rec->editIndex = editIndex;
    rec->r = r;

    recurseApprox(rec, com);

}