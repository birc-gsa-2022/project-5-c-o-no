#include <malloc.h>
#include <stdio.h>
#include "approx.h"
#include "rotater.h"
#define reset rec->r->start = rStart; rec->r->end = rEnd; rec->k=k; rec->patIndex=patIndex;
#define forAlphabet(code) for(int sym=1; sym<5; sym++) {code(sym, rec, com); reset;} //This is weird, I love it!

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
    rec->editString[rec->editIndex++] = 'M';
    if(sym != rec->patChar) ++rec->k;
    ++rec->patIndex;
    limitRangeByChar(sym, rec->r, com->C, com->O);
    if(rec->r->start == rec->r->end) return;
    recurseApprox(rec, com);
}

void recurseI(struct Recur* rec, struct CommenRec* com) {
    rec->editString[rec->editIndex++] = 'I';
    ++rec->k;
    ++rec->patIndex;
    //Do not limit range, since insert whatever we need without looking at x
    recurseApprox(rec, com);
}

void recurseD(int sym, struct Recur* rec, struct CommenRec* com) {
    rec->editString[rec->editIndex++] = 'D';
    ++rec->k;
    limitRangeByChar(sym, rec->r, com->C, com->O);
    if(rec->r->start == rec->r->end) return;
    recurseApprox(rec, com);
}

//TODO make not recursive
void recurseApprox(struct Recur* rec, struct CommenRec* com) {
    //TODO recurse with values instead of structs

    if(com->D[rec->patIndex] < rec->k) return;
    if(rec->patIndex == com->m) {
        //TODO report
        //TODO will need to copy editString if it is not printed directly
    }
    int patChar = com->pattern[rec->patIndex];
    rec->patChar = patChar;
    int rStart = rec->r->start;
    int rEnd = rec->r->end;
    int k = rec->k;
    int patIndex = rec->patIndex;

    recurseI(rec, com);
    reset;
    if(rec->patIndex) { //Removes initial deletions
        forAlphabet(recurseD);
    }
    reset;
    forAlphabet(recurseM);
    reset;
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