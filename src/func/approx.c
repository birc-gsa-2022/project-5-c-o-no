#include <malloc.h>
#include <stdio.h>
#include "approx.h"
#include "rotater.h"
#include <string.h>
#define reset rec->r->start = rStart; rec->r->end = rEnd; rec->k=k; rec->patIndex=patIndex;
#define forAlphabet(code) for(int sym=1; sym<5; sym++) {code(sym, rec, com); reset;} //This is weird, I love it!

void makeD(int* D, int* C, int** RO, const int* pattern, int n, int m, struct Range* r) {
    int jumpChar;
    r->start = 0;
    r->end = n;
    for(int i=0; i<m; i++) {
        jumpChar = pattern[i];
        limitRangeByChar(jumpChar, r, C, RO);
        D[i] = i ? D[i-1] : 0;
        if(r->start >= r->end) {
            D[i]++;
            r->start = 0;
            r->end = n;
        }
        printf("D[%d]=%d\n",i,D[i]);
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
        struct ApproxMatch* match = malloc(sizeof *match);
        match->rStart = rec->r->start;
        match->rEnd = rec->r->end;
        match->editStringLen = rec->k;
        rec->editString[rec->k] = '\0';
        match->editString = strdup(rec->editString);
        if(com->appCont->amount>= com->appCont->listSize) {
            com->appCont->listSize <<= 1;
            com->appCont->AMs = realloc(com->appCont->AMs, com->appCont->listSize*sizeof *(com->appCont->AMs));
        }
        com->appCont->AMs[(com->appCont->amount)++] = match;
        return;
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
    struct ApproxMatchContainer* appCont = malloc(sizeof *appCont);
    appCont->listSize = 8;
    appCont->AMs = malloc(8*sizeof *appCont->AMs);
    appCont->amount = 0;

    struct CommenRec* com = malloc(sizeof *com);
    com->pattern = pattern;
    com->n = n;
    com->m = m;
    com->D = D;
    com->C = C;
    com->O = O;
    com->appCont = appCont;
    struct Recur* rec = malloc(sizeof *rec);
    rec->patIndex = patIndex;
    rec->k = k;
    rec->editString = editString;
    rec->editIndex = editIndex;
    rec->r = r;

    recurseApprox(rec, com);
}