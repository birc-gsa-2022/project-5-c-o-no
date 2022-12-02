#include "approx.h"
#include "rotater.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "helper.h"
#define forAlphabet(code) for(int sym=1; sym<5; sym++) code(sym, rec, com); //This is weird, I love it!

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
    }
}


void recurseM(int sym, struct Recur rec, struct CommenRec* com) {
    rec.editString[rec.editIndex++] = 'M';
    if(sym != rec.patChar) ++rec.editAmount;
    --rec.patIndex;
    com->r->start = rec.rStart;
    com->r->end = rec.rEnd;
    //TODO make limit take ints instead of range
    limitRangeByChar(sym, com->r, com->C, com->O);
    rec.rStart = com->r->start;
    rec.rEnd = com->r->end;
    if(rec.rStart == rec.rEnd) return;
    recurseApprox(rec, com);
}

void recurseI(struct Recur rec, struct CommenRec* com) {
    rec.editString[rec.editIndex++] = 'I';
    --rec.patIndex;
    ++rec.editAmount;
    //Do not limit range, since insert whatever we need without looking at x
    recurseApprox(rec, com);
}

void recurseD(int sym, struct Recur rec, struct CommenRec* com) {
    rec.editString[rec.editIndex++] = 'D';
    ++rec.editAmount;
    com->r->start = rec.rStart;
    com->r->end = rec.rEnd;
    limitRangeByChar(sym, com->r, com->C, com->O);
    rec.rStart = com->r->start;
    rec.rEnd = com->r->end;
    if(rec.rStart == rec.rEnd) return;
    recurseApprox(rec, com);
}


void recurseApprox(struct Recur rec, struct CommenRec* com) {
    if(rec.editAmount > com->allowedEdits) {
        return;
    }

    if(rec.patIndex == -1) {
        //Report
        //Save in struct for easy debug, instead of printing
        struct ApproxMatch* match = malloc(sizeof *match); // Todo: get from pool?
        match->rStart = rec.rStart;
        match->rEnd = rec.rEnd;
        match->editStringLen = rec.editIndex;
        rec.editString[rec.editIndex] = '\0';
        match->editString = strdup(rec.editString);
        if(com->appCont->amount>= com->appCont->listSize) {
            com->appCont->listSize <<= 1;
            com->appCont->AMs = realloc(com->appCont->AMs, com->appCont->listSize*sizeof *(com->appCont->AMs));
        }
        com->appCont->AMs[(com->appCont->amount)++] = match;
        return;
    }

    if(rec.patIndex && ((com->allowedEdits)-(rec.editAmount) < com->D[rec.patIndex-1])) {
        return;
    }

    int patChar = com->pattern[rec.patIndex];
    rec.patChar = patChar;

    if(rec.editAmount < com->allowedEdits) {
        recurseI(rec, com);
    }

    if((rec.patIndex != com->m-1) && (rec.editAmount < com->allowedEdits)) { //Removes initial deletions
        forAlphabet(recurseD);
    }

    forAlphabet(recurseM);
}


struct ApproxMatchContainer* runApprox(int* pattern, int n, int m, int* D, int* C, int** O, int allowedEdits, char* editString, struct Range* r) {
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
    com->allowedEdits = allowedEdits;
    com->r = r;
    struct Recur* rec = malloc(sizeof *rec);
    rec->patIndex = m-1;
    rec->editString = editString;
    rec->editIndex = 0;
    rec->rStart = 0;
    rec->rEnd = n;
    rec->editAmount = 0;

    recurseApprox(*rec, com);
    return appCont;
}

void freeApproxMatch(struct ApproxMatch* approxMatch) {
    free(approxMatch->editString);
}

void freeApproxMatchContainer(struct ApproxMatchContainer* approxMatchContainer) {
    for(int i=0; i<approxMatchContainer->amount; i++) {
        freeApproxMatch(approxMatchContainer->AMs[i]);
        free(approxMatchContainer->AMs[i]);
    }
}


