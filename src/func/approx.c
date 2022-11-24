#include <malloc.h>
#include <stdio.h>
#include "approx.h"
#include "rotater.h"

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

//TODO make not recursive
void recurseApprox(int patIndex, int k, char* editString, int editIndex, struct CommenRec* com) {
    if(com->D[patIndex] < k) return;
    if(patIndex == com->m) {
        //TODO report
        //TODO will need to copy editString if it is not printed directly
    }

    int patChar = com->pattern[patIndex];
    int rStart = com->r->start;
    int rEnd = com->r->end;

    //Insert
    editString[editIndex] = 'I';


    //TODO



    //Delete
    if(patIndex) { //Removes initial deletions
        editString[editIndex] = 'D';
        //TODO
    }

    //Substitution / mismatch
    editString[editIndex] = 'M';
    for(int sym=1; sym<5; sym++) {


    }
    //TODO
}

void runApprox(int* pattern, int patIndex, int n, int m, int* D, int* C, int** O, int k, char* editString, int editIndex, struct Range* r) {
    struct CommenRec* com = malloc(sizeof *com);
    com->pattern = pattern;
    com->n = n;
    com->m = m;
    com->D = D;
    com->C = C;
    com->O = O;
    com->r = r;

    recurseApprox(patIndex, k, editString, editIndex, com);

}