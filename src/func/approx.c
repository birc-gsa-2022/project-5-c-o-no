#include "approx.h"
#include "rotater.h"
#include "rotater.h"

void makeD(int* D, int* C, int** RO, int* pattern, int n, int m, struct Range* r) {
    int jumpChar;
    r->start = 0;
    r->end = n+1;
    for(int i=0; i<m; i++) {
        jumpChar = pattern[i];
        r->start = jump(r->start, jumpChar, C, RO);
        r->end = jump(r->end, jumpChar, C, RO);
        D[i] = i ? D[i-1] : 0
        if(r->start>=r->end) {
            D[i]++;
            r->start = 0;
            r->end = n+1;
        }
    }
}

void runApprox(int* pattern, int patIndex, int n, int m, int* D, int* C, int** O, int k, char* editString, int editIndex, struct Range* r) {

    //TODO make not recursive
    void recurseApprox(int patIndex, int k, char* editString, int editIndex, struct Range* r) {
        if(D[patIndex] < k) return;
        if(patIndex == m) {
            //TODO report
            //TODO will need to copy editString if it is not printed directly
        }

        int patChar = pattern[patIndex];
        int rStart = r->start;
        int rEnd = r->end;

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
        //TODO
    }

    recurseApprox(patIndex, k, editString, editIndex, r);


}