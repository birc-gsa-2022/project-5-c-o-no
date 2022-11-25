#include "rotater.h"

#ifndef APPROX
#define APPROX

struct ApproxMatch {
    char* editString;
    int editStringLen;
    int index;
};

struct ApproxMatchContainer {
    struct ApproxMatch** AMs;
    int amount;
};

struct Recur {
    int patIndex;
    int patChar;
    int k;
    char* editString;
    int editIndex;
    struct Range* r;
};

struct CommenRec {
    int* pattern;
    int n;
    int m;
    int* D;
    int* C;
    int** O;
};

void makeD(int* D, int* C, int** RO, int* pattern, int n, int m, struct Range* r);
void recurseM(int sym, struct Recur* rec, struct CommenRec* com);
void recurseI(struct Recur* rec, struct CommenRec* com);
void recurseD(int sym, struct Recur* rec, struct CommenRec* com);
void recurseApprox(struct Recur* rec, struct CommenRec* com);
void runApprox(int* pattern, int patIndex, int n, int m, int* D, int* C, int** O, int k, char* editString, int editIndex, struct Range* r);

#endif