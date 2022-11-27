#include "rotater.h"

#ifndef APPROX
#define APPROX

struct ApproxMatch {
    char* editString;
    int editStringLen;
    int rStart;
    int rEnd;
};

struct ApproxMatchContainer {
    struct ApproxMatch** AMs;
    int amount;
    int listSize;
};

struct Recur {
    int patIndex;
    int patChar;
    char* editString;
    int editIndex;
    int editAmount;
    struct Range* r;
};

struct CommenRec {
    int* pattern;
    int n;
    int m;
    int allowedEdits;
    int* D;
    int* C;
    int** O;
    struct ApproxMatchContainer* appCont;
};

void makeD(int* D, int* C, int** RO, const int* pattern, int n, int m, struct Range* r);
void recurseM(int sym, struct Recur* rec, struct CommenRec* com);
void recurseI(struct Recur* rec, struct CommenRec* com);
void recurseD(int sym, struct Recur* rec, struct CommenRec* com);
void recurseApprox(struct Recur* rec, struct CommenRec* com);
struct ApproxMatchContainer* runApprox(int* pattern, int n, int m, int* D, int* C, int** O, int allowedEdits, char* editString, struct Range* r);
#endif