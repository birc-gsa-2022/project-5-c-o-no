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
    int k;
    char* editString;
    int editIndex;
    int rangeStart;
    int rangeEnd;
};

struct CommenRec {
    int* pattern;
    int n;
    int m;
    int* D;
    int* C;
    int** O;
    struct Range* r;
};

void makeD(int* D, int* C, int** RO, int* pattern, int n, int m, struct Range* r);
void runApprox(int* pattern, int patIndex, int n, int m, int* D, int* C, int** O, int k, char* editString, int editIndex, struct Range* r);

#endif