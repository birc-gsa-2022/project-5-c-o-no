#ifndef ROTATER_H
#define ROTATER_H

struct Range {
    int start;
    int end;
};

int oLookUp(int** o, int searchChar, int i);
int jump(int bwtIndex, int jumpChar, int* c, int** o);
void limitRangeByChar(int jumpChar, struct Range* r, int* C, int** O);
void rotateString(const int* string, int stringLen, int* c, int** o, int btwLen, struct Range* r);
void makeOandC(const int* string, int n, int** O, int* C, int alphabetSize);

#endif