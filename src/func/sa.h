#ifndef SA_h
#define SA_h
#include "parsers/simple-fasta-parser.h"

struct Interval {
    int start;
    int mid;
    int end;
};

int *constructSA(struct Fasta fasta, int reverse);
int **constructMultipleSA(struct FastaContainer *fastaContainer);
int **constructMultipleRevSA(struct FastaContainer *fastaContainer);
int *constructSARadix(struct Fasta fasta, int reverse);
int **constructMultipleSARadix(struct FastaContainer *fastaContainer);
int **constructMultipleRevSARadix(struct FastaContainer *fastaContainer);
struct Interval binarySearch(const char* x, const int* sa, char patchar, int parIndex, struct Interval interval, int mode);
struct Interval searchPatternInSA(struct Fasta fasta, const char* pattern, int* sa, int m);

#endif // SA_h