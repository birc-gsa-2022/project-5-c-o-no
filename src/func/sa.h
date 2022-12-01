#ifndef SA_h
#define SA_h
#include "parsers/simple-fasta-parser.h"
#include <stdint.h>

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
int *constructSAPrefixDoubling(struct Fasta fasta, int reverse);
int **constructMultipleSAPrefixDoubling(struct FastaContainer *fastaContainer);
int **constructMultipleRevSAPrefixDoubling(struct FastaContainer *fastaContainer);

struct Interval binarySearch(const char* x, const int* sa, char patchar, int parIndex, struct Interval interval, int mode);
struct Interval searchPatternInSA(struct Fasta fasta, const char* pattern, int* sa, int m);
int radixSort64Interval(int start, int end, int* sa, const uint64_t* keys);
uint64_t getByte(uint64_t key, int index);

#endif // SA_h