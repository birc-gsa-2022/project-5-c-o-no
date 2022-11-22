#ifndef SA_h
#define SA_h
#include "parsers/simple-fasta-parser.h"
//#include "parsers/simple-fastq-parser.h"

int *constructSARadix(struct Fasta fasta);
int **constructMultipleSARadix(struct FastaContainer *fastaContainer);
struct Interval {
    int start;
    int mid;
    int end;
};
struct Interval binarySearch(const char* x, const int* sa, char patchar, int parIndex, struct Interval interval, int mode);
struct Interval searchPatternInSA(struct Fasta fasta, const char* pattern, int* sa, int m);

#endif // SA_h