#include "sa.h"
#include <string.h>
#include <malloc.h>


int *constructSARadix(struct Fasta fasta, int reverse) {
    int n = fasta.fasta_len;
    char* x = malloc(n*sizeof *x);
    x[n-1] = 0;
    for(int i=0; i<n-1; i++) {
        x[i] = fasta.fasta_sequence[reverse ? n-i-2 : i];
    }
    //char *x = fasta.fasta_sequence;
    int *saReader = malloc(n * sizeof *saReader);
    int *saWriter = malloc(n * sizeof *saReader);

    for (int i = 0; i<n; i++) {
        saReader[i] = i;
    }
    int *bucketsIndices = malloc(fasta.alphabet.size*sizeof *bucketsIndices);
    int accumSum = 0;
    for(int i=0; i<fasta.alphabet.size; i++) {
        bucketsIndices[i] = accumSum;
        accumSum += fasta.alphabet.sightings[i];
    }
    int *buckets = malloc(fasta.alphabet.size*sizeof *buckets);

    for(int i=n-1; i>=0; i--) {
        memset(buckets, 0, fasta.alphabet.size * sizeof *buckets);
        for(int j=0; j<n; j++) {
            int charIndex = (saReader[j] + i) % n;
            char c = x[charIndex];
            int elemInBucket = buckets[c]++;
            int saIndex = bucketsIndices[c] + elemInBucket;
            saWriter[saIndex] = saReader[j];
        }
        int *temp = saReader;
        saReader = saWriter;
        saWriter = temp;
    }

    free(saWriter);
    free(buckets);
    free(x);

    return saReader;
}

int **constructMultipleSARadix(struct FastaContainer *fastaContainer) {
    int **SAs = malloc(fastaContainer->numberOfFastas*sizeof *SAs);

    for(int i=0; i<fastaContainer->numberOfFastas; i++) {
        int *sa = constructSARadix(*(fastaContainer->fastas)[i], 0);
        SAs[i] = sa;
    }
    return SAs;
}

//TODO combine with constructMultipleSARadix
int **constructMultipleRevSARadix(struct FastaContainer *fastaContainer) {
    int **SAs = malloc(fastaContainer->numberOfFastas*sizeof *SAs);

    for(int i=0; i<fastaContainer->numberOfFastas; i++) {
        int *sa = constructSARadix(*(fastaContainer->fastas)[i], 1);
        SAs[i] = sa;
    }
    return SAs;
}



struct Interval binarySearch(const char* x, const int* sa, char patchar, int parIndex, struct Interval interval, int mode) {
    char xchar;
    while (interval.start != interval.end) {
        interval.mid = (interval.start+interval.end)/2;
        xchar = x[sa[interval.mid]+parIndex];
        if(patchar == xchar) {
            if(!mode) return interval;
            if(mode==1) interval.start = interval.mid + 1;
            else interval.end = interval.mid;
        }
        else {
            if (xchar > patchar) {
                interval.end = interval.mid;
            }
            else interval.start = interval.mid + 1;
        }
    }
    return interval;
}

struct Interval searchPatternInSA(struct Fasta fasta, const char* pattern, int* sa, int m) {
    char * x = fasta.fasta_sequence;
    struct Interval *interval = malloc(sizeof *interval);
    struct Interval *intervalSaver = malloc(sizeof *interval);
    interval->start = 0;
    interval->end = fasta.fasta_len;

    for(int i=0; i<m; i++) {
        char patchar = (char) fasta.alphabet.symbols[pattern[i]];
        *interval = binarySearch(x, sa, patchar, i, *interval, 0);
        if(interval->start == interval->end) return *interval;
        intervalSaver->mid = interval->mid;
        intervalSaver->end = interval->end;

        interval->end = interval->mid;
        *interval = binarySearch(x, sa, patchar, i, *interval, -1);
        intervalSaver->start = interval->start;

        interval->start = intervalSaver->mid;
        interval->end = intervalSaver->end;

        *interval = binarySearch(x, sa, patchar, i, *interval, 1);
        interval->start = intervalSaver->start;
    }

    free(intervalSaver);
    return *interval;
}



