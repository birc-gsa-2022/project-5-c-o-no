#include "sa.h"
#include <string.h>
#include <malloc.h>


int *constructSARadix(struct Fasta fasta) {
    char *x = fasta.fasta_sequence;
    int n = fasta.fasta_len;
    int *sa = malloc(n * sizeof *sa);
    int *saCopy = malloc(n * sizeof *sa);

    for (int i = 0; i<n; i++) {
        sa[i] = i;
    }
    int *bucketsIndices = fasta.alphabet.sightings;
    int accumSum = 0;
    for(int i=0; i<fasta.alphabet.size; i++) {
        int sighting = bucketsIndices[i];
        bucketsIndices[i] = accumSum;
        accumSum += sighting;
    }
    int *buckets = malloc(fasta.alphabet.size*sizeof *buckets);

    for(int i=n-1; i>=0; i--) {
        memset(buckets, 0, fasta.alphabet.size * sizeof *buckets);
        for(int j=0; j<n; j++) {
            int charIndex = (sa[j] + i) % n;
            char c = x[charIndex];
            int elemInBucket = buckets[c]++;
            int saIndex = bucketsIndices[c] + elemInBucket;
            saCopy[saIndex] = sa[j];
        }
        int *temp = sa;
        sa = saCopy;
        saCopy = temp;
    }

    free(saCopy);
    free(buckets);

    return sa;
}

int **constructMultipleSARadix(struct FastaContainer *fastaContainer) {
    int **SAs = malloc(fastaContainer->numberOfFastas*sizeof *SAs);

    for(int i=0; i<fastaContainer->numberOfFastas; i++) {
        int *sa = constructSARadix(*(fastaContainer->fastas)[i]);
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



