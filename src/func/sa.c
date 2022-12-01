#include "sa.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "constants.h"
#define BYTEPOSSIBILITIES 256

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

uint64_t getByte(uint64_t key, int index) {
    uint32_t maskLow = 255; //11111111
    uint64_t maskHigh = (uint64_t)maskLow << index*8;
    uint64_t rightShift = (key & maskHigh) >> index*8;
    return (rightShift & maskLow);
}

int radixSort64Interval(int start, int end, int* sa, const uint64_t* keys) {
    //Sort the interval
    int intervalSize = end-start+1;
    int* saWriter = malloc(intervalSize*sizeof *saWriter);

    //Sort 8 bytes
    for(int byteIndex=0; byteIndex<8; byteIndex++) {
        int* bucketIndices = calloc(BYTEPOSSIBILITIES, sizeof *bucketIndices);
        int* buckets = calloc(BYTEPOSSIBILITIES, sizeof *buckets);

        //Count number of each combination
        for(int j=start; j<end+1; j++) {
            uint64_t byte = getByte(keys[j], byteIndex);
            bucketIndices[byte]++;
        }

        int accum = 0;
        for(int j=0; j<BYTEPOSSIBILITIES; j++) {
            int sight = bucketIndices[j];
            bucketIndices[j] = accum;
            accum += sight;
        }

        //Sort
        for(int j=start; j<=end; j++) {
            uint64_t byte = getByte(keys[j], byteIndex);
            int saIndex = bucketIndices[byte] + (buckets[byte]++);
            saWriter[saIndex-start] = sa[j];
        }

        //Write to sa
        for(int j=start; j<=end; j++) {
            sa[j] = saWriter[j-start];
        }
        //memcpy(&sa[start], saWriter, intervalSize*sizeof *saWriter);

    }
    free(saWriter);
}

int radixSort64(int* sa, const uint64_t* keys, uint32_t* rank, const int n) {
    //Return 1 if sorted

    uint64_t prevKey = keys[0];
    int start = 0;
    int end = 0;
    int sorted = 1;
    int rankVal = 1;

    for(int i=1; i<n; i++) {
        //Find interval to sort
        rank[i] = rankVal;
        if(keys[i]==prevKey) end++;
        else{
            prevKey = keys[i];
            ++rankVal;

            if(start!=end) {
                sorted = 0;
                radixSort64Interval(start, end, sa, keys);
            }

            start = i;
            end = i;
        }

    }

    return sorted;
}


int* constructSAPrefixDoubling(struct Fasta fasta, int reverse) {
    int n = fasta.fasta_len;
    //int* x = malloc(n*sizeof *x);
    int* sa = malloc(n*sizeof *sa);
    uint32_t* rank = malloc(n*sizeof *rank);
    uint64_t* key = malloc(n*sizeof *key);


    //Initial sorting
    sa[n-1] = 11;
    rank[n-1] = 0;
    key[n-1] = 0;

    for(int i=0; i<n-1; i++) {
        sa[i] = i;
        rank[i] = (uint32_t)fasta.fastaSeqVal[reverse ? n-i-2 : i];
    }

    for(int i=0; i<n-2; i++) {
        key[i] = (uint64_t)rank[i]<<32 | rank[i+1];
    }
    key[n-2] = (uint64_t)rank[n-2]<<32;
    radixSort64Interval(0, n-1, sa, key);


    /*sa[0] = n-1;
    rank[0] = 0;
    key[0] = 0;

    int AIndex = 1;
    int CIndex = fasta.alphabet.sightings[1]+1;
    int GIndex = CIndex + fasta.alphabet.sightings[2];
    int TIndex = GIndex + fasta.alphabet.sightings[3];
    int bucketIndecies[4] = {AIndex, CIndex, GIndex, TIndex};
    int* buckets = calloc(NUMDIFFCHARS, sizeof *buckets);

    for(int i=0; i<n-1; i++) {
        uint32_t curRank = (uint32_t)fasta.fastaSeqVal[reverse ? n-i-2 : i];
        int index = (buckets[curRank-1]++) + bucketIndecies[curRank-1];
        rank[index] = curRank;
        sa[index] = i;
    }

    for(int i=0; i<n-1; i++) {
        key
    }*/


    for(int k=1; k<(int)log2(n); k <<= 1) {
        for(int i=0; i<n; i++) {
            int saVal = sa[i];
            uint64_t lower = saVal+k<n ? rank[saVal+k] : 0;
            key[i] = (uint64_t)rank[i]<<32 | lower;
            printf("key is %lld made from %d -> %lld and %lld\n", key[i], rank[i], (uint64_t)rank[i]<<32, lower);
        }


        if(radixSort64(sa, key, rank, n)) break;
    }


    return sa;
}




int **constructMultipleSAPrefixDoubling(struct FastaContainer *fastaContainer) {
    int **SAs = malloc(fastaContainer->numberOfFastas*sizeof *SAs);

    for(int i=0; i<fastaContainer->numberOfFastas; i++) {
        int *sa = constructSAPrefixDoubling(*(fastaContainer->fastas)[i], 0);
        SAs[i] = sa;
    }
    return SAs;
}

int **constructMultipleRevSAPrefixDoubling(struct FastaContainer *fastaContainer) {
    int **SAs = malloc(fastaContainer->numberOfFastas*sizeof *SAs);

    for(int i=0; i<fastaContainer->numberOfFastas; i++) {
        int *sa = constructSAPrefixDoubling(*(fastaContainer->fastas)[i], 1);
        SAs[i] = sa;
    }
    return SAs;
}








int *constructSA(struct Fasta fasta, int reverse) {
    return constructSARadix(fasta, reverse);
}

int **constructMultipleSA(struct FastaContainer *fastaContainer) {
    return constructMultipleSARadix(fastaContainer);
}

int **constructMultipleRevSA(struct FastaContainer *fastaContainer) {
    return constructMultipleRevSARadix(fastaContainer);
}
