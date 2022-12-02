#include <stdio.h>
#include <stdlib.h>
#include "../src/func/parsers/simple-fasta-parser.h"
#include "../src/func/helper.h"
#include "../src/func/sa.h"
#include "../src/func/rotater.h"
#include <sys/time.h>
#define FPA 32
#define LPA 126
#define headerBufferSizeFasta 6
#define headerBufferSizeFastq 5


//Change for different tests:
#define maxN 10000
#define minN 1000
#define stepN 1000
#define minM 100
#define stepM 100
#define alphabetSize 1
#define TIMEFILENAMEPROCESS "timingProcess1.csv"
#define TIMEFILENAMESEARCH "timingSearch1.csv"
#define TIMEFILENAMEEDITS "timingEdits1.csv"


char * getPattern(int m, char* p) {
    for(int i=0; i<m; i++) {
        p[i] = 'a';
    }
    p[m] = '\0';
    return p;
}

void genLetters(int length, char* str, char startingChar) {
    str[0] = startingChar;
    str[1] = 'a';
    str[2] = '\n';
    for(int i=3; i<length+3; i+= alphabetSize) {
        for(int j=0; j<alphabetSize; j++) {
            if(i+j<length+3) {
                char curChar = 'a';
                curChar += (curChar >= '>') + (curChar >= '@');
                str[i+j] = curChar;
            }
            else break;
        }
    }
    str[length+3] = '\n';
    str[length+4] = '\0';
}


int getTime(struct timeval start, struct timeval stop) {
    // Source: https://stackoverflow.com/questions/10192903/time-in-milliseconds-in-c
    return (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
}

struct FastaContainer * timeParseFasta(FILE* file, char *fasta_str) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);

    struct FastaContainer *con = parse_fasta(fasta_str);

    gettimeofday(&stop, NULL);
    fprintf(file, "%d,", getTime(start, stop));

    return con;
}

int** timeConstructSAsRadix(FILE* file, struct FastaContainer* fastaContainer, int revserse) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);

    int** SAs = revserse ? constructMultipleRevSA(fastaContainer) : constructMultipleSA(fastaContainer);

    gettimeofday(&stop, NULL);
    fprintf(file, "%d,", getTime(start, stop));

    return SAs;
}

void timeProcessFasta(FILE* res_file, FILE* process_file, struct FastaContainer* fastaContainer, int** SAs, int** revSAs) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);

    processFastas(process_file, fastaContainer, SAs, revSAs);

    gettimeofday(&stop, NULL);
    fprintf(res_file, "%d,", getTime(start, stop));
}

void timeReadProcessedFile(char* readString, FILE* resFile, int edits) {
    // Do the whole thing

    char* processString = read_file("processed.txt");

    struct timeval stop, start;
    gettimeofday(&start, NULL);

    readFromProcessed(processString, readString, edits);

    gettimeofday(&stop, NULL);
    fprintf(resFile, "%d,", getTime(start, stop));
}

int main() {
    FILE* resfileProcess = fopen(TIMEFILENAMEPROCESS, "w+");
    FILE* processFileProcess = fopen("processed.txt", "w+");
    fprintf(resfileProcess,"n,parseTime,SATime,revSA,ProcessTime,σ\n");




    char *x = malloc((sizeof *x)*(maxN + headerBufferSizeFasta));
    for(int n=minN; n<maxN+1; n+=stepN) {
        fprintf(resfileProcess, "%d,", n);
        genLetters(n, x, '>');
        struct FastaContainer* fc = timeParseFasta(resfileProcess, x);
        int** SAs = timeConstructSAsRadix(resfileProcess, fc, 0);
        int** revSAs = timeConstructSAsRadix(resfileProcess, fc, 1);
        timeProcessFasta(resfileProcess, processFileProcess, fc, SAs, revSAs);
        fprintf(resfileProcess, "%d\n", alphabetSize);
    }
    fclose(processFileProcess);
    fclose(resfileProcess);

    FILE* resfileSearch = fopen(TIMEFILENAMESEARCH, "w+");
    fprintf(resfileSearch,"m,processTime,n,step,σ\n");

    char *readString = malloc((sizeof *readString)*(maxN+headerBufferSizeFastq));

    for(int m=minM; m<minN+1; m+=stepM) {
        fprintf(resfileSearch, "%d,", m);
        genLetters(m, readString, '@');
        timeReadProcessedFile(readString, resfileSearch, 0);
        fprintf(resfileSearch, "%d,%d,%d\n", maxN,stepN, alphabetSize);
    }
    fclose(resfileSearch);

    FILE* resfileEdit = fopen(TIMEFILENAMEEDITS, "w+");
    fprintf(resfileEdit,"edits,time,m, n,σ\n");

    for(int e=0; e<3; e++) {
        genLetters(100, readString, '@');
        fprintf(resfileEdit, "%d", e);
        timeReadProcessedFile(readString, resfileEdit, e);
        fprintf(resfileEdit, "%d,%d,%d\n", 100,stepN,alphabetSize);
    }
    fclose(resfileEdit);

}