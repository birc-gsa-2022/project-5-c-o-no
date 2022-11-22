#include <stdio.h>
#include <malloc.h>
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
#define maxN 100000
#define minN 1000
#define stepN 1000
#define minM 1000
#define stepM 1000
#define alphabetSize 90
#define TIMEFILENAMEPROCESS "timingProcess4.csv"
#define TIMEFILENAMESEARCH "timingSearch4.csv"


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
                char curChar = FPA+j;
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

int** timeConstructSAsRadix(FILE* file, struct FastaContainer* fastaContainer) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);

    int** SAs = constructMultipleSARadix(fastaContainer);

    gettimeofday(&stop, NULL);
    fprintf(file, "%d,", getTime(start, stop));

    return SAs;
}

void timeProcessFasta(FILE* res_file, FILE* process_file, struct FastaContainer* fastaContainer, int** SAs) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);

    processFastas(process_file, fastaContainer, SAs);

    gettimeofday(&stop, NULL);
    fprintf(res_file, "%d,", getTime(start, stop));
}

void timeReadProcessedFile(char* readString, FILE* resFile) {
    // Do the whole thing
    char* processString = read_file("processed.txt");

    struct timeval stop, start;
    gettimeofday(&start, NULL);

    readFromProcessed(processString, readString);

    gettimeofday(&stop, NULL);
    fprintf(resFile, "%d,", getTime(start, stop));
}

/*
 * Generates a BWT for a string of length n, if we know that
 * the string is a repeated pattern - i.e:
 * abcdeabcdeabcdeab
*/
int *generateBwt(int n) {
    int *bwt = malloc((n+1) * sizeof *bwt);
    int numOfBlocks = n / alphabetSize;
    int lastBlockSize = n%alphabetSize;
    if(numOfBlocks>0){
        bwt[0] = lastBlockSize ? lastBlockSize : alphabetSize; //Char before sentinel
        int endOfFirstRepetion = numOfBlocks + (n % alphabetSize > 0);
        for (int i = 1; i < endOfFirstRepetion; i++) {
            bwt[i] = alphabetSize;
        }
        bwt[endOfFirstRepetion] = 0;

        int repetionStart = endOfFirstRepetion+1;
        int currChar = 1;
        while(repetionStart<n+1){
            //Fill
            int inLastBlock = (currChar<=lastBlockSize) - (currChar == bwt[0]);
            for(int i=0; i < numOfBlocks + inLastBlock; i++) {
                bwt[repetionStart+i] = currChar;
            }
            //Jump based on size of repetion
            repetionStart += numOfBlocks + inLastBlock;
            currChar++;
        }
    }
    else {
        bwt[0] = lastBlockSize;
        bwt[1] = 0;
        for(int i=2; i<n+1; i++) {
            bwt[i] = i-1;
        }
    }

    return bwt;
}

void timeRotation(int m, const char* pattern, int n, FILE *resFile) {
    int* alphabet = calloc(128, sizeof *alphabet);
    for(int j=0; j<alphabetSize; j++) {
        char curChar = FPA+j;
        alphabet[curChar+(curChar >= '>')+(curChar >= '@')] = j+1;
    }
    struct Range* saRange = malloc(sizeof *saRange);

    int **O = malloc((n+1)*sizeof (*O));
    int* C = calloc(alphabetSize+1, sizeof *C);

    int* bwt = generateBwt(n);

    makeOandC(bwt, n+1, O, C, alphabetSize+1);

    struct timeval stop, start;
    gettimeofday(&start, NULL);

    int *patternConvert = malloc((m+1)*sizeof *patternConvert);
    for(int i=0; i<m; i++) {
        patternConvert[i] = alphabet[pattern[i+3]];
    }
    free(alphabet);
    patternConvert[m] = 0;

    rotateString(patternConvert, m, C, O, n, saRange);
    free(saRange);
    free(C);
    for(int i=0; i<n; i++) {
        free(O[i]);
    }
    free(O);
    free(patternConvert);
    free(bwt);

    gettimeofday(&stop, NULL);
    fprintf(resFile, "%d,", getTime(start, stop));

}


int main() {
    FILE* resfileProcess = fopen(TIMEFILENAMEPROCESS, "w+");
    FILE* processFileProcess = fopen("processed.txt", "w+");
    fprintf(resfileProcess,"n,parseTime,SATime,ProcessTime,σ\n");

    char *x = malloc((sizeof *x)*(maxN + headerBufferSizeFasta));
    for(int n=minN; n<maxN+1; n+=stepN) {
        fprintf(resfileProcess, "%d,", n);
        genLetters(n, x, '>');
        struct FastaContainer* fc = timeParseFasta(resfileProcess, x);
        int** SAs = timeConstructSAsRadix(resfileProcess, fc);
        timeProcessFasta(resfileProcess, processFileProcess, fc, SAs);
        fprintf(resfileProcess, "%d\n", alphabetSize);
    }
    fclose(processFileProcess);

    FILE* resfileSearch = fopen(TIMEFILENAMESEARCH, "w+");
    fprintf(resfileSearch,"m,processTime,SearchTime,maxn,minn,step,σ\n");

    char *readString = malloc((sizeof *readString)*(maxN+headerBufferSizeFastq));

    for(int m=minM; m<maxN+1; m+=stepM) {
        fprintf(resfileSearch, "%d,", m);
        genLetters(m, readString, '@');
        timeReadProcessedFile(readString, resfileSearch);
        timeRotation(m, readString, maxN, resfileSearch);
        fprintf(resfileSearch, "%d,%d,%d,%d\n", maxN,minN,stepN, alphabetSize);
    }
}