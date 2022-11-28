#include "helper.h"
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include "parsers/simple-fastq-parser.h"
#include "parsers/simple-fasta-parser.h"
#include "rotater.h"

void printIntArray(int * a, int len) {
    printf("[");
    for(int i=0; i<len; i++) {
        printf("%d,", a[i]);
    }
    printf("]\n");
}

void printString(char * a, int len) {
    for(int i=0; i<len; i++) {
        printf("%c", a[i]);
    }
    printf("\n");
}

void printO(int** o, int numRows) {
    printf("----------------\n");
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%d, ", o[i][j]);
        }
        printf("\n");
    }
    printf("----------------\n");
}


char *read_file(const char *file_name) {
    FILE * fp = fopen(file_name, "rb");

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);  /* same as rewind(f); */

    char *string = malloc(fsize + 1);
    fread(string, fsize, 1, fp);

    string[fsize] = '\0'; // terminate with zero
    fclose(fp);
    return string;
}

FILE* get_file(const char* file_name) {
    return fopen(file_name, "w+");
}

char* get_file_name_by_fa(const char* faName) {
    size_t nameLen = strlen(faName)-2;
    char* processingFileName = malloc(nameLen+4);
    strcpy(processingFileName, faName);
    strcpy(processingFileName+nameLen, "txt");
    return processingFileName;
}

FILE* get_file_by_fa(const char* faName) {
    char* fileName = get_file_name_by_fa(faName);
    FILE* file = get_file(fileName);
    free(fileName);
    return file;
}

void write_to_file(FILE* fpt, char* content, int len) {
    for(int i=0; i<len; i++) {
        fprintf(fpt, "%c", *(content++));
    }
}

struct ReadContainer* makeReadContainer(char* readString) {
    struct ReadContainer* rc = malloc(sizeof *rc);

    int listSize = 8;
    int* patLens = malloc(listSize*sizeof *patLens);
    char** heads = malloc(listSize*sizeof *heads);
    char** patterns = malloc(listSize*sizeof *patterns);

    int count = 0;
    while (*readString != '\0') {
        if(count>=listSize) {
            listSize <<= 1;
            patLens = realloc(patLens, listSize*sizeof *patLens);
            heads = realloc(heads, listSize*sizeof *heads);
            patterns = realloc(patterns, listSize*sizeof *patterns);
        }
        heads[count] = read_fastq_head(&readString);
        patterns[count] = read_fastq_pattern(&readString);
        patLens[count] = (int) strlen(patterns[count]); // TODO in parsing
        count++;
    }

    rc->count=count;
    rc->patterns=patterns;
    rc->heads=heads;
    rc->patLens=patLens;
    return rc;
}


void processFastas(FILE* processFile, struct FastaContainer* fastaContainer, int** SAs) {
    //TODO make binary file instead of txt
    for(int i=0; i<fastaContainer->numberOfFastas; i++) {
        struct Fasta* fasta = fastaContainer->fastas[i];
        fprintf(processFile, "%s\n", fasta->fasta_head); //Save head
        fprintf(processFile, "%d\n", fasta->fasta_len); //Save length
        fprintf(processFile, "%d\n", fasta->alphabet.size); //Save alphabetsize
        for(int j=0; j<fasta->fasta_len; j++) { //Save bwt
            //TODO We can change this to compress, but not necessary
            fprintf(processFile, "%d,", SAs[i][j] ? fasta->fasta_sequence[SAs[i][j]-1] : 0);
        }
        fprintf(processFile, "\n");
        for(int j=0; j<fasta->fasta_len; j++) { //Save sa
            fprintf(processFile, "%d,", SAs[i][j]);
        }
        fprintf(processFile, "\n");

        for(int j=0; j<128; j++) {
            //TODO do in parser
            if(fasta->alphabet.symbols[j]) fprintf(processFile, "%c", j);
        }
        fprintf(processFile, "\n");

        free(SAs[i]);
    }
}


void readFromProcessed(char *processString, char* readString) {
    struct Range* saRange = malloc(sizeof *saRange);
    struct ReadContainer *read_container = makeReadContainer(readString);

    while(*processString != '\0') {
        char *fastaHead = processString;
        while (*(++processString) != '\n') {}
        if(*(processString-1) == '\r') {
            *(processString-1) = '\0';
        }
        else {
            *(processString) = '\0';
        }

        processString++;
        int n = atoi(processString); //atoi stops at first non-int
        while(*(processString++) != '\n') {}
        int alphabetSize = atoi(processString);
        while(*(processString++) != '\n') {}

        //bwt
        //TODO could make each symbol distinguishable without ,
        //TODO Don't make bwt, go directly to O and C

        int *bwt = malloc(n * sizeof *bwt);
        for(int i=0; i<n; i++) {
            bwt[i] = atoi(processString);
            while(*(processString++) != ',') {}
        }
        int **O = malloc(n*sizeof (*O));
        int* C = calloc(alphabetSize, sizeof *C);
        makeOandC(bwt, n, O, C, alphabetSize);

        while(*(processString++) != '\n') {}

        //sa
        int *sa = malloc(n * sizeof *sa);
        for(int i=0; i<n; i++) {
            sa[i] = atoi(processString);
            while(*(processString++) != ',') {}
        }

        while(*(processString++) != '\n') {}

        int *alphabet = calloc(128, sizeof *alphabet);
        for(int i=0; i<alphabetSize-1; i++) {
            alphabet[*(processString++)] = i+1;
        }



        for (int j = 0; j < read_container->count; j++) {
            char *readHead = read_container->heads[j];
            char *pattern = read_container->patterns[j];
            int pattern_len = read_container->patLens[j];

            int *patternConvert = malloc(pattern_len*sizeof *patternConvert);
            for(int i=0; i<pattern_len; i++) {
                patternConvert[i] = alphabet[pattern[i]];
            }

            rotateString(patternConvert, pattern_len, C, O, n, saRange);
            free(patternConvert);

            int start = saRange->start;
            int end = saRange->end;
            for(int i=start; i<end; i++) {
                printf("%s\t%s\t%d\t%dM\t%s\n", readHead, fastaHead, sa[i]+1, pattern_len, pattern);
            }
        }

        while(*(processString++) != '\n') {}

        for(int i=0; i<n; i++) {
            free(O[i]);
        }
        free(O);
        free(C);
        free(sa);
    }
}
