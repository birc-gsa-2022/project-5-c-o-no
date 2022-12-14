#include "helper.h"
#include "approx.h"
#include "parsers/simple-fasta-parser.h"
#include "parsers/simple-fastq-parser.h"
#include "rotater.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debugger.h"


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
    struct Fastq** reads = malloc(listSize*sizeof *reads);

    int count = 0;
    while (*readString != '\0') {
        if(count>=listSize) {
            listSize <<= 1;
            reads = realloc(reads,listSize*sizeof *reads);
        }
        struct Fastq* read = parseFastq(&readString);
        reads[count] = read;
        count++;
    }

    rc->count=count;
    rc->reads=reads;
    return rc;
}


void processFastas(FILE* processFile, struct FastaContainer* fastaContainer, int** SAs, int** revSAs) {
    //TODO make binary file instead of txt
    for(int i=0; i<fastaContainer->numberOfFastas; i++) {
        struct Fasta* fasta = fastaContainer->fastas[i];
        fprintf(processFile, "%s\n", fasta->fasta_head); //Save head
        fprintf(processFile, "%d\n", fasta->fasta_len); //Save length
        //fprintf(processFile, "%d\n", fasta->alphabet.size); //Save alphabetsize
        for(int j=0; j<fasta->fasta_len; j++) { //Save bwt
            //TODO We can change this to compress, but not necessary
            fprintf(processFile, "%d,", SAs[i][j] ? fasta->fastaSeqVal[SAs[i][j]-1] : 0);
        }
        fprintf(processFile, "\n");
        for(int j=0; j<fasta->fasta_len; j++) { //Save sa
            fprintf(processFile, "%d,", SAs[i][j]);
        }
        fprintf(processFile, "\n");
        for(int j=0; j<fasta->fasta_len; j++) { //Save revbwt
            //TODO We can change this to compress, but not necessary
            fprintf(processFile, "%d,", revSAs[i][j] ? fasta->fastaSeqVal[fasta->fasta_len-revSAs[i][j]-1] : 0);
        }
        fprintf(processFile, "\n");
        /*for(int j=0; j<128; j++) {
            //TODO do in parser
            if(fasta->alphabet.symbols[j]) fprintf(processFile, "%c", j);
        }*/
        fprintf(processFile, "\n");

        free(SAs[i]);
    }
}

void printEditString(struct ApproxMatch* am) {
    char* editString = am->editString;
    int len = am->editStringLen;
    char curChar = editString[len-1];
    int curCount = 1;
    for(int i=len-2; i>=0; i--) {
        if(editString[i] != curChar) {
            printf("%d%c", curCount, curChar);
            curChar = editString[i];
            curCount = 1;
        } else ++curCount;
    }
    printf("%d%c", curCount, curChar);
}

void readFromProcessed(char *processString, char* readString, int allowedEdits) {
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
        /*int alphabetSize = atoi(processString);
        while(*(processString++) != '\n') {}*/

        //bwt
        //TODO could make each symbol distinguishable without ,
        //TODO Don't make bwt, go directly to O and C

        int *bwt = malloc(n * sizeof *bwt);
        for(int i=0; i<n; i++) {
            bwt[i] = atoi(processString);
            while(*(processString++) != ',') {}
        }
        int **O = malloc(n*sizeof (*O));
        int* C = calloc(5, sizeof *C);

        makeOandC(bwt, n, O, C);


        while(*(processString++) != '\n') {}


        //sa
        int *sa = malloc(n * sizeof *sa);
        for(int i=0; i<n; i++) {
            sa[i] = atoi(processString);
            while(*(processString++) != ',') {}
        }


        while(*(processString++) != '\n') {}

        int *revbwt = malloc(n * sizeof *revbwt);
        for(int i=0; i<n; i++) {
            revbwt[i] = atoi(processString);
            while(*(processString++) != ',') {}
        }


        int **RO = malloc(n*sizeof (*RO));
        //TODO don't remake C
        int* rC = calloc(5, sizeof *rC);
        makeOandC(revbwt, n, RO, rC);



        while(*(processString++) != '\n') {}

        //int *alphabet = calloc(128, sizeof *alphabet);
        //for(int i=0; i<alphabetSize-1; i++) {
        //    alphabet[*(processString++)] = i+1;
        //}


        for (int j = 0; j < read_container->count; j++) {
            char *readHead = read_container->reads[j]->head;
            char *pattern = read_container->reads[j]->seq;
            int* patVal = read_container->reads[j]->seqVal;
            int m = read_container->reads[j]->length;

            int* D = malloc(m*sizeof *D);
            struct Range* Drange = malloc(sizeof *Drange);
            makeD(D, rC, RO, patVal, n, m, Drange);

            char* editString = malloc(allowedEdits+m+1);

            struct ApproxMatchContainer* matches = runApprox(patVal, n, m, D, C, O, allowedEdits, editString, saRange);

            for (int i = 0; i < matches->amount; ++i) {
                int start = matches->AMs[i]->rStart;
                int end = matches->AMs[i]->rEnd;

                for(int saIndex=start; saIndex<end; saIndex++) {
                    printf("%s\t%s\t%d\t", readHead, fastaHead, sa[saIndex]+1);
                    printEditString(matches->AMs[i]);
                    printf("\t%s\n", pattern);
                }
            }
            free(editString);
            freeApproxMatchContainer(matches);
            free(matches);
        }

        while(*(processString++) != '\n') {}

        for(int i=0; i<n; i++) {
            free(O[i]);
            free(RO[i]);
        }

        free(O);
        free(C);
        free(RO);
        free(rC);
        free(sa);
    }
    for(int i=0; i<read_container->count; i++) {
        freeFastq(read_container->reads[i]);
    }
}
