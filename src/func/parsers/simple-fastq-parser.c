#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "simple-fastq-parser.h"

char *read_fastq_head(char **strptr) {
    char *string = *strptr;
    // We assume first character is '@', so skip it
    string++;
    // trim leading spaces
    while (isspace(*string)) {
        string++;
    }
    // Our header starts from this point
    char* header_start = string;

    // Move forward until a newline is encountered
    int i = 0;
    int containR = 0;
    while (string[i] != '\n') {
        if(string[i] == '\r') {
            containR = 1;
            break;
        }
        i++;
    }

    // Go back until first non-space character and thus trim trailing spaces and the \n
    int j = i;
    while (isspace(header_start[j])) {
        header_start[j] = '\0';
        j--;
    }

    // Move the pointer of the string to start of header plus count of all characters in header plus 1 or 2
    *strptr = header_start+i+1+containR;
    // Now the start of the string points to the first sequence character, and we can return the FASTA header

    return header_start;
}

char *read_fastq_pattern(char **strptr) {
    char *string = *strptr;
    int i = 0;
    int shift = 0;
    while (1) {
        // Increase shift for each \n we encounter
        while (string[i+shift] == '\n' || string[i+shift] == '\r') {
            shift++;
        }
        // If we encounter a \0 or >, then we are done
        if (string[i+shift] == '\0' || string[i+shift] == '@') {
            break;
        }
        // Add shift
        string[i] = string[i+shift];
        i++;
    }
    // Terminate with \0 and return the pointer to the first instance of >
    string[i] = '\0';
    *strptr = *strptr+i+shift;
    return string;
}

void updateFastqBySeq(char **strptr, struct Fastq* read) {
    char *string = *strptr;
    int i = 0;
    int shift = 0;
    while (1) {
        // Increase shift for each \n we encounter
        while (string[i+shift] == '\n' || string[i+shift] == '\r') {
            shift++;
        }
        // If we encounter a \0 or >, then we are done
        if (string[i+shift] == '\0' || string[i+shift] == '@') {
            break;
        }
        // Add shift
        string[i] = string[i+shift];
        i++;
    }
    // Terminate with \0 and return the pointer to the first instance of >
    string[i] = '\0';
    *strptr = *strptr+i+shift;
    read->seq = string;
    read->length = i-1;
    printf("i-1=%d\n", read->length);
    int* seqVal = malloc(i*sizeof *seqVal);
    for(int j=0; j<i; j++) {
        switch (string[j]) {
            case 'A':
                seqVal[j] = 1;
                break;
            case 'C':
                seqVal[j] = 2;
                break;
            case 'G':
                seqVal[j] = 3;
                break;
            case 'T':
                seqVal[j] = 4;
                break;
        }
    }
    read->seqVal = seqVal;

}

struct Fastq* parseFastq(char **strptr) {
    struct Fastq* read = malloc(sizeof *read);
    read->head = read_fastq_head(strptr);
    updateFastqBySeq(strptr, read);
    return read;
}

void freeFastq(struct Fastq* fastq) {
    free(fastq->seqVal);
    free(fastq->head);
    free(fastq->seq);
}