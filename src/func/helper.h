#ifndef HELPER_H
#define HELPER_H
#include <stdio.h>
#include "parsers/simple-fasta-parser.h"

struct ReadContainer {
    int count;
    char** heads;
    char** patterns;
    int* patLens;
};

void printIntArray(int * a, int len);
void printString(char * a, int len);
char *read_file(const char *file_name);
void write_to_file(FILE* fpt, char* content, int len);
char* get_file_name_by_fa(const char* faName);
FILE* get_file(const char* file_name);
struct ReadContainer* makeReadContainer(char* readString);
void processFastas(FILE* processFile, struct FastaContainer* fastaContainer, int** SAs, int** revSAs);
void readFromProcessed(char *processString, char* readString, int allowedEdits);

#endif //HELPER_H
