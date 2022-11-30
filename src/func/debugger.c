#include <stdio.h>
#include "debugger.h"

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
