#ifndef SA_TESTHELPER_H
#define SA_TESTHELPER_H
#include <malloc.h>
#include <string.h>

char *test_data_path;

char *mis_fa_file;
char *mis_fastq_file;

void lazy_alloc() {
    test_data_path = malloc(1000);
    mis_fa_file = malloc(1000);
    mis_fastq_file = malloc(1000);
}

void init(char *root) {
    lazy_alloc();

    test_data_path = root;
    strcpy(mis_fa_file, test_data_path);
    strcpy(mis_fastq_file, test_data_path);

    // Append test-data/ to all files
    strcat(mis_fa_file, "test-data/");
    strcat(mis_fastq_file, "test-data/");

    strcat(mis_fa_file, "mis.fa");
    strcat(mis_fastq_file, "mis.fastq");
}

#endif //SA_TESTHELPER_H
