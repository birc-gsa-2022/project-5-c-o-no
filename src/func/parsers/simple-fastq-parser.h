#ifndef SIMPLE_FASTQ_PARSER_H
#define SIMPLE_FASTQ_PARSER_H

struct Fastq {
    char* head;
    char* seq;
    int* seqVal;
    int length;
};

char *read_fastq_head(char **strptr);
char *read_fastq_pattern(char **strptr);
struct Fastq* parseFastq(char **strptr);
void freeFastq(struct Fastq* fastq);

#endif //SIMPLE_FASTQ_PARSER_H
