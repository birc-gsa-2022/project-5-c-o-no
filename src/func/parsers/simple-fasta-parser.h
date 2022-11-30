#ifndef SIMPLE_FASTA_PARSER_H
#define SIMPLE_FASTA_PARSER_H

struct Alphabet {
    int size;
    int *symbols;
    int *sightings;
};

struct Fasta {
    char *fasta_head;
    char *fasta_sequence;
    //char *fasta_sequence_debugger;
    int fasta_len;
    struct Alphabet alphabet;
};

struct FastaContainer {
    struct Fasta ** fastas;
    int numberOfFastas;
};

int *alloc_sightings(int *bigAlphabet, int alphabetSize);
char *read_fasta_head(char **strptr);
void update_fasta_by_sequence(char **strptr, struct Fasta *f);
struct FastaContainer *parse_fasta(char *fasta_str);
void freeFasta(struct Fasta *fasta);
void free_fastas(struct Fasta **fastas, int count);
void free_fasta_container(struct FastaContainer *fc);

#endif //SIMPLE_FASTA_PARSER_H
