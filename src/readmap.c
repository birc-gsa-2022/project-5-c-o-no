#include "func/helper.h"
#include "func/parsers/simple-fasta-parser.h"
#include "func/sa.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(const char *progname)
{
    fprintf(stderr,
            "Usage: %s -p genome\n       %s -d dist genome reads\n",
            progname, progname);
    exit(1);
}

int main(int argc, char const *argv[])
{
    if ((argc == 3) && strcmp("-p", argv[1]) == 0)
    {
        // preprocessing
        //printf("Preprocessing genome %s.\n", argv[2]);

        // preprocessing
        char* fastaStr = read_file(argv[2]);
        char* processFileName = get_file_name_by_fa(argv[2]);
        FILE* processFile = get_file(processFileName);
        free(processFileName);
        struct FastaContainer* fastaContainer = parse_fasta(fastaStr);
        int** SAs = constructMultipleSA(fastaContainer);
        int** revSAs = constructMultipleRevSA(fastaContainer);
        processFastas(processFile, fastaContainer, SAs, revSAs);
        fclose(processFile);
        free_fasta_container(fastaContainer);
        free(SAs);
    }
    else if ((argc == 5) && strcmp("-d", argv[1]) == 0)
    {
        int d = atoi(argv[2]); // you can do better than atoi, but I'm feeling lazy
        const char *genome = argv[3];
        const char *reads = argv[4];
        //printf("Mapping in genome %s for reads in %s within distance %d.\n", genome, reads, d);

        char* processFileName = get_file_name_by_fa(genome);
        char* processString = read_file(processFileName);
        free(processFileName);
        char* readString = read_file(reads);
        readFromProcessed(processString, readString, d);
        free(processString);
        free(readString);
    }
    else
    {
        print_usage(argv[0]);
    }

    return 0;
}
