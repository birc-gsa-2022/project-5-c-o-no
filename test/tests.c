#include <stdbool.h>
#include "minunit.h"
#include "../src/func/parsers/simple-fasta-parser.h"
#include "../src/func/helper.h"
#include "testHelper.h"


void test_setup(void) {
    /* Nothing */
}

void test_teardown(void) {
    /* Nothing */
}

//int (*test_fun)(char*, char*, int, int);


MU_TEST(test_parser_abc) {
    char * mal = malloc(sizeof(int)*3);
    char * sequence = "ABC";
    strcpy(mal, sequence);
    struct Fasta * f = malloc(sizeof *f);
    update_fasta_by_sequence(&mal, f);

    mu_assert_int_eq(4, f->fasta_len);
    mu_assert_int_eq(4, f->alphabet.size);
    for(int i = 0; i<4; i++) mu_assert_int_eq(1, f->alphabet.sightings[i]);

    mu_assert_int_eq(1, f->fasta_sequence[0]);
    mu_assert_int_eq(2, f->fasta_sequence[1]);
    mu_assert_int_eq(3, f->fasta_sequence[2]);
    mu_assert_int_eq(0, f->fasta_sequence[3]);

    mu_assert_int_eq('A', f->fasta_sequence_debugger[0]);
    mu_assert_int_eq('B', f->fasta_sequence_debugger[1]);
    mu_assert_int_eq('C', f->fasta_sequence_debugger[2]);
    mu_assert_int_eq('\0', f->fasta_sequence_debugger[3]);

    mu_assert_int_eq(1, f->alphabet.symbols['A']);
    mu_assert_int_eq(2, f->alphabet.symbols['B']);
    mu_assert_int_eq(3, f->alphabet.symbols['C']);

}


MU_TEST(test_parser_aaaa) {
    char * mal = malloc(sizeof(int)*4);
    char * sequence = "AAAA";
    strcpy(mal, sequence);
    struct Fasta * f = malloc(sizeof *f);
    update_fasta_by_sequence(&mal, f);

    mu_assert_int_eq(5, f->fasta_len);
    mu_assert_int_eq(2, f->alphabet.size);
    mu_assert_int_eq(1, f->alphabet.sightings[0]);
    mu_assert_int_eq(4, f->alphabet.sightings[1]);
}

MU_TEST(test_parser_aLong) {
    char * mal = malloc(sizeof(int)*100);
    char * sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    strcpy(mal, sequence);
    struct Fasta * f = malloc(sizeof *f);
    update_fasta_by_sequence(&mal, f);

    mu_assert_int_eq(101, f->fasta_len);
    mu_assert_int_eq(2, f->alphabet.size);
    mu_assert_int_eq(1, f->alphabet.sightings[0]);
    mu_assert_int_eq(100, f->alphabet.sightings[1]);
}

MU_TEST(test_parser_mis) {
    char *file = read_file(mis_fa_file);
    struct FastaContainer *fastasContainer = parse_fasta(file);
    struct Fasta **fastas = fastasContainer->fastas;

    mu_assert_string_eq("chr1", (*fastas)->fasta_head);
    // I = 1, M = 2, P = 3, S = 4
    int chr1[11] = {2,1,4,4,1,4,4,1,3,3,1};
    mu_assert_int_arr_eq(chr1, (*fastas)->fasta_sequence);
    mu_assert_int_eq(12, (*fastas)->fasta_len);

    // Go to next fasta
    fastas++;
    mu_assert_string_eq("chr2", (*fastas)->fasta_head);
    // I = 1, M = 2, P = 3, S = 4
    int chr2[22] = {2,1,4,4,1,4,4,1,3,3,1,2,1,4,4,1,4,4,1,3,3,1};
    mu_assert_int_arr_eq(chr2, (*fastas)->fasta_sequence);
    mu_assert_int_eq(23, (*fastas)->fasta_len);

    free_fasta_container(fastasContainer);
}

void run_all_tests() {
    MU_RUN_TEST(test_parser_abc);
    MU_RUN_TEST(test_parser_aaaa);
    MU_RUN_TEST(test_parser_aLong);
    MU_RUN_TEST(test_parser_mis);
}

MU_TEST_SUITE(fasta_parser_test_suite) {
    run_all_tests();
}


int main(int argc, char *argv[]) {
    // Path to the folder where test-data is located should be given.
    // (typically ../ if running from the debugger (project_root/cmake-build-debug/....
    // or ./ if running directly from project root)
    init(argv[1]);
    MU_RUN_SUITE(fasta_parser_test_suite);
    MU_REPORT();
    return MU_EXIT_CODE;
}