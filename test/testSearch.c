#include "minunit.h"
#include "../src/sa.h"
#include "../src/helper.h"
#include "testHelper.h"


void test_setup(void) {
    /* Nothing */
}

void test_teardown(void) {
    /* Nothing */
}

MU_TEST(test_search_abc) {
    char * seq = malloc(sizeof(* seq)*6);
    char * sequence = "ABCABC";
    strcpy(seq, sequence);
    char * pattern = malloc(sizeof(* pattern)*3);
    char * pat = "ABC";
    strcpy(pattern, pat);
    int * sa = malloc(sizeof(* sa)*7);
    sa[0] = 6;
    sa[1] = 3;
    sa[2] = 0;
    sa[3] = 4;
    sa[4] = 1;
    sa[5] = 5;
    sa[6] = 2;
    struct Fasta * fasta = malloc(sizeof *fasta);
    update_fasta_by_sequence(&seq, fasta);

    struct Interval inter = searchPatternInSA(*fasta, pattern, sa, 3);
    mu_assert_int_eq(3, inter.end);
    mu_assert_int_eq(1, inter.start);
}


MU_TEST(test_search_multiple_sa_constr) {
    char *file = read_file(mis_fa_file);
    struct FastaContainer *fastasContainer = parse_fasta(file);
    struct Fasta **fastas = fastasContainer->fastas;

    mu_assert_int_eq(2, fastasContainer->numberOfFastas);
    int** sa = constructMultipleSARadix(fastasContainer);

    int *sa1 = malloc(12*sizeof *sa1);
    sa1[0] = 11;
    sa1[1] = 10;
    sa1[2] = 10-3;
    sa1[3] = 10-6;
    sa1[4] = 10-9;
    sa1[5] = 0;
    sa1[6] = 9;
    sa1[7] = 8;
    sa1[8] = 6;
    sa1[9] = 3;
    sa1[10] = 5;
    sa1[11] = 2;

    int *sa2 = malloc(23*sizeof *sa2);
    sa2[0] = 22;
    sa2[1] = 21;
    sa2[2] = 10;
    sa2[3] = 18;
    sa2[4] = 7;
    sa2[5] = 15;
    sa2[6] = 4;
    sa2[7] = 12;
    sa2[8] = 1;
    sa2[9] = 11;
    sa2[10] = 0;
    sa2[11] = 20;
    sa2[12] = 9;
    sa2[13] = 19;
    sa2[14] = 8;
    sa2[15] = 17;
    sa2[16] = 6;
    sa2[17] = 14;
    sa2[18] = 3;
    sa2[19] = 16;
    sa2[20] = 5;
    sa2[21] = 13;
    sa2[22] = 2;

    mu_assert_int_arr_eq_n(sa1, *sa, 11);
    sa++;
    mu_assert_int_arr_eq_n(sa2, *sa, 22);
}

MU_TEST(test_binary_search) {
    char * seq = "mississippi";
    int sa[12] = {11,10,7,4,1,0,9,8,6,3,5,2};
    char patChar = 'i';
    struct Interval *interval = malloc(sizeof *interval);

    // Find 'i'
    interval->start = 0;
    interval->end = 12;

    struct Interval res = binarySearch(seq, sa, patChar, 0, *interval, 0);
    mu_assert_int_eq(0, res.start);
    mu_assert_int_eq(3, res.mid);
    mu_assert_int_eq(6, res.end);

    interval->end = 3;
    res = binarySearch(seq, sa, patChar, 0, *interval, -1);
    mu_assert_int_eq(1, res.start);

    interval->start = 3;
    interval->end = 6;
    res = binarySearch(seq, sa, patChar, 0, *interval, 1);
    mu_assert_int_eq(5, res.end);

    // Find 's'
    interval->start = 1;
    interval->end = 5;
    patChar = 's';
    res = binarySearch(seq, sa, patChar, 1, *interval, 0);
    mu_assert_int_eq(1, res.start);
    mu_assert_int_eq(3, res.mid);
    mu_assert_int_eq(5, res.end);


}


void run_all_tests() {
    MU_RUN_TEST(test_binary_search);
    MU_RUN_TEST(test_search_abc);
    MU_RUN_TEST(test_search_multiple_sa_constr);
}

MU_TEST_SUITE(fasta_parser_test_suite) {
    run_all_tests();
}

int main(int argc, char *argv[]) {
    init(argv[1]);
    MU_RUN_SUITE(fasta_parser_test_suite);
    MU_REPORT();
    return MU_EXIT_CODE;
}