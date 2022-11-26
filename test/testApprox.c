#include <stdbool.h>
#include "minunit.h"
#include "testHelper.h"
#include "../src/func/approx.h"



MU_TEST(test_makeDeq) {
    int* D = malloc(4*sizeof *D);
    int* C = calloc(5, sizeof *C);
    int** RO = malloc(sizeof *RO);
    int nX = 5;
    int bwt[5] = {1,2,2,1,0}; //bwt of 21210
    makeOandC(bwt, nX, RO, C, 5);

    int pattern[4] = {1,2,1,2};
    int m = 4;
    struct Range* r = malloc(sizeof *r);
    makeD(D, C, RO, pattern, nX, m, r);
    int exp[4] = {0,0,0,0};
    mu_assert_int_eq(0, *D);
    mu_assert_int_eq(0, D[1]);
    mu_assert_int_eq(0, D[2]);
    mu_assert_int_eq(0, D[3]);
    mu_assert_int_arr_eq(exp, D);
}


void run_all_tests() {
    MU_RUN_TEST(test_makeDeq);
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