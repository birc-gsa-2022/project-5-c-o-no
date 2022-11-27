#include <stdbool.h>
#include "minunit.h"
#include "testHelper.h"
#include "../src/func/approx.h"



MU_TEST(test_makeDeq) {
    int* C = calloc(5, sizeof *C);
    int** RO = malloc(5*sizeof *RO);
    int nX = 5;
    int bwt[5] = {1,2,2,1,0}; //bwt of 21210
    makeOandC(bwt, nX, RO, C, 5);
    int pattern[4] = {1,2,1,2};
    int m = 4;
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(4*sizeof *D);
    makeD(D, C, RO, pattern, nX, m, r);
    int exp[4] = {0,0,0,0};
    mu_assert_int_arr_eq(exp, D);
}

MU_TEST(test_makeDnotEq) {
    int* C = calloc(5, sizeof *C);
    int** RO = malloc(5*sizeof *RO);
    int nX = 5;
    int bwt[5] = {1,2,2,1,0}; //bwt of 21210
    makeOandC(bwt, nX, RO, C, 5);

    int pattern[4] = {3,3,3,3};
    int m = 4;
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(4*sizeof *D);
    makeD(D, C, RO, pattern, nX, m, r);
    int exp[4] = {1,2,3,4};
    mu_assert_int_arr_eq(exp, D);
}

MU_TEST(test_makeD1Edit) {
    int* C = calloc(7, sizeof *C);
    int** RO = malloc(7*sizeof *RO);
    int nX = 7;
    int bwt[7] = {1,2,2,1,2,0,1}; //bwt of baabba0
    makeOandC(bwt, nX, RO, C, 5);
    int pattern[5] = {1,2,1,2,1};
    int m = 5;
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(5*sizeof *D);
    makeD(D, C, RO, pattern, nX, m, r);
    int exp[5] = {0,0,1,1,1};
    mu_assert_int_arr_eq(exp, D);
}

MU_TEST(test_runApproxExactEqSmall) {
    int pattern[1] = {1};
    int n = 3;
    int m = 1;

    //D
    int* C = calloc(3, sizeof *C);
    int** RO = malloc(4*sizeof *RO);
    int rbwt[3] = {1,1,0}; //bwt of 110
    makeOandC(rbwt, n, RO, C, 3);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(1*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, 3);
    int** O = malloc(4*sizeof *RO);
    int bwt[3] = {1,1,0}; //bwt of 110
    makeOandC(bwt, n, O, C, 3);

    char* editString = malloc(m+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, 0, editString, r);

    mu_assert_int_eq(0, *D);
    mu_assert_int_eq(1, amc->amount);
    mu_assert_int_eq(1, amc->AMs[0]->rStart);
    mu_assert_int_eq(3, amc->AMs[0]->rEnd);
    mu_assert_int_eq(1, amc->AMs[0]->editStringLen);
    mu_assert_string_eq("M", amc->AMs[0]->editString);
}

MU_TEST(test_runApproxExactEq) {
    int pattern[2] = {1,2};
    int n = 5;
    int m = 2;

    //D
    int* C = calloc(5, sizeof *C);
    int** RO = malloc(6*sizeof *RO);
    int rbwt[5] = {1,2,2,1,0}; //bwt of 21210
    makeOandC(rbwt, n, RO, C, 5);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(2*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, 5);
    int** O = malloc(6*sizeof *RO);
    int bwt[5] = {2,2,0,1,1}; //bwt of 12120
    makeOandC(bwt, n, O, C, 5);

    char* editString = malloc(m+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, 0, editString, r);

    mu_assert_int_eq(1, amc->amount);
    mu_assert_int_eq(1, amc->AMs[0]->rStart);
    mu_assert_int_eq(3, amc->AMs[0]->rEnd);
    mu_assert_int_eq(2, amc->AMs[0]->editStringLen);
    mu_assert_string_eq("MM", amc->AMs[0]->editString);
}

MU_TEST(test_runApproxDiff) {
    int pattern[1] = {1};
    int n = 3;
    int m = 1;

    //D
    int* C = calloc(5, sizeof *C);
    int** RO = malloc(4*sizeof *RO);
    int rbwt[3] = {2,2,0}; //bwt of 220
    makeOandC(rbwt, n, RO, C, 3);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(1*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, 5*sizeof (*C));
    int** O = malloc(4*sizeof *RO);
    int bwt[3] = {2,2,0}; //bwt of 220
    makeOandC(bwt, n, O, C, 3);

    int k=1;
    char* editString = malloc(m+k+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, k, editString, r);

    mu_assert_int_eq(1, *D);
    mu_assert_int_eq(2, amc->amount);
    mu_assert_int_eq(1, amc->AMs[0]->editStringLen);
    mu_assert_int_eq(1, amc->AMs[1]->editStringLen);
    mu_assert_string_eq("I", amc->AMs[0]->editString);
    mu_assert_string_eq("M", amc->AMs[1]->editString);
    mu_assert_int_eq(0, amc->AMs[0]->rStart);
    mu_assert_int_eq(3, amc->AMs[0]->rEnd);
    mu_assert_int_eq(1, amc->AMs[1]->rStart);
    mu_assert_int_eq(3, amc->AMs[1]->rEnd);
}


void run_all_tests() {
    MU_RUN_TEST(test_makeDeq);
    MU_RUN_TEST(test_makeDnotEq);
    MU_RUN_TEST(test_makeD1Edit);
    MU_RUN_TEST(test_runApproxExactEqSmall);
    MU_RUN_TEST(test_runApproxExactEq);
    MU_RUN_TEST(test_runApproxDiff);
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