#include "minunit.h"
#include "../src/func/helper.h"
#include "testHelper.h"
#include "../src/func/rotater.h"

int lenMis = 12;

void test_setup(void) {

}

void test_teardown(void) {
    //Nothing
}


MU_TEST(lookup) {
    int alphabetSize = 3;
    int row1[3] = {1,0,0};
    int row2[3] = {2,0,0};
    int row3[3] = {2,1,0};
    int row4[3] = {2,1,1};
    int row5[3] = {3,1,2};

    //int a[6] = {1,2,2,2,2,3};
    //int b[6] = {0,0,1,1,1,1};
    //int c[6] = {0,0,0,1,2,2};
    int** o = malloc(5*sizeof(*o));
    o[0] = row1;
    o[1] = row2;
    o[2] = row3;
    o[3] = row4;
    o[4] = row5;

    mu_assert_int_eq(1, oLookUp(o, 3, 4)); //c4

    int i[12] = {1,1,1,1,1,1,1,2,2,2,3,4};
    int m[12] = {0,0,0,0,1,1,1,1,1,1,1,1};
    int p[12] = {0,1,1,1,1,1,2,2,2,2,2,2};
    int s[12] = {0,0,1,2,2,2,2,2,3,4,4,4};
    int **oMis = malloc(12*sizeof (*oMis));
    for(int j=0; j<12; j++) {
        int* row = malloc(4*sizeof (*row));
        row[0] = i[j];
        row[1] = m[j];
        row[2] = p[j];
        row[3] = s[j];
        oMis[j] = row;
    }
    mu_assert_int_eq(1, oLookUp(oMis, 3, 6));
    mu_assert_int_eq(2, oLookUp(oMis, 3, 7));
    mu_assert_int_eq(0, oLookUp(oMis, 2, 4));
    mu_assert_int_eq(1, oLookUp(oMis, 2, 5));
    mu_assert_int_eq(0, oLookUp(oMis, 0, 12));
    mu_assert_int_eq(0, oLookUp(oMis, 4, 0));
    mu_assert_int_eq(0, oLookUp(oMis, 4, 1));
}

MU_TEST(test_jump) {
    int cMis[5] = {0, 1, 5,6,8};
    int i[12] = {1,1,1,1,1,1,1,2,2,2,3,4};
    int m[12] = {0,0,0,0,1,1,1,1,1,1,1,1};
    int p[12] = {0,1,1,1,1,1,2,2,2,2,2,2};
    int s[12] = {0,0,1,2,2,2,2,2,3,4,4,4};
    int **oMis = malloc(12*sizeof (*oMis));
    for(int j=0; j<12; j++) {
        int* row = malloc(4*sizeof (*row));
        row[0] = i[j];
        row[1] = m[j];
        row[2] = p[j];
        row[3] = s[j];
        oMis[j] = row;
    }
    int jumpPoint = jump(0, 1, cMis, oMis);
    mu_assert_int_eq(1, jumpPoint);
}

MU_TEST(test_rotateString) {
    int i[12] = {1,1,1,1,1,1,1,2,2,2,3,4};
    int m[12] = {0,0,0,0,1,1,1,1,1,1,1,1};
    int p[12] = {0,1,1,1,1,1,2,2,2,2,2,2};
    int s[12] = {0,0,1,2,2,2,2,2,3,4,4,4};
    int** oMis = malloc(12*sizeof (*oMis));
    for(int j=0; j<12; j++) {
        int* row = malloc(4*sizeof (*row));
        row[0] = i[j];
        row[1] = m[j];
        row[2] = p[j];
        row[3] = s[j];
        oMis[j] = row;
    }


    int cMis[5] = {0, 1, 5,6,8};

    int ssi[3] = {4, 4, 1};
    struct Range* rotation = malloc(sizeof *rotation);
    rotateString(ssi, 3, cMis, oMis,  lenMis, rotation);
    mu_assert_int_eq(10, rotation->start);
    mu_assert_int_eq(12, rotation->end);

    int isi[3] = {1,4,1};
    rotateString(isi, 3, cMis, oMis, lenMis, rotation);
    mu_assert_int_eq(rotation->start, rotation->end);
    mu_assert_int_eq(3, rotation->start);
    mu_assert_int_eq(3, rotation->end);
}

void run_all_tests() {
    MU_RUN_TEST(test_jump);
    MU_RUN_TEST(lookup);
    MU_RUN_TEST(test_rotateString);
}

MU_TEST_SUITE(fasta_parser_test_suite) {
    MU_SUITE_CONFIGURE(&test_setup, &test_teardown);
    run_all_tests();
}

int main(int argc, char *argv[]) {
    init(argv[1]);
    MU_RUN_SUITE(fasta_parser_test_suite);
    MU_REPORT();
    return MU_EXIT_CODE;
}