#include "minunit.h"
#include "testHelper.h"
#include "../src/func/approx.h"
#include "../src/func/sa.h"
#include "../src/func/helper.h"
#include "../src/func/debugger.h"

MU_TEST(test_makeDeq) {
    int* C = calloc(5, sizeof *C);
    int** RO = malloc(5*sizeof *RO);
    int nX = 5;
    int bwt[5] = {1,2,2,1,0}; //bwt of 21210
    makeOandC(bwt, nX, RO, C);
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
    makeOandC(bwt, nX, RO, C);

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
    makeOandC(bwt, nX, RO, C);
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
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(1*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, 3);
    int** O = malloc(4*sizeof *RO);
    int bwt[3] = {1,1,0}; //bwt of 110
    makeOandC(bwt, n, O, C);

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
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(2*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, 5);
    int** O = malloc(6*sizeof *RO);
    int bwt[5] = {2,2,0,1,1}; //bwt of 12120
    makeOandC(bwt, n, O, C);

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
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(1*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, 5*sizeof (*C));
    int** O = malloc(4*sizeof *RO);
    int bwt[3] = {2,2,0}; //bwt of 220
    makeOandC(bwt, n, O, C);

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

MU_TEST(test_runApproxForceD) {
    int pattern[4] = {1,1,1,1};
    int n = 6; // 11211
    int m = 4; // 1111

    //D
    int* C = calloc(n+1, sizeof *C);
    int** RO = malloc(n*sizeof *RO);
    int rbwt[6] = {1,1,2,0,1,1}; //bwt of 112110
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(m*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, n*sizeof (*C));
    int** O = malloc(n*sizeof *RO);
    int bwt[6] = {1,1,2,0,1,1}; //bwt of 112110
    makeOandC(bwt, n, O, C);

    int allowedEdits=1;
    char* editString = malloc(m+allowedEdits+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);

    mu_assert_int_eq(3, amc->amount);
    mu_assert_int_eq(5, amc->AMs[0]->editStringLen);
    mu_assert_string_eq("MMDMM", amc->AMs[0]->editString);
}

MU_TEST(test_runApprox2EditsSmall) {
    int pattern[2] = {1,1};
    int n = 3;
    int m = 2;

    //D
    int* C = calloc(n+1, sizeof *C);
    int** RO = malloc(n*sizeof *RO);
    int rbwt[3] = {2,2,0}; //bwt of 220
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(m*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, n*sizeof (*C));
    int** O = malloc(n*sizeof *RO);
    int bwt[3] = {2,2,0}; //bwt of 220
    makeOandC(bwt, n, O, C);

    int allowedEdits=2;
    char* editString = malloc(m+allowedEdits+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);

    mu_assert_int_eq(4, amc->amount);
    mu_assert_string_eq("II", amc->AMs[0]->editString);
    mu_assert_string_eq("IM", amc->AMs[1]->editString);
    mu_assert_string_eq("MI", amc->AMs[2]->editString);
    mu_assert_string_eq("MM", amc->AMs[3]->editString);

    mu_assert_int_eq(3, amc->AMs[0]->rEnd-amc->AMs[0]->rStart);
    mu_assert_int_eq(2, amc->AMs[1]->rEnd-amc->AMs[1]->rStart);
    mu_assert_int_eq(2, amc->AMs[2]->rEnd-amc->AMs[2]->rStart);
    mu_assert_int_eq(1, amc->AMs[3]->rEnd-amc->AMs[3]->rStart);
}

MU_TEST(test_runApproxAsym) {
    int pattern[3] = {1,2,3};
    int n = 4;
    int m = 3;

    //D
    int* C = calloc(n+1, sizeof *C);
    int** RO = malloc(n*sizeof *RO);
    int rbwt[4] = {1,2,4,0}; //bwt of 4210
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(m*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, n*sizeof (*C));
    int** O = malloc(n*sizeof *RO);
    int bwt[4] = {4,0,1,2}; //bwt of 1240
    makeOandC(bwt, n, O, C);

    int allowedEdits=1;
    char* editString = malloc(m+allowedEdits+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);

    mu_assert_int_eq(2, amc->amount);
    mu_assert_string_eq("IMM", amc->AMs[0]->editString);
    mu_assert_string_eq("MMM", amc->AMs[1]->editString);

    mu_assert_int_eq(1, amc->AMs[0]->rEnd-amc->AMs[0]->rStart);
    mu_assert_int_eq(1, amc->AMs[1]->rEnd-amc->AMs[1]->rStart);
}

MU_TEST(test_runApproxAsym2) {
    int pattern[6] = {1,2,3,4,1,2};
    int n = 8;
    int m = 6;

    //D
    int* C = calloc(n+1, sizeof *C);
    int** RO = malloc(n*sizeof *RO);
    int rbwt[8] = {1,2,3,2,1,0,4,1}; //bwt of 21431210
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(m*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, n*sizeof (*C));
    int** O = malloc(n*sizeof *RO);
    int bwt[8] = {2,4,0,2,1,1,1,3}; //bwt of 12134120
    makeOandC(bwt, n, O, C);

    int allowedEdits=1;
    char* editString = malloc(m+allowedEdits+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);

    mu_assert_int_eq(2, amc->amount);
    mu_assert_string_eq("MMMMIM", amc->AMs[0]->editString);
    mu_assert_string_eq("MMMMDMM", amc->AMs[1]->editString);

    mu_assert_int_eq(1, amc->AMs[0]->rEnd-amc->AMs[0]->rStart);
    mu_assert_int_eq(1, amc->AMs[1]->rEnd-amc->AMs[1]->rStart);
}

MU_TEST(test_runApprox2Edits) {
    int pattern[4] = {3,1,2,1};
    int n = 7;
    int m = 4;

    //D
    int* C = calloc(n+1, sizeof *C);
    int** RO = malloc(n*sizeof *RO);
    int rbwt[7] = {3,4,0,4,1,2,1}; //bwt of 1424130
    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(m*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    //O and C
    memset(C, 0, n*sizeof (*C));
    int** O = malloc(n*sizeof *O);
    int bwt[7] = {1,4,3,4,0,2,1}; //bwt of 3142410
    makeOandC(bwt, n, O, C);

    int allowedEdits=2;
    char* editString = malloc(m+allowedEdits+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);

    mu_assert_int_eq(8, amc->amount);
    mu_assert_string_eq("IIMM", amc->AMs[0]->editString);
    mu_assert_string_eq("IMDMM", amc->AMs[1]->editString);
    mu_assert_string_eq("IMMM", amc->AMs[2]->editString);
    mu_assert_string_eq("MIIM", amc->AMs[3]->editString);
    mu_assert_string_eq("MDMDMM", amc->AMs[4]->editString);
    mu_assert_string_eq("MMMM", amc->AMs[5]->editString);
    mu_assert_string_eq("MIMM", amc->AMs[6]->editString);
    mu_assert_string_eq("MMDMM", amc->AMs[7]->editString);
}

MU_TEST(test_saAlgSmall) {
    struct Fasta* fasta = malloc(sizeof *fasta);
    fasta->fasta_head = "head";

    char * mal = malloc(sizeof(*mal)*4);
    char* seq = "acg";
    strcpy(mal, seq);
    update_fasta_by_sequence(&mal, fasta);
    mu_assert_int_eq(4, fasta->fasta_len);
    mu_assert_int_eq(4, fasta->alphabet.size);
    mu_assert_int_eq(1, *fasta->alphabet.sightings);
    mu_assert_int_eq(1, fasta->alphabet.sightings[1]);
    mu_assert_int_eq(1, fasta->alphabet.sightings[2]);
    mu_assert_int_eq(1, fasta->alphabet.sightings[3]);

    int* sa = constructSA(*fasta, 0);
    int* revSa = constructSA(*fasta, 1);

    struct Fasta* fasta2 = malloc(sizeof *fasta);
    fasta->fasta_head = "head";

    char * mal2 = malloc(sizeof(*mal2)*4);
    char* seq2 = "gca";
    strcpy(mal2, seq2);
    update_fasta_by_sequence(&mal2, fasta2);
    mu_assert_int_eq(4, fasta->fasta_len);
    mu_assert_int_eq(4, fasta->alphabet.size);
    mu_assert_int_eq(1, *fasta2->alphabet.sightings);
    mu_assert_int_eq(1, fasta2->alphabet.sightings[1]);
    mu_assert_int_eq(1, fasta2->alphabet.sightings[2]);
    mu_assert_int_eq(1, fasta2->alphabet.sightings[3]);

    int revSaExp[4] = {3,2,1,0};

    int* sa2 = constructSA(*fasta2, 0);
    mu_assert_int_arr_eq(revSaExp, sa2);

    int saExp[4] = {3,0,1,2};
    //012345678
    //214213210
    mu_assert_int_arr_eq(saExp, sa);
    mu_assert_int_arr_eq(revSaExp, revSa);
    //free(mal);
    free(sa);
    //free(revSa);
}


MU_TEST(test_saAlg) {
    struct Fasta* fasta = malloc(sizeof *fasta);
    fasta->fasta_head = "head";
    //char* seq = "12312412";

    char * mal = malloc(sizeof(*mal)*9);
    char* seq = "abcabdab";
    strcpy(mal, seq);
    update_fasta_by_sequence(&mal, fasta);
    mu_assert_int_eq(9, fasta->fasta_len);
    mu_assert_int_eq(5, fasta->alphabet.size);
    mu_assert_int_eq(1, *fasta->alphabet.sightings);
    mu_assert_int_eq(3, fasta->alphabet.sightings[1]);
    mu_assert_int_eq(3, fasta->alphabet.sightings[2]);
    mu_assert_int_eq(1, fasta->alphabet.sightings[3]);
    mu_assert_int_eq(1, fasta->alphabet.sightings[4]);

    int* sa = constructSA(*fasta, 0);
    int* revSa = constructSA(*fasta, 1);

    //012345678
    //123124120
    int saExp[9] = {8, 6, 0,3,7,1,4,2,5};
    //012345678
    //214213210
    int revSaExp[9] = {8,7,4,1,6,3,0,5,2};
    mu_assert_int_arr_eq(saExp, sa);
    mu_assert_int_arr_eq(revSaExp, revSa);
    //free(mal);
    free(sa);
    free(revSa);
}

MU_TEST(test_runApproxExactMis) {
    int pattern[3] = {2,1,4};
    int n = 12; // 214414413310
    int m = 3; // 214

    //D
    int* C = calloc(n+1, sizeof *C);
    int** RO = malloc(n*sizeof *RO);

    struct Fasta* fasta = malloc(sizeof *fasta);
    fasta->fasta_head = "head";

    char * mal = malloc(sizeof(*mal)*n);
    char* seq = "21441441331";
    strcpy(mal, seq);
    update_fasta_by_sequence(&mal, fasta);
    int* sa = constructSA(*fasta, 0);
    int* revSa = constructSA(*fasta, 1);

    int* rbwt = malloc(n*sizeof *rbwt); //bwt of 133144144120
    int* bwt = malloc(n*sizeof *bwt); //bwt of 214414413310

    for(int i=0; i<n; i++) {
        rbwt[i] = revSa[i] ? fasta->fasta_sequence[n-revSa[i] - 1] : 0;
        bwt[i] = sa[i] ? fasta->fasta_sequence[sa[i] - 1] : 0;
    }
   /*
    * 012345678901
    * 133144144120
    * 0 -> 2
    * 1 :
    *   120      -> 4
    *   1331...  -> 0
    *   144120   -> 4
    *   144144120-> 3
    * 2 -> 1
    * 3 :
    *   31... -> 3
    *   33... -> 1
    * 4:
    *   4120    -> 4
    *   4144120 -> 4
    *   44120   -> 1
    *   44144120-> 1
    **/
    int expRevSa[12] = {11, 9, 0, 6,3,10,2,1,8,5,7,4};
    int expRbwt[12] = {2,4,0,4,3, 1, 3,1,4,4,1,1}; // 133144144120
    mu_assert_int_arr_eq(expRevSa, revSa);
    mu_assert_int_eq(2, *rbwt);
    mu_assert_int_arr_eq(expRbwt, rbwt);

    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(m*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);
    mu_assert_int_eq(0, D[0]);
    mu_assert_int_eq(0, D[1]);
    mu_assert_int_eq(0, D[2]);

    //O and C
    memset(C, 0, n*sizeof (*C));
    int** O = malloc(n*sizeof *O);
    makeOandC(bwt, n, O, C);

    int allowedEdits=0;
    char* editString = malloc(m+allowedEdits+1);


    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);

    mu_assert_int_eq(1, amc->amount);
    mu_assert_string_eq("MMM", amc->AMs[0]->editString);
}


MU_TEST(test_nonConfirmed) {
    /* This is a test written based on an old version of the code
     * Since we assume that version is correct, we assume these tests can be used as regression tests
    */
    int pattern[10] = {3,4,1,1,2,4,1,3,2,1};
    int n = 100;
    int m = 10;

    //D
    int* C = calloc(n+1, sizeof *C);
    int** RO = malloc(n*sizeof *RO);
    struct Fasta* fasta = malloc(sizeof *fasta);
    char * mal = malloc(sizeof(*mal)*n);
    char* seq = "231212334124123412312234114212234123421342114233132222342134211324443124141141432112424432121322111";
    strcpy(mal, seq);
    update_fasta_by_sequence(&mal, fasta);

    int* sa = constructSA(*fasta, 0);
    int* revSa = constructSA(*fasta, 1);

    int* rbwt = malloc(n*sizeof *rbwt);
    int* bwt = malloc(n*sizeof *bwt);
    for(int i=0; i<n; i++) {
        rbwt[i] = revSa[i] ? fasta->fasta_sequence[n-revSa[i] - 1] : 0;
        bwt[i] = sa[i] ? fasta->fasta_sequence[sa[i] - 1] : 0;
    }

    int expSa[100] = {99,98,97,96,81,61,74,24,42,2,90,19,28,16,4,12,33,9,69,82,92,48,62,57,38,72,75,25,43,77,95,80,60,41,89,27,3,91,56,37,94,50,51,20,29,52,0,17,45,5,21,13,30,53,34,10,70,83,85,64,1,18,68,47,79,88,93,49,63,46,6,22,14,
            31,7,58,39,54,35,73,23,15,11,32,8,71,76,59,40,26,55,36,44,84,67,78,87,66,86,65};
    int expRevSa[100] = {99,0,1,16,55,36,23,73,6,94,2,17,8,70,56,37,60,41,96,79,50,29,21,24,26,86,82,65,74,89,98,15,93,7,69,95,78,28,85,81,64,88,68,77,45,46,3,47,34,4,48,18,9,53,71,13,57,38,61,42,35,5,59,40,49,97,92,84,80,63,67,76,44,52,91,51,19,10,30,54,22
            ,72,20,25,14,27,87,33,12,58,39,83,62,66,75,43,90,32,11,31};
    int expbwt[100] = {1,1,1,2,2,2,4,4,2,3,2,3,2,4,2,4,4,4,3,1,2,3,1,2,2,4,1,1,1,4,2,3,4,4,3,4,1,1,4,4,3,3,2,1,1,2,0,1,4,1,2,1,2,2,1,1,1,1
            ,4,3,2,2,4,3,4,4,1,1,1,2,2,2,2,2,3,1,1,2,2,1,3,3,2,3,3,2,1,3,3,1,3,3,1,2,4,1,4,4,2,2};
    int expRbwt[100] = {2,0,1,2,4,3,4,4,3,2,1,1,2,2,1,1,3,3,2,2,3,2,4,1,4,2,2,2,1,2,3,4,3,1,2,1,2,4,3,3,3,4,3,3,3,2,1,2,4,2,2,1,1,3,1,4,1,
                        1,1,1,2,2,4,4,2,1,3,4,1,4,4,4,4,3,4,1,2,2,1,2,1,2,3,1,2,1,1,4,4,2,2,1,2,1,1,2,1,4,3,3};

    mu_assert_int_arr_eq(expSa, sa);
    mu_assert_int_arr_eq(expRevSa, revSa);
    mu_assert_int_arr_eq(expbwt, bwt);
    mu_assert_int_arr_eq(expRbwt, rbwt);


    makeOandC(rbwt, n, RO, C);
    struct Range* r = malloc(sizeof *r);
    int* D = malloc(m*sizeof *D);
    makeD(D, C, RO, pattern, n, m, r);

    int expD[10] = {0,0,0,0,1,1,1,2,2,2};
    mu_assert_int_arr_eq(expD, D);

    //O and C
    memset(C, 0, n*sizeof (*C));
    int** O = malloc(n*sizeof *O);
    makeOandC(bwt, n, O, C);

    int allowedEdits=0;
    char* editString = malloc(m+allowedEdits+1);

    struct ApproxMatchContainer* amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);

    mu_assert_int_eq(0, amc->amount);

    allowedEdits = 1;
    amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);
    mu_assert_int_eq(0, amc->amount);
    allowedEdits = 2;
    amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);
    mu_assert_int_eq(0, amc->amount);

    allowedEdits = 3;
    amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);
    mu_assert_int_eq(9, amc->amount);

    int editStringLen3[9] = {10,10,10,11,11,11,11,10,10};
    char* editStrings3[9] = {"IMIMMMIMMM","IMIMMMMIMM","MMIIMIMMMM","MMMDMMIMMMI","MMMDMMIMMMM","MMMDMMMIMMM","MMMDMMMMIMM","MMIMMMIMMM","MMIMMMMIMM"};
    int rStart3[9] = {74,74,71,79,25,74,74,74,74};
    int rEnd3[9] = {75,75,72,80,26,75,75,75,75};



    for(int i=0; i<9; i++) {
        mu_assert_int_eq(editStringLen3[i], amc->AMs[i]->editStringLen);
        mu_assert_int_eq(rStart3[i], amc->AMs[i]->rStart);
        mu_assert_int_eq(rEnd3[i], amc->AMs[i]->rEnd);
        mu_assert_string_eq(editStrings3[i], amc->AMs[i]->editString);
    }



    allowedEdits = 4;
    amc = runApprox(pattern, n, m, D, C, O, allowedEdits, editString, r);
    mu_assert_int_eq(190, amc->amount);
    int editStringLen4[190] = {10,10
            ,10,10,11,11,11,11,10,10,11,11,10,10,10,10,10,11,10,10,11,10,10,10,10,10,10,11,11,10,10,10,11,11,10,10,10,10,10,10,10,11
            ,11,10,10,10,12,12,10,10,10,10,10,10,10,12,12,11,11,11,12,12,12,12,11,11,10,10,10,10,10,10,10,11,11,10,10,10,10,10,10,10
            ,11,11,11,11,11,11,11,11,13,13,12,10,10,10,10,11,11,11,11,12,11,11,11,11,11,12,12,11,11,11,11,11,12,12,12,11,11,12,10,10
            ,12,10,10,12,10,10,10,10,10,10,10,11,10,10,10,10,10,11,11,11,11,12,11,11,12,11,11,11,11,11,10,10,11,11,10,10,10,10,10,10
            ,10,10,10,10,10,11,10,10,11,10,10,10,10,10,11,11,10,10,10,11,11,11,11,10,10,10,10,10};
    char* editStrings4[190] = {"IIIMMMIMMM","IIIMMMMIMM","IIMMMMIMMM","IIMMMMMIMM","IIMDMMMIMMM","IIMDMMMMIMM","IDMIMMMIMMM","IDMIMMMMIMM","IMIIMIMMMM","IMIIMMMMMM","IMIMIMDMMMM","IMIMDMIMMMM","IMIMMMMMMM","IMIMMMIMMM","IMIMMMMIMM","IMIMMIMMMM","IMIMMMIMMI","IMIMMMIMMDM","IMIMMMIMMM","IMIMMMMIMI","IMIMMMMIMDM","IMIMMMMIMM","IMIMMMMMIM","IMIMMMMMMM","IMMIMMIMMM","IMMIMMMIMM","IMMMMIMMMM","IMMMIMDMMMM","IMMMDMIMMMM","IMMMMMMMMM","IMMMMMIMMM","IMMMMMMIMM","IMMDMMIMMMI","IMMDMMIMMMM","IMMMIIMMMM","IMMMIMIMMM","IMMMIMMIMM","IMMMIMMMIM","IMMMMMIIMM","IMMMMMIMMM","IMMMMMMIMM","IMMDMMMIMMM","IMMDMMMMIMM","MIIIMMIMMM","MIIIMMMIMM","MIIMMIMMMM","MIDMDMMMIMMM","MIDMDMMMMIMM","MIMIMIMMMM","MIMIIMIMMM","MIMIIMMIMM","MIMMMMIMIM","MIMMMMMIIM","MIMMMMMMMI","MIMMMMMMMM","MDMMDMMIMMMI","MDMMDMMIMMMM","MDMMMIMIMMM","MDMMMIMMIMM","MDMMMIMMMIM","MDIMDMMMIMMM","MDIMDMMMMIMM","MDDMIMMMIMMM","MDDMIMMMMIMM","MDMMMMMIMMM","MDMMMMMMIMM","MMMMMMIMIM","MMMMMMMIIM","MMMMMMMMMI","MMMMMMMMMM","MMIIIMMMMM","MMIIMIMMMI","MMIIMIMMMM","MMIIMDMIMMM","MMIIMDMMIMM","MMIIMMIMMM","MMIIMMMIMM","MMIIMMMMMM","MMIIMMMMMM","MMIMIIMMMM","MMIMMMIMMM","MMIMMMMIMM","MMDMMIMIMMM","MMDMMIMMIMM","MMDMMIMMMIM","MMDMMMIMMMI","MMDMMMIMMMM","MMDMIIMIMMM","MMDMIIMMIMM","MMDMMIMIIMM","MMDMMDMDMIMMM","MMDMMDMDMMIMM","MMDMMDMMMMMM","MMMIMIIMMI","MMMIMIIMMM","MMMIMIMIMI","MMMIMIMIMM","MMMIMDMMMMI","MMMIMDMMMMM","MMMIMMDMMMI","MMMIMMDMMMM","MMMDDMMMMMMM","MMMDMIMMMMI","MMMDMIMMMMM","MMMDMMIMMII","MMMDMMIMMIM","MMMDMMIMMMI","MMMDMMIMMMDI","MMMDMMIMMMDM","MMMDMMIMMMM","MMMDMMMIMMI","MMMDMMMIMMM","MMMDMMMMIMI","MMMDMMMMIMM","MMMDMMMMDMMI","MMMDMMMMDMMM","MMMDMDMMMMMM","MMMMDMIMMMI","MMMMDMIMMMM","MMMMDDMMMMMM","MMMMMMMMMI","MMMMMMMMMM","MMMMMDDMMMMM","MMMMMMIIMI","MMMMMMIIMM","MMMMMMDDMMMM","MMMIIIMMMM","MMMIMMIMMM","MMMIMMMIMM","MMMMIIIMMM","MMMMIIMIMM","MMMMIMIMMM","MMMMIMMIMM","MMMMMIMMDMM","MMMMMMMMMM","MMIMMMIMIM","MMIMMMMIIM","MMIMMMMMMI","MMIMMMMMMM","MMDMMMMIMMM","MMDMMMMMIMM","MMMDMMIMMMM","MMMDMMMIMMI","MMMDMMMIMMDM","MMMDMMMIMMM","MMMDMMMMIMI","MMMDMMMMIMDM","MMMDMMMMIMM","MMMDMMMMMIM","MMMDMMMMMMM","MMMMDMMIMMM","MMMMDMMMIMM","MIIMMMIMMM","MIIMMMMIMM","MMIMIMDMMMM","MMIMDMIMMMM","MMIMMMMMMM","MMIMMMIMMM","MMIMMMMIMM","MMMMMIMMMM","MMMMIMIMMM","MMMMIMMIMM","MMMMIMMMIM","MIMMMMIMMM","MIMMMMMIMM","MMIMMIMMMM","MMIMMMIMMI","MMIMMMIMMDM","MMIMMMIMMM","MMIMMMMIMI","MMIMMMMIMDM","MMIMMMMIMM","MMIMMMMMIM","MMIMMMMMMM","MMMIMMIMMM","MMMIMMMIMM","MMMMIMDMMMM","MMMMDMIMMMM","MMMMMMMMMM","MMMMMMIMMM","MMMMMMMIMM","MIMDMMMIMMM","MIMDMMMMIMM","MDMIMMMIMMM","MDMIMMMMIMM","MMIIMMMMMM","MMMMIIMMMM","MMMMMMIIMM","MMMMMMIMMM","MMMMMMMIMM"
    };
    int rStart4[190] = {74,74
            ,74,74,74,74,74,74,71,64,71,71,71,72,72,74,84,70,74,84,70,74,74,70,74,74,71,71,71,71,72,72,79,25,75,65,65,65,75,74,74,74
            ,74,74,74,71,74,74,71,72,72,62,62,94,97,79,25,65,65,65,74,74,74,74,74,74,62,62,94,97,71,80,71,73,73,71,71,77,73,71,73,73
            ,65,65,65,79,25,73,73,77,73,73,73,86,26,86,26,79,25,79,25,64,79,25,6,79,79,25,85,25,79,25,79,25,85,56,64,79,25,64,79,25,
                        64,93,57,64,71,73,73,71,71,73,73,76,76,62,62,94,97,74,74,74,84,70,74,84,70,74,74,70,74,74,74,74,71,71,71,72,72,71,65,65,
                        65,74,74,74,84,70,74,84,70,74,74,70,74,74,71,71,71,72,72,74,74,74,74,64,75,75,74,74};
    int rEnd4[190] = {75,75
            ,75,75,75,75,75,75,72,65,72,72,72,73,73,75,85,71,75,85,71,75,75,71,75,75,72,72,72,72,73,73,80,26,76,66,66,66,76,75,75,75
            ,75,75,75,72,75,75,72,73,73,63,63,95,98,80,26,66,66,66,75,75,75,75,75,75,63,63,95,98,72,81,72,74,74,72,72,79,74,72,74,74
            ,66,66,66,80,26,74,74,79,74,74,74,87,27,87,27,80,26,80,26,65,80,26,7,80,80,26,86,26,80,26,80,26,86,57,65,80,26,65,80,26,
                      65,94,58,65,72,74,74,72,72,74,74,77,77,63,63,95,98,75,75,75,85,71,75,85,71,75,75,71,75,75,75,75,72,72,72,73,73,72,66,66,
                      66,75,75,75,85,71,75,85,71,75,75,71,75,75,72,72,72,73,73,75,75,75,75,65,76,76,75,75};

    for(int i=0; i<190; i++) {
        mu_assert_int_eq(editStringLen4[i], amc->AMs[i]->editStringLen);
        mu_assert_int_eq(rStart4[i], amc->AMs[i]->rStart);
        mu_assert_int_eq(rEnd4[i], amc->AMs[i]->rEnd);
        mu_assert_string_eq(editStrings4[i], amc->AMs[i]->editString);
    }
}


void compareSAExpected(int* expected, struct Fasta fasta, int reverse) {
    int* saRadix = constructSARadix(fasta, reverse);
    int* saPrefixDoubling = constructSAPrefixDoubling(fasta, reverse);
    int* saMain = constructSA(fasta, reverse);
    mu_assert_int_arr_eq(expected, saRadix);
    mu_assert_int_arr_eq(expected, saPrefixDoubling);
    mu_assert_int_arr_eq(expected, saMain);

}

void compareSA(struct Fasta fasta, int reverse) {
    int* saRadix = constructSARadix(fasta, reverse);
    int* saPrefixDoubling = constructSAPrefixDoubling(fasta, reverse);
    int* saMain = constructSA(fasta, reverse);
    mu_assert_int_arr_eq(saRadix, saPrefixDoubling);
    mu_assert_int_arr_eq(saRadix, saMain);
}

struct Fasta* makeFasta(char* seq) {
    char** seqp = malloc(sizeof *seqp);
    seqp[0] = seq;
    //strcpy(seqp, seq);
    struct Fasta* fasta = malloc(sizeof *fasta);
    update_fasta_by_sequence(seqp, fasta);
    return fasta;
}


MU_TEST(test_saConstructionKnow) {
    struct Fasta* fasta = malloc(sizeof *fasta);
/*
    char * malsmall = malloc(sizeof(*malsmall)*3);
    char* seqsmall = "AA";
    strcpy(malsmall, seqsmall);
    update_fasta_by_sequence(&malsmall, fasta);
    int expsmall[3] = {2,1,0};
    int expsmallRev[3] = {2,1,0};
    compareSAExpected(expsmall, *fasta, 0);
    compareSAExpected(expsmallRev, *fasta, 1);*/

    char * mal0 = malloc(sizeof(*mal0)*5);
    char* seq0 = "ACAC";
    strcpy(mal0, seq0);
    update_fasta_by_sequence(&mal0, fasta);
    int exp0[5] = {4,2,0,3,1};
    int exp0Rev[5] = {4,3,1,2,0};
    compareSAExpected(exp0, *fasta, 0);
    compareSAExpected(exp0Rev, *fasta, 1);

/*
    char * mal1 = malloc(sizeof(*mal1)*9);
    char* seq1 = "ACGTACGT";
    strcpy(mal1, seq1);
    update_fasta_by_sequence(&mal1, fasta);
    int exp1[9] = {8,4,0,5,1,6,2,7,3};
    int exp1Rev[9] = {8,7,3,6,2,5,1,4,0};
    compareSAExpected(exp1, *fasta, 0);
    compareSAExpected(exp1Rev, *fasta, 1);

    char * mal2 = malloc(sizeof(*mal2)*12);
    char* seq2 = "CATTATTAGGA";
    strcpy(mal2, seq2);
    update_fasta_by_sequence(&mal2, fasta);
    int exp2[12] = {11,10,7,4,1,0,9,8,6,3,5,2};
    compareSAExpected(exp2, *fasta, 0);
    freeFasta(fasta);*/

}

MU_TEST(test_saConstructionRandomSeed) {
    srand(1);
    int n = 1000;
    char** seqp = malloc(sizeof *seqp);
    char* seq = malloc((n+1)*sizeof *seq);
    seqp[0] = seq;
    struct Fasta* fasta = malloc(sizeof *fasta);
    for(int len=1; len<=n; len++) {
        for(int i=0; i<len; i++) {
            seq[i] = (rand() % 4) + 'A';
        }
        seq[len] = '\0';
        update_fasta_by_sequence(seqp, fasta);

        compareSA(*fasta, 0);
        compareSA(*fasta, 1);
    }
    freeFasta(fasta);
}

MU_TEST(test_saConstructionRandom) {
    int n = 1000;
    char** seqp = malloc(sizeof *seqp);
    char* seq = malloc((n+1)*sizeof *seq);
    struct Fasta* fasta = malloc(sizeof *fasta);
    seqp[0] = seq;
    for(int len=1; len<=n; len++) {
        for(int i=0; i<len; i++) {
            seq[i] = (rand() % 4)+'A';
        }
        seq[len] = '\0';
        update_fasta_by_sequence(seqp, fasta);
        compareSA(*fasta, 0);
        compareSA(*fasta, 1);
    }
    freeFasta(fasta);
}

MU_TEST(test_radixGetByte) {
    uint64_t* res = malloc(8*sizeof *res);
    for(int i=0; i<8; i++) {
        res[i] = getByte((uint64_t)2<<32 | 1, i);
    }
    uint64_t exp[8] = {1, 0,0,0, 2, 0,0,0};
    mu_assert_int_arr_eq(exp, res);

    for(int i=0; i<8; i++) {
        res[i] = getByte((uint64_t)4<<32 | 4, i);
    }
    uint64_t exp2[8] = {4, 0,0,0, 4, 0,0,0};
    mu_assert_int_arr_eq(exp2, res);

}

uint64_t com(uint32_t a, uint32_t b) {
    return (uint64_t)a<<32 | b;
}

MU_TEST(test_radixSort64Interval) {

    int sa[2] = {0,1};
    // Keys: 14, 23
    uint64_t keys[2] = { (uint64_t)1<<32|4, (uint64_t)2<<32 | 3};
    radixSort64Interval(0,1, sa, keys);
    int exp[2] = {0,1};
    mu_assert_int_arr_eq(exp, sa);


    int sa2[5] = {0,1,2,3,4};
    // Keys: 21, 14,44,41,14
    uint64_t keys2[5] = {(uint64_t)2<<32 | 1, (uint64_t)1<<32|4, (uint64_t)4<<32|4, (uint64_t)4<<32|1, (uint64_t)1<<32|4};
    radixSort64Interval(0,4, sa2, keys2);
    int exp2[5] = {1,4, 0,3,2 };
    mu_assert_int_arr_eq(exp2, sa2);


    int sa3[14] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    // 25, 95, 75, 22, 41, 16, 17, 51, 42, 94, 64, 82, 16, 41
    uint64_t keys3[14] = {com(2,5), com(9,5), com(7,5), com(2,2), com(4,1), com(1,6), com(1,7), com(5,1), com(4,2), com(9,4), com(6,4), com(8,2), com(1,6), com(4,1)};
    radixSort64Interval(0,13, sa3, keys3);
    // 25, 95, 75, 22, 41, 16, 17, 51, 42, 94, 64, 82, 16, 41
    //  0   1   2  3    4   5   6   7   8   9  10  11  12  13
    //  4  13  10  3    5   0   2   8   7   12 9   11   1   6
    int exp3[14] = {5,12,6,3,0,4,13,8,7,10, 2, 11,9,1};
    mu_assert_int_arr_eq(exp3, sa3);


    int sa4[5] = {0,1,2,3,4};
    // 25, 95, 75, 22, 41
    uint64_t keys4[5] = {com(2,5), com(9,5), com(7,5), com(2,2), com(4,1)};
    radixSort64Interval(1,3, sa4, keys4);
    int exp4[5] = {0,3,2,1,4};
    mu_assert_int_arr_eq(exp4, sa4);


    int sa5[14] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    // 25, 95, 75, 22, 41, 16, 17, 51, 42, 94, 64, 82, 16, 41
    uint64_t keys5[14] = {com(2,5), com(9,5), com(7,5), com(2,2), com(4,1), com(1,6), com(1,7), com(5,1), com(4,2), com(9,4), com(6,4), com(8,2), com(1,6), com(4,1)};
    radixSort64Interval(5,10, sa5, keys5);
    // 16, 17, 51, 42, 94, 64
    //  5   6   7   8   9  10
    //  5   6   8   7   10  9
    int exp5[14] = {0,1,2,3,4,5,6,8,7,10,9,11,12,13};
    mu_assert_int_arr_eq(exp5, sa5);



    int sa6[5] = {0,1,2,3,4};
    // ACAC
    // 12, 21, 12, 20, 00
    uint64_t keys6[5] = {com(1,2), com(2,1), com(1,2), com(2,0), com(0,0)};
    radixSort64Interval(0,4, sa6, keys6);
    int exp6[5] = {4,0,2,3,1};
    uint64_t exp6Key[5] = {com(0,0), com(1,2), com(1,2), com(2,0), com(2,1)};
    mu_assert_int_arr_eq(exp6, sa6);

    mu_assert_int_arr_eq(exp6Key, keys6);

}


MU_TEST_SUITE(fasta_parser_test_suite) {
    /*MU_RUN_TEST(test_makeDeq);
    MU_RUN_TEST(test_makeDnotEq);
    MU_RUN_TEST(test_makeD1Edit);
    MU_RUN_TEST(test_runApproxExactEqSmall);
    MU_RUN_TEST(test_runApproxExactEq);
    MU_RUN_TEST(test_runApproxDiff);
    MU_RUN_TEST(test_runApproxForceD);
    MU_RUN_TEST(test_runApprox2EditsSmall);
    MU_RUN_TEST(test_runApproxAsym);
    MU_RUN_TEST(test_runApproxAsym2);
    MU_RUN_TEST(test_runApprox2Edits);
    MU_RUN_TEST(test_saAlgSmall);
    MU_RUN_TEST(test_saAlg);
    MU_RUN_TEST(test_runApproxExactMis);
    MU_RUN_TEST(test_nonConfirmed);*/
    MU_RUN_TEST(test_saConstructionKnow);
    //MU_RUN_TEST(test_saConstructionRandom);
    //MU_RUN_TEST(test_saConstructionRandomSeed);
    MU_RUN_TEST(test_radixSort64Interval);
    MU_RUN_TEST(test_radixGetByte);
}


int main(int argc, char *argv[]) {
    // Path to the folder where test-data is located should be given.
    // (typically ../ if running from the debugger (project_root/cmake-build-debug/....
    // or ./ if running directly from project root)
    MU_RUN_SUITE(fasta_parser_test_suite);
    MU_REPORT();
    return MU_EXIT_CODE;
}