#ifndef TEST_RM_HPP_
#define TEST_RM_HPP_
#include "Encoder/Encoder.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Encoder/Encoder_RM_CSS.hpp"
#include "Decoder/Decoder.hpp"
#include "Decoder/Decoder_RM_SC.hpp"

using namespace std;
// #define CHN_AWGN
// #define IS_VERBOSE
// #define VANILLA_LIST_DECODING
#define DEGENERACY_LIST_DECODING
// #define CHECK
// generate codeword or information bits
void generate_random(int N, int *Y_N);

void generate_symmetric_noise(int m, vector<int>& noise, int level);
void generate_bent(int m, vector<int>& noise);
void generate_bent_guess(int m, vector<int>& noise);
void generate_bent_third_order(int m, vector<int>& noise);
void generate_bent_m6_r2_order(vector<int>& noise);
void generate_all_codewords(int m, int r, vector<vector<int>>& codewords);
void generate_all_equiv_classes(int m, int rx, int rz, vector<vector<int>>& codewords, vector<vector<int>>& equiv_classes);

// verify the encoded and decoded info bits are the same
inline bool verify(int K, int *U_K_1, int *U_K_2) {
    bool is_same = true;
    for (int i = 0; i < K; i++) {
        if (U_K_1[i] != U_K_2[i]) {
            is_same = false;
            break;
        }
    }
    return is_same;
}

// check RM is_codeword
void verify_RM_is_codeword();

void verify_parity_check();

void test_is_logical();

void test_RM_SCL_symmetry();

void test_RM_syndrome_SC(int m, int r);

void test_linearity_xor(int m, int r);

inline void xor_vec(int N, int* a, const int* b, int* c) {
    for (int i = 0; i < N; i++)
        c[i] = a[i] ^ b[i];
}

inline int count_flip(int N, const int* a, const int* b) {
    int cnt = 0;
    for (int i = 0; i < N; i++) 
        if (a[i] != b[i]) cnt++;
    return cnt;
}

inline int count_weight(const vector<int>& a) {
    int cnt = 0;
    for (int i : a)
        if (i) cnt++;
    return cnt;
}

inline double db2val(double x) {
  return exp(log(10.0) * x / 10.0);
}

inline void decimal2bianry(const int& n, vector<int>& b)
{
    int x = n;
    for (int i = 0; x > 0; i++) {
        b[i] = x % 2;
        x >>= 1;
    }
}

inline int binary2decimal(const int *b, int size)
{
    int n = 0;
    for (int i = 0; i < size; i++) {
        if (b[i] == 1)
            n += (1 << i);
    }
    return n;
}
#endif