#ifndef TEST_RM_HPP_
#define TEST_RM_HPP_
#include "Encoder/Encoder.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Encoder/Encoder_RM_CSS.hpp"
#include "Decoder/Decoder.hpp"
#include "Decoder/Decoder_RM_SC.hpp"

using namespace std;

// generate codeword or information bits
void generate_random(int N, int *Y_N);

// verify the encoded and decoded info bits are the same
bool verify(int K, int *U_K_1, int *U_K_2);

// check RM is_codeword
void verify_RM_is_codeword();

void verify_parity_check();

void test_is_logical();

void test_RM_SCL_symmetry();

void test_RM_syndrome_SC();

int simulation_RM_SCL();

int simulation_RM_CSS(int m, int rx, int rz, int list_size);

inline void xor_vec(int N, const int* a, const int* b, int* c) {
    for (int i = 0; i < N; i++)
        c[i] = a[i] ^ b[i];
}

inline int count_flip(int N, const int* a, const int* b) {
    int cnt = 0;
    for (int i = 0; i < N; i++) 
        if (a[i] != b[i]) cnt++;
    return cnt;
}
#endif