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

void test_is_logical();

int simulation_RM_SCL();

int simulation_RM_CSS();

#endif