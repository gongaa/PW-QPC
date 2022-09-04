#ifndef TEST_RM_HPP_
#define TEST_RM_HPP_
#include "Encoder/Encoder.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder.hpp"
#include "Decoder/Decoder_RM_SC.hpp"

using namespace std;

// generate codeword or information bits
void generate_random(int N, int *Y_N);

// verify the encoded and decoded info bits are the same
bool verify(int K, int *U_K_1, int *U_K_2);

#endif