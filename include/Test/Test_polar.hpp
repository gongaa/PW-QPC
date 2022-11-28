#ifndef TEST_POLAR_HPP_
#define TEST_POLAR_HPP_
#include "Encoder/Encoder.hpp"
#include "Encoder/Encoder_polar.hpp"
#include "Decoder/Decoder.hpp"
#include "Decoder/Decoder_polar.hpp"
#include "Util/Util.hpp"
using namespace std;

void test_crc();

void test_polar_stabilizer();

void test_encoding_inverse();

void test_syndrome_extraction(int N, int K, bool print=false);

#endif // TEST_POLAR_HPP_