#ifndef TEST_POLAR_HPP_
#define TEST_POLAR_HPP_
#include "Encoder/Encoder.hpp"
#include "Encoder/Encoder_polar.hpp"
#include "Decoder/Decoder.hpp"
#include "Decoder/Decoder_polar.hpp"
#include "Util/Util.hpp"
using namespace std;

void test_crc();

void test_polar_stabilizer(int N, int K, double p);

void test_encoding_inverse();

void test_syndrome_extraction(int N, int K, bool print);

void test_noise_dist(int N, int num_total, int seed);

void test_generate_exact_flip(int N, double p, int num_total, int seed);

#endif // TEST_POLAR_HPP_