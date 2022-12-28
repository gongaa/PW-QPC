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

void test_syndrome_extraction(int N, int K, bool print=false);

void test_direct_syndrome_decoding(int N, int K, int list_size, double p);

void test_noise_dist(int N, int num_total=10000, int seed=42);

void test_generate_exact_flip(int N, double p, int num_total=10000, int seed=42);

void print_polar_con(int N, int K, CONSTRUCTION con);

#endif // TEST_POLAR_HPP_