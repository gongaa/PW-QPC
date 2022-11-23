#ifndef TEST_POLAR_HPP_
#define TEST_POLAR_HPP_
#include "Encoder/Encoder.hpp"
#include "Encoder/Encoder_polar.hpp"
#include "Decoder/Decoder.hpp"
#include "Decoder/Decoder_polar.hpp"
#include "Util/Util.hpp"
using namespace std;

void test_polar(int N, int K, int L, double p, double db, double design_snr); 

void test_crc();

void test_polar_stabilizer();

#endif // TEST_POLAR_HPP_