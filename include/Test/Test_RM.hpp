#ifndef TEST_RM_HPP_
#define TEST_RM_HPP_
#include "Encoder/Encoder.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Encoder/Encoder_RM_CSS.hpp"
#include "Decoder/Decoder.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Util/Util.hpp"

using namespace std;
// #define CHN_AWGN
// #define IS_VERBOSE
// #define VANILLA_LIST_DECODING
#define DEGENERACY_LIST_DECODING
// #define CHECK

void generate_symmetric_noise(int m, vector<int>& noise, int level);
void generate_bent(int m, vector<int>& noise);
void generate_bent_guess(int m, vector<int>& noise);
void generate_bent_third_order(int m, vector<int>& noise);
void generate_bent_m6_r2_order(vector<int>& noise);
void generate_all_codewords(int m, int r, vector<vector<int>>& codewords);
void generate_all_equiv_classes(int m, int rx, int rz, vector<vector<int>>& codewords, vector<vector<int>>& equiv_classes);

void test_encoder_decode();
// check RM is_codeword
void verify_RM_is_codeword();

void verify_parity_check();

void test_is_logical();

void test_RM_SCL_symmetry();

void test_RM_syndrome_SC(int m, int r);

void test_linearity_xor(int m, int r);

#endif