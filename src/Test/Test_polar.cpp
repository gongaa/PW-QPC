#include <iostream>
#include <random>
#include <cassert>
#include "Test/Test_polar.hpp"
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_polar.hpp"
#include "Decoder/Decoder_polar.hpp"
#include "Channel/Channel.hpp"

using namespace std;

/*
  from Arikan's paper: A Performance Comparison of Polar Codes and Reed-Muller Codes
  generate the information bit set
  compute the vector z[n] = (z[N][1],...,z[N][N]) through recursion
  z[1][1]  = 0.5
  z[2k][j] = 2*z[k][j] - z[k][j]**2     for 1<=j<=k
           = z[k][j-k]**2               for k+1<=j<=2k
  then choose the positions corresponding to the K highest ones
*/

void calculate_frozen_bits_positions(int N, int K, vector<bool>& frozen_bits) {
    vector<double> z(N);
    vector<double> z_temp(N);
    z[0] = 0.5;
    for (int k = 1; k < N; k <<= 1) 
        for (int j = 0; j < k; j++) {
            z_temp[j] = 2 * z[j] - z[j] * z[j];
        for (int j = k; j < 2*k; j++) 
            z_temp[j] = z[j-k] * z[j-k];
        for (int j = 0; j < 2*k; j++) 
            z[j] = z_temp[j];
    }
    vector<std::tuple<int, double>> idx_z_pair(N);
    for (int i = 0; i < N; i++) idx_z_pair[i] = make_tuple(i, z[i]);
    std::sort(idx_z_pair.begin(), idx_z_pair.end(),
        [](std::tuple<int, double> x, std::tuple<int, double> y) {
            return std::get<1>(x) > std::get<1>(y);
        });
    for (int i = 0; i < N; i++) frozen_bits[i] = 0;
    for (int i = 0; i < N-K; i++) frozen_bits[std::get<0>(idx_z_pair[i])] = 1;
}

void test_polar() {
    int K = 2048, N = 4096, L=4;
    vector<bool> frozen_bits(N, 0);
    calculate_frozen_bits_positions(N, K, frozen_bits);
    // cerr << "frozen bits mask: ";
    // for (int i : frozen_bits) cerr << i;
    // cerr << endl;
    double p = 0.1;
    // vector<bool> frozen_bits = {1,0,1,0,1,0,0,0}; // info bits at 8,4,6,7,2
    Encoder_polar* encoder = new Encoder_polar(K, N, frozen_bits);
    Decoder_polar_SCL* decoder = new Decoder_polar_SCL(K, N, L, frozen_bits);
    Channel_BSC* chn_bsc = new Channel_BSC(N, p, 42);
    vector<int> info_bits(K);
    vector<int> codeword(N);
    vector<int> noisy_codeword(N);
    vector<double> llr_noisy_codeword(N, 0);
    vector<int> denoised_codeword(N);
    for (int k = 1; k < 1000; k++) {
        generate_random(K, info_bits.data());
        // cerr << "info bits: ";
        // for (int i : info_bits) cerr << i;
        // cerr << endl;
        encoder->encode(info_bits.data(), codeword.data(), 0);
        // cerr << "codeword: ";
        // for (int i : codeword) cerr << i;
        // cerr << endl << "is codeword? " << encoder->is_codeword(codeword.data()) << endl;
        assert(encoder->is_codeword(codeword.data()));
        // int ml_flips = chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
        // cerr << endl << ml_flips << " flips, noisy codeword: ";
        // for (int i : noisy_codeword) cerr << i;
        // cerr << endl;
        for (int i = 0; i < N; i++)
            llr_noisy_codeword[i] = codeword[i] ? -1.0 : 1.0;
            // llr_noisy_codeword[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
        decoder->decode(llr_noisy_codeword.data(), denoised_codeword.data(), 1);
        if (verify(N, codeword.data(), denoised_codeword.data()))
            cerr << k << " decode successfully" << endl;
        else 
            cerr << "decode fail, num_flips: " << count_flip(N, codeword.data(), denoised_codeword.data()) << endl;
    }
}