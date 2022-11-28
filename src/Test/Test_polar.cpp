#include <iostream>
#include <random>
#include <cassert>
#include "Test/Test_polar.hpp"
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_polar.hpp"
#include "Decoder/Decoder_polar.hpp"
#include "Channel/Channel.hpp"
#include "Util/CRC_polynomial.hpp"

using namespace std;

void test_polar_stabilizer() 
{
    // int K = 1536, N = 2048;
    int K = 50, N = 64, n = log2(N); 
    vector<bool> frozen_bits(N, 0);
    frozen_bits_generator_PW(N, K, frozen_bits);
    // assume X and Z are using the same code
    // expect K - (N-K) = 1024 stabilizers
    vector<int> Z_code(K, 0);
    vector<vector<int>> Z_code_bin(K, vector<int>(n, 0));
    int Z_cnt = 0;
    vector<int> X_stab(N-K, 0);
    vector<vector<int>> X_stab_bin(N-K, vector<int>(n, 0));
    int X_stab_cnt = 0;
    for (int i = 0; i < N; i++) {
        if (frozen_bits[i] == 0) { // is info bit, definite in Z_code
            Z_code[Z_cnt] = i;
            decimal2binary(i, Z_code_bin[Z_cnt]);
            Z_cnt++;
            if (frozen_bits[N-1-i] == 1) { // dual is also an info bit
                X_stab[X_stab_cnt] = i;
                decimal2binary(i, X_stab_bin[X_stab_cnt]);
                X_stab_cnt++;
            }
        }
    }
    assert (Z_cnt == K);
    assert (X_stab_cnt == (N-K));
    for (int i : Z_code) cerr << i << " ";
    cerr << endl;
    for (int i : X_stab) cerr << i << " ";
    cerr << endl;

    for (auto& cz : Z_code_bin) {
        for (auto& sx : X_stab_bin) {
            int i = 0;
            for (; i < n; i++) {
                if (cz[i] && sx[i]) break;
            }
            if (i == n) {
                cerr << "violation at cz=";
                for (int k : cz) cerr << k;
                cerr << " , sx=";
                for (int k : sx) cerr << k;
                cerr << endl;
            }
        }
    }

}

void test_crc() 
{
    int K = 34, crc_size = 8;
    // vector<int> info_bits(K);
    vector<int> info_bits_crc(K + crc_size);
    CRC_polynomial<int>* crc = new CRC_polynomial<int>(K, "8-WCDMA", crc_size);
    for (int idx = 0; idx < 100; idx++) {
        generate_random(K, info_bits_crc.data());
        crc->build(info_bits_crc.data(), info_bits_crc.data(), 0);
        cerr << "CRC remainder: ";
        for (int i = 0; i < crc_size; i++) cerr << info_bits_crc[K+i];
        cerr << endl;
        assert (crc->check(info_bits_crc.data(), 0));
    }
}

void test_encoding_inverse()
{
    // the inverse of the circuit is the circuit itself.
    int N = 1024;
    vector<int> input(N);
    vector<int> output(N);
    vector<int> outoutput(N);
    vector<bool> frozen_bits(N, 0);
    Encoder_polar* encoder = new Encoder_polar(N, N, frozen_bits);
    for (int idx = 0; idx < 10000; idx++) {
        generate_random(N, input.data());
        output = input;
        encoder->light_encode(output.data());
        outoutput = output;
        encoder->light_encode(outoutput.data());
        assert (verify(N, input.data(), outoutput.data()));
    }
}


void test_syndrome_extraction(int N, int K, bool print)
{
    // use the checks (the N-K rows of best channels) to create N-K syndromes
    // the i^th row is the codeword obtained by setting only the i^th position
    // as info and others as frozen.
    vector<bool> frozen_bits(N, 0);
    vector<bool> stab_frozen_bits(N, 1);
    vector<int>  syndromes(N-K, 0);
    vector<int>  syndromes_at_frozen(N-K, 0);
    vector<int>  input(N, 0);
    vector<int>  output(N, 0);
    vector<int>  noise(N, 0);
    vector<vector<int>> parity_checks(N-K, vector<int>(N,0));
    frozen_bits_generator_PW(N, K, frozen_bits);
    for (int i = 0; i < N; i++)
        if (frozen_bits[i] == 0 && frozen_bits[N-1-i] == 1)
            stab_frozen_bits[i] = 0;
    Encoder_polar* encoder = new Encoder_polar(N, N, frozen_bits);
    int j = 0;
    for (int i = 0; i < N; i++) {
        if (!stab_frozen_bits[i]) {
            input[i] = 1;
            output = input;
            encoder->light_encode(output.data());
            parity_checks[j++] = output;
            input[i] = 0;
            // cerr << "parity check at i=" << i << " is ";
            // for (int k : output) cerr << k;
            // cerr << endl;
        }
    }
    for (int idx = 0; idx < 10000; idx++) {
        generate_random(N, noise.data());
        j = 0;
        for (int i = 0; i < N-K; i++) {
            syndromes[i] = dot_product(N, parity_checks[i], noise);
        }
        if (print) {
            cerr << "syndromes                ";
            for (int k : syndromes) cerr << k;
        }

        // compared with bits at frozen after encoding the noise
        output = noise;
        // bit reversal
        bit_reversal(output); 
        encoder->light_encode(output.data());
        j = 0;
        for (int i = 0; i < N; i++) {
            if (frozen_bits[i]) {
                syndromes_at_frozen[j++] = output[i];
            }
        }
        // bit reversal again
        bit_reversal(syndromes_at_frozen);
        if (print) {
            cerr << endl << "syndromes at frozen bits ";
            for (int k : syndromes_at_frozen) cerr << k;
            cerr << endl;
        }

        assert (verify(N-K, syndromes.data(), syndromes_at_frozen.data()));
    }


}

/*
  from Arikan's paper: A Performance Comparison of Polar Codes and Reed-Muller Codes
  generate the information bit set
  compute the vector z[n] = (z[N][1],...,z[N][N]) through recursion
  z[1][1]  = 0.5
  z[2k][j] = 2*z[k][j] - z[k][j]**2     for 1<=j<=k
           = z[k][j-k]**2               for k+1<=j<=2k
  then choose the positions corresponding to the K highest ones
*/
/*
void calculate_frozen_bits_positions(int N, int K, double p, vector<bool>& frozen_bits) {
    // for BEC and BSC only
    vector<double> z(N, 0);
    vector<double> z_temp(N, 0);
    z[0] = p;
    for (int k = 1; k < N; k <<= 1) {
        for (int j = 0; j < k; j++) 
            z_temp[j] = 2 * z[j] - z[j] * z[j];
        for (int j = k; j < 2*k; j++) 
            z_temp[j] = z[j-k] * z[j-k];
        for (int j = 0; j < 2*k; j++) 
            z[j] = z_temp[j];
    }
    // for (auto i : z) cerr << i << " ";
    // cerr << endl;
    vector<std::tuple<int, double>> idx_z_pair(N);
    for (int i = 0; i < N; i++) idx_z_pair[i] = make_tuple(i, z[i]);
    std::sort(idx_z_pair.begin(), idx_z_pair.end(),
        [](std::tuple<int, double> x, std::tuple<int, double> y) {
            return std::get<1>(x) > std::get<1>(y);
        });
    for (auto p : idx_z_pair) cerr << std::get<0>(p)+1 << " ";
    cerr << endl;
    for (auto p : idx_z_pair) cerr << std::get<1>(p) << " ";
    cerr << endl;
    for (int i = 0; i < N; i++) frozen_bits[i] = 0;
    for (int i = 0; i < N-K; i++) frozen_bits[std::get<0>(idx_z_pair[i])] = 1;
}
*/

