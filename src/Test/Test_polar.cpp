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

void print_polar_con(int N, int K, CONSTRUCTION con)
{
    cerr << "N=" << N << ", K=" << K << endl;
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> X_code_frozen_bits(N, 0);
    vector<bool> X_stab_frozen_bits(N, 1);
    vector<int>  X_stab_info_indices; // expect to have size 2*K-N
    construct_frozen_bits(con, N, K, Z_code_frozen_bits);
    for (int i = 0; i < N; i++) X_code_frozen_bits[i] = Z_code_frozen_bits[N-1-i];
    for (int i = 0; i < N; i++) 
        if (Z_code_frozen_bits[i] == 0 && Z_code_frozen_bits[N-1-i] == 1) 
            X_stab_frozen_bits[i] = 0;

    for (int i = 0; i < N; i++)
        if (X_stab_frozen_bits[i] && !Z_code_frozen_bits[i])
            X_stab_info_indices.push_back(i);

    cerr << "info bits at ";
    for (int i : X_stab_info_indices) cerr << i << " ";
    cerr << endl << "logical operators" << endl;

    Encoder_polar* encoder_Z = new Encoder_polar(K, N, Z_code_frozen_bits);
    vector<int> output(N, 0);
    for (int i : X_stab_info_indices) {
        fill(output.begin(), output.end(), 0);
        output[i] = 1;
        encoder_Z->light_encode(output.data());
        cerr << "weight=" << count_weight(output) << endl;
    }
    cerr << "Z code frozen bits mask" << endl;
    for (auto i : Z_code_frozen_bits) cerr << i;
    cerr << endl;
}

void test_polar_stabilizer(int N, int K, double p) 
{
    // int K = 1536, N = 2048;
    int n = log2(N); 
    vector<bool> frozen_bits(N, 0);
    // frozen_bits_generator_PW(N, K, frozen_bits);
    frozen_bits_generator_BEC(N, K, p, frozen_bits);
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
        assert (verify(N, input, outoutput));
    }
}


void test_syndrome_extraction(int N, int K, bool print)
{
    // use the checks (the N-K rows of best channels) to create N-K syndromes
    // the i^th row is the codeword obtained by setting only the i^th position
    // as info and others as frozen.
    vector<bool> frozen_bits(N, 0);
    vector<bool> reversed_frozen_bits(N, 0);
    vector<bool> stab_frozen_bits(N, 1);
    vector<int>  syndromes(N-K, 0);
    vector<int>  syndromes_at_frozen(N-K, 0);
    vector<int>  input(N, 0);
    vector<int>  output(N, 0);
    vector<int>  noise(N, 0);
    vector<vector<int>> parity_checks(N-K, vector<int>(N,0));
    frozen_bits_generator_PW(N, K, frozen_bits);
    for (int i = 0; i < N; i++) reversed_frozen_bits[i] = frozen_bits[N-1-i];
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
        }

        assert (verify(N-K, syndromes, syndromes_at_frozen));

        // from syndrome back to a noisy codeword
        bit_reversal(syndromes_at_frozen);
        j = 0;
        for (int i = 0; i < N; i++) {
            if (frozen_bits[i]) {
                output[i] = syndromes_at_frozen[j++];
            } else output[i] = 0;
        }
        encoder->light_encode(output.data());
        bit_reversal(output);
        xor_vec(N, output, noise, output);
        assert (encoder->is_codeword(output.data()));

        // use the transpose encoding circuit
        output = noise;
        encoder->transpose_encode(output.data());
        j = 0;
        for (int i = 0; i < N; i++) {
            if (reversed_frozen_bits[i]) {
                syndromes_at_frozen[j++] = output[i];
            }
        }
        if (print) {
            cerr << endl << "transposed circuit       ";
            for (int k : syndromes_at_frozen) cerr << k;
            cerr << endl;
        }
        assert (verify(N-K, syndromes, syndromes_at_frozen));

        // from syndrome back to a noisy codeword
        j = 0;
        for (int i = 0; i < N; i++) {
            if (reversed_frozen_bits[i]) {
                output[i] = syndromes_at_frozen[j++];
            } else output[i] = 0;
        }
        encoder->transpose_encode(output.data());
        xor_vec(N, output, noise, output);
        assert (encoder->is_codeword(output.data()));
    }


}

void test_noise_dist(int N, int num_total, int seed)
{
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    int num_flips;
    vector<int> flips(N+1, 0);
    vector<double> pzs = {0.08};
    // vector<double> pzs = {0.1, 0.12, 0.14, 0.16, 0.18, 0.2};
    for (auto pz : pzs) {
        std::fill(flips.begin(), flips.end(), 0);
        chn_bsc_q->set_prob(0, pz);
        for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
            chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
            num_flips = count_weight(noise_Z);
            flips[num_flips]++;
        }
        for (int k : flips) cerr << k << " ";
        cerr << endl;
    }

}

void test_generate_exact_flip(int N, double p, int num_total, int seed)
{
    int cnt = 0;
    int floor = N * p, ceil = floor + 1;
    cerr << "floor " << floor << ", ceiling " << ceil << endl; 
    int num_flips;
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(0, p);
    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        do {
            chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
            num_flips = count_weight(noise_Z);
            cnt++;
        } while (num_flips != floor && num_flips != ceil);
    }
    cerr << "generate " << num_total << " noise that is with exact flip using " << cnt << " samples" << endl; 

}

void test_direct_syndrome_decoding(int N, int K, int list_size, double pz)
{
    vector<bool> frozen_bits(N, 0);
    vector<int>  frozen_values(N, 0);
    vector<bool> stab_frozen_bits(N, 1);
    vector<int>  syndromes(N-K, 0);
    vector<int>  syndromes_at_frozen(N-K, 0);
    vector<int>  input(N, 0);
    vector<int>  output(N, 0);
    vector<int>  noise_Z(N, 0);
    vector<int>  noise_X(N, 0);
    vector<int>  noisy_codeword_Z(N, 0);
    vector<double> llr_noisy_codeword_Z(N);
    vector<int>  SCL_denoised_codeword_Z(N);
    double pm_best;
    vector<vector<int>> parity_checks(N-K, vector<int>(N,0));
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, 42);
    chn_bsc_q->set_prob(0, pz);
    frozen_bits_generator_PW(N, K, frozen_bits);
    for (int i = 0; i < N; i++)
        if (frozen_bits[i] == 0 && frozen_bits[N-1-i] == 1)
            stab_frozen_bits[i] = 0;
    Encoder_polar* encoder = new Encoder_polar(N, N, frozen_bits);
    Decoder_polar_SCL* decoder = new Decoder_polar_SCL(K, N, list_size, frozen_bits);
    int j = 0;
    for (int i = 0; i < N; i++) {
        if (!stab_frozen_bits[i]) {
            input[i] = 1;
            output = input;
            encoder->light_encode(output.data());
            parity_checks[j++] = output;
            input[i] = 0;
        }
    }
    for (int idx = 0; idx < 100; idx++) {
        chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
        cerr << "noise number of flips: " << count_weight(noise_Z) << endl;
        cerr << "noise: ";
        for (auto i : noise_Z) cerr << i;
        cerr << endl;
        for (int i = 0; i < N-K; i++) 
            syndromes[i] = dot_product(N, parity_checks[i], noise_Z);
        // put reversed syndrome into frozen values
        j = N-K-1;
        for (int i = 0; i < N; i++) 
            if (frozen_bits[i]) frozen_values[i] = syndromes[j--];
        decoder->set_frozen_values(frozen_values);
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-pz)/pz) : log((1-pz)/pz); // 0 -> 1.0; 1 -> -1.0
        pm_best = decoder->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        cerr << "reSCL: ";
        bit_reversal(SCL_denoised_codeword_Z);
        for (auto i : SCL_denoised_codeword_Z) cerr << i;
        cerr << endl;
        cerr << "reSCL number of flips: " << count_weight(SCL_denoised_codeword_Z) << endl;

    }
}