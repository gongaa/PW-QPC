#include <iostream>
#include <random>
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Channel/Channel.hpp"

using namespace std;

// generate codeword of information bits
// or int K, int *U_K
void generate_random(int N, int *Y_N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution d(0.5);
 
    for(int i = 0; i < N; i++) {
        Y_N[i] = d(gen);
    }
}

// verify the encoded and decoded info bits are the same
bool verify(int K, int *U_K_1, int *U_K_2) {
    bool is_same = true;
    for (int i = 0; i < K; i++) {
        if (U_K_1[i] != U_K_2[i])
            is_same = false;
    }
    return is_same;
}

// check with no added noise, decoder can correctly decode


// check RM(m, 1) is working, both dumer and the LLR version


// check triangle decode RM(m, 1) perform not worse than to RM(m, 0) 


// check both dumer SC and LLR SC give the same result under no noise
void test_dumer_llr(bool no_noise = true) {
    int m = 20, r = 5;
    Encoder* encoder = new Encoder_RM(m, r);
    Decoder* decoder = new Decoder_RM_SC(m ,r, 1);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    generate_random(K, info_bits.data());
    vector<double> noisy_codeword(N, 0);
    vector<int> decoded(K, 0);
    encoder->encode(info_bits.data(), codeword.data(), 1);
    // cerr << "codeword is ";
    // for (int i : codeword)
    //   cerr << i << " ";
    // cerr << endl;
    if (no_noise)
        for (int i = 0; i < N; i++)
            noisy_codeword[i] = codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
    decoder->decode(noisy_codeword.data(), decoded.data(), 1);
    // cerr << "decoded result: ";
    // for (int i : decoded)
    //   cerr << i << " ";
    // cerr << endl;
    if (verify(K, info_bits.data(), decoded.data()))
        cerr << "decode successfully" << endl;
    else 
        cerr << "decoding failed" << endl;
}

// test copy_until of Decoder_RM_SCL
void test_copy_until() {
    int m = 20, r = 5;
    Decoder_RM_SCL* decoder = new Decoder_RM_SCL(m ,r, 2);
    decoder->test_copy_until();
}



int simulation() {
    // each time a decoding failure occurred, 
    // we checked whether the decoded codeword was more likely than the transmitted codeword.
    // If so, then the optimal ML decoder would surely misdecode y as well.
    int m = 10, r = 3;
    Decoder_RM_SCL* decoder = new Decoder_RM_SCL(m ,r, 10);
    Encoder* encoder = new Encoder_RM(m, r);
    // Decoder* decoder = new Decoder_RM_SC(m ,r, 1);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
    Channel_BSC* chn_bsc = new Channel_BSC(N, 1e-1, 42);
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    vector<int> noisy_codeword(N, 0);
    vector<double> llr_noisy_codeword(N, 0);
    vector<int> denoised_codeword(N, 0);
    vector<int> decoded(K, 0);
    int num_total = 1000, num_err = 0, num_ml_failed = 0, num_flips = 0, ml_flips=0;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
        // cerr << "codeword is " << endl;
        // for (int i : codeword)
        //     cerr << i << " ";
        // cerr << endl;
        num_flips = chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
        // cerr << "noisy codeword is" << endl;
        // for (int i : noisy_codeword)
        // cerr << i << " ";
        // cerr << endl;
        for (int i = 0; i < N; i++)
            llr_noisy_codeword[i] = noisy_codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
        decoder->decode(llr_noisy_codeword.data(), denoised_codeword.data(), 1);
        // cerr << "denoised codeword result: " << endl;
        // for (int i : denoised_codeword)
        //   cerr << i << " ";
        // cerr << endl;
        if (!verify(N, codeword.data(), denoised_codeword.data())) {
            num_err++;
            // test whether ml decoding will fail
            ml_flips = 0;
            for (int i = 0; i < N; i++)
                if (codeword[i] != denoised_codeword[i]) 
                    ml_flips++;
            if (num_flips < ml_flips) {
                // cerr << "ml failed" << endl;
                num_ml_failed++;
            }
        }
    }
    cerr << "num_err: " << num_err << endl;
    cerr << "Frame Error Rate: " << (double)num_err / num_total << endl;
    cerr << "ML decoding failed rate: " << (double)num_ml_failed / num_total << endl;
    return 0;
}