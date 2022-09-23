#include <iostream>
#include <random>
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Channel/Channel.hpp"
// #define CHN_AWGN
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
        if (U_K_1[i] != U_K_2[i]) {
            is_same = false;
            break;
        }
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

static double db2val(double x) {
  return exp(log(10.0) * x / 10.0);
}

int simulation() {
    // each time a decoding failure occurred, 
    // we checked whether the decoded codeword was more likely than the transmitted codeword.
    // If so, then the optimal ML decoder would surely misdecode y as well.
    int m = 7, r = 3, list_size = 16;
    Encoder* encoder = new Encoder_RM(m, r);
    // Decoder_RM_SC* SC_decoder = new Decoder_RM_SC(m ,r, 1);
    Decoder_RM_SCL* SCL_decoder = new Decoder_RM_SCL(m ,r, list_size);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
    cerr << "List size=" << list_size << endl;
#ifdef CHN_AWGN
    double db = 4;
    double code_rate = (double)K / N;
    double sigma = 1 / sqrt(2 * code_rate * db2val(db));
    cerr << "sigma=" << sigma << ", db=" << db << endl;
    Channel_AWGN* chn_awgn = new Channel_AWGN(N, sigma, 42);
#else
    double p = 0.07;
    cerr << "p=" << p << endl;
    Channel_BSC* chn_bsc = new Channel_BSC(N, p, 42);
#endif // USE_AWGN
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    vector<int> noisy_codeword(N, 0);
    vector<double> llr_noisy_codeword(N, 0);
    vector<double> dumer_noisy_codeword(N, 0);
    vector<int> SC_denoised_codeword(N, 0);
    vector<int> SCL_denoised_codeword(N, 0);
    vector<int> decoded(K, 0);
    int num_total = 10000, SC_num_err = 0, SCL_num_err = 0, num_ml_failed = 0;
    int SC_num_flips = 0, SCL_num_flips = 0, ml_flips=0;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
#ifdef CHN_AWGN
        chn_awgn->add_noise(codeword.data(), llr_noisy_codeword.data(), 0);
#else
        ml_flips = chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
        for (int i = 0; i < N; i++) {
            llr_noisy_codeword[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            // llr_noisy_codeword[i] = noisy_codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
            // dumer_noisy_codeword[i] = noisy_codeword[i] ? (2*p-1) : (1-2*p);
        }
#endif // USE_AWGN
        // SC_decoder->decode(llr_noisy_codeword.data(), SC_denoised_codeword.data(), 1);
        // SC_decoder->decode(dumer_noisy_codeword.data(), SC_denoised_codeword.data(), 1);
        SCL_decoder->decode(llr_noisy_codeword.data(), SCL_denoised_codeword.data(), 1);
        /*
        if (!verify(N, codeword.data(), SC_denoised_codeword.data())) {
            // cerr << "codeword is " << endl;
            // for (int i : codeword)
            //     cerr << i << " ";
            // cerr << endl;
            // cerr << "denoised codeword result: " << endl;
            // for (int i : SC_denoised_codeword)
            //     cerr << i << " ";
            // cerr << endl;
            SC_num_err++;
            // test whether ml decoding will fail
            SC_num_flips = 0;
            for (int i = 0; i < N; i++)
                if (codeword[i] != SC_denoised_codeword[i]) 
                    SC_num_flips++;
            if (SC_num_flips < ml_flips) {
                // cerr << "ML decoding failed for SC, SC_num_flips=" << SC_num_flips << ", ml_flips=" << ml_flips << endl;
                num_ml_failed++;
            }
        }
        */
        if (!verify(N, codeword.data(), SCL_denoised_codeword.data())) {
            SCL_num_err++;
            // test whether ml decoding will fail
            /*
            SCL_num_flips = 0;
            for (int i = 0; i < N; i++)
                if (codeword[i] != SCL_denoised_codeword[i]) 
                    SCL_num_flips++;
            if (SCL_num_flips > SC_num_flips) {
                cerr << "SCL is worse than SC, SC_num_flips=" << SC_num_flips << ", SCL_num_flips=" << SCL_num_flips << endl;
                int path_idx = SCL_decoder->is_codeword_in_list(SC_denoised_codeword.data());
                cerr << "path_idx=" << path_idx << endl;
            }
            */
        }
    }
    cerr << "SC_num_err: " << SC_num_err << ". SCL_num_err: " << SCL_num_err << endl;
    cerr << "SC Frame Error Rate: " << (double)SC_num_err / num_total << endl;
    cerr << "SCL Frame Error Rate: " << (double)SCL_num_err / num_total << endl;
    cerr << "ML decoding failed rate: " << (double)num_ml_failed / num_total << endl;
    return 0;
}