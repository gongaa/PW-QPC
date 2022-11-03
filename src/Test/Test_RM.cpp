#include <iostream>
#include <random>
#include <cassert>
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Decoder/Decoder_RM_syndrome_SC.hpp"
#include "Channel/Channel.hpp"
#include "Util/CRC_polynomial.hpp"

using namespace std;

void test_linearity_xor(int m, int r) {
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N();
    vector<int> info_bits1(K, 0);
    vector<int> info_bits2(K, 0);
    vector<int> info_bits_dest(K, 0);
    vector<int> codeword1(N, 0);
    vector<int> codeword2(N, 0);
    vector<int> codeword_dest(N, 0);
    vector<int> codeword_dest_xor(N, 0);
    int x1, x2, x_dest;
    int max = 1 << K;
    for (int i = 0; i < 100; i++) {
        x1 = rand() % max;
        x2 = rand() % max;
        decimal2bianry(x1, info_bits1);
        decimal2bianry(x2, info_bits2);
        for (int l = 0; l < K; l++) info_bits_dest[l] = info_bits1[l] ^ info_bits2[l];
        encoder->encode(info_bits1.data(), codeword1.data(), 1);
        encoder->encode(info_bits2.data(), codeword2.data(), 1);
        encoder->encode(info_bits_dest.data(), codeword_dest.data(), 1);
        for (int l = 0; l < K; l++) codeword_dest_xor[l] = codeword1[l] ^ codeword2[l];
        assert (verify(N, codeword_dest.data(), codeword_dest_xor.data()));
    }
    cerr << "all tests completed, RM is linear under xor of info bits" << endl;
}

void generate_symmetric_noise(int m, vector<int>& noise, int level) {
    vector<int> binary_repr(m, 0);
    for (int i = 0; i < (1 << m); i++) {
        decimal2bianry(i, binary_repr);
        if (count_weight(binary_repr) <= level) noise[i] = 1;
    }
}

void generate_bent(int m, vector<int>& noise) {
    vector<int> binary_repr(m, 0);
    int temp;
    for (int i = 0; i < (1 << m); i++) {
        decimal2bianry(i, binary_repr);
        temp = 0;
        for (int k = 0; k < m; k += 2) {
            temp ^= binary_repr[k] * binary_repr[k+1];
        }
        if (temp) noise[i] = 1;
    }
}

void generate_bent_guess(int m, vector<int>& noise) {
    vector<int> binary_repr(m, 0);
    int temp;
    for (int i = 0; i < (1 << m); i++) {
        decimal2bianry(i, binary_repr);
        temp = 0;
        // for (int k = 0; k < m; k += 1) {
        //     temp ^= binary_repr[k] * binary_repr[(k+1) % m];
        // }
        for (int k = 0; k < (m-2); k += 2) {
            temp ^= binary_repr[k] * binary_repr[k+1];
        }
        if (m % 2 == 1) temp ^= binary_repr[m-1];
        else temp ^= binary_repr[m-2] * binary_repr[m-1];
        if (temp) noise[i] = 1;
    }
}

void generate_bent_third_order(int m, vector<int>& noise) {
    vector<int> binary_repr(m, 0);
    int temp;
    for (int i = 0; i < (1 << m); i++) {
        decimal2bianry(i, binary_repr);
        temp = 0;
        for (int k = 0; k < m; k += 3) {
            temp ^= binary_repr[k] * binary_repr[k+1] * binary_repr[k+2];
        }
        if (temp) noise[i] = 1;
    }
}

void generate_bent_m6_r2_order(vector<int>& noise) {
    vector<int> b(6, 0);
    int temp;
    for (int i = 0; i < (1 << 6); i++) {
        decimal2bianry(i, b);
        temp = 0;
        temp = b[0]*b[1]*b[2] ^ b[0]*b[3]*b[4] ^ b[1]*b[3]*b[5] ^ b[2]*b[4]*b[5] ^ b[3]*b[4]*b[5];
        if (temp) noise[i] = 1;
    }
}

void generate_all_codewords(int m, int r, vector<vector<int>>& codewords) {
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N();
    vector<int> info_bits(K, 0);
    vector<int> codeword(N, 0);
    for (int i = 0; i < (1 << K); i++) {
        decimal2bianry(i, info_bits);
        encoder->encode(info_bits.data(), codeword.data(), 1);
        codewords[i] = codeword;
        // if (Encoder_RM::is_codeword(codeword.data(), m, r-1)) {
        //     for (int k : info_bits) cerr << k;
        //     cerr << endl;
        // }
    }
}

void generate_all_equiv_classes(int m, int rx, int rz, vector<vector<int>>& codewords, vector<vector<int>>& equiv_class) {
    assert (m-rx-1 < rz);
    Encoder_RM_CSS* encoder = new Encoder_RM_CSS(m, rx, rz);
    Encoder_RM* encoder_Z = new Encoder_RM(m, rz);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "m=" << m << ", m-rx-1=" << m-rx-1 << ", rz=" << rz
    << ", K=" << K << ", N=" << N << endl;
    int num_equiv_class = 1 << K;
    cerr << "number of equiv classes: " << num_equiv_class << endl;
    generate_all_codewords(m, rz, codewords);
    int num_codewords = codewords.size();
    cerr << "have generated all " << num_codewords << " codewords of RM(m,rz)" << endl;
    vector<int> noise_X_diff(N, 0);
    equiv_class.push_back({0});
    bool is_in_one_class;
    for (int i = 1; i < num_codewords; i++) {
        is_in_one_class = false;
        for (auto& ec : equiv_class) {
            xor_vec(N, codewords[ec[0]].data(), codewords[i].data(), noise_X_diff.data());
            if (encoder->is_X_stabilizer(noise_X_diff.data())) {
                ec.push_back(i);
                is_in_one_class = true;
                break;
            }
        }
        if (!is_in_one_class)
            equiv_class.push_back({i});
    }
    cerr << "have generated all " << equiv_class.size() << " equivalence classes" << endl;
    // for (auto ec : equiv_class) cerr << ec.size() << " ";
}

// check RM is_codeword
void verify_RM_is_codeword() {
  int m = 5, r = 2;
  Encoder* encoder = new Encoder_RM(m, r);
  int K = encoder->get_K(), N = encoder->get_N();
  cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
  vector<int> codeword(N, 0);
  cerr << Encoder_RM::is_codeword(codeword.data(), m, r) << endl;
  for (int i = 1; i < 9; i++) codeword[i] = 1;
  cerr << Encoder_RM::is_codeword(codeword.data(), m, r) << endl;
  /*
  Channel_BSC* chn_bsc = new Channel_BSC(N, 1e-1, 42);
  vector<int> info_bits(K, 1);
  vector<int> codeword(N, 0);
  vector<int> noisy_codeword(N, 0);
  generate_random(K, info_bits.data());
  encoder->encode(info_bits.data(), codeword.data(), 1);
  assert (Encoder_RM::is_codeword(codeword.data(), m, r));
  assert (Encoder_RM::is_codeword(codeword.data(), m, r+1));
  chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
  assert (!Encoder_RM::is_codeword(noisy_codeword.data(), m, r));
  */
}

// check with no added noise, decoder can correctly decode
void test_encoder_decode() {
    int m = 20, r = 10;
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N();
    vector<int> info_bits(K), decoded_info_bits(K);
    vector<int> codeword(N);
    for (int idx = 0; idx < 100; idx++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
        encoder->decode(codeword.data(), decoded_info_bits.data());
        assert (verify(K, info_bits.data(), decoded_info_bits.data()));
    } 
}

// check RM(m, 1) is working, both dumer and the LLR version


// check triangle decode RM(m, 1) perform not worse than to RM(m, 0) 


// check both dumer SC and LLR SC give the same result under no noise
void test_dumer_llr() {
    int m = 20, r = 5;
    Encoder* encoder = new Encoder_RM(m, r);
    Decoder_RM_SC* decoder = new Decoder_RM_SC(m ,r, 1);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
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
    generate_random(K, info_bits.data());
    vector<int> noisy_codeword(N, 0);
    vector<double> llr_noisy_codeword(N, 0);
    vector<double> dumer_noisy_codeword(N, 0);
    vector<int> decoded(K, 0);
    vector<int> denoised_codeword(N, 0);
    int num_total = 10000, SC_num_err = 0, num_ml_failed = 0;
    int SC_num_flips = 0, SCL_num_flips = 0, ml_flips=0;
    encoder->encode(info_bits.data(), codeword.data(), 1);
#ifdef CHN_AWGN
    chn_awgn->add_noise(codeword.data(), llr_noisy_codeword.data(), 0);
#else
    ml_flips = chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
    if (p == 0.0)
        for (int i = 0; i < N; i++)
            noisy_codeword[i] = codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
    else 
        for (int i = 0; i < N; i++) {
            llr_noisy_codeword[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            dumer_noisy_codeword[i] = noisy_codeword[i] ? (2*p-1) : (1-2*p);
        }
#endif // USE_AWGN
    decoder->decode(llr_noisy_codeword.data(), denoised_codeword.data(), 1);
    if (verify(N, codeword.data(), denoised_codeword.data()))
        cerr << "decode successfully" << endl;
    else 
        cerr << "decoding failed" << endl;
}

void test_RM_syndrome_SC(int m, int r) {
    // int m = 7, r = 3;
    Encoder_RM* encoder = new Encoder_RM(m, r);
    Decoder_RM_SC* decoder = new Decoder_RM_SC(m ,r, 1);
    Decoder_RM_syndrome_SC* syndrome_decoder = new Decoder_RM_syndrome_SC(m, m-r-1, 1);
    int K = encoder->get_K(), N = encoder->get_N(), dual_K = syndrome_decoder->K;
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << 
    ", dual_K=" << dual_K << ", N=" << N << endl;
    double p = 0.04;
    cerr << "p=" << p << endl;
    Channel_BSC* chn_bsc = new Channel_BSC(N, p, 42);
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    vector<int> noise(N, 0);
    vector<int> noisy_codeword(N, 0);
    vector<int> syndrome(dual_K, 0);
    vector<double> llr_noisy_codeword(N, 0);
    vector<double> syndrome_llr(N, 0);
    vector<int> codeword_error(N, 0);
    vector<int> syndrome_error(N, 0);
    vector<int> syndrome_of_syndrome_error(dual_K);
    vector<int> denoised_codeword(N, 0);
    int num_total = 100, SC_num_err = 0, syndrome_SC_num_err = 0, num_ml_failed = 0;
    int SC_num_flips = 0, SCL_num_flips = 0, ml_flips=0;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
        ml_flips = chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
        xor_vec(N, codeword.data(), noisy_codeword.data(), noise.data());
        encoder->parity_check(noisy_codeword.data(), syndrome.data());
        if (p == 0.0)
            for (int i = 0; i < N; i++) {
                syndrome_llr[i] = 1.0; 
                llr_noisy_codeword[i] = noisy_codeword[i] ? -1 : 1; // 0 -> 1.0; 1 -> -1.0
            }
        else 
            for (int i = 0; i < N; i++) {
                syndrome_llr[i] = log((1-p)/p);
                llr_noisy_codeword[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            }
            
        decoder->decode(llr_noisy_codeword.data(), denoised_codeword.data(), 1);
        xor_vec(N, denoised_codeword.data(), noisy_codeword.data(), codeword_error.data());
        syndrome_decoder->decode(syndrome_llr.data(), syndrome.data(), syndrome_error.data(), 1);
        encoder->parity_check(syndrome_error.data(), syndrome_of_syndrome_error.data());
        assert (verify(dual_K, syndrome_of_syndrome_error.data(), syndrome.data()));
        cerr << endl << "syndrome error: ";
        for (int k : syndrome_error) cerr << k;
        cerr << endl << "codeword error: ";
        for (int k : codeword_error) cerr << k;
        cerr << endl << "real noise    : ";
        for (int k : noise) cerr << k;
        cerr << endl;
        if (verify(N, syndrome_error.data(), codeword_error.data()))
            cerr << "syndrome decoding gave the same result as codeword decoding" << endl;
        else 
            cerr << "not the same" << endl;
        if (verify(N, syndrome_error.data(), noise.data())) {
            cerr << "syndrome decoding succeeds" << endl;
        } else {
            cerr << "syndrome decoding fails" << endl;
            syndrome_SC_num_err++;
        }
        if (!verify(N, codeword_error.data(), noise.data())) {
            SC_num_err++;
        }
    }
    cerr << "SC_num_err: " << SC_num_err << ". syndrome_SC_num_err: " << syndrome_SC_num_err << endl;
    cerr << "SC Frame Error Rate: " << (double)SC_num_err / num_total << endl;
    cerr << "syndrome SC Frame Error Rate: " << (double)syndrome_SC_num_err / num_total << endl;
}

// test copy_until of Decoder_RM_SCL
void test_copy_until() {
    int m = 20, r = 5;
    Decoder_RM_SCL* decoder = new Decoder_RM_SCL(m ,r, 2);
    decoder->test_copy_until();
}

void verify_parity_check() {
    int m = 8, r = 5;
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N();
    int dual_K = Encoder_RM::calculate_K(m, m-r-1);
    cerr << "m=" << m << ", r=" << r << ", K=" << K << ", N=" << N << ", dual_K=" << dual_K << endl;
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    vector<int> syndrome(dual_K, 0);

    int num_total = 1000;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
        encoder->parity_check(codeword.data(), syndrome.data());
        for (int j = 0; j < dual_K; j++) {
            if (syndrome[j] != 0)
                cerr << "fail at j=" << j << endl;
        }
        //     cerr << syndrome[j] << " ";
        // cerr << endl;
    }
}

void test_is_logical() {
    int m = 8, r = 5;
    Encoder* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N();
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    int num_total = 10000;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
        if (Encoder_RM::is_logical(codeword.data(), m, r, r-1)) {
            for (int i : codeword) cerr << i << " ";
            cerr << endl;
        }
    }
}

// check SCL decoder behaves symmetrically around every codeword under BSC
void test_RM_SCL_symmetry() {
    int m = 7, r = 3, list_size = 128;
    Encoder* encoder = new Encoder_RM(m, r);
    // Decoder_RM_SC* SC_decoder = new Decoder_RM_SC(m ,r, 1);
    Decoder_RM_SCL* SCL_decoder = new Decoder_RM_SCL(m ,r, list_size);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
    cerr << "List size=" << list_size << endl;

    double p = 0.01;
    cerr << "p=" << p << endl;
    Channel_BSC* chn_bsc = new Channel_BSC(N, p, 42);
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    vector<int> all_zero(N, 0);
    vector<int> noise(N, 0); // all-zero + noise
    vector<int> noisy_codeword(N, 0);
    vector<int> diff_to_noise(N, 0);
    vector<double> llr_noisy_codeword_1(N, 0);
    vector<double> llr_noisy_codeword_2(N, 0);
    vector<int> SCL_denoised_codeword_1(N, 0);
    vector<int> SCL_denoised_codeword_2(N, 0);
    vector<vector<int>> X_list_1(list_size, vector<int>(N, 0));
    vector<double> pm_X_list_1(list_size, 0.0);
    vector<vector<int>> X_list_2(list_size, vector<int>(N, 0));
    vector<double> pm_X_list_2(list_size, 0.0);
    int num_total = 3, SCL_num_err = 0;
    int SCL_num_flips_1 = 0, SCL_num_flips_2 = 0, bsc_flips = 0;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
        bsc_flips = chn_bsc->add_noise(all_zero.data(), noise.data(), 0);
        cerr << "BSC #flips=" << bsc_flips << endl;
        xor_vec(N, noise.data(), codeword.data(), noisy_codeword.data());
        for (int i = 0; i < N; i++) {
            llr_noisy_codeword_1[i] = noise[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            llr_noisy_codeword_2[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            // llr_noisy_codeword[i] = noisy_codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
        }
        SCL_decoder->decode(llr_noisy_codeword_1.data(), SCL_denoised_codeword_1.data(), 0);
        SCL_decoder->copy_codeword_list(X_list_1, pm_X_list_1);
        SCL_decoder->decode(llr_noisy_codeword_2.data(), SCL_denoised_codeword_2.data(), 0);
        SCL_decoder->copy_codeword_list(X_list_2, pm_X_list_2);
        for (int i = 0; i < list_size; i++) {
            SCL_num_flips_1 = count_flip(N, noise.data(), X_list_1[i].data());
            cerr << "1: pm=" << pm_X_list_1[i] << "\t";
            cerr << "#flips=" << SCL_num_flips_1 << "\t";
            for (int j : X_list_1[i]) cerr << j;
            cerr << endl;
            SCL_num_flips_2 = count_flip(N, noisy_codeword.data(), X_list_2[i].data());
            cerr << "2: pm=" << pm_X_list_2[i] << "\t";
            cerr << "#flips=" << SCL_num_flips_2 << "\t";
            for (int j = 0; j < N; j++) X_list_2[i][j] ^= codeword[j];
            for (int j : X_list_2[i]) cerr << j;
            cerr << endl;
        }
    }
}