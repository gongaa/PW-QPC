#include <iostream>
#include <random>
#include <cassert>
#include <bits/stdc++.h>
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Decoder/Decoder_RM_syndrome_SC.hpp"
#include "Channel/Channel.hpp"
// #define CHN_AWGN
// #define IS_VERBOSE
// #define VANILLA_LIST_DECODING
#define DEGENERACY_LIST_DECODING
// #define CHECK
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

// check RM is_codeword
void verify_RM_is_codeword() {
  int m = 13, r = 7;
  Encoder* encoder = new Encoder_RM(m, r);
  int K = encoder->get_K(), N = encoder->get_N();
  cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
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
}

// check with no added noise, decoder can correctly decode


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

void test_RM_syndrome_SC() {
    int m = 7, r = 3;
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

static double db2val(double x) {
  return exp(log(10.0) * x / 10.0);
}

int simulation_RM_SCL() {
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
        }
#endif // USE_AWGN
        // SC_decoder->decode(llr_noisy_codeword.data(), SC_denoised_codeword.data(), 0);
        SCL_decoder->decode(llr_noisy_codeword.data(), SCL_denoised_codeword.data(), 0);
#ifdef IS_VERBOSE
        if (!verify(N, codeword.data(), SC_denoised_codeword.data())) {
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
#endif // IS_VERBOSE
        if (!verify(N, codeword.data(), SCL_denoised_codeword.data())) {
            SCL_num_err++;
#ifdef IS_VERBOSE
            // test whether ml decoding will fail
            SCL_num_flips = 0;
            for (int i = 0; i < N; i++)
                if (codeword[i] != SCL_denoised_codeword[i]) 
                    SCL_num_flips++;
            if (SCL_num_flips > SC_num_flips) {
                cerr << "SCL is worse than SC, SC_num_flips=" << SC_num_flips << ", SCL_num_flips=" << SCL_num_flips << endl;
                int path_idx = SCL_decoder->is_codeword_in_list(SC_denoised_codeword.data());
                cerr << "path_idx=" << path_idx << endl;
            }
#endif // IS_VERBOSE
        }
    }
    cerr << "SC_num_err: " << SC_num_err << ". SCL_num_err: " << SCL_num_err << endl;
    cerr << "SC Frame Error Rate: " << (double)SC_num_err / num_total << endl;
    cerr << "SCL Frame Error Rate: " << (double)SCL_num_err / num_total << endl;
    cerr << "ML decoding failed rate: " << (double)num_ml_failed / num_total << endl;
    return 0;
}

int simulation_RM_CSS(int m, int rx, int rz, int list_size) {
    // int m = 5, rx = 3, rz = 2, list_size = 16384;
    Encoder_RM_CSS* encoder = new Encoder_RM_CSS(m, rx, rz);
    // if two X-type errors differ by a codeword of RM(m, rz)
    // they will have the same syndrome
    Decoder_RM_SCL* SCL_decoder_X = new Decoder_RM_SCL(m, rz, list_size);
    Decoder_RM_SCL* SCL_decoder_Z = new Decoder_RM_SCL(m, rx, list_size);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", rx="<< rx << ", rz=" << rz
         << ", K=" << K << ", N=" << N << endl;
    cerr << "List size=" << list_size << endl;
    double px = 0.1, pz = 0;
    cerr << "px=" << px << ", pz=" << pz << endl;
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, px, pz, 41);
    // ideally, the noise_X should be decoded to the all-zero codeword
    // but it also can be any codeword of RM(m, rz)
    vector<int> noise_X(N, 0);
    vector<int> noise_X_diff(N, 0);
    vector<int> desired_X(N, 0);
    vector<int> noise_Z(N, 0);
    vector<double> llr_noisy_codeword_X(N, 0);
    vector<double> llr_noisy_codeword_Z(N, 0);
    vector<int> SCL_denoised_codeword_X(N, 0);
    vector<vector<int>> X_list(list_size, vector<int>(N, 0));
    vector<double> pm_X_list(list_size, 0.0);
    // vector<int> SCL_denoised_codeword_Z(N, 0);
    // vector<vector<int>> Z_list(list_size, vector<int>(N, 0));
    // vector<double> pm_Z_list(list_size, 0.0);
    vector<vector<int>> equiv_class;
    bool is_in_one_class; int largest_class_size; int *largest_class;
    int num_total = 100, SCL_num_X_err = 0, SCL_num_Z_err = 0;
    int SCL_num_X_err_deg = 0, SCL_num_X_err_deg_list = 0;
    int num_Z_list_err = 0, SCL_num_Z_err_deg = 0, SCL_num_Z_err_deg_list = 0;
    double pm_best;
    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
        for (int i = 0; i < N; i++) {
            llr_noisy_codeword_X[i] = noise_X[i] ? -log((1-px)/px) : log((1-px)/px); // 0 -> 1.0; 1 -> -1.0
            // llr_noisy_codeword_Z[i] = noise_Z[i] ? -log((1-pz)/px) : log((1-px)/px); // 0 -> 1.0; 1 -> -1.0
        }
        pm_best = SCL_decoder_X->decode(llr_noisy_codeword_X.data(), SCL_denoised_codeword_X.data(), 0);
        SCL_decoder_X->copy_codeword_list(X_list, pm_X_list);
        // SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        // SCL_decoder_Z->copy_codeword_list(Z_list, pm_Z_list);
// #ifdef VANILLA_LIST_DECODING
        cerr << "******* turn " << turn_idx << endl;
        if (!std::all_of(SCL_denoised_codeword_X.begin(), SCL_denoised_codeword_X.end(), [](int i) { return i==0; })) {
            cerr << "SCL was wrong, the codeword it gave has pm=" << pm_best << " and content" << endl;
            for (int k : SCL_denoised_codeword_X) cerr << k;
            cerr << ", pm=" << pm_best << ", #flips=" << 
            count_flip(N, SCL_denoised_codeword_X.data(), noise_X.data()) << endl;
            SCL_num_X_err++;
            for (int i = 0; i < list_size; i++) {
                if (encoder->is_X_stabilizer(X_list[i].data())) {
                    cerr << "wrong but idx=" << i << " differs by only a stabilizer" << endl;
                    for (int k : X_list[i]) cerr << k;
                    cerr << ", pm=" << pm_X_list[i] << ", #flips=" <<
                    count_flip(N, X_list[i].data(), noise_X.data()) << endl;
                }
            }
        }
        if (!encoder->is_X_stabilizer(SCL_denoised_codeword_X.data()))
            SCL_num_X_err_deg++;
// #elif defined DEGENERACY_LIST_DECODING
        equiv_class.clear();
        equiv_class.push_back({0});
        for (int i = 1; i < list_size; i++) {
            is_in_one_class = false;
            for (auto& ec : equiv_class) {
                xor_vec(N, X_list[ec[0]].data(), X_list[i].data(), noise_X_diff.data());
                if (encoder->is_X_stabilizer(noise_X_diff.data())) {
                    ec.push_back(i);
                    is_in_one_class = true;
                    break;
                }
            }
            if (!is_in_one_class)
                equiv_class.push_back({i});
        }
        // cerr << "there are " << equiv_class.size() << " equiv classes" << endl;
#ifdef CHECK
        for (int i = 0; i < list_size; i++)
            assert (Encoder_RM::is_codeword(X_list[i].data(), m, rz));
        int num_equiv_class = equiv_class.size();
        for (int i = 0; i < num_equiv_class; i++) {
            for (int j = 0; j < equiv_class[i].size(); j++) 
                cerr << equiv_class[i][j] << " ";
            cerr << endl;
        }

        for (int i = 0; i < num_equiv_class; i++) {
            for (int j = i + 1; j < num_equiv_class; j++) {
                xor_vec(N, X_list[equiv_class[i][0]].data(), X_list[equiv_class[j][0]].data(), noise_X_diff.data());
                assert (!encoder->is_X_stabilizer(noise_X_diff.data()));
                assert (encoder->is_logical_X(noise_X_diff.data()));
            }
            int ec_size = equiv_class[i].size();
            if (ec_size > 1) {
                for (int j = 0; j < ec_size; j++) {
                    for (int k = j + 1; k < ec_size; k++) {
                        xor_vec(N, X_list[equiv_class[i][j]].data(), X_list[equiv_class[i][k]].data(), noise_X_diff.data());
                        // assert (Encoder_RM::is_codeword(noise_X_diff.data(), m, rz));
                        assert (!encoder->is_logical_X(noise_X_diff.data()));
                    }
                }
            }
        }
#endif // CHECK
        largest_class_size = 0;
        cerr << "classes sizes are: ";
        for (auto& ec : equiv_class) {
            if (encoder->is_X_stabilizer(X_list[ec[0]].data())) cerr << "stab_";
            cerr << ec.size() << " ";
            if (ec.size() > largest_class_size) {
                largest_class_size = ec.size();
                largest_class = X_list[ec[0]].data();
            }
        }
        cerr << endl;
        // cerr << "larget equiv class size " << largest_class_size << endl;
        if (!encoder->is_X_stabilizer(largest_class))
            SCL_num_X_err_deg_list++;
        // else
            // cerr << "turn " << turn_idx << " 00....0 is correctly decoded considering degeneracy, largest equiv class size " << largest_class_size << endl;
// endif
    }
    cerr << "SCL_num_err: " << SCL_num_X_err << endl;
    // cerr << "SCL_num_X_list_0_err: " << SCL_num_X_list_0_err << endl;
    cerr << "SCL_num_err_deg: " << SCL_num_X_err_deg << endl;
    cerr << "SCL_num_err_deg_list: " << SCL_num_X_err_deg_list << endl;
    cerr << "SCL Frame Error Rate: " << (double)SCL_num_X_err / num_total << endl;
    return 0;
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