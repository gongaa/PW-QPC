#include "Decoder/Decoder_RM_syndrome_SC.hpp"
#include "Encoder/Encoder_RM.hpp"
#include <vector>
#include <iostream>
using namespace std;

Decoder_RM_syndrome_SC::Decoder_RM_syndrome_SC(const int& m, const int& r, const int& L)
: Decoder(Encoder_RM::calculate_K(m, r), 1 << m), m(m), r(r), L(L)
{
}

double Decoder_RM_syndrome_SC::decode(const double *Y_N, const int *S_K, int *E_N, const size_t frame_id)
{  // recursive decoding down to (1, m) and (m, m)
    int S_K_len = K;
    this->_decode_llr(Y_N, S_K, E_N, m, r, N, S_K_len);
    return 0.0;
}

int Decoder_RM_syndrome_SC::_decode_llr(const double *Y_N, const int *S_K, int *X_dec_N, const int& m, const int& r, const int& N, int& S_K_len)
{   // return 0 if succeed, return 1 if fail (e.g. in FHT decoder)
    // | G_{m-1, r} | G_{m-1, r}   | | E_L | = | G_{m-1, r} (E_L + E_R) |   | S_L | 
    // | 0          | G_{m-1, r-1} | | E_R | = | G_{m-1, r-1} E_R       | = | S_R |
    int i = 0;
    vector<int> tmp(N); 
    if (r == m) {
        S_K_len = N;
        Encoder_RM::encode_mm_code(S_K, X_dec_N, N);
        return 0; 
    }

    if (r == 0) {
        S_K_len = 1; // one information bit.
        for (i = 0; i < N; i++)
            X_dec_N[i] = 0;
        if (S_K[0] != 0) {
            // double min_llr = std::abs(Y_N[0]);
            double min_llr = Y_N[0];
            int min_idx = 0;
            double temp;
            for (i = 1; i < N; i++) {
                // temp = abs(Y_N[i]);
                temp = Y_N[i];
                if (temp < min_llr) {
                    min_idx = i;
                    min_llr = temp;
                }
            }
            X_dec_N[min_idx] = 1;
        } 
        return 0;
    }

    int N_half = N / 2;
    // first solve G_{m-1, r} (E_L + E_R) = S_L
    const double *Y_fst = Y_N, *Y_snd = Y_N + N_half;
    int *X_dec_fst = X_dec_N, *X_dec_snd = X_dec_N + N_half;
    vector<double> LLRs_buffer(N_half); // store LLR of E_L + E_R
    f_plus(Y_fst, Y_snd, N_half, LLRs_buffer.data());
    int curr_S_K_len = 0;
    _decode_llr(LLRs_buffer.data(), S_K, X_dec_fst, m-1, r, N_half, curr_S_K_len); // decode E_L + E_R 
    S_K_len = curr_S_K_len;

    // then solve  G_{m-1, r-1} E_R       = S_R
    f_minus(Y_fst, Y_snd, X_dec_fst, N_half, LLRs_buffer.data());
    _decode_llr(LLRs_buffer.data(), S_K + curr_S_K_len, X_dec_snd, m-1, r-1, N_half, curr_S_K_len);
    S_K_len += curr_S_K_len;

    // E_L = (E_L + E_R) + E_R
    for (i = 0; i < N_half; i++)
        X_dec_fst[i] ^= X_dec_snd[i]; 

    return 0;
}