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
#ifdef HAMMING
    if (r == 1) {
        S_K_len = m + 1;
        for (i = 0; i < N; i++)
            X_dec_N[i] = 0;
        bool all_zero = std::all_of(S_K, S_K + m, [](int i) { return i == 0; });
        if (S_K[m] == 0 && all_zero) {}
        else if (S_K[m] == 1 && all_zero) {
            X_dec_N[0] = 1; 
        } else if (S_K[m] == 1) {
            X_dec_N[binary2decimal(S_K, m)] = 1;
        } else { // two errors
            X_dec_N[0] = 1;
            X_dec_N[binary2decimal(S_K, m)] = 1;
        }
        return 0;
    }
    if (r == m) {
        S_K_len = N;
        copy(S_K, S_K + N, X_dec_N);
        return 0;
    }
#else
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
            double min_llr = Y_N[0];
            int min_idx = 0;
            double temp;
            for (i = 1; i < N; i++) {
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
#endif // HAMMING
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
    // _decode_llr(Y_snd, S_K + curr_S_K_len, X_dec_snd, m-1, r-1, N_half, curr_S_K_len);
    S_K_len += curr_S_K_len;

    // E_L = (E_L + E_R) + E_R
    for (i = 0; i < N_half; i++)
        X_dec_fst[i] ^= X_dec_snd[i]; 

    return 0;
}