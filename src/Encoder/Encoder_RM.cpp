#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>

#include "Encoder/Encoder_RM.hpp"
using namespace std;

Encoder_RM::Encoder_RM(const int& m, const int& r)
: Encoder(calculate_K(m, r), 1 << m), m(m), r(r), X_dec_N(1 << m, 0)
{
}

int Encoder_RM::calculate_K(const int& m, const int& r) 
{
    vector<int> s(r+1, 1);
    for (int i = 1; i <= m; i++) 
        for (int j = r; j >= 1; j--) 
            s[j] = s[j] + s[j-1];
    // recursion formula: {m\choose \leq r} = {(m-1)\choose \leq r} + {(m-1)\choose \leq (r-1)} 
    return s[r];

}

int Encoder_RM::encode_mm_code(const int *U_K, int *X_N, int N) // U_K: information vector, X_N: output codeword
{   // Encode (m, m) code for BSC (X_N \in {0,1}).
    // Or decode, since the inverse matrix is the same.
    // recursive_encode_mm_code(U_K, X_N, N);
    // or light encode
    for (int i = 0; i < N; i++)
        X_N[i] = U_K[i];
    for (auto k = (N >> 1); k > 0; k >>= 1)
        for (auto j = 0; j < N; j += 2 * k)
            for (auto i = 0; i < k; i++)
                X_N[j + i] = X_N[j + i] ^ X_N[k + j + i];
    // the above code has the same effect as:
    // recursive_encode_mm_code(U_K, X_N, N);
    return N;

}

void Encoder_RM::recursive_encode_mm_code(const int* U_K, int *X_N, int N)
{   // the generator matrix of RM(m, m) code is {{1,1},{0,1}}^{\otimes m}
    if (N == 1) { // m = 0
        X_N[0] = U_K[0];
        return;
    }
    int N_half = N >> 1;
    recursive_encode_mm_code(U_K, X_N, N_half);
    recursive_encode_mm_code(U_K + N_half, X_N + N_half, N_half);
    for (int i = 0; i < N_half; i++)
        X_N[i] ^= X_N[i + N_half];
}

void Encoder_RM::encode(const int *U_K, int *X_N, const size_t frame_id)
{
    _encode(U_K, X_N, this->m, this->r);
}

int Encoder_RM::_encode(const int *U_K, int *X_N, int m, int r) 
{   // perhaps terrible performance
    int N = 1 << m;
    if (r == 0) { // K = 1
        for (int i = 0; i < N; i++)
            X_N[i] = U_K[0];
        return 1;
    } 
    if (r == m) {
        return encode_mm_code(U_K, X_N, N);
    } 
    int N_half = N >> 1;
    int K_v = _encode(U_K, X_N, m-1, r-1);
    int K_u = _encode(U_K + K_v, X_N + N_half, m-1, r);
    for (int i = 0; i < N_half; i++)
        X_N[i] ^= X_N[i + N_half]; // (u+v, u)
    return K_v + K_u;
}

void Encoder_RM::decode(const int *X_N, int *V_K)
{
   copy(X_N, X_N + N, X_dec_N.data());
   _decode(X_dec_N.data(), V_K, this->m, this->r);
}

int Encoder_RM::_decode(int *X_N, int *V_K, const int& m, const int& r)
{
   if (r == 0) {
      V_K[0] = X_N[0];
      return 1;
   }
   if (r == m) {
      int N = 1 << m;
      encode_mm_code(X_N, V_K, N);
      return N;
   }
   // (u+v, u)
   int N_half = 1 << (m-1);
   for (int i = 0; i < N_half; i++) X_N[i] ^= X_N[N_half + i];
   int K_v = _decode(X_N, V_K, m-1, r-1);
   int K_u = _decode(X_N + N_half, V_K + K_v, m-1, r);
   return K_v + K_u;
}

bool Encoder_RM::is_codeword(const int *X_N, int m, int r)
{   // (u+v, u) where u\in RM(m-1, r), v\in RM(m-1, r-1)
    if (m == r) return true;
    int N = 1 << m;
    if (r == 0) {// return std::equal(X_N + 1, X_N + N, X_N);
        for (int i = 1; i < N; i++) {
            if (X_N[i] != X_N[0]) return false;
        }
        return true;
    }
    int N_half = N >> 1;
    if (!is_codeword(X_N + N_half, m-1, r)) return false;
    vector<int> X_N_tmp(X_N, X_N + N_half);
    for (int i = 0; i < N_half; i++)
        X_N_tmp[i] = X_N[i] ^ X_N[i + N_half];
    // std::transform(X_N, X_N + N_half, X_N + N_half,
    //            X_N_tmp.begin(), std::bit_xor<bool>());
    return is_codeword(X_N_tmp.data(), m-1, r-1);
}

bool Encoder_RM::is_logical(const int *X_N, int m, int r1, int r2)
{
    assert (r1 > r2);
    int N = 1 << m;
    vector<int> X_N_tmp(X_N, X_N + N);
    if (!is_codeword(X_N_tmp.data(), m, r1)) return false;
    copy(X_N, X_N + N, X_N_tmp.data());
    return !is_codeword(X_N_tmp.data(), m, r2);
}

void Encoder_RM::parity_check(const int *X_N, int *S_K)
{
    this->_parity_check(X_N, m, m-r-1, S_K);
}

int Encoder_RM::_parity_check(const int *X_N, int m, int r, int *S_K)
{   // X_N: (u+v, u) a codeword of RM(m, m-r-1), u\in RM(m-1, m-r-1), v\in RM(m-1, m-r-2)
    //                                           u\perp RM(m-1, r-1), v\perp RM(m-1, r)
    // parity check matrix size: {m\choose <= r} * N 
    // | G_{m-1, r} | G_{m-1, r}   | | X_L | = | G_{m-1, r} (X_L + X_R) |   | G_{m-1, r}   v |
    // | 0          | G_{m-1, r-1} | | X_R | = | G_{m-1, r-1} X_R       | = | G_{m-1, r-1} u |
    // cerr << "in parity check m=" << m << ", r=" << r << endl;
    int N = 1 << m;
#ifdef HAMMING
    if (r == 1) { // Hamming parity check
        for (int i = 0; i < m+1; i++) S_K[i] = 0;
        vector<int> temp(m, 0);
        for (int i = 0; i < N; i++) {
            if (X_N[i] != 0) {
                S_K[m] ^= 1;
                decimal2bianry(i, temp);
                for (int l = 0; l < m; l++)
                    S_K[l] ^= temp[l];
            }
        }
        return (m+1);
    }
    if (r == m) {
        copy(X_N, X_N + N, S_K);
        return N;
    }
#else
    if (r == m) {
        return encode_mm_code(X_N, S_K, N);
    }
    if (r == 0) {
        S_K[0] = X_N[0];
        for (int i = 1; i < N; i++)
            S_K[0] ^= X_N[i];
        return 1;
    }
#endif // HAMMING
    int N_half = N >> 1; 
    vector<int> X_fst(N_half, 0);
    for (int i = 0; i < N_half; i++)
        X_fst[i] = X_N[i] ^ X_N[i + N_half]; // v
    int K_v = _parity_check(X_fst.data(), m-1, r, S_K);
    int K_u = _parity_check(X_N + N_half, m-1, r-1, S_K + K_v);
    return K_v + K_u;
}