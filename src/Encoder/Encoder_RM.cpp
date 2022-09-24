#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

#include "Encoder/Encoder_RM.hpp"
using namespace std;

Encoder_RM::Encoder_RM(const int& m, const int& r)
: Encoder(calculate_K(m, r), 1 << m), m(m), r(r)
{
}

int Encoder_RM::calculate_K(const int& m, const int& r) 
{
    vector<int> s(r+1, 1);
    for (int i = 1; i <= m; i++) 
        for (int j = r; j >= 0; j--) 
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
{   // the generator matrix of RM(m, m) code is {{1,1,},{0,1}}^{\otimes m}
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
    // WARGNING: the recursion formula for the generator matrix is subject to bit reversal
    int K_v = _encode(U_K, X_N, m-1, r-1);
    int K_u = _encode(U_K + K_v, X_N + N_half, m-1, r);
    for (int i = 0; i < N_half; i++)
        X_N[i] ^= X_N[i + N_half]; // after bit reversal, get Plotkin (u, u+v)
    return K_v + K_u;
}

bool Encoder_RM::is_codeword(const int *X_N, int m, int r)
{   // (u+v, u) where u\in RM(m-1, r), v\in RM(m-1, r-1)
    int N = 1 << m;
    vector<int> X_N_tmp(X_N, X_N + N);
    return _is_codeword(X_N_tmp.data(), m, r);
}

bool Encoder_RM::_is_codeword(int *X_N, int m, int r)
{
    cerr << "m=" << m << ", r=" << r << endl;
    if (m == r) return true;
    int N = 1 << m, N_half = N >> 1;
    if (r == 0) return std::equal(X_N + 1, X_N + N, X_N);
    if (!_is_codeword(X_N + N_half, m-1, r)) return false;
    cerr << "before folding" << endl;
    for (int i = 0; i < N; i++) cerr << X_N[i] << " ";
    std::transform(X_N, X_N + N_half, X_N + N_half,
               X_N, std::bit_xor<bool>());
    cerr << endl << "after folding" << endl;
    for (int i = 0; i < N; i++) cerr << X_N[i] << " ";
    return _is_codeword(X_N, m-1, r-1); 
}

bool Encoder_RM::is_logical(const int *X_N, int m, int r1, int r2)
{
    assert (r1 > r2);
    int N = 1 << m;
    vector<int> X_N_tmp(X_N, X_N + N);
    if (!_is_codeword(X_N_tmp.data(), m, r1)) return false;
    copy(X_N, X_N + N, X_N_tmp.data());
    return !_is_codeword(X_N_tmp.data(), m, r2);
}