#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>

#include "Encoder_RM.hpp"
using namespace std;

Encoder_RM::Encoder_RM(const int& m, const int& r)
: m(m), r(r), Encoder(calculate_K(m, r), 1 << m)
{
}

int Encoder_RM::calculate_K(const int& m, const int& r) 
{
    // vector<vector<int>> s(m, vector<int>(r, 0));
    vector<int> s(r+1, 1);
    for (int i = 1; i <= m; i++) 
        for (int j = r; j >= 0; j--) 
            s[j] = s[j] + s[j-1];
    // recursion formula: {m\choose \leq r} = {(m-1)\choose \leq r} + {(m-1)\choose \leq (r-1)} 
    return s[r];

}

int Encoder_RM::encode_mm_code(const int* U_K, int *X_N, int N) // U_K: information vector, X_N: output codeword
{   // Encode (m, m) code for BSC (X_N \in {0,1}).
    // Or decode, since the inverse matrix is the same.
    recursive_encode_mm_code(U_K, X_N, N);
    // or light encode
    // for (int i = 0; i < N; i++)
    //     X_N[i] = U_K[i];
    // for (auto k = (N >> 1); k > 0; k >>= 1)
    //     for (auto j = 0; j < N; j += 2 * k)
    //         for (auto i = 0; i < k; i++)
    //             X_N[j + i] = X_N[j + i] ^ X_N[k + j + i];
    // TODO: test that they are equal 
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