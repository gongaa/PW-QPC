#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>

#include "Encoder_polar.hpp"
using namespace std;

Encoder_polar::Encoder_polar(const int& K, const int& N, const vector<bool>& frozen_bits)
: Encoder(K, N), m((int)log2(N)), frozen_bits(frozen_bits), X_N_tmp(this->N) 
{
}

void Encoder_polar::_encode(const int *U_K, int *X_N, const size_t frame_id)
{
    this->convert(U_K, X_N);
    this->light_encode(X_N);
}

void Encoder_polar::light_encode(int *bits)
{
    for (auto k = (this->N >> 1); k > 0; k >>= 1)
        for (auto j = 0; j < this->N; j += 2 * k)
            for (auto i = 0; i < k; i++)
                bits[j + i] = bits [j + i] ^ bits[k + j + i];
}
