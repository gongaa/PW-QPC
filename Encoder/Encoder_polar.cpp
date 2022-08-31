#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <vector>

#include "Encoder_polar.hpp"
using namespace std;

Encoder_polar::Encoder_polar(const int& K, const int& N, const vector<bool>& frozen_bits)
: Encoder(K, N), m((int)log2(N)), frozen_bits(frozen_bits), X_N_tmp(this->N), info_bits_pos(K) 
{
}

void Encoder_polar::encode(const int *U_K, int *X_N, const size_t frame_id)
{
    this->convert(U_K, X_N);
    this->light_encode(X_N);
}

void Encoder_polar::convert(const int *U_K, int *U_N)
{
    if (U_K == U_N) {
        vector<int> U_K_tmp(this->K);
        copy(U_K, U_K + this->K, U_K_tmp.begin());

        auto j = 0;
        for (unsigned i = 0; i < frozen_bits.size(); i++) 
            U_N[i] = (frozen_bits[i]) ? 0 : U_K_tmp[j++];
    } else {
        auto j = 0;
        for (unsigned i = 0; i < frozen_bits.size(); i++)
            U_N[i] = (frozen_bits[i]) ? 0 : U_K[j++];
    }

}

void Encoder_polar::light_encode(int *bits)
{
    for (auto k = (this->N >> 1); k > 0; k >>= 1)
        for (auto j = 0; j < this->N; j += 2 * k)
            for (auto i = 0; i < k; i++)
                bits[j + i] = bits [j + i] ^ bits[k + j + i];
}

bool Encoder_polar::is_codeword(const int *X_N)
{
    copy(X_N, X_N + this->N, this->X_N_tmp.data());

    for (auto k = (this->N >> 1); k > 0; k >>= 1) {
        for (auto j = 0; j < this->N; j += 2 * k) {
            for (auto i = 0; i < k; i++)
                this->X_N_tmp[j+i] = this->X_N_tmp[j+i] ^ this->X_N_tmp[k+j+i];

            if (this->frozen_bits[j+k-1] && this->X_N_tmp[j+k-1])
                return false;
        }
    }
    return true;
}


void Encoder_polar::set_frozen_bits(const vector<bool>& frozen_bits)
{
    copy(frozen_bits.begin(), frozen_bits.end(), this->frozen_bits.begin());
    auto k = 0;
    for (auto n = 0; n < this->N; n++)
        if (!this->frozen_bits[n])
            this->info_bits_pos[k++] = n;
}