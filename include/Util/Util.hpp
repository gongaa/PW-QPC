#ifndef UTIL_HPP_
#define UTIL_HPP_
#include <random>
#include <cmath>
#include <vector>
using namespace std;
inline bool verify(int K, int *U_K_1, int *U_K_2) {
    bool is_same = true;
    for (int i = 0; i < K; i++) {
        if (U_K_1[i] != U_K_2[i]) {
            is_same = false;
            break;
        }
    }
    return is_same;
}

inline void xor_vec(int N, int* a, const int* b, int* c) {
    for (int i = 0; i < N; i++)
        c[i] = a[i] ^ b[i];
}

inline int count_flip(int N, const int* a, const int* b) {
    int cnt = 0;
    for (int i = 0; i < N; i++) 
        if (a[i] != b[i]) cnt++;
    return cnt;
}

inline int count_weight(const vector<int>& a) {
    int cnt = 0;
    for (int i : a)
        if (i) cnt++;
    return cnt;
}

inline double db2val(double x) {
  return exp(log(10.0) * x / 10.0);
}

inline void decimal2bianry(const int& n, vector<int>& b)
{
    int x = n;
    for (int i = 0; x > 0; i++) {
        b[i] = x % 2;
        x >>= 1;
    }
}

inline int binary2decimal(const int *b, int size)
{
    int n = 0;
    for (int i = 0; i < size; i++) {
        if (b[i] == 1)
            n += (1 << i);
    }
    return n;
}

inline void generate_random(int N, int *Y_N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution d(0.5);
 
    for(int i = 0; i < N; i++) {
        Y_N[i] = d(gen);
    }
}

void frozen_bits_generator_BEC(int N, int K, double p, vector<bool>& frozen_bits); 

void frozen_bits_generator_AWGN(int N, int K, double db, vector<bool>& frozen_bits);

void frozen_bits_generator_AWGN_SC(int N, int K, double db, vector<bool>& frozen_bits);

double phi_inv(double t);

double phi(double t);

double square_plus_DE(const double& zl, const double& zr);

#endif // UTIL_HPP_