#ifndef UTIL_HPP_
#define UTIL_HPP_
#include <random>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <string>
using namespace std;

enum CONSTRUCTION { BEC = 0, RM = 1, PW = 2, HPW = 3 };
static std::map<std::string, CONSTRUCTION> construction_map = {
    {"BEC", BEC},
    {"RM" , RM },
    {"PW" , PW },
    {"HPW", HPW}
};

inline bool verify(int K, vector<int>& U_K_1, vector<int>& U_K_2) {
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

inline int dot_product(int N, vector<int>& a, vector<int>& b) {
    int sum = 0;
    for(int i = 0; i < N; i++)
        if (a[i] && b[i]) sum++;
    return (sum % 2);
}

inline void bit_reversal(vector<int>& a) {
    int temp;
    int N = a.size();
    for (int i = 0; i < N/2; i++) {
        temp = a[i];
        a[i] = a[N-1-i];
        a[N-1-i] = temp;
    }
}

inline int count_flip(int N, vector<int>& a, vector<int>& b) {
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

inline void decimal2binary(const int& n, vector<int>& b)
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
 
    for (int i = 0; i < N; i++) {
        Y_N[i] = d(gen);
    }
}

inline void print_wt_dist(vector<int>& wt) {
    sort(wt.begin(), wt.end());
    int i = wt[0];
    int cnt = 1;
    for (int k = 1; k < wt.size(); k++) {
        if (wt[k] == i) cnt++;
        else {
            cerr << i << ":" << cnt << "  ";
            i = wt[k]; cnt = 1;
        }
    }
    cerr << i << ":" << cnt << endl;
}

inline void cal_wt_dist_prob(vector<int>& wt, double& p, const int& offset, const double& weight) {
    int i = wt[0];
    int cnt = 1;
    for (int k = 1; k < wt.size(); k++) {
        if (wt[k] == i) cnt++;
        else {
            p += pow(weight, i-offset) * cnt;
            i = wt[k]; cnt = 1;
        }
    }
    p += pow(weight, i-offset) * cnt;
}

void frozen_bits_generator_BEC(const int& N, const int& K, const double& p, vector<bool>& frozen_bits); 

void frozen_bits_generator_AWGN(const int& N, const int& K, const double& db, vector<bool>& frozen_bits);

void frozen_bits_generator_AWGN_SC(const int& N, const int& K, const double& db, vector<bool>& frozen_bits);

void frozen_bits_generator_BSC_SC(const int& N, const int& K, const double& p, vector<bool>& frozen_bits);

void frozen_bits_generator_PW(const int& N, const int& K, vector<bool>& frozen_bits);

void frozen_bits_generator_RM(const int& N, const int& K, vector<bool>& frozen_bits);

void frozen_bits_generator_HPW(const int& N, const int& K, vector<bool>& frozen_bits);

bool construct_frozen_bits(CONSTRUCTION con, int& N, int& K, vector<bool>& frozen_bits);

double phi_inv(double t);

double phi(double t);

double square_plus_DE(const double& zl, const double& zr);

#endif // UTIL_HPP_