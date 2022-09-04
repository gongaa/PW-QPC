#include <iostream>
#include <random>
#include "Test/Test_RM.hpp"

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
        if (U_K_1[i] != U_K_2[i])
            is_same = false;
    }
    return is_same;
}

// check with no added noise, decoder can correctly decode


// check RM(m, 1) is working, both dumer and the LLR version


// check triangle decode RM(m, 1) perform not worse than to RM(m, 0) 


// check both dumer SC and LLR SC give the same result



