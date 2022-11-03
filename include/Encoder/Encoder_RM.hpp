#ifndef ENCODER_RM_HPP_
#define ENCODER_RM_HPP_

#include "Encoder.hpp"
#include "Test/Test_RM.hpp"

using namespace std;
// #define HAMMING
class Encoder_RM : public Encoder
{
protected:
    const int m, r;  // N = 2^m
    int K;           // K = m\choose \leq r
    vector<int> X_dec_N;
  
public:
    Encoder_RM(const int& m, const int& r);
    virtual ~Encoder_RM() = default;
    static int calculate_K(const int& m, const int& r);
    void encode(const int *U_K, int *X_N, const size_t frame_id);
    void decode(const int* X_N, int *V_K);
    void parity_check(const int *X_N, int *S_K);
    static int encode_mm_code(const int* U_K, int *X_N, int N);
    static void recursive_encode_mm_code(const int* U_K, int *X_N, int N);
    static int _encode(const int *U_K, int *X_N, int m, int r); // internal function called by encode()
    static int _decode(int *X_N, int *V_K, const int& m, const int& r);
    static int _parity_check(const int *X_N, int m, int r, int *S_K);

    static bool is_codeword(const int *X_N, int m, int r);
    static bool is_logical(const int *X_N, int m, int r1, int r2);
};
      
#endif // ENCODER_RM_HPP_

