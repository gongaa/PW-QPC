#ifndef ENCODER_RM_HPP_
#define ENCODER_RM_HPP_

#include "Encoder.hpp"
using namespace std;

class Encoder_RM : public Encoder
{
protected:
    const int m, r;  // N = 2^m
    int K;           // K = m\choose \leq r
  
public:
    Encoder_RM(const int& m, const int& r);
    virtual ~Encoder_RM() = default;
    static int calculate_K(const int& m, const int& r);
    void encode(const int *U_K, int *X_N, const size_t frame_id);
    static int encode_mm_code(const int* U_K, int *X_N, int N);
    static void recursive_encode_mm_code(const int* U_K, int *X_N, int N);
    static int _encode(const int *U_K, int *X_N, int m, int r); // internal function called by encode()

    static bool is_codeword(const int *X_N, int m, int r);
    static bool _is_codeword(int *X_N, int m, int r);
    static bool is_logical(const int *X_N, int m, int r1, int r2);
};
      
#endif

