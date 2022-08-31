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

    void encode(const int *U_K, int *X_N, const size_t frame_id);
    int encode_mm_code(const int* U_K, int *X_N, int N);

    bool is_codeword(const int *X_N);

private:
    int calculate_K(const int& m, const int& r);
    int _encode(const int *U_K, int *X_N, int m, int r); // internal function called by encode()
    // two ways to encode an RM(m, m) code
    void recursive_encode_mm_code(const int* U_K, int *X_N, int N);
    void light_encode_mm(int *bits);



};
      
#endif

