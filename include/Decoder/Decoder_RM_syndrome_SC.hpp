#ifndef DECODER_RM_SYNDROME_SC_HPP_
#define DECODER_RM_SYNDROME_SC_HPP_
#include "Decoder.hpp"
using namespace std;

class Decoder_RM_syndrome_SC : public Decoder
{
protected:
    const int m, r;
    const int L;

    
public: 
    Decoder_RM_syndrome_SC(const int& m, const int& r, const int& L);
    virtual ~Decoder_RM_syndrome_SC() = default;
    double decode(const double *Y_N, const int *S_K, int *E_N, const size_t frame_id);

private:
    int _decode_llr(const double *Y_N, const int *S_K, int *X_dec_N, const int& m, const int& r, const int& N, int& S_K_len);
};


#endif /* DECOODER_RM_SYNDROME_SC_HPP_ */