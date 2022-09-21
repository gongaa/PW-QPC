#ifndef DECODER_RM_SC_HPP_
#define DECODER_RM_SC_HPP_
#include "Decoder.hpp"
using namespace std;

class Decoder_RM_SC : public Decoder
{
protected:
    const int m, r;
    const int L;

    
public: 
    Decoder_RM_SC(const int& m, const int& r, const int& L);
    virtual ~Decoder_RM_SC() = default;
    virtual int decode(const double *Y_N, int *V_K, const size_t frame_id);

private:
    // int FHT_decode(); // r = 1, ML decoder based on Fast Hadarmard Transform.
    // int ML_decode();  // r = m, ML decoder.
    int _decode_dumer(const double *Y_N, double *Y_dec_N, int *V_K, const int& m, const int& r, const int& N, int& V_K_len, double *buf);
    int _decode_m1(const double *Y_N, int *V_K, const int& m, const int& N);
    int _decode_llr(const double *Y_N, int *Y_dec_N, int *V_K, const int& m, const int& r, const int& N, int& V_K_len);
    int _decode_m1_llr(const double *Y_N, int *V_K, const int& m, const int& N);
};


#endif /* DECOODER_RM_SC_HPP_ */