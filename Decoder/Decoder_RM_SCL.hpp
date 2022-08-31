#ifndef DECODER_RM_SCL_HPP_
#define DECODER_RM_SCL_HPP_
#include "Decoder.hpp"
#include <vector>
#include "Decoder.hpp"
using namespace std;

class Decoder_RM_SCL : Decoder
{
protected:
    const int m, r;
    const int L;
    
public: 
    Decoder_RM_SCL(const int& m, const int& r, const int& L);
    virtual ~Decoder_RM_SCL();
    virtual int decode(const float *Y_N, int *V_K, const size_t frame_id);

private:
    void FHT_decode(); // r = 1, ML decoder based on Fast Hadarmard Transform.
    void ML_decode();  // r = m, ML decoder.
};


#endif /* DECOODER_RM_SCL_HPP_ */