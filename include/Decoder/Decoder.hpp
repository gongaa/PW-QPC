#ifndef DECODER_HPP_
#define DECODER_HPP_

#include <vector>
#include <cmath>
using namespace std;
 
class Decoder
{
public:
    const int K, N;
    // vector<double> Y_N;
    // static: let different decoders share the lambdas
    // lambdas are LLR update rules
    static vector<function<double(const vector<double> &LLRs, const vector<int> &bits)>> lambdas;

public:
    explicit Decoder(const int K, const int N);
    virtual int decode(const double *Y_N, int *V_K, const size_t frame_id) = 0;
    static double phi(const double& mu, const double& lambda, const int& u); // path metric update function
};

#endif /* DECODER_HPP_ */