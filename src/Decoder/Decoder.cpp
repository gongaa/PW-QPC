#include "Decoder/Decoder.hpp"
#include <cassert>
using namespace std;
// #define USE_APPROXIMATION

vector<function<double(const vector<double> &LLRs, const vector<int> &bits)>> my_lambdas =  {
    [](const vector<double> &LLRs, const vector<int> &bits) -> double
    {   // the hardware-efficient f- function
        auto sign = std::signbit(LLRs[0]) ^ std::signbit(LLRs[1]);
        auto abs0 = std::abs(LLRs[0]);
        auto abs1 = std::abs(LLRs[1]);
        auto min = std::min(abs0, abs1);
        return sign ? -min : min;
    },
    [](const vector<double> &LLRs, const vector<int> &bits) -> double
    {   // the f+ function
        return ((bits[0] == 0) ? LLRs[0] : -LLRs[0]) + LLRs[1];
    }
};

// auto Decoder::lambdas = vector<function<double(const vector<double> &LLRs, const vector<int> &bits)>>(my_lambdas);

void Decoder::f_plus(const double* LLR_fst, const double* LLR_snd, const int size, double* LLR_new)
{
    for (int i = 0; i < size; i++) {
#ifdef USE_APPROXIMATION
        auto sign = signbit(LLR_fst[i]) ^ signbit(LLR_snd[i]);
        auto abs0 = abs(LLR_fst[i]);
        auto abs1 = abs(LLR_snd[i]);
        auto min = std::min(abs0, abs1);
        LLR_new[i] = sign ? -min : min;
#else
        LLR_new[i] = log(exp(LLR_fst[i] + LLR_snd[i]) + 1) - log(exp(LLR_fst[i]) + exp(LLR_snd[i]));
#endif // USE_APPROXIMATION
    }
}

void Decoder::f_minus(const double* LLR_fst, const double* LLR_snd, const int* bits, const int size, double* LLR_new)
{
    for (int i = 0; i < size; i++) {
        LLR_new[i] = ((bits[i] == 0) ? LLR_fst[i] : -LLR_fst[i]) + LLR_snd[i];
    }
}

Decoder::Decoder(const int K, const int N) : K(K), N(N)
{
}

double Decoder::phi(const double& mu, const double& lambda, const int& u)
{   // path metric update function
    assert(!isnan(mu));
    assert(!isnan(lambda));
    double new_mu;
#ifdef USE_APPROXIMATION
    if (u == 0 && lambda < 0)
        new_mu = mu - lambda;
    else if (u != 0 && lambda > 0)
        new_mu = mu + lambda;
    else // if u = [1-sign(lambda)]/2 correct prediction
        new_mu = mu;
#else
    new_mu = mu + ((u == 0) ? log(1 + exp(-lambda)) : log(1 + exp(lambda)));
    assert(!isnan(new_mu));
#endif // USE_APPROXIMATION
    return new_mu;
}