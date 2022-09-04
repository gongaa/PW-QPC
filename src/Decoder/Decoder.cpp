#include "Decoder/Decoder.hpp"
using namespace std;

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

auto Decoder::lambdas = vector<function<double(const vector<double> &LLRs, const vector<int> &bits)>>(my_lambdas);

Decoder::Decoder(const int K, const int N) : K(K), N(N)
{
}

double Decoder::phi(const double& mu, const double& lambda, const int& u)
{   // path metric update functino
    double new_mu;
    if (u == 0 && lambda < 0)
        new_mu = mu - lambda;
    else if (u != 0 && lambda > 0)
        new_mu = mu + lambda;
    else // if u = [1-sign(lambda)]/2 correct prediction
        new_mu = mu;

    return new_mu;
}