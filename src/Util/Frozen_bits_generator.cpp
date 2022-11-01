#include "Util/Util.hpp"
#include "Decoder/Decoder_polar.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>
using namespace std;

void frozen_bits_generator_BEC(int N, int K, double p, vector<bool>& frozen_bits) 
{    // for BEC and BSC only
    for (int i = 0; i < N; i++) frozen_bits[i] = 1;

    vector<uint32_t> best_channels(N);
    std::iota(best_channels.begin(), best_channels.end(), 0);

    vector<double> z(N, 0);
	z[0] = std::log(static_cast<double>(p));
    int m = std::log2(N);
	for (auto l = 1; l <= m; l++) {
		auto b = exp2(l);
		auto stride = N / b;
		for (auto j = 0; j < b / 2; j++)
		{
			auto t = z[j * 2 * stride ];
			z[j * 2 * stride         ] = std::log(2) + t + std::log1p(- std::exp(2 * t - (std::log(2) + t)));
			z[j * 2 * stride + stride] = 2 * t;
		}
	}

	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] < z[i2]; });
    for (int i = 0; i < K; i++) frozen_bits[best_channels[i]] = 0;
}

static constexpr double alpha = -0.4527;
static constexpr double beta  =  0.0218;
static constexpr double gamma =  0.8600;

static constexpr double a =  1.0  / alpha;
static constexpr double b = -beta / alpha;
static constexpr double c =  1.0  / gamma;
static constexpr double phi_pivot     = 0.867861;
static constexpr double phi_inv_pivot = 0.6845772418;


double phi_inv(double t)
{
	if (t > phi_inv_pivot)
		return 4.304964539 * (1 - sqrt(1 + 0.9567131408 * std::log(t)));
	else
		return std::pow(a * std::log(t) + b, c);
}

double phi(double t)
{
	if (t < phi_pivot)
		return std::exp(0.0564 * t * t - 0.48560 * t);
	else // if(t >= phi_pivot)
		return std::exp(alpha * std::pow(t, gamma) + beta);
}

double square_plus_DE(const double& zl, const double& zr)
{
	auto z = phi_inv(1.0 - ((1.0 - phi(zl)) * (1.0 - phi(zr))));
	return (z == HUGE_VAL) ? zl + M_LN2 / (alpha * gamma) : z;
	// return z;
}

void frozen_bits_generator_AWGN(int N, int K, double db, vector<bool>& frozen_bits) 
{
    for (int i = 0; i < N; i++) frozen_bits[i] = 1;
    vector<uint32_t> best_channels(N);
    std::iota(best_channels.begin(), best_channels.end(), 0);
    vector<double> z(N, 0);
    int m = std::log2(N);

    for (auto i = 0; i < N; i++)
		z[i] = 2.0 / std::pow(db, 2.0);

	for (auto l = 1; l <= m; l++)
	{
		auto o1 = (int)std::exp2(m - l + 1);
		auto o2 = (int)std::exp2(m - l   );

		for (auto t = 0; t < (int)std::exp2(l -1); t++)
		{
			double T = z[t * o1];

			z[t * o1] = phi_inv(1.0 - std::pow(1.0 - phi(T), 2.0));
			if (z[t * o1] == HUGE_VAL)
				z[t * o1] = T + M_LN2 / (alpha * gamma);

			z[t * o1 + o2] = 2.0 * T;
		}
	}

	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < K; i++) frozen_bits[best_channels[i]] = 0;
}

void frozen_bits_generator_AWGN_SC(int N, int K, double db, vector<bool>& frozen_bits) 
{
    vector<double> z(N, 0);
    std::fill(z.begin(), z.end(), 2.0 / std::pow(db, 2.0));
    vector<bool> fake_frozen_bits(N);
	std::fill(fake_frozen_bits.begin(),     fake_frozen_bits.begin() + K, 0);
	std::fill(fake_frozen_bits.begin() + K, fake_frozen_bits.end(),       1);
    Decoder_polar_SCL* decoder = new Decoder_polar_SCL(K, N, 1, fake_frozen_bits);
    vector<int> denoised_codeword(N);
    decoder->decode_SC(z.data(), denoised_codeword.data(), 0);
    decoder->get_llr_for_frozen_bits(z.data());
    for (int i = 0; i < N; i++) frozen_bits[i] = 1;
    vector<uint32_t> best_channels(N);
	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < K; i++) frozen_bits[best_channels[i]] = 0;
}