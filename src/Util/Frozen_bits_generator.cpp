#include "Util/Util.hpp"
#include "Decoder/Decoder_polar.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cassert>
using namespace std;

void frozen_bits_generator_BEC(const int& N, const int& K, const double& p, vector<bool>& frozen_bits) 
{    // for BEC and BSC only
    cerr << "BEC construction" << endl;
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
	// for (int i = 0; i < K; i++) assert (best_channels[i] == (N-1-best_channels[N-1-i]));
}

double my_alpha = -0.4527;
double my_beta  =  0.0218;
double my_gamma =  0.8600;

double a =  1.0  / my_alpha;
double b = -my_beta / my_alpha;
double c =  1.0  / my_gamma;
double phi_pivot     = 0.867861;
double phi_inv_pivot = 0.6845772418;

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
		return std::exp(my_alpha * std::pow(t, my_gamma) + my_beta);
}

double square_plus_DE(const double& zl, const double& zr)
{
	auto z = phi_inv(1.0 - ((1.0 - phi(zl)) * (1.0 - phi(zr))));
	return (z == HUGE_VAL) ? zl + M_LN2 / (my_alpha * my_gamma) : z;
	// return z;
}

void frozen_bits_generator_AWGN(const int& N, const int& K, const double& db, vector<bool>& frozen_bits) 
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
				z[t * o1] = T + M_LN2 / (my_alpha * my_gamma);

			z[t * o1 + o2] = 2.0 * T;
		}
	}

	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < K; i++) frozen_bits[best_channels[i]] = 0;
}

void frozen_bits_generator_AWGN_SC(const int& N, const int& K, const double& db, vector<bool>& frozen_bits) 
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
    std::iota(best_channels.begin(), best_channels.end(), 0);
	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < K; i++) frozen_bits[best_channels[i]] = 0;
}

void frozen_bits_generator_BSC_SC(const int& N, const int& K, const double& p, vector<bool>& frozen_bits) 
{
    vector<double> z(N, 0);
	std::fill(z.begin(), z.end(), log((1-p)/p));
    vector<bool> fake_frozen_bits(N);
	std::fill(fake_frozen_bits.begin(),     fake_frozen_bits.begin() + K, 0);
	std::fill(fake_frozen_bits.begin() + K, fake_frozen_bits.end(),       1);
    Decoder_polar_SCL* decoder = new Decoder_polar_SCL(K, N, 1, fake_frozen_bits);
    vector<int> denoised_codeword(N);
    decoder->decode_SC(z.data(), denoised_codeword.data(), 0);
    decoder->get_llr_for_frozen_bits(z.data());
    for (int i = 0; i < N; i++) frozen_bits[i] = 1;
    vector<uint32_t> best_channels(N);
    std::iota(best_channels.begin(), best_channels.end(), 0);
	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < K; i++) frozen_bits[best_channels[i]] = 0;
}

void frozen_bits_generator_PW(const int& N, const int& Kz, const int& Kx, vector<bool>& Z_code_frozen_bits, vector<bool>& X_stab_info_bits, double beta)
{
    cerr << "PW construction, beta = " << beta << endl;
	vector<double> z(N, 0);
	int n = log2(N);
	vector<int> bin(n, 0);
	double temp = 0.0;
	for (int i = 0; i < N; i++) {
		decimal2binary(i, bin);
		temp = 0.0;
		for (int k = 0; k < n; k++) 
			if (bin[k]) temp += pow(beta, k);
		z[i] = temp;
	}
	std::fill(Z_code_frozen_bits.begin(), Z_code_frozen_bits.end(), 1);
	std::fill(X_stab_info_bits.begin(), X_stab_info_bits.end(), 0);
    vector<uint32_t> best_channels(N);
    std::iota(best_channels.begin(), best_channels.end(), 0);
	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < Kz; i++) Z_code_frozen_bits[best_channels[i]] = 0;
	for (int i = 0; i < Kz; i++) assert (best_channels[i] == (N-1-best_channels[N-1-i])); // CSS check
	for (int i = 0; i < N-Kx; i++) X_stab_info_bits[best_channels[i]] = 1;
}

void frozen_bits_generator_RM(const int& N, const int& Kz, const int& Kx, vector<bool>& Z_code_frozen_bits, vector<bool>& X_stab_info_bits)
{
    cerr << "RM construction" << endl;
	vector<double> z(N, 0);
	int n = log2(N);
	vector<int> bin(n, 0);
	int weight;
	vector<int> cnt(n+1, 0);
	double eps = 1.0 / N;
	for (int i = 0; i < N; i++) {
		decimal2binary(i, bin);
		weight = count_weight(bin);
		z[i] = weight + cnt[weight] * eps;
		cnt[weight]++;
	}
	std::fill(Z_code_frozen_bits.begin(), Z_code_frozen_bits.end(), 1);
	std::fill(X_stab_info_bits.begin(), X_stab_info_bits.end(), 0);
    vector<uint32_t> best_channels(N);
    std::iota(best_channels.begin(), best_channels.end(), 0);
	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < Kz; i++) Z_code_frozen_bits[best_channels[i]] = 0;
	for (int i = 0; i < Kz; i++) assert (best_channels[i] == (N-1-best_channels[N-1-i]));
	for (int i = 0; i < N-Kx; i++) X_stab_info_bits[best_channels[i]] = 1;
}

void frozen_bits_generator_HPW(const int& N, const int& Kz, const int& Kx, vector<bool>& Z_code_frozen_bits, vector<bool>& X_stab_info_bits)
{
    cerr << "HPW construction" << endl;
	vector<double> z(N, 0);
	int n = log2(N);
	vector<int> bin(n, 0);
	double temp = 0.0;
	double beta1 = pow(2, 0.25);
	double beta2 = pow(2, 0.25*0.25);
	for (int i = 0; i < N; i++) {
		decimal2binary(i, bin);
		temp = 0.0;
		for (int k = 0; k < n; k++) 
			if (bin[k]) temp += pow(beta1, k) + 0.25 * pow(beta2, k);
		z[i] = temp;
	}
	std::fill(Z_code_frozen_bits.begin(), Z_code_frozen_bits.end(), 1);
	std::fill(X_stab_info_bits.begin(), X_stab_info_bits.end(), 0);
    vector<uint32_t> best_channels(N);
    std::iota(best_channels.begin(), best_channels.end(), 0);
	std::sort(best_channels.begin(), best_channels.end(), [z](int i1, int i2) { return z[i1] > z[i2]; });
    for (int i = 0; i < Kz; i++) Z_code_frozen_bits[best_channels[i]] = 0;
	for (int i = 0; i < Kz; i++) assert (best_channels[i] == (N-1-best_channels[N-1-i]));
	for (int i = 0; i < N-Kx; i++) X_stab_info_bits[best_channels[i]] = 1;
}

bool construct_frozen_bits(CONSTRUCTION con, const int& N, const int& Kz, const int& Kx, vector<bool>& Z_code_frozen_bits, vector<bool>& X_stab_info_bits, double beta) {
	// the best Kz positions will be chosen as the info bits for Z code (C_Z, X-type codewords)
	// the best N-Kx positions will be chosen as the info bits for X-type stabilizers (C_X^\perp)
    switch (con) {
        case BEC:
            cerr << "Do not use BEC for CSS code" << endl;
            return false;
        case RM:
            frozen_bits_generator_RM(N, Kz, Kx, Z_code_frozen_bits, X_stab_info_bits);
            break;
        case PW:
            frozen_bits_generator_PW(N, Kz, Kx, Z_code_frozen_bits, X_stab_info_bits, beta);
            break;
        case HPW:
            frozen_bits_generator_HPW(N, Kz, Kx, Z_code_frozen_bits, X_stab_info_bits);
            break;
		case Q1: // Kz+Kx=N+1, so N-Kx=Kz+1
			std::fill(Z_code_frozen_bits.begin(), Z_code_frozen_bits.end(), 1);
			std::fill(X_stab_info_bits.begin(), X_stab_info_bits.end(), 0);
			for (int i = 0; i < Kz; i++) Z_code_frozen_bits[N-i-1] = 0;
			for (int i = 0; i < N-Kx; i++) X_stab_info_bits[N-i-1] = 1;
			break;
        default:
            cerr << "Currently not supported" << endl;
            return false;
    }
    return true;
}

void print_mixing_factor(vector<bool>& frozen_bits) {
	// the mixing factor is the number of information bits that appear before the last frozen bit
	int last_frozen_bit = -1;
	int mixing_factor = 0;
	for (int i = 0; i < frozen_bits.size(); i++)
		if (frozen_bits[i]) last_frozen_bit = i;
	for (int i = 0; i < last_frozen_bit; i++)
		if (!frozen_bits[i]) mixing_factor++;
	cerr << "mixing factor is " << mixing_factor << endl;
}