#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>
#include <limits>
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Decoder/Decoder_RM_syndrome_SC.hpp"
#include "Channel/Channel.hpp"
#include "Encoder/Encoder_polar.hpp"
#include "Decoder/Decoder_polar.hpp"
#include "Util/CRC_polynomial.hpp"
#include "Util/Util.hpp"

// #define CHN_AWGN

// Reed-Muller
int simulation_RM_SCL(int m, int r, int list_size, double p, double db, double design_snr, int num_total=1000); 

int simulation_RM_CSS(int m, int rx, int rz, int list_size);

int simulation_RM_CSS_weighted_degeneracy(int m, int rx, int rz, int list_size, double px, int num_total, bool use_crc);

int simulation_RM_degeneracy(int m, int rx, int rz, double px, double pz);

int simulation_RM_closest_codewords(int m, int rx, int rz); 

int simulation_RM_d_star(int m, int r);

int simulation_symmetric_noise(int m, int r);

int test_RM_d_star();

int compare_equiv_classes();

// Polar
void simulation_polar_SCL(int N, int K, int L, double p, double db, double design_snr, CONSTRUCTION con=CONSTRUCTION::PW, int num_total=1000); 

void simulation_polar_CSS(int N, int K, int L, double p, int num_total, int seed);

void simulation_polar_syndrome(int N, int K, int L, double p, int num_total, CONSTRUCTION con, int exact_t=0, int seed=42);

void simulation_polar_syndrome_fast(int N, int K, int list_size, double pz, int num_total, CONSTRUCTION con, int exact_t=0, int seed=42);

void simulation_stab_MW_codewords(int N, int K);

#endif