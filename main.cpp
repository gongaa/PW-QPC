#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Test/Test_RM.hpp"
#include "Test/Test_polar.hpp"
#include "Simulation/Simulation.hpp"
#include "Channel/Channel.hpp"

// #define _OPENMP
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num () { return 0; }
inline int omp_get_num_threads() { return 1; }
#endif

using namespace std;

static void show_usage(string name)
{
  cerr << "Usage: " << name << " <options> SOURCES"
       << "Options:\n\t"
       << left << setw(20) << "-h,--help" << "Show this help message\n\t"
       << left << setw(20) << "-N" << "Specify the blocklength (a power of two)\n\t"
       << left << setw(20) << "-Kz" << "Specify the constituent Z code rate\n\t"
       << left << setw(20) << "-Kx" << "Specify the constituent X code rate\n\t"
       << left << setw(20) << "-px" << "Specify the bit-flip error rate\n\t"
       << left << setw(20) << "-n,--num_samples" << "Specify the number of samples to perform simulation on\n\t"
       << left << setw(20) << "-l,--list_size" << "Specify the list size used in the SCL decoder\n\t"
       << left << setw(20) << "-seed" << "Specify the random seed\n\t"
       << left << setw(20) << "-con" << "Specify the construction, choose from PW, HPW, RM, Q1, BEC\n\t"
       << left << setw(20) << "-version" << "Specify the decoding method, choose from 0(codeword decoder, default), 1(syndrome decoder)\n\t"
       << left << setw(20) << "-interval" << "Specify the printing interval, default is 1000\n\t"
       << left << setw(20) << "-beta" << "Specify the beta used in the PW construction, default is 2^{1/4}"
       << endl;
}

int main(int argc, char** argv)
{
  int m, rx, rz; // for Reed-Muller together with the Dumer's list decoder 
  int list_size = 32, n = 1000;
  int N = 1024, Kz = 513, Kx = 513;
  double px; 
  int version = 0;
  int seed = 42;
  int print_interval = 1000;
	double beta = pow(2, 0.25);
  bool is_depolarize = false;
  string con_str;
  CONSTRUCTION con = CONSTRUCTION::PW;
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      return 0;
    } 
    std::istringstream iss( argv[++i] );
    if (arg == "-m") {
      iss >> m;
    } else if (arg == "-rx") {
      iss >> rx;
    } else if (arg == "-rz") {
      iss >> rz;
    } else if (arg == "-N")  {
      iss >> N;
    } else if (arg == "-Kz") {
      iss >> Kz;
    } else if (arg == "-Kx") {
      iss >> Kx;
    } else if (arg == "-px") {
      iss >> px;  
    } else if ((arg == "-n") || (arg == "--num_samples")) {
      iss >> n;
    } else if ((arg == "-l") || (arg == "--list_size")) {
      iss >> list_size;
    } else if (arg == "-seed") {
      iss >> seed;
    } else if (arg == "-con") {
      iss >> con_str;
      con = construction_map[con_str];
    } else if (arg == "-version") {
      iss >> version;
    } else if (arg == "-interval") {
      iss >> print_interval;
    } else if (arg == "-beta") {
      iss >> beta;
    } else if ((arg == "-dep") || (arg == "--depolarize")) {
      iss >> is_depolarize;
    } else {
      cerr << "Argument not supported." << endl;
      return 1;
    }
  }

  if (N & (N-1)) {
    cerr << "N is not a power of two, abort." << endl;
    return 1;
  } else if (Kx+Kz <= N) {
    cerr << "Encoding zero qubit, abort." << endl;
    return 1;
  }
  // int db = 0;
  // int design_snr = 1.0; // 1.0dB is the best for Gaussian Approximation Construction
  // simulation_polar_SCL(N, K, list_size, pz, db, design_snr, CONSTRUCTION::RM, n);
  // simulation_RM_SCL(m, rx, list_size, pz, db, design_snr, n);
  auto start = std::chrono::high_resolution_clock::now();
  
  if (is_depolarize)
    simulation_polar_depolarize(N, Kz, list_size, px, n, con, beta, seed, print_interval); // only codeword decoding is done for depolarizing channel
  else if (version == 0)
    simulation_polar_codeword(N, Kz, Kx, list_size, px, n, con, beta, seed, print_interval);
  else if (version == 1)
    simulation_polar_syndrome_direct(N, Kz, Kx, list_size, px, n, con, beta, seed, print_interval);
  else
    simulation_polar_syndrome(N, Kz, Kx, list_size, px, n, con, beta, seed, print_interval);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  cerr << "Finish in " << duration.count() << " s" << endl;

  return 0;
}