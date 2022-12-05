#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
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
       << "Options:\n"
       << "\t-h,--help\t\tShow this help message\n"
       << "\t-d,--destination DESTINATION\tSpecify the destination path"
       << endl;
}

int main(int argc, char** argv)
{
  // write a simple parser for selecting number of threads
  // which code to use
  // which channel to use and what's their parameters.
  int m, rx, rz, list_size, n = 100;
  int N = 1024, K = 513;
  double px, pz; 
  int p_min, p_max; // in percentage
  bool use_crc = false, use_exact = false;
  int seed;
  string con_str;
  CONSTRUCTION con;
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
    } else if (arg == "-K")  {
      iss >> K;
    } else if (arg == "-px") {
      iss >> px;  
    } else if (arg == "-pz") {
      iss >> pz;  
    } else if (arg == "-pmin") {
      iss >> p_min;
    } else if (arg == "-pmax") {
      iss >> p_max;
    } else if ((arg == "-n") || (arg == "--num_samples")) {
      iss >> n;
    } else if ((arg == "-l") || (arg == "--list_size")) {
      iss >> list_size;
    } else if (arg == "-seed") {
      iss >> seed;
    } else if (arg == "-con") {
      iss >> con_str;
      con = construction_map[con_str];
    } else if (arg == "-crc") {
      iss >> use_crc;
    } else if (arg == "-exact") {
      iss >> use_exact;
    } else {
      cerr << "argument not supported" << endl;
      return 1;
    }
  }
  // simulation_RM_CSS(m, rx, rz, list_size);
  // cerr << "Use CRC: " << (use_crc ? "True" : "False") << endl;
  // simulation_RM_CSS_weighted_degeneracy(m ,rx, rz, list_size, px, n, use_crc);
  // int N = 256, K = 128, L = 1024;
  // int N = 32, K = 24, L = 1024;
  // m = 11; rx = 5;
  // double db = 3, design_snr = 3;
  // simulation_polar_CSS(N, K, list_size, pz, n);
  // simulation_polar_syndrome(N, K, list_size, pz, n, con, use_exact, seed);

  /*
  vector<bool> frozen_bits_HPW(N, 0);
  vector<bool> frozen_bits_PW(N, 0);
  frozen_bits_generator_HPW(N, K, frozen_bits_HPW);
  frozen_bits_generator_PW(N, K, frozen_bits_PW);
  for (int i = 0; i < N; i++)
    if (frozen_bits_HPW[i] != frozen_bits_PW[i])
      cerr << "differ at i=" << i << endl;
  */
  int db = 0;
  int design_snr = 1.0; // 1.0dB is the best for Gaussian Approximation Construction
  auto start = std::chrono::high_resolution_clock::now();
  simulation_polar_SCL(N, K, list_size, pz, db, design_snr);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  cerr << "polar takes " << duration.count() << " s" << endl;
  start = std::chrono::high_resolution_clock::now();
  simulation_RM_SCL(m, rx, list_size, pz, db, design_snr);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  cerr << "RM takes " << duration.count() << " s" << endl;
  // Now RM takes twice the time of Polar to decode
  // TODO: optimize the copy in RM SCL decoder
  return 0;
}