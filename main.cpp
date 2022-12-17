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
  int m, rx, rz, list_size, n = 1000;
  int N = 1024, K = 513;
  double px, pz; 
  int p_min, p_max; // in percentage
  bool use_crc = false, use_fast = true;
  int exact_t = 1;
  int seed = 0;
  int print_interval = 1000;
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
      iss >> exact_t;
    } else if (arg == "-fast") {
      iss >> use_fast;
    } else if (arg == "-interval") {
      iss >> print_interval;
    } else {
      cerr << "argument not supported" << endl;
      return 1;
    }
  }

  // int db = 0;
  // int design_snr = 1.0; // 1.0dB is the best for Gaussian Approximation Construction
  // simulation_polar_SCL(N, K, list_size, pz, db, design_snr, CONSTRUCTION::RM, n);
  // simulation_RM_SCL(m, rx, list_size, pz, db, design_snr, n);
  // print_polar_con(N, K, con);
  // test_polar_stabilizer(N, K, con);
  auto start = std::chrono::high_resolution_clock::now();
  if (!use_fast)
    simulation_polar_syndrome(N, K, list_size, pz, n, con, exact_t, seed, print_interval);
  else
    simulation_polar_syndrome_fast(N, K, list_size, pz, n, con, exact_t, seed, print_interval);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  cerr << "Finish in " << duration.count() << " s" << endl;

  return 0;
}