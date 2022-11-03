#include <iostream>
#include <string>
#include <sstream>
#include <vector>
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
  double px, pz; 
  int p_min, p_max; // in percentage
  bool use_crc = false;
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
    } else if (arg == "--crc") {
      iss >> use_crc;
    } else {
      cerr << "argument not supported" << endl;
      return 1;
    }
  }
  // verify_RM_is_codeword();
  // simulation_RM_CSS(m, rx, rz, list_size);
  // test_RM_SCL_symmetry();
  // verify_parity_check();
  // test_RM_syndrome_SC(m, rx);
  // simulation_RM_degeneracy(m, rx, rz, px, pz);
  // simulation_RM_closest_codewords(m, rx, rz);
  // test_linearity_xor(6, 3);
  // simulation_RM_d_star(m, rx);
  // test_RM_d_star();
  // compare_equiv_classes();
  // simulation_symmetric_noise(m, rx);
  cerr << "Use CRC: " << (use_crc ? "True" : "False") << endl;
  simulation_RM_CSS_weighted_degeneracy(m ,rx, rz, list_size, p_min, p_max, n, use_crc);
  // test_polar();
  // test_crc();
  // test_encoder_decode();
  // simulation_RM_SCL();
  return 0;
}