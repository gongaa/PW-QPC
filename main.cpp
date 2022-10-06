#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Test/Test_RM.hpp"
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
  int m, rx, rz, list_size;
  double px, pz;
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
    } else if ((arg == "-l") || (arg == "--list_size")) {
      iss >> list_size;
    } else {
      cerr << "argument not supported" << endl;
      return 1;
    }
  }
  // simulation_RM_CSS(m, rx, rz, list_size);
  // test_RM_SCL_symmetry();
  // verify_parity_check();
  // test_RM_syndrome_SC();
  simulation_RM_degeneracy(m, rx, rz, px, pz);
  return 0;
}