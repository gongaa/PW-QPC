#include <iostream>
#include <string>
#include <vector>
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Test/Test_RM.hpp"
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
  vector <string> sources;
  string destination;
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      return 0;
    } else if ((arg == "-d") || (arg == "--destination")) {
      if (i + 1 < argc) {
        destination = argv[i++];
      } else {
        cerr << "--destination option requires one argument." << endl;
        return 1;
      }
    } else {
      sources.push_back(argv[i]);
    }
  }
  simulation();
  /*
  int m = 10, r = 3;
  Decoder_RM_SCL* decoder = new Decoder_RM_SCL(m ,r, 10);
  // decoder->test_copy_until();
  // decoder->test_assign_path_idx();
  Encoder* encoder = new Encoder_RM(m, r);
  // Decoder* decoder = new Decoder_RM_SC(m ,r, 1);
  int K = encoder->get_K(), N = encoder->get_N();
  cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
  Channel_BSC* chn_bsc = new Channel_BSC(N, 1e-2, 42);
  vector<int> info_bits(K, 1);
  vector<int> codeword(N, 0);
  generate_random(K, info_bits.data());
  vector<double> noisy_codeword(N, 0);
  vector<int> denoised_codeword(N, 0);
  vector<int> decoded(K, 0);
  encoder->encode(info_bits.data(), codeword.data(), 1);
  // cerr << "codeword is " << endl;
  // for (int i : codeword)
  //   cerr << i << " ";
  // cerr << endl;
  chn_bsc->add_noise(codeword.data(), 0);
  // cerr << "flipped codeword is" << endl;
  // for (int i : flipped_codeword)
  //   cerr << i << " ";
  // cerr << endl;
  for (int i = 0; i < N; i++)
    noisy_codeword[i] = codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
  // decoder->decode(noisy_codeword.data(), decoded.data(), 1);
  decoder->decode(noisy_codeword.data(), denoised_codeword.data(), 1);
  // cerr << "denoised codeword result: " << endl;
  // for (int i : denoised_codeword)
  //   cerr << i << " ";
  // cerr << endl;
  // if (verify(K, info_bits.data(), decoded.data()))
  if (verify(N, codeword.data(), denoised_codeword.data()))
    cerr << "decode successfully" << endl;
  else 
    cerr << "decoding failed" << endl;
  */
  return 0;
}