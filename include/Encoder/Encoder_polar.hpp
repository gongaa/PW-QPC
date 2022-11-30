#ifndef ENCODER_POLAR_HPP_
#define ENCODER_POLAR_HPP_

#include <vector>
#include "Encoder.hpp"
using namespace std;

class Encoder_polar : public Encoder
{
protected:
  const int    m;
  vector<bool> frozen_bits; // a mask of size N, true if current position is frozen
  vector<int>  X_N_tmp;
  vector<uint32_t> info_bits_pos; /*!< Positions of the information bits in the codeword */

public:
  Encoder_polar(const int& K, const int& N, const vector<bool>& frozen_bits);
  virtual ~Encoder_polar() = default;
  
  virtual void encode(const int *U_K, int *X_N, const size_t frame_id);
  void light_encode(int *bits);
  void transpose_encode(int *bits);

  bool is_codeword(const int *X_N);

  virtual const vector<bool> &get_frozen_bits() const { return frozen_bits; }
  virtual void set_frozen_bits(const vector<bool>& frozen_bits);

protected:
  void convert(const int *U_K, int *U_N);
};

#endif
