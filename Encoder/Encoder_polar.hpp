#ifndef ENCODER_POLAR_HPP_
#define ENCODER_POLAR_HPP_

#include <vector>
using namespace std;

class Encoder_polar : public Encoder
{
protected:
  const int    m;
  vector<bool> frozen_bits;
  vector<int>  X_N_tmp;

public:
  Encoder_polar(const int& K, const int& N, const vector<bool>& frozen_bits);
  virtual ~Encoder_polar() = default;
  
  virtual Encoder_polar* clone() const;

  void light_encode(int *bits);

  bool is_codeword(const int *X_N);

  virtual const vector<bool>& get_grozen_bits() const;
  virtual void set_frozen_bits(const vector<bool>& frozen_bits);

protected:
  virtual void _encode(const int *U_K, B *X_N, const size_t frame_id);
  void convert(const int *U_K, const *U_N);

};

#endif
