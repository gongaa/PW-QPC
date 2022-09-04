#ifndef ENCODER_HPP_
#define ENCODER_HPP_

#include <string>
using namespace std;

class Encoder
{
public:
  int K = 0;
  int N = 0;
  
  explicit Encoder();
  Encoder(const int K, const int N) : K(K), N(N) {}
  
  virtual ~Encoder() = default;

  virtual void encode(const int *U_K, int *X_N, const size_t frame_id) = 0;

  virtual int get_K() { return K; }
  virtual int get_N() { return N; }

  Encoder* build() const;

protected:
  Encoder(const string &n, const string &p);

};

#endif
