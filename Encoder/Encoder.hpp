#ifndef ENCODER_HPP
#define ENCODER_HPP

#include <string>
using namespace std;

class Encoder
{
public:
  int K = 0;
  int N = 0;
  
  explicit Encoder();
  Encoder(const int K, const int T);
  
  virtual ~Encoder() = default;
  virtual Encoder* clone() const;

  Encoder* build() const;

protected:
  Encoder(const string &n, const string &p);

};

#endif
