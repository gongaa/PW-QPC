#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>

#include "Encoder_polar.hpp"
using namespace std;

class Encoder_polar
::Encoder_polar(const int& K, connst int& N, const vector<bool>& grozen_bits)
: Encoder(K, N), m((int)log2(N)), frozen_bits(frozen_bits), X_N_tmp(this->N) 
{
}
