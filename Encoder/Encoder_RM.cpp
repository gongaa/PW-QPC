#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>

#include "Encoder_RM.hpp"
using namespace std;

Encoder_RM::Encoder_RM()
: Encoder(K, N)
{

}

// calculate dimension of RM(r, m) code
int calc_rm_k(int m, int r) {
   if (r == 0) return 1;
   if (m == r) return 1 << m;
   return calc_rm_k(m - 1, r - 1) + calc_rm_k(m - 1, r);
}

// Encode (m, m) code (for BSC).
// (Or decode, since the inverse matrix is the same.)
// Modified generator matrix.
int mrm_enc_mm_bsc(
   int m,
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   int n = 1 << m;
   int n2;
   int i;

   if (m == 0) {
      y[0] = x[0];
      return 1;
   }

   n2 = n >> 1;

   mrm_enc_mm_bsc(m - 1, x, y);
   mrm_enc_mm_bsc(m - 1, x + n2, y + n2);
   
   for (i = 0; i < n2; i++) y[i] ^= y[i + n2];

   return n;
}

// Encode RM code (for BSC).
// Modified generator matrix.
int mrm_enc_bsc(
   int m, // RM m parameter.
   int r, // RM r parameter.
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   int n = 1 << m;
   int n2, kv, ku;
   int i;

   if (r == 0) {
      for (i = 0; i < n; i++) y[i] = x[0];
      return 1;
   }

   n2 = n >> 1;

   if (r == m) {
      kv = mrm_enc_bsc(m - 1, m - 1, x, y);
      ku = mrm_enc_bsc(m - 1, m - 1, x + kv, y + n2);
   }
   else {
      kv = mrm_enc_bsc(m - 1, r - 1, x, y);
      ku = mrm_enc_bsc(m - 1, r, x + kv, y + n2);
   }
   for (i = 0; i < n2; i++) y[i] ^= y[i + n2];
   return kv + ku;
}