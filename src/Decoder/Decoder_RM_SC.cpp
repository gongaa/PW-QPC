#include "Decoder/Decoder_RM_SC.hpp"
#include "Encoder/Encoder_RM.hpp"
#include <vector>
using namespace std;

Decoder_RM_SC::Decoder_RM_SC(const int& m, const int& r, const int& L)
: Decoder(Encoder_RM::calculate_K(m, r), 1 << m), m(m), r(r), L(L)
{
}

int Decoder_RM_SC::decode(const double *Y_N, int *V_K, const size_t frame_id)
{  // recursive decoding down to (1, m) and (m, m)
   vector<double> Y_dec_N(N, 0); // TODO: move it to protected attribute
   vector<double> buf(N, 0);
   int V_K_len = K;
   // return this->_decode_dumer(Y_N, Y_dec_N.data(), V_K, m, r, N, V_K_len, buf.data());
   return this->_decode_llr(Y_N, Y_dec_N.data(), V_K, m, r, N, V_K_len);
   // TODO: test whether they are equal
}

int Decoder_RM_SC::_decode_dumer(const double *Y_N, double *Y_dec_N, int *V_K, const int& m, const int& r, const int& N, int& V_K_len, double *aux_buf)
{  // return 0 if succeed, return 1 if fail (e.g. in FHT decoder)
   int i = 0;
   vector<int> tmp(N); 
   if (r == m) {
      V_K_len = N;
      for (i = 0; i < N; i++) {
         if (Y_N[i] > 0) { 
            tmp[i] = 0;
            Y_dec_N[i] = 1.0;
         } else {          
            tmp[i] = 1;
            Y_dec_N[i] = -1.0;
         }
      }
      Encoder_RM::encode_mm_code(tmp.data(), V_K, N);
      return 0; 
   }

   // if (r == 0) {
   //    V_K_len = 1; // one information bit.
   //    double p = 0.0;
   //    for (i = 0; i < N; i++) {
   //       p += log((1.0 + Y_N[i]) / (1.0 - Y_N[i]));
   //    }
   //    if (p > 0.0) {
   //       V_K[0] = 0;
   //       for (i = 0; i < N; i++) 
   //          Y_dec_N[i] = 1.0;
   //    }
   //    else {
   //       V_K[0] = 1;
   //       for (i = 0; i < N; i++) 
   //          Y_dec_N[i] = -1.0;
   //    }
   //    return 0;
   // }

   if (r == 1) {
      V_K_len = m + 1; // {m\choose 1} + {m\choose 0}
      // TODO: use FHT decoder
      _decode_m1(Y_N, V_K, m, N);
      Encoder_RM::_encode(V_K, tmp.data(), m, r);
      for (i = 0; i < N; i++)
         Y_dec_N[i] = (double)(1 - (tmp[i] << 1));
      return 0;
   }

   int N_half = N / 2;
   int curr_V_K_len = 0;
   const double *Y_fst = Y_N, *Y_snd = Y_N + N_half;
   double *Y_dec_fst = Y_dec_N, *Y_dec_snd = Y_dec_N + N_half;
   // decode (u+v) + u = v
   for (i = 0; i < N_half; i++)
      aux_buf[i] = Y_fst[i] * Y_snd[i];
   _decode_dumer(aux_buf, Y_dec_fst, V_K, m-1, r-1, N_half, curr_V_K_len, aux_buf + N_half);
   V_K_len = curr_V_K_len;
   
   // decode u
   double Y_tmp;
   for (i = 0; i < N_half; i++) {
      Y_tmp = Y_fst[i] * Y_dec_fst[i]; // \hat{u} = (u+v) + \hat{v}
      aux_buf[i] = (Y_tmp + Y_snd[i]) / (1.0 + Y_tmp * Y_snd[i]);
   }
   _decode_dumer(aux_buf, Y_dec_snd, V_K + curr_V_K_len, m-1, r, N_half, curr_V_K_len, aux_buf + N_half);
   V_K_len += curr_V_K_len;

   // (u+v, u)
   for (i = 0; i < N_half; i++)
      Y_dec_fst[i] *= Y_dec_snd[i];
   
   return 0;
}

int Decoder_RM_SC::_decode_m1(const double *Y_N, int *V_K, const int& m, const int& N)
{
   int i = 0, j = 0, N_, N_half = N / 2; // counters
   vector<double> S_list(2 * N, 0);      // TODO: why size 2*N instead of N?
   double s1, s2, y1;
   vector<double> Y_list((m+1) * N, 0); 
   for (i = 0; i < N; i++)
      Y_list[i] = Y_N[i];
   double *Y_1 = Y_list.data(), *Y_2, *Y_3, *Y_4;
   for (N_ = N; N_ >= 2; N_ >>= 1) {      // m loops 
      N_half = N_ / 2;
      for (i = 0; i < N; i += N_, Y_1 += N_) { // TODO: Y_1 += N_ or Y_1 += N???
         Y_2 = Y_1 + N_half;
         Y_3 = Y_1 + N; 
         Y_4 = Y_2 + N;
         s1 = 0.0; s2 = 0.0;
         for (j = 0; j < N_half; j++) {
            y1 = Y_1[j] * Y_2[j];
            if (y1 > 0.0) {
               s1 += log(1.0 + y1);          // y = 2q-1 e.g. BSC(1-q)
               s2 += log(1.0 - y1);
            } else {
               s1 += log(1.0 + y1); // TODO: either merge it, or is it the reverse way here?
               s2 += log(1.0 - y1);
            }
            Y_3[j] = (Y_1[j] + Y_2[j]) / (1.0 + Y_1[j] * Y_2[j]);
            Y_4[j] = (-Y_1[j] + Y_2[j]) / (1.0 - Y_1[j] * Y_2[j]);
         }
         S_list[i + N_half] = S_list[i] + s2; // likelihood corresponds to 1
         S_list[i] += s1;                     // likelihood corresponds to 0
      }
   }

   // last end node (0,0)
   s2 = S_list[0] + log(1.0 + Y_1[0]);
   j = 0;
   for (i = 0; i < N; i++) {
      s1 = S_list[i] + log(1.0 + Y_1[0]);
      if (s1 > s2) {
         s2 = s1;
         j = i << 1;
      }
      s1 = S_list[i] + log(1.0 - Y_1[0]);
      if (s1 > s2) {
         s2 = s1;
         j = (i << 1) + 1;
      }
   }

   // the information bits V_K is the bit representation of j
   for (i = m; i >= 0; i--, j >>= 1) 
      V_K[i] = j & 1;

   return 0;
}

int Decoder_RM_SC::_decode_llr(const double *Y_N, double *Y_dec_N, int *V_K, const int& m, const int& r, const int& N, int& V_K_len)
{  // return 0 if succeed, return 1 if fail (e.g. in FHT decoder)
   int i = 0;
   vector<int> tmp(N); 
   if (r == m) {
      V_K_len = N;
      for (i = 0; i < N; i++) {
         if (Y_N[i] > 0) { // LLR > 0, more likely to be 0
            tmp[i] = 0;
            Y_dec_N[i] = 1.0;
         } else {          // LLR <= 0, more likely to be 1
            tmp[i] = 1;
            Y_dec_N[i] = -1.0;
         }
      }
      Encoder_RM::encode_mm_code(tmp.data(), V_K, N);
      return 0; 
   }

   if (r == 0) {
      V_K_len = 1; // one information bit.
      double sum_llr = 0.0;
      for (i = 0; i < N; i++) {
         sum_llr += Y_N[i];
      }
      if (sum_llr > 0.0) { // 1 + (-1)
         V_K[0] = 0;
         for (i = 0; i < N; i++)
            Y_dec_N[i] = 1.0;
      }
      else {
         V_K[0] = 1;
         for (i = 0; i < N; i++)
            Y_dec_N[i] = -1.0;
      }
      return 0;
   }

   // if (r == 1) {
   //    V_K_len = m + 1; // {m\choose 1} + {m\choose 0}
   //    // 2^{m+1} = 2N codewords
   //    // TODO: use FHT decoder
   //    return 0;
   // }

   int N_half = N / 2;
   const double *Y_fst = Y_N, *Y_snd = Y_N + N_half;
   double *Y_dec_fst = Y_dec_N, *Y_dec_snd = Y_dec_N + N_half;
   vector<double> LLRs_buffer(N_half); // store LLR of (u+v) + u = v 
   vector<double> LLRs(2);
   vector<int> bits(1);
   for (i = 0; i < N_half; i++) {
      LLRs = {Y_fst[i], Y_snd[i]};
      LLRs_buffer[i] = lambdas[0](LLRs, bits);
   }
   int curr_V_K_len = 0;
   _decode_llr(LLRs_buffer.data(), Y_dec_fst, V_K, m-1, r-1, N_half, curr_V_K_len);
   V_K_len = curr_V_K_len;

   // decode u
   for (i = 0; i < N_half; i++) {
      bits[0] = (Y_dec_fst[i] > 0) ? 0 : 1;
      LLRs = {Y_fst[i], Y_snd[i]};
      LLRs_buffer[i] = lambdas[1](LLRs, bits); // \hat{u} = (u+v) + \hat{v}
   } 
   _decode_llr(LLRs_buffer.data(), Y_dec_snd, V_K + curr_V_K_len, m-1, r, N_half, curr_V_K_len);
   V_K_len += curr_V_K_len;

   // first half = \hat{u} + \hat{v}
   for (i = 0; i < N_half; i++) {
      LLRs = {Y_dec_fst[i], Y_dec_snd[i]};
      Y_dec_fst[i] = lambdas[0](LLRs, bits);
   }

   return 0;
}

int Decoder_RM_SC::_decode_m1_llr(const double *Y_N, int *V_K, const int& m, const int& N)
{
   return 0;
}