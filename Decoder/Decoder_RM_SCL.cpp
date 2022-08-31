
// Recursive decoding on triangle down to (1, m) and (m, m).
int
tri1_dec(
   double *y_in,
   int r,
   int m,
   int n,
   double *y_dec,
   int *x_dec,
   int *x_dec_len,
   double *aux_buf
)
{
   int i, n2, cur_x_dec_len;
   int tmp[4096];
   double *y1, *y2, *y_dec1, *y_dec2;

   // RM(m, m) case.
   if (r == m) {

      (*x_dec_len) = n; // n information bits.

      for (i = 0; i < n; i++) {
         if (y_in[i] > 0) {
            tmp[i] = 0;
            y_dec[i] = 1.0;
         }
         else {
            tmp[i] = 1;
            y_dec[i] = -1.0;
         }
      }
      mrm_enc_mm_bsc(m, tmp, x_dec);

      return 0;

   } // if (r == m).

   // RM(1, m) case.
   if (r == 1) {
      (*x_dec_len) = m + 1; // # of information bits.
#ifdef FORMAT_EPS4_H
      if (rm1_dec_sh_dotp(m, y_in, x_dec)) return 1;
#else
      if (rm1_dec_sh_finp(m, y_in, x_dec)) return 1;
#endif
      mrm_enc_bsc(m, r, x_dec, tmp);
      for (i = 0; i < n; i++) y_dec[i] = (double)(1 - (tmp[i] << 1));
      return 0;
   } // if (r == 1).

   // Further reduction.

   n2 = n / 2;
   y1 = y_in;
   y2 = y_in + n2;
   y_dec1 = y_dec;
   y_dec2 = y_dec + n2;

   // Decode u1.
   // [x_dec1, y_dec1] = dtrm01_she(r - 1, m - 1, y1 .* y2);
   for (i = 0; i < n2; i++) aux_buf[i] = y1[i] * y2[i];
   tri1_dec(
      aux_buf, r - 1, m - 1, n2,
      y_dec1, x_dec, &cur_x_dec_len, aux_buf + n2
   );
   (*x_dec_len) = cur_x_dec_len;

   // Decode u2.
   // [x_dec2, y_dec2] = dtrm01_she(r, m - 1, add_est_eps(y1 .* y_dec1, y2, n2));
   for (i = 0; i < n2; i++) y1[i] *= y_dec1[i];
   VADD_EST(y1, y2, aux_buf, n2);
   // ADD_EST_EPS(y1, y2, aux_buf, n2);
   tri1_dec(
      aux_buf, r, m - 1, n2,
      y_dec2, x_dec + cur_x_dec_len, &cur_x_dec_len, aux_buf + n2
   );
   (*x_dec_len) += cur_x_dec_len;

   // y_dec = [y_dec1 .* y_dec2, y_dec2];
   for (i = 0; i < n2; i++) y_dec1[i] *= y_dec2[i];

   return 0;
}

//-----------------------------------------------------------------------------
// codec interface functions.

int
cdc_init(
   char param_str[],
   void **cdc
)
{
   cdc_inst_type *dd = NULL;
   char *str1 = NULL, *token;
   int c_n = 0, c_k = 0, c_m = 0, c_r = 0;

   // Allocate decoder instance.
   dd = (cdc_inst_type *)malloc(sizeof(cdc_inst_type));
   if (dd == NULL) {
      show_msg("cdc_init: Short of memory!");
      goto ret_err;
   }
   memset(dd, 0, sizeof(cdc_inst_type));

   // Allocate temporary string for parsing.
   str1 = (char *)malloc(SP_STR_MAX);
   if (str1 == NULL) {
      err_msg("cdc_init: Short of memory!");
      goto ret_err;
   }
   str1[0] = 0;

   // ----- Parse parameters.
   strcpy(str1, param_str);

   token = strtok(str1, tk_seps_prepared);
   while (token != NULL) {
      TRYGET_INT_TOKEN(token, "RM_m", c_m);
      TRYGET_INT_TOKEN(token, "RM_r", c_r);
      SPF_SKIP_UNKNOWN_PARAMETER(token);
   }

   free(str1); str1 = NULL;

   dd->c_m = c_m;
   dd->c_r = c_r;

   // Check if all necessary parameters are specified.
   if ((c_m == 0) || (c_r == 0)) {
      err_msg("cdc_init: RM_m or RM_r are not specified.");
      goto ret_err;
   }

   // Complete settings.
   dd->c_n = c_n = 1 << c_m;
   dd->c_k = c_k = calc_rm_k(c_m, c_r);

#ifdef FORMAT_Q
   if (dec_tables_init()) {
      err_msg("cdc_init: error in dec_tables_init().");
      goto ret_err;
   }
#endif // FORMAT_Q

   (*cdc) = dd;

   return 0;

ret_err:

   if (str1 != NULL) free(str1);
   if (dd != NULL) {
      free(dd);
      (*cdc) = NULL;
   }
   return 1;
}

void
cdc_close(
   void *cdc
)
{
   cdc_inst_type *dd;

   if (cdc == NULL) return;
   dd = (cdc_inst_type *)cdc;
   free(dd);
#ifdef FORMAT_Q
   dec_tables_free();
#endif // FORMAT_Q
}

int cdc_get_n(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->c_n;
}

int cdc_get_k(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->c_k;
}

void
cdc_set_sg(
   void *cdc,
   double noise_sg
)
{
   cdc_inst_type *dd;

   dd = (cdc_inst_type *)cdc;
   dd->sg22 = 2.0 / (noise_sg * noise_sg);
}

int
enc_bpsk(
   void *cdc,
   int x[],
   double y[]
)
{
   cdc_inst_type *dd;

   dd = (cdc_inst_type *)cdc;
   return mrm_enc_bpsk(dd->c_m, dd->c_r, x, y);
}

int
dec_bpsk(
   void *cdc,
   double c_out[],
   int x_dec[]
)
{
   cdc_inst_type *dd;
   double *dec_input, *y_dec, *aux_buf;
   // int max_peak_lsiz;
   void *mem_buf;
   int i;

   dd = (cdc_inst_type *)cdc;

   mem_buf = (double *)malloc(dd->c_n * 3 * sizeof(double));
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }

   dec_input = (double *)mem_buf;
   y_dec = dec_input + dd->c_n;
   aux_buf = y_dec + dd->c_n;
   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   VY_TO_FORMAT(dec_input, dec_input, dd->c_n);

   tri1_dec(dec_input, dd->c_r, dd->c_m, dd->c_n, y_dec, x_dec, &i, aux_buf);

   free(mem_buf);

   return 0;
}