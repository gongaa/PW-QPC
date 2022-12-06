#include "Simulation/Simulation.hpp"

int simulation_RM_SCL(int m, int r, int list_size, double p, double db, double design_snr, int num_total) {
    // each time a decoding failure occurred, 
    // we checked whether the decoded codeword was more likely than the transmitted codeword.
    // If so, then the optimal ML decoder would surely misdecode y as well.
    // int m = 11, r = 5, list_size = 4;
    Encoder_RM* encoder = new Encoder_RM(m, r);
    Decoder_RM_SC* SC_decoder = new Decoder_RM_SC(m ,r, 1);
    Decoder_RM_SCL* SCL_decoder = new Decoder_RM_SCL(m ,r, list_size);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
    cerr << "List size=" << list_size << endl;
#ifdef CHN_AWGN
    // double db = 4;
    double code_rate = (double)K / N;
    double sigma = 1 / sqrt(2 * code_rate * db2val(db));
    cerr << "sigma=" << sigma << ", db=" << db << endl;
    Channel_AWGN* chn_awgn = new Channel_AWGN(N, sigma, 42);
#else
    // double p = 0.1;
    cerr << "p=" << p << endl;
    Channel_BSC* chn_bsc = new Channel_BSC(N, p, 42);
#endif // USE_AWGN
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    vector<int> noisy_codeword(N, 0);
    vector<double> llr_noisy_codeword(N, 0);
    vector<double> dumer_noisy_codeword(N, 0);
    vector<int> SC_denoised_codeword(N, 0);
    vector<int> SCL_denoised_codeword(N, 0);
    vector<int> decoded(K, 0);
    int SC_num_err = 0, SCL_num_err = 0, ML_err = 0;
    int SC_num_flips = 0, SCL_num_flips = 0, num_flips=0;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
#ifdef CHN_AWGN
        chn_awgn->add_noise(codeword.data(), llr_noisy_codeword.data(), 0);
#else
        num_flips = chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
        for (int i = 0; i < N; i++) {
            llr_noisy_codeword[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            // llr_noisy_codeword[i] = noisy_codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
        }
#endif // USE_AWGN
        SC_decoder->decode(llr_noisy_codeword.data(), SC_denoised_codeword.data(), 0);
        SCL_decoder->decode(llr_noisy_codeword.data(), SCL_denoised_codeword.data(), 0);
        if (!verify(N, codeword, SC_denoised_codeword)) SC_num_err++;
        if (!verify(N, codeword, SCL_denoised_codeword)) {
            SCL_num_err++;
            if (num_flips >= count_flip(N, noisy_codeword, SCL_denoised_codeword))
                ML_err++;
        }
    }
    cerr << "RM SC FER: " << (double)SC_num_err / num_total << endl;
    cerr << "RM SCL FER: " << (double)SCL_num_err / num_total << endl;
    cerr << "RM ML FER: " << (double)ML_err / num_total << endl;
    return 0;
}
