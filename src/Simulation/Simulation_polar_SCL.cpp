#include "Simulation/Simulation.hpp"
using namespace std;
void simulation_polar_SCL(int N, int K, int L, double p, double db, double design_snr) {
    // int K = 92, N = 256, L = 64;
    // int K = 16, N = 32, L = 1;
    // int K = 1024, N = 2048, L = 4;
    // int K = 1723, N = 2048, L = 4;
    vector<int> info_bits(K);
    vector<int> codeword(N);
    vector<int> noisy_codeword(N);
    vector<double> llr_noisy_codeword(N, 0);
    vector<int> denoised_codeword(N);
    vector<bool> frozen_bits(N, 0);
#ifdef CHN_AWGN      // flag in Simulation.hpp
    // double db = 1;
    // double design_snr = 1;
    double code_rate = (double)K / N;
    double sigma = 1 / sqrt(2 * code_rate * db2val(db));
    double design_sigma = 1 / sqrt(2 * code_rate * db2val(design_snr));
    cerr << "sigma=" << sigma << ", db=" << db << endl;
    cerr << "design sigma=" << design_sigma << ", design snr db=" << design_snr << endl;
    Channel_AWGN* chn_awgn = new Channel_AWGN(N, sigma, design_sigma, 42);
    // frozen_bits_generator_AWGN(N, K, design_snr, frozen_bits);
    // frozen_bits_generator_AWGN_SC(N, K, design_snr, frozen_bits);
    frozen_bits_generator_PW(N, K, frozen_bits);
#else
    // double p = 0.1;
    Channel_BSC* chn_bsc = new Channel_BSC(N, p, 42);
    // frozen_bits_generator_BEC(N, K, p, frozen_bits);
    // frozen_bits_generator_BSC_SC(N, K, p, frozen_bits);
    frozen_bits_generator_PW(N, K, frozen_bits);
    // frozen_bits_generator_RM(N, K, frozen_bits);
#endif
    int cnt = 0;
    // for (auto i : frozen_bits) if (i) cnt++;
    // cerr << "frozen bits mask contains " << cnt << " ones : ";
    // for (auto i : frozen_bits) cerr << i;
    cerr << endl;
    // vector<bool> frozen_bits = {1,0,1,0,1,0,0,0}; // info bits at 8,4,6,7,2
    Encoder_polar* encoder = new Encoder_polar(K, N, frozen_bits);
    Decoder_polar_SCL* decoder = new Decoder_polar_SCL(K, N, L, frozen_bits);
    int num_total = 1000, SCL_num_err = 0;
    for (int k = 0; k < num_total; k++) {
        generate_random(K, info_bits.data());
        // cerr << "info bits: ";
        // for (int i : info_bits) cerr << i;
        // cerr << endl;
        encoder->encode(info_bits.data(), codeword.data(), 0);
        // cerr << "codeword: ";
        // for (int i : codeword) cerr << i;
        // cerr << endl;
        // cerr << endl << "is codeword? " << encoder->is_codeword(codeword.data()) << endl;
        assert(encoder->is_codeword(codeword.data()));
#ifdef CHN_AWGN
        chn_awgn->add_noise(codeword.data(), llr_noisy_codeword.data(), 0);
#else
        int ml_flips = chn_bsc->add_noise(codeword.data(), noisy_codeword.data(), 0);
        // cerr << endl << ml_flips << " flips, noisy codeword: ";
        // for (int i : noisy_codeword) cerr << i;
        // cerr << endl;
        for (int i = 0; i < N; i++)
            // llr_noisy_codeword[i] = codeword[i] ? -1.0 : 1.0;
            llr_noisy_codeword[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
#endif
        decoder->decode(llr_noisy_codeword.data(), denoised_codeword.data(), 1);
        // decoder->decode_SC(llr_noisy_codeword.data(), denoised_codeword.data(), 1);
        assert(encoder->is_codeword(denoised_codeword.data()));
        if (!verify(N, codeword.data(), denoised_codeword.data())) {
            // cerr << "decode fail, num_flips: " << count_flip(N, noisy_codeword.data(), denoised_codeword.data()) << " , ml_flips: " << ml_flips << endl;
            SCL_num_err++;
        } 
        // else cerr << k << " decode successfully" << endl;
        // else cerr << k << " decode successfully, ml_flips: " << ml_flips << endl;
    }
    cerr << "Polar SCL FER: " << (double)SCL_num_err / num_total << endl;
}

void simulation_stab_MW_codewords(int N, int K)
{

    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> X_stab_frozen_bits(N, 1);
    vector<bool> worst_frozen_bits(N, 1);
    vector<int> codeword_Z(N);
    vector<int> stab_X(N);
    vector<int> codeword(N);
    frozen_bits_generator_PW(N, K, Z_code_frozen_bits);
    for (int i = 0; i < N; i++) 
        if (Z_code_frozen_bits[i] == 0 && Z_code_frozen_bits[N-1-i] == 1) 
            X_stab_frozen_bits[i] = 0;

    for (int i = 0; i < N; i++)
        worst_frozen_bits[i] = Z_code_frozen_bits[i] ? 0 : 1;
    Encoder_polar* Z_code = new Encoder_polar(K, N, Z_code_frozen_bits);
    Encoder_polar* X_stab = new Encoder_polar(N-K, N, X_stab_frozen_bits);
    Encoder_polar* stab_worst = new Encoder_polar(N-K, N, worst_frozen_bits);

    vector<int> info_bits_Z(K);
    vector<int> info_bits_X_stab(N-K);

    int num_total = 10000;
    vector<int> wt(num_total);
    vector<int> wt_worst(num_total);

    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        generate_random(K, info_bits_Z.data());
        Z_code->encode(info_bits_Z.data(), codeword_Z.data(), 0);
        generate_random(N-K, info_bits_X_stab.data());
        X_stab->encode(info_bits_X_stab.data(), stab_X.data(), 0);
        stab_worst->encode(info_bits_X_stab.data(), codeword.data(), 0);

        int cnt = 0;
        for (int i = 0; i < N; i++)
            if (codeword_Z[i] && stab_X[i]) cnt++;
        
        assert (cnt % 2 == 0);

        wt[turn_idx] = count_weight(stab_X);
        wt_worst[turn_idx] = count_weight(codeword);
    }
    cerr << "best  N-K wd ";
    print_wt_dist(wt);
    cerr << "worst N-K wd ";
    print_wt_dist(wt_worst);


}