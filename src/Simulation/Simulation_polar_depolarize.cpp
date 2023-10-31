#include "Simulation/Simulation.hpp"
using namespace std;

void simulation_polar_depolarize(int N, int K, int L, double p, int min_samples, CONSTRUCTION con, double beta, int seed, int print_interval) {
    int Kq = 2*K - N;
    vector<int> info_bits(K);
    vector<int> codeword_X(N), codeword_Z(N);
    vector<int> noisy_codeword_X(N), noisy_codeword_Z(N);
    vector<int> noise_X(N), noise_Z(N);
    vector<vector<double>> p_noisy_codeword(N, vector<double>(4,0));
    vector<vector<int>> denoised_codeword(N, vector<int>(2,0));
    vector<int> denoised_X(N, 0), denoised_Z(N, 0);
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> X_stab_info_bits(N, 0);
    cerr << "N=" << N << ", K=" << K << ", L=" << L << ", p=" << p << endl;
    Channel_depolarize_q* chn_depolarize = new Channel_depolarize_q(N, p, seed);
    construct_frozen_bits(con, N, K, K, Z_code_frozen_bits, X_stab_info_bits, beta);
    print_mixing_factor(Z_code_frozen_bits);

    vector<int> info_indices; 
    for (int i = 0; i < N; i++)
        if (!X_stab_info_bits[i] && !Z_code_frozen_bits[i]) 
            info_indices.push_back(i);
    cerr << "info indices: ";
    for (auto i : info_indices) cerr << i << " ";
    cerr << endl;

    vector<int> info_indices_K; // in V_K
    int j = 0;
    for (int i = 0; i < N; i++)
        if (X_stab_info_bits[i] && !Z_code_frozen_bits[i]) ++j;
        else if (!X_stab_info_bits[i] && !Z_code_frozen_bits[i]) info_indices_K.push_back(j++);
    cerr << "info indices in V_K: ";
    for (auto i : info_indices_K) cerr << i << " ";
    cerr << endl;

    vector<int> info_X(Kq, 0), info_Z(Kq, 0);
    vector<int> info_X_hat(Kq, 0), info_Z_hat(Kq, 0);

    Encoder_polar* encoder = new Encoder_polar(K, N, Z_code_frozen_bits);
    Decoder_depolarize* decoder = new Decoder_depolarize(K, N, L, Z_code_frozen_bits);
    int SCL_num_err = 0, ML_err = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(seed);
    std::bernoulli_distribution d(0.5);

    int num_total = 0;
    while (1) {
        num_total++;
        for (int i = 0; i < K; i++) info_bits[i] = d(gen);
        for (int i = 0; i < Kq; i++) info_X[i] = info_bits[info_indices_K[i]];
        encoder->encode(info_bits.data(), codeword_X.data(), 0);
        for (int i = 0; i < K; i++) info_bits[i] = d(gen);
        for (int i = 0; i < Kq; i++) info_Z[i] = info_bits[info_indices_K[i]];
        encoder->encode(info_bits.data(), codeword_Z.data(), 0);
        bit_reversal(codeword_Z);
        chn_depolarize->add_noise(noise_X.data(), noise_Z.data(), 0);
        xor_vec(N, codeword_X, noise_X, noisy_codeword_X);
        xor_vec(N, codeword_Z, noise_Z, noisy_codeword_Z);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 4; j++) p_noisy_codeword[i][j] = p/3;
            p_noisy_codeword[i][noisy_codeword_Z[i]*2 + noisy_codeword_X[i]] = 1-p; 
        }
        decoder->decode(p_noisy_codeword, denoised_codeword);
        // decoder->decode_SC(p_noisy_codeword, denoised_codeword);
        for (int i = 0; i < N; i++) denoised_X[i] = denoised_codeword[i][0];
        for (int i = 0; i < N; i++) denoised_Z[i] = denoised_codeword[i][1];
        bit_reversal(denoised_Z);
        encoder->light_encode(denoised_X.data());
        encoder->light_encode(denoised_Z.data());
        for (int i = 0; i < Kq; i++) info_X_hat[i] = denoised_X[info_indices[i]];
        for (int i = 0; i < Kq; i++) info_Z_hat[i] = denoised_Z[info_indices[i]];
        if (!verify(Kq, info_X, info_X_hat) || !verify(Kq, info_Z, info_Z_hat)) SCL_num_err++;
        if (num_total >= min_samples && SCL_num_err >= 100 && num_total % print_interval == 0) break;
        if (num_total % print_interval == 0) cerr << "#err: " << SCL_num_err << " / " << num_total << endl;
    }
    cerr << "#err: " << SCL_num_err << " / " << num_total << endl;
    cerr << "Polar depolarize SCL FER: " << (double)SCL_num_err / num_total << endl;
    // cerr << "Polar ML FER: " << (double)ML_err / num_total << endl;
}
