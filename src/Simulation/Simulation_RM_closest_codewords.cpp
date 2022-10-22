#include "Simulation/Simulation.hpp"
#include <cstddef>
#include <cstring>
int simulation_RM_closest_codewords(int m, int rx, int rz) {
    Encoder_RM_CSS* encoder = new Encoder_RM_CSS(m, rx, rz);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", rx="<< rx << ", rz=" << rz
         << ", K=" << K << ", N=" << N << endl;
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_X_diff(N, 0);
    vector<int> desired_X(N, 0);
    vector<vector<int>> codewords;
    vector<vector<int>> stabilizers;
    generate_all_codewords(m, rz, codewords);
    generate_all_codewords(m, m-rx-1, stabilizers);
    int num_codewords = codewords.size(), num_stab = stabilizers.size();
    cerr << "have generated all " << num_codewords << " codewords and " << num_stab << " stabilizers." << endl;
    int num_total = 2000;
    vector<int> min_flip_indices;
    int min_num_flips, num_flips;
    vector<int> flips(N + 1, 0);

    vector<double> pxs;
    for (int i = 9; i < 40; i++) pxs.push_back((double)(i + 5) / 100.0);
    for (double px : pxs) {
        cerr << "px = " << px << endl;
        Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, px, 0.0, 41);
        for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
            chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0); // only use noise_X for now
            // ML_err is just take the codeword with smallest flip, 
            // see whether it is in the all-zero's equiv class (stabilizers)
            min_num_flips = N; min_flip_indices.clear();
            for (int i = 0; i < num_codewords; i++) {
                num_flips = count_flip(N, noise_X.data(), codewords[i].data());
                if (num_flips < min_num_flips) {
                    min_flip_indices.clear();
                    min_flip_indices.push_back(i);
                    min_num_flips = num_flips;
                } else if (num_flips == min_num_flips) {
                    min_flip_indices.push_back(i);
                }
            }
            // cerr << "there are " << min_flip_indices.size() << " codewords that have min flips" << endl;
            if (min_flip_indices.size() > 1) {
                cerr << "min_num_flips=" << min_num_flips <<
                ", min_flip_indices size=" << min_flip_indices.size() << endl;
                cerr << "noise is: " << endl;;
                for (int n : noise_X) cerr << n;
                // vector<int> one_indices;
                // for (int i = 0; i < N; i++) {
                //     if (noise_X[i] == 1) one_indices.push_back(i);
                // }
                cerr << " , neighboring codewords are: " << endl;
                for (int i : min_flip_indices) {
                    // for (int n : one_indices) cerr << codewords[i][n];
                    for (int n : codewords[i]) cerr << n;
                    cerr << endl;                
                    for (int k = 0; k < N + 1; k++) flips[k] = 0;
                    for (auto& s : stabilizers) {
                        xor_vec(N, codewords[i].data(), s.data(), noise_X_diff.data());
                        num_flips = count_flip(N, noise_X.data(), noise_X_diff.data());
                        flips[num_flips]++;
                    }
                    for (int n : flips) cerr << n << " ";
                    cerr << endl;
                }
                cerr << "finish";

        
                cerr << endl;
            }

        }
    }
    return 0;
}

int simulation_RM_d_star(int m, int r) {
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r=" << r << ", K=" << K << ", N=" << N << ", d=" << (1 << (m-r)) << endl;
    vector<int> all_zero(N, 0);
    vector<int> noise(N, 0);
    vector<vector<int>> codewords(1<<K, vector<int>(N));
    generate_all_codewords(m, r, codewords);
    int num_codewords = codewords.size();
    cerr << "have generated all " << num_codewords << " codewords" << endl;
    int num_total = 2000;
    vector<int> min_flip_indices;
    int min_num_flips, num_flips;

    // **********************
    int d_star = 0;
    int d_star_star = 0;
    // **********************
 
    vector<double> pxs;
    for (int i = 2; i <= 20; i+=1) pxs.push_back((double)(i) / 100.0);
    for (double px : pxs) {
        cerr << "px = " << px << endl;
        Channel_BSC* chn_bsc = new Channel_BSC(N, px, 42);
        for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
            chn_bsc->add_noise(all_zero.data(), noise.data(), 0); 
            // ML_err is just take the codeword with smallest flip, 
            min_num_flips = N; min_flip_indices.clear();
            for (int i = 0; i < num_codewords; i++) {
                num_flips = count_flip(N, noise.data(), codewords[i].data());
                if (num_flips < min_num_flips) {
                    min_flip_indices.clear();
                    min_flip_indices.push_back(i);
                    min_num_flips = num_flips;
                } else if (num_flips == min_num_flips) {
                    min_flip_indices.push_back(i);
                }
            }
            cerr << "noise ";
            for (int l : noise) cerr << l;
            cerr << endl << "min_num_flips = " << min_num_flips << ", size = " << min_flip_indices.size() << endl;
            // cerr << "there are " << min_flip_indices.size() << " codewords that have min flips" << endl;
            if (min_flip_indices.size() > 1) {
                if (min_num_flips > d_star_star) {
                    d_star_star = min_num_flips;
                    cerr << "updated d_star_star=" << d_star_star << ", noise ";
                    for (int l : noise) cerr << l;
                    cerr << endl;
                }
            } else {
                if (min_num_flips > d_star) {
                    d_star = min_num_flips;
                    cerr << "updated d_star=" << d_star << ", noise ";
                    for (int l : noise) cerr << l;
                    cerr << endl;
                }
            }
        }
    }
    return 0;
}

int test_RM_d_star() {
    int m = 6, r = 2;
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N(), N_half = N >> 1;
    cerr << "For m=" << m << ", r=" << r << ", K=" << K << ", N=" << N << ", d=" << (1 << (m-r)) << endl;
    vector<int> all_zero(N, 0);
    vector<int> noise(N, 0);
    vector<vector<int>> codewords(1<<K, vector<int>(N));
    generate_all_codewords(m, r, codewords);
    int num_codewords = codewords.size();
    cerr << "have generated all " << num_codewords << " codewords" << endl;
    vector<int> min_flip_indices;
    int min_num_flips, num_flips;

    string noise_str = "0000000010001000000000000000000100110010000111111000001000101000";
    assert (noise_str.length() == N);
    for (int i = 0; i < N; i++) {
        noise[i] = int(noise_str[i]) - int('0');
    }

    // string u_str = "10000000001000110100000000001000";
    // string v_str = "00110010000111111000001000101000";
    // assert (u_str.length() == N_half);
    // for (int i = 0; i < N_half; i++) {
    //     noise[i] = int(u_str[i]) - int('0');
    //     noise[i + N_half] = int(v_str[i]) - int ('0');
    //     noise[i + N_half] ^= noise[i];
    // }

    for (int i : noise) cerr << i;
    cerr << endl;
    min_num_flips = N; min_flip_indices.clear();
    for (int i = 0; i < num_codewords; i++) {
        num_flips = count_flip(N, noise.data(), codewords[i].data());
        if (num_flips < min_num_flips) {
            min_flip_indices.clear();
            min_flip_indices.push_back(i);
            min_num_flips = num_flips;
        } else if (num_flips == min_num_flips) {
            min_flip_indices.push_back(i);
        }
    }
    cerr << "min_flip_indices size " << min_flip_indices.size() << endl;
    if (min_flip_indices.size() > 1) {
        cerr << "d_star_star=" << min_num_flips << endl;
        for (int idx : min_flip_indices) {
            for (int l : codewords[idx]) cerr << l;
            cerr << endl;
        }
    } else {
        cerr << "d_star=" << min_num_flips << endl;
        for (int l : codewords[min_flip_indices[0]]) cerr << l;
        cerr << endl;
        cerr << "diffenrece to noise" << endl;
        vector<int> noise_diff(N, 0);
        xor_vec(N, codewords[min_flip_indices[0]].data(), noise.data(), noise_diff.data());
        for (int l : noise_diff) cerr << l;
        cerr << endl << "binary representation of the one's" << endl;
        vector<int> binary_repr(m, 0);
        for (int i = 0; i < N; i++) 
            if (noise_diff[i]) {
                decimal2bianry(i, binary_repr);
                for (int l : binary_repr) cerr << l;
                cerr << endl;
            }

    }
    
    return 0;
}

int compare_equiv_classes() {
    int m = 5, r = 2;
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N(), N_half = N >> 1;
    cerr << "For m=" << m << ", r=" << r << ", K=" << K << ", N=" << N << ", d=" << (1 << (m-r)) << endl;
    vector<int> all_zero(N, 0);
    vector<int> noise(N, 0);
    vector<int> repr(N, 0);
    vector<int> temp(N, 0);
    vector<vector<int>> stabilizers(1<<K, vector<int>(N));
    generate_all_codewords(m, r, stabilizers);
    int num_stab = stabilizers.size();
    cerr << "have generated all " << num_stab << " stabilizers" << endl;
    vector<int> min_flip_indices;
    int num_flips; 
    string noise_str = "00000000000000000000000000000001";
    string repr_str  = "00000000000000000000000000110011";
    assert (noise_str.length() == N);
    for (int i = 0; i < N; i++) {
        noise[i] = int(noise_str[i]) - int('0');
        repr[i] = int(repr_str[i]) - int('0');
    }
    assert(Encoder_RM::is_codeword(repr.data(), m, r+1));

    vector<int> flips1(N + 1, 0); // for \bar{0}
    vector<int> flips2(N + 1, 0); // for \bar{repr}
    
    for (int i = 0; i < num_stab; i++) {
        num_flips = count_flip(N, noise.data(), stabilizers[i].data());
        flips1[num_flips]++;
        xor_vec(N, repr.data(), stabilizers[i].data(), temp.data());
        num_flips = count_flip(N, noise.data(), temp.data());
        flips2[num_flips]++;
    }
    cerr << "flips1: ";
    for (int l : flips1) cerr << l << " ";
    cerr << endl << "flips2: ";
    for (int l : flips2) cerr << l << " ";
    cerr << endl;

    vector<double> flips_err_prob(N + 1, 0);
    double px = 0.5;
    // double log_err = N * log(1-px), log_err_diff = log(px) - log(1-px);
    // flips_err_prob[0] = exp(log_err);
    // for (int i = 1; i < N + 1; i++) {
    //     log_err += log_err_diff;
    //     flips_err_prob[i] = exp(log_err);
    //     // c++ double exponent min is -1022, may underflow
    // }
    for (int i = 0; i < N + 1; i++) {
        flips_err_prob[i] = 1 << (32-i);
    }
    double ec_err_prob = 0.0; 
    for (int k = 0; k < N+1; k++)
        ec_err_prob += flips1[k] * flips_err_prob[k];
    cerr << "flips1 error prob: " << ec_err_prob << endl;
    ec_err_prob = 0.0; 
    for (int k = 0; k < N+1; k++)
        ec_err_prob += flips2[k] * flips_err_prob[k];
    cerr << "flips2 error prob: " << ec_err_prob << endl;

    return 0;
}

int simulation_symmetric_noise(int m, int r) {
    Encoder_RM* encoder = new Encoder_RM(m, r);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r=" << r << ", K=" << K << ", N=" << N << ", d=" << (1 << (m-r)) << endl;
    vector<int> all_zero(N, 0);
    vector<int> noise(N, 0);
    vector<vector<int>> codewords(1<<K, vector<int>(N));
    generate_all_codewords(m, r, codewords);
    int num_codewords = codewords.size();
    cerr << "have generated all " << num_codewords << " codewords" << endl;
    vector<int> min_flip_indices;
    int min_num_flips, num_flips;
    vector<int> flips(N + 1, 0);

    for (int level = 1; level <= (m/2); level++) {
        for (int i = 0; i < flips.size(); i++) flips[i] = 0;
        cerr << "level = " << level << endl;
        generate_symmetric_noise(m, noise, level);
        cerr << "noise has weight " << count_weight(noise) << endl;
        min_num_flips = N; min_flip_indices.clear();
        for (int i = 0; i < num_codewords; i++) {
            num_flips = count_flip(N, noise.data(), codewords[i].data());
            flips[num_flips]++;
            if (num_flips < min_num_flips) {
                min_flip_indices.clear();
                min_flip_indices.push_back(i);
                min_num_flips = num_flips;
            } else if (num_flips == min_num_flips) {
                min_flip_indices.push_back(i);
            }
        }
        cerr << "min_num_flips=" << min_num_flips << endl;
        for (int f : flips) cerr << f << " ";
        cerr << endl;
    }

    // generate_symmetric_noise(m, noise, 1);
    // int noise_weight = count_weight(noise);
    // cerr << "noise has weight " << noise_weight << endl;
    // double p = (double)noise_weight / N;
    // Decoder_RM_SCL* SCL_decoder = new Decoder_RM_SCL(m, r+1, 1024);
    // vector<double> llr_noisy_codeword(N, 0);
    // for (int i = 0; i < N; i++) {
    //     llr_noisy_codeword[i] = noise[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
    // }
    // vector<int> SCL_denoised_codeword(N, 0);
    // SCL_decoder->decode(llr_noisy_codeword.data(), SCL_denoised_codeword.data(), 0);
    // vector<vector<int>> c_list(1024, vector<int>(N));
    // vector<double> pm_list(1024);
    // SCL_decoder->copy_codeword_list(c_list, pm_list);
    // for (int l : SCL_denoised_codeword) cerr << l;
    // cerr << endl;
    // cerr << "SCL decoded codeword is distance " << count_flip(N, noise.data(), SCL_denoised_codeword.data()) << " from noise" << endl;
    // for (auto& c : c_list) cerr << count_flip(N, c.data(), noise.data()) << " ";
    // cerr << endl;
    return 0; 
}
// want to perform RM(6,4)/RM(6,3), need a clever way
// use list decoder to obtain the codewords that are within distance 1.5d to the noise
// only add the minimum weight stabilizers to them (need a generator for this)

