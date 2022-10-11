#include "Simulation/Simulation.hpp"
int simulation_RM_closest_codewords(int m, int rx, int rz) {
    Encoder_RM_CSS* encoder = new Encoder_RM_CSS(m, rx, rz);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", rx="<< rx << ", rz=" << rz
         << ", K=" << K << ", N=" << N << endl;
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_X_diff(N, 0);
    vector<int> desired_X(N, 0);
    vector<vector<int>> codewords;
    generate_all_codewords(m, rz, codewords);
    int num_codewords = codewords.size();
    cerr << "have generated all " << num_codewords << " codewords" << endl;
    int num_total = 2000;
    vector<int> min_flip_indices;
    int min_num_flips, num_flips;
 
    vector<double> pxs;
    for (int i = 9; i < 40; i++) pxs.push_back((double)(i + 1) / 100.0);
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
                }
                cerr << "finish";
        
                cerr << endl;
            }

        }
    }
    return 0;
}