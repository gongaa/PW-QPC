#include "Simulation/Simulation.hpp"
void simulation_polar_CSS(int N, int K, int list_size, double pz)
{
    vector<int> info_bits_Z(K);
    vector<int> desired_Z(N);
    vector<int> noisy_codeword_Z(N);
    vector<double> llr_noisy_codeword_Z(N);
    vector<int> SCL_denoised_codeword_Z(N);
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> X_stab_frozen_bits(N, 1);

    frozen_bits_generator_PW(N, K, Z_code_frozen_bits);
    for (int i = 0; i < N; i++) 
        if (Z_code_frozen_bits[i] == 0 && Z_code_frozen_bits[N-1-i] == 1) 
            X_stab_frozen_bits[i] = 0;

    Encoder_polar* encoder_Z = new Encoder_polar(K, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(K, N, list_size, Z_code_frozen_bits);
    Encoder_polar* X_stab = new Encoder_polar(N-K, N, X_stab_frozen_bits);

    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, 42);
    chn_bsc_q->set_prob(0, pz);
    vector<double> weight = {1.0, 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] will be set to the real p 
    weight[0] = pz / (1.0-pz);
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_Z_diff(N, 0);
    vector<vector<int>> Z_list(list_size, vector<int>(N, 0));
    vector<double> pm_Z_list(list_size, 0.0);

    int total_flips = 0, num_flips = 0, SCL_num_flips = 0;
    double pm_best;

    vector<vector<int>> equiv_class;
    bool is_in_one_class; int largest_class_size; int* largest_class;
    int desired_class_idx;

    int SCL_num_X_err = 0, SCL_num_Z_err = 0;
    int SCL_num_Z_err_deg = 0, SCL_num_Z_err_equal_weight = 0;
    int SCL_equal_weight_guess_correct = 0;
    bool is_SCL_wrong = false, is_SCL_deg_wrong = false, is_SCL_weighted_deg_wrong = false;

    int num_total = 100;
    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        is_SCL_wrong = false; is_SCL_deg_wrong = false; is_SCL_weighted_deg_wrong = false;
        generate_random(K, info_bits_Z.data());
        encoder_Z->encode(info_bits_Z.data(), desired_Z.data(), 0);
        chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
        num_flips = count_weight(noise_Z);
        total_flips += num_flips;
        xor_vec(N, desired_Z.data(), noise_Z.data(), noisy_codeword_Z.data());
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-pz)/pz) : log((1-pz)/pz); // 0 -> 1.0; 1 -> -1.0
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        SCL_decoder_Z->copy_codeword_list(Z_list, pm_Z_list);
        // cerr << "begin of copied list" << endl;
        // for (auto& c : Z_list) {
        //     for (int k : c) cerr << k;
        //     cerr << endl;
        // }
        // cerr << "end of copied list" << endl;
        if (!verify(N, SCL_denoised_codeword_Z.data(), desired_Z.data())) { is_SCL_wrong = true; SCL_num_Z_err++; }
        xor_vec(N, SCL_denoised_codeword_Z.data(), desired_Z.data(), SCL_denoised_codeword_Z.data());
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z.data(), noise_Z.data());
        if (!X_stab->is_codeword(SCL_denoised_codeword_Z.data())) { is_SCL_deg_wrong = true; SCL_num_Z_err_deg++; }
        if (is_SCL_wrong) cerr << "num_flips: " << num_flips << " , SCL_num_flips: " << SCL_num_flips << endl;
        for (auto& dz : Z_list) xor_vec(N, dz.data(), desired_Z.data(), dz.data());

        if (turn_idx % 10 == 9) {
            cerr << "SCL #err : " << SCL_num_Z_err << " / " << (turn_idx+1) << endl;
            cerr << "SCL #err deg : " << SCL_num_Z_err_deg << " / " << (turn_idx+1) << endl;
        }

        equiv_class.clear();
        equiv_class.push_back({0});
        for (int i = 1; i < list_size; i++) {
            is_in_one_class = false;
            for (auto& ec : equiv_class) {
                xor_vec(N, Z_list[ec[0]].data(), Z_list[i].data(), noise_Z_diff.data());
                if (X_stab->is_codeword(noise_Z_diff.data())) {
                    ec.push_back(i);
                    is_in_one_class = true;
                    break;
                }
            }
            if (!is_in_one_class)
                equiv_class.push_back({i});
        }
        cerr << "there are " << equiv_class.size() << " equiv classes" << endl;
        largest_class_size = 0;
        for (auto& ec : equiv_class) {
            if (ec.size() > largest_class_size) {
                largest_class_size = ec.size();
                largest_class = ec.data();
            }
        }
        cerr << "largest class size: " << largest_class_size << ". Weight distribution ";
        vector<int> wt(largest_class_size);
        for (int i = 0; i < largest_class_size; i++) 
            wt[i] = count_weight(Z_list[largest_class[i]]);
        print_wt_dist(wt);

        for (auto ec : equiv_class) {
            if (X_stab->is_codeword(Z_list[ec[0]].data())) {
                cerr << "stabilizer class size: " << ec.size() << ". Weight distribution ";
                vector<int> wt(ec.size());
                for (int i = 0; i < ec.size(); i++) 
                    wt[i] = count_weight(Z_list[ec[i]]);
                print_wt_dist(wt);
                break;
            }
        }
        // for (int i = 0; i < largest_class_size; i++) {
        //     for (int k : Z_list[largest_class[i]]) cerr << k;
        //     cerr << endl;
        // }


    }
    cerr << "average #flips: " << (double)total_flips / num_total << endl;
    cerr << "SCL #err    : " << SCL_num_Z_err << ", SCL Frame Error Rate: " << (double)SCL_num_Z_err / num_total << endl;
    cerr << "SCL #err considering degeneracy: " << SCL_num_Z_err_deg << endl;
}