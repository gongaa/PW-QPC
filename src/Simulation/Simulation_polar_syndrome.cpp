#include "Simulation/Simulation.hpp"
using namespace std;
void simulation_polar_syndrome(int N, int K, int list_size, double pz, int num_total, CONSTRUCTION con, int exact_t, int seed, int interval)
{
    int Kx = K, Kz= K;
    cerr << "Simulation Polar syndrome decoding N=" << N << ", Kx=" << Kx << ", Kz=" << Kz << ", l=" << list_size << ", pz=" << pz 
         << ", #samples=" << num_total << ", seed=" << seed << endl;
    vector<int>  desired_Z(N, 1);
    vector<int>  noisy_codeword_Z(N);
    vector<double> llr_noisy_codeword_Z(N);
    vector<int>  SCL_denoised_codeword_Z(N);
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> Z_stab_info_bits(N, 0);
    vector<bool> X_code_frozen_bits(N, 0);
    vector<bool> X_stab_info_bits(N, 0);
    // Z type noisy codeword has the syndrome at the N-Kz frozen positions.

    int num_stab = N-Kz;
    vector<int>  X_stab_syndromes(num_stab, 0);
    vector<vector<int>> X_stab(num_stab, vector<int>(N,0));

    vector<int>  input(N, 0);
    vector<int>  output(N, 0);

    construct_frozen_bits(con, N, Kz, Kx, Z_code_frozen_bits, Z_stab_info_bits);

    // N-Kx of Z_stab_info_bits are 1. N-Kx of X_code_frozen_bits are 1.
    // N-Kz of Z_code_frozen_bits are 1. N-Kz of X_stab_info_bits are 1.
    for (int i = 0; i < N; i++) X_code_frozen_bits[N-1-i] = Z_stab_info_bits[i];
    for (int i = 0; i < N; i++) X_stab_info_bits[N-1-i] = Z_code_frozen_bits[i];
            
    Encoder_polar* encoder_Z         = new Encoder_polar(K, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(K, N, list_size, Z_code_frozen_bits);

    int j = 0;
    for (int i = 0; i < N; i++) {
        if (X_stab_info_bits[i]) {
            input[i] = 1;
            output = input;
            encoder_Z->light_encode(output.data());
            X_stab[j++] = output;
            input[i] = 0;
        }
    }
    
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(0, pz);
    int floor_Z = N * pz, ceil_Z = floor_Z + 1;
    if (exact_t)
        cerr << "Noise #flips in range [" << floor_Z-exact_t << ", " << ceil_Z+exact_t << "]" << endl;
    vector<double> weight = {pz/(1.0-pz), 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] will be set to the real p 
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_Z_diff(N, 0);
    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    vector<int>  Z_class_info_indices; // expect to have size Kz+Kx-N
    for (int i = 0; i < N; i++)
        if (!Z_stab_info_bits[i] && !Z_code_frozen_bits[i])
            Z_class_info_indices.push_back(i);

    cerr << "Z class info indices: ";
    for (auto i : Z_class_info_indices) cerr << i << " ";
    cerr << endl; 

    int info_size = Kz + Kx - N; 
    assert (Z_class_info_indices.size() == info_size);
    map<int, vector<int>> equiv_class; // expect there to be (at most) Kz+Kx-N classes
    map<int, vector<int>> flips;
    vector<int> path_info(info_size); // store the binary representation of desired_class_idx
    bool is_in_one_class;

    int equal_flips_err = 0, SCL_smaller = 0;
    int SCL_num_Z_err_deg = 0, SCL_num_Z_err_equal_weight = 0;
    int SCL_equal_weight_guess_correct = 0;
    vector<int> SCL_equal_weight_deg_guess_correct(weight.size(), 0);
    vector<int> SCL_num_Z_err_deg_list(weight.size(), 0);
    vector<int> degeneracy_helps(weight.size(), 0);
    vector<int> degeneracy_worse(weight.size(), 0);
    vector<double> max_class_prob(weight.size());
    vector<vector<int>> max_class_idx(weight.size());
    int desired_class_idx, SCL_best_class_idx;
    double temp_prob;
    bool is_SCL_deg_wrong = false, is_SCL_weighted_deg_wrong = false;

    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        // reset flags
        is_SCL_deg_wrong = false; is_SCL_weighted_deg_wrong = false;
        // generate noise
        if (!exact_t) {
            chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
            num_flips = count_weight(noise_Z);
        } else {
            do {
                chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
                num_flips = count_weight(noise_Z);
            } while (num_flips < (floor_Z - exact_t) || num_flips > (ceil_Z + exact_t));
        }
        total_flips += num_flips;
        // obtain syndrome
        for (int i = 0; i < num_stab; i++) {
            X_stab_syndromes[i] = dot_product(N, X_stab[i], noise_Z);
        }
        // from syndrome back to a noisy codeword
        j = 0;
        for (int i = 0; i < N; i++) {
            if (X_code_frozen_bits[N-i-1]) {
                noisy_codeword_Z[i] = X_stab_syndromes[j++];
            } else noisy_codeword_Z[i] = 0; // TODO: make it random
        }
        encoder_Z->transpose_encode(noisy_codeword_Z.data());
        // SCL decode
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-pz)/pz) : log((1-pz)/pz); // 0 -> 1.0; 1 -> -1.0
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z);
        xor_vec(N, noise_Z, noisy_codeword_Z, desired_Z);

        // obtain the partition and number of flips
        equiv_class.clear(); flips.clear();
        SCL_decoder_Z->partition(Z_class_info_indices, equiv_class, flips, noisy_codeword_Z, SCL_best_class_idx);
        encoder_Z->light_encode(desired_Z.data()); // find pre-image
        for (int i = 0; i < info_size; i++)
            path_info[i] = desired_Z[Z_class_info_indices[i]];
        desired_class_idx = binary2decimal(path_info, info_size);
        if (SCL_best_class_idx != desired_class_idx) { is_SCL_deg_wrong = true; SCL_num_Z_err_deg++; }
        if (is_SCL_deg_wrong) {
            if (num_flips == SCL_num_flips) equal_flips_err++;
            if (SCL_num_flips < num_flips)  SCL_smaller++;
        }
        // determine max class with respect to weighted degeneracy
        std::fill(max_class_prob.begin(), max_class_prob.end(), 0.0);
        for (int i = 0; i < max_class_idx.size(); i++) max_class_idx[i].clear();
        min_num_flips = (SCL_num_flips > num_flips) ? num_flips : SCL_num_flips;

        for (auto& x : equiv_class) {
            auto& ec = x.second;
            int k = x.first;
            if (ec.size() == 0) continue;
            auto& wt = flips[k];
            sort(wt.begin(), wt.end());
            // print_wt_dist(wt); // sort is included
            for (int i = 0; i < weight.size(); i++) {
                temp_prob = 0.0;
                cal_wt_dist_prob(wt, temp_prob, min_num_flips, weight[i]);
                if (temp_prob > max_class_prob[i]) {
                    max_class_prob[i] = temp_prob;
                    max_class_idx[i].clear();
                    max_class_idx[i].push_back(k);
                } else if (temp_prob > (max_class_prob[i] - std::numeric_limits<double>::epsilon())) {
                    max_class_idx[i].push_back(k);
                }
            }
        }

        for (int i = 0; i < weight.size(); i++) {
            is_SCL_weighted_deg_wrong = false;
            if (max_class_idx[i].size() > 1) {
                cerr << "turn idx " << turn_idx << " SCL-C multiple guess" << endl;
            }
            if (max_class_idx[i][0] != desired_class_idx) {
                SCL_num_Z_err_deg_list[i]++;
                is_SCL_weighted_deg_wrong = true;
            }
            if (!is_SCL_weighted_deg_wrong && is_SCL_deg_wrong && (max_class_idx[i].size() == 1)) {
                // cerr << "turn " << turn_idx << " degeneracy helps by not random guessing" << endl;
                degeneracy_helps[i]++;
            }
            if (is_SCL_weighted_deg_wrong && !is_SCL_deg_wrong) degeneracy_worse[i]++;
        }

        if (turn_idx % interval == (interval-1)) {
            cerr << "SCL #err deg : " << SCL_num_Z_err_deg << " / " << (turn_idx+1) << endl;
            cerr << "SCL #err weighted degeneracy : ";
            for (int i : SCL_num_Z_err_deg_list) cerr << i << " ";
            cerr << endl << "degeneracy helps : ";
            for (int i : degeneracy_helps) cerr << i << " ";
            cerr << endl << "degeneracy worse : ";
            for (int i : degeneracy_worse) cerr << i << " ";
            cerr << endl;
        }


    }
    cerr << "average #flips: " << (double)total_flips / num_total << endl;
    cerr << "SCL #err considering degeneracy: " << SCL_num_Z_err_deg << endl;
    cerr << "error due to equal " << equal_flips_err << endl;
    cerr << "error due to SCL smaller " << SCL_smaller << endl;
}

void simulation_polar_syndrome_direct(int N, int K, int list_size, double pz, int num_total, CONSTRUCTION con, int exact_t, int seed, int interval)
{
    int Kx = K, Kz = K;
    cerr << "Simulation Polar direct syndrome decoding N=" << N << ", Kx=" << Kx << ", Kz=" << Kz << ", l=" << list_size << ", pz=" << pz 
         << ", #samples=" << num_total << ", seed=" << seed << endl;
    vector<int>  desired_Z(N, 1);
    vector<int>  noisy_codeword_Z(N, 0);
    vector<double> llr_noisy_codeword_Z(N, log((1-pz)/pz));
    vector<int>  SCL_denoised_codeword_Z(N);
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> Z_stab_info_bits(N, 0);
    vector<int>  Z_code_frozen_values(N, 0);
    vector<bool> X_stab_info_bits(N, 1);

    int num_stab = N-K;
    vector<int>  X_stab_syndromes(num_stab, 0);
    vector<vector<int>> X_stab(num_stab, vector<int>(N,0));

    vector<int>  input(N, 0);
    vector<int>  output(N, 0);

    construct_frozen_bits(con, N, Kz, Kx, Z_code_frozen_bits, Z_stab_info_bits);

    for (int i = 0; i < N; i++) X_stab_info_bits[N-1-i] = Z_code_frozen_bits[i];

    Encoder_polar* encoder_Z         = new Encoder_polar(K, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(K, N, list_size, Z_code_frozen_bits);

    int j = 0;
    for (int i = 0; i < N; i++) {
        if (X_stab_info_bits[i]) {
            input[i] = 1;
            output = input;
            encoder_Z->light_encode(output.data());
            X_stab[j++] = output;
            input[i] = 0;
        }
    } 

    vector<int>  Z_class_info_indices; // expect to have size Kz+Kx-N
    for (int i = 0; i < N; i++)
        if (!Z_stab_info_bits[i] && !Z_code_frozen_bits[i])
            Z_class_info_indices.push_back(i);

    cerr << "Z class info indices: ";
    for (auto i : Z_class_info_indices) cerr << i << " ";
    cerr << endl; 

    int info_size = Kz + Kx - N;
    assert (Z_class_info_indices.size() == info_size);
    
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(0, pz);
    int floor_Z = N * pz, ceil_Z = floor_Z + 1;
    if (exact_t)
        cerr << "Noise #flips in range [" << floor_Z-exact_t << ", " << ceil_Z+exact_t << "]" << endl;
    vector<double> weight = {pz/(1.0-pz), 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] will be set to the real p 
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_Z_diff(N, 0);

    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    map<int, vector<int>> equiv_class; // expect there to be 2*K-N classes
    map<int, vector<int>> flips;
    vector<int> path_info(info_size); // store the binary representation of desired_class_idx
    bool is_in_one_class;

    int equal_flips_err = 0, SCL_smaller = 0;
    int SCL_num_Z_err_deg = 0, SCL_num_Z_err_equal_weight = 0;
    int SCL_equal_weight_guess_correct = 0;
    vector<int> SCL_equal_weight_deg_guess_correct(weight.size(), 0);
    vector<int> SCL_num_Z_err_deg_list(weight.size(), 0);
    vector<int> degeneracy_helps(weight.size(), 0);
    vector<int> degeneracy_worse(weight.size(), 0);
    vector<double> max_class_prob(weight.size());
    vector<vector<int>> max_class_idx(weight.size());
    int desired_class_idx, SCL_best_class_idx;
    double temp_prob;
    bool is_SCL_deg_wrong = false, is_SCL_weighted_deg_wrong = false;
    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        // reset flags
        is_SCL_deg_wrong = false; is_SCL_weighted_deg_wrong = false;
        // generate noise
        if (!exact_t) {
            chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
            num_flips = count_weight(noise_Z);
        } else {
            do {
                chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
                num_flips = count_weight(noise_Z);
            } while (num_flips < (floor_Z - exact_t) || num_flips > (ceil_Z + exact_t));
        }
        total_flips += num_flips;
        // obtain syndrome
        for (int i = 0; i < num_stab; i++) {
            X_stab_syndromes[i] = dot_product(N, X_stab[i], noise_Z);
        }
        // put reversed syndrome into frozen values
        j = num_stab - 1;
        for (int i = 0; i < N; i++) 
            if (Z_code_frozen_bits[i]) Z_code_frozen_values[i] = X_stab_syndromes[j--];
        SCL_decoder_Z->set_frozen_values(Z_code_frozen_values);
        // SCL decode
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z);
        // reverse noise
        bit_reversal(noise_Z);

        // obtain the partition and number of flips
        equiv_class.clear(); flips.clear();
        SCL_decoder_Z->partition(Z_class_info_indices, equiv_class, flips, noisy_codeword_Z, SCL_best_class_idx);

        xor_vec(N, noise_Z, noisy_codeword_Z, desired_Z);
        encoder_Z->light_encode(desired_Z.data()); // find pre-image
        for (int i = 0; i < info_size; i++)
            path_info[i] = desired_Z[Z_class_info_indices[i]];
        desired_class_idx = binary2decimal(path_info, info_size);
        if (SCL_best_class_idx != desired_class_idx) { is_SCL_deg_wrong = true; SCL_num_Z_err_deg++; }
        if (is_SCL_deg_wrong) {
            if (num_flips == SCL_num_flips) equal_flips_err++;
            if (SCL_num_flips < num_flips)  SCL_smaller++;
        }

        // cerr << "num_flips: " << num_flips << " , SCL_num_flips: " << SCL_num_flips << endl;

        // determine max class with respect to weighted degeneracy
        std::fill(max_class_prob.begin(), max_class_prob.end(), 0.0);
        for (int i = 0; i < max_class_idx.size(); i++) max_class_idx[i].clear();
        min_num_flips = (SCL_num_flips > num_flips) ? num_flips : SCL_num_flips;

        for (auto& x : equiv_class) {
            auto& ec = x.second;
            int k = x.first;
            if (ec.size() == 0) continue;

            auto& wt = flips[k];
            sort(wt.begin(), wt.end());
            // print_wt_dist(wt); // sort is included
            for (int i = 0; i < weight.size(); i++) {
                temp_prob = 0.0;
                cal_wt_dist_prob(wt, temp_prob, min_num_flips, weight[i]);
                if (temp_prob > max_class_prob[i]) {
                    max_class_prob[i] = temp_prob;
                    max_class_idx[i].clear();
                    max_class_idx[i].push_back(k);
                } else if (temp_prob > (max_class_prob[i] - std::numeric_limits<double>::epsilon())) {
                    max_class_idx[i].push_back(k);
                }
            }
        }
       
        for (int i = 0; i < weight.size(); i++) {
            // if (max_class_idx[i].size() > 1) cerr << "multiple max class" << endl;
            is_SCL_weighted_deg_wrong = false;
            if (max_class_idx[i][0] != desired_class_idx) {
                SCL_num_Z_err_deg_list[i]++;
                is_SCL_weighted_deg_wrong = true;
            }
            if (!is_SCL_weighted_deg_wrong && is_SCL_deg_wrong && (max_class_idx[i].size() == 1)) {
                // cerr << "turn " << turn_idx << " degeneracy helps by not random guessing" << endl;
                degeneracy_helps[i]++;
            }
            if (is_SCL_weighted_deg_wrong && !is_SCL_deg_wrong) degeneracy_worse[i]++;
        }

        if (turn_idx % interval == (interval-1)) {
            cerr << "SCL #err deg : " << SCL_num_Z_err_deg << " / " << (turn_idx+1) << endl;
            cerr << "SCL #err weighted degeneracy : ";
            for (int i : SCL_num_Z_err_deg_list) cerr << i << " ";
            cerr << endl << "degeneracy helps : ";
            for (int i : degeneracy_helps) cerr << i << " ";
            cerr << endl << "degeneracy worse : ";
            for (int i : degeneracy_worse) cerr << i << " ";
            cerr << endl;
        }

    }
    cerr << "average #flips: " << (double)total_flips / num_total << endl;
    cerr << "SCL #err considering degeneracy: " << SCL_num_Z_err_deg << endl;
    cerr << "error due to equal " << equal_flips_err << endl;
    cerr << "error due to SCL smaller " << SCL_smaller << endl;
}