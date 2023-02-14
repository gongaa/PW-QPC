#include "Simulation/Simulation.hpp"
using namespace std;
void simulation_polar_codeword(int N, int Kz, int Kx, int list_size, double pz, int num_total,  CONSTRUCTION con, int exact_t, int seed, int print_interval)
{
    cerr << "Simulation Polar codeword decoding N=" << N << ", Kx=" << Kx << ", Kz=" << Kz << ", l=" << list_size << ", pz=" << pz 
         << ", #samples=" << num_total << ", seed=" << seed << endl;
    vector<int> info_bits_Z(Kz);
    vector<int> desired_Z(N, 1);
    vector<int> noisy_codeword_Z(N);
    vector<double> llr_noisy_codeword_Z(N);
    vector<int> SCL_denoised_codeword_Z(N);
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> Z_stab_info_bits(N, 1);
    // Z-type error is detected via the N-Kz X-type stabilizers.
    // Errors that differ by a Z-type stabilizer are deemed the same
    // Z equivalence classes (size Kz+Kx-N) = Z code (size Kz) quotient Z stablizers (size N-Kx)

    construct_frozen_bits(con, N, Kz, Kx, Z_code_frozen_bits, Z_stab_info_bits);
    print_mixing_factor(Z_code_frozen_bits);

    Encoder_polar* encoder_Z = new Encoder_polar(Kz, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(Kz, N, list_size, Z_code_frozen_bits);

    vector<int>  Z_class_info_indices; // expect to have size Kz+Kx-N
    for (int i = 0; i < N; i++)
        if (!Z_stab_info_bits[i] && !Z_code_frozen_bits[i]) // Z stab frozen and Z code info
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

    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(0, pz);
    int floor_Z = N * pz, ceil_Z = floor_Z + 1;
    if (exact_t)
        cerr << "Noise #flips in range [" << floor_Z-exact_t << ", " << ceil_Z+exact_t << "]" << endl;
    vector<double> weight = {pz/(1.0-pz), 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] is set to the real p 
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_Z_diff(N, 0);

    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    int SCL_num_Z_err = 0, SCL_num_Z_err_deg = 0, SCL_num_Z_err_equal_weight = 0;
    int equal_flips_err = 0, SCL_smaller = 0; // out of the SCL_num_Z_err, these are the non-ML errors
    vector<int> SCL_equal_weight_deg_guess_correct(weight.size(), 0);
    vector<int> SCL_num_Z_err_deg_list(weight.size(), 0);
    vector<int> degeneracy_helps(weight.size(), 0);
    vector<int> degeneracy_worse(weight.size(), 0);
    vector<double> max_class_prob(weight.size());
    vector<vector<int>> max_class_idx(weight.size());
    int desired_class_idx, SCL_best_class_idx;
    double temp_prob;
    bool is_SCL_wrong = false, is_SCL_deg_wrong = false, is_SCL_weighted_deg_wrong = false;

    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(seed);
    std::bernoulli_distribution d(0.5);

    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        is_SCL_wrong = false; is_SCL_deg_wrong = false; is_SCL_weighted_deg_wrong = false;
        // generate_random(K, info_bits_Z.data()); // otherwise if all-zero is among the closest to the noise, SCL will always guess correctly
        // generate_random is causing extra randomness other than the BSC channel
        for (int i = 0; i < Kz; i++) info_bits_Z[i] = d(gen);
        encoder_Z->encode(info_bits_Z.data(), desired_Z.data(), 0);
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
        xor_vec(N, desired_Z, noise_Z, noisy_codeword_Z); // noisy_codeword_Z = desired_Z + noise_Z
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-pz)/pz) : log((1-pz)/pz); // 0 -> 1.0; 1 -> -1.0
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        if (!verify(N, SCL_denoised_codeword_Z, desired_Z)) { is_SCL_wrong = true; SCL_num_Z_err++; }
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z);

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
                cal_wt_dist_prob(wt, temp_prob, min_num_flips, weight[i]); // normalize by the min_num_flips, then add probabilities up
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
            // if (max_class_idx[i].size() > 1) cerr << "turn idx " << turn_idx << " SCL-C multiple guess max_class_idx" << endl;
            auto& c = max_class_idx[i];
            int to_compare = c[0];
            // SCL-C and SCL-E should guess the same if the SCL result is among the closest to the noise
            if (c.size() > 1 && (std::find(c.begin(), c.end(), SCL_best_class_idx) != c.end()))
                to_compare = SCL_best_class_idx;
            if (to_compare != desired_class_idx) {
                SCL_num_Z_err_deg_list[i]++;
                is_SCL_weighted_deg_wrong = true;
            }
            if (!is_SCL_weighted_deg_wrong && is_SCL_deg_wrong && (max_class_idx[i].size() == 1)) {
                // cerr << "turn " << turn_idx << " degeneracy helps by not random guessing" << endl;
                degeneracy_helps[i]++;
            }
            if (is_SCL_weighted_deg_wrong && !is_SCL_deg_wrong) degeneracy_worse[i]++;
        }

        if (turn_idx % print_interval == (print_interval-1)) {
            cerr << "SCL #err : " << SCL_num_Z_err << " / " << (turn_idx+1) << endl;
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
    cerr << "SCL #err    : " << SCL_num_Z_err << ", SCL Frame Error Rate: " << (double)SCL_num_Z_err / num_total << endl;
    cerr << "SCL #err considering degeneracy: " << SCL_num_Z_err_deg << ", deg FER: " << (double)SCL_num_Z_err_deg / num_total << endl;
    cerr << "error due to equal " << equal_flips_err << endl;
    cerr << "error due to SCL smaller " << SCL_smaller << endl;
}