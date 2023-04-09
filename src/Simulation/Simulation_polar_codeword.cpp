#include "Simulation/Simulation.hpp"
using namespace std;
void simulation_polar_codeword(int N, int Kz, int Kx, int list_size, double px, int num_total, CONSTRUCTION con, double beta, int seed, int print_interval)
{
    cerr << "Simulation Polar codeword decoding N=" << N << ", Kx=" << Kx << ", Kz=" << Kz << ", l=" << list_size << ", px=" << px 
         << ", #samples=" << num_total << ", seed=" << seed << endl;
    // the suffix Z means in the Z basis, the type is specified in the comment
    // I choose this convention so that X and Z won't be mixed up in correlated X and Z decoding in the future
    vector<int> info_bits_Z(Kz); // a random infomation vector used to generate a random X-type codeword
    vector<int> desired_Z(N, 1); // the random X-type codeword generated in each iteration for the Steane EC
    vector<int> noisy_codeword_Z(N); // what individual Z-type measurements on ancilla bits (Stean EC) give you
    vector<double> llr_noisy_codeword_Z(N); // log-likelihood-ratio of the received noisy codeword y
    vector<int> SCL_denoised_codeword_Z(N); // answer of the codeword SCL decoder to y
    vector<bool> Z_code_frozen_bits(N);     // frozen bit mask for C_Z (X-type codewords), i.e., the frozen set \mathcal{A}_Z^c
    vector<bool> X_stab_info_bits(N);       // info bit mask for C_X^{\perp} (X-type stabilizers)
    // X-type error is detected via the N-Kz Z-type stabilizers (columns \mathcal{A}_Z^c of the encoding matrix E).
    // Errors that differ by a X-type stabilizer are dgenerate (deemed the same)
    // X equivalence classes (size Kz+Kx-N) = X-type codewords (C_Z, size Kz) quotient X-type stablizers (C_X^{\perp}, size N-Kx)

    construct_frozen_bits(con, N, Kz, Kx, Z_code_frozen_bits, X_stab_info_bits, beta);
    print_mixing_factor(Z_code_frozen_bits);

    Encoder_polar* encoder_Z = new Encoder_polar(Kz, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(Kz, N, list_size, Z_code_frozen_bits);

    vector<int> info_indices; // \mathcal{A}_X \bigcap \mathcal{A}_Z, expect to have size Kz+Kx-N, carry quantum logical state
    for (int i = 0; i < N; i++)
        if (!X_stab_info_bits[i] && !Z_code_frozen_bits[i]) // C_Z / C_X^{\perp}
            info_indices.push_back(i);

    cerr << "info indices: "; // you can use this to calculate the distance of the quantum code
    for (auto i : info_indices) cerr << i << " "; // row i of E (a logical X operator generator) has weight 2^{wt(bin(i))}
    cerr << endl; // linear code, take the minimum of those weight to obtain X distance
    // logical Z operators are formed by row i of E^T, which is the bit-reversal of row N-1-i of E, weight 2^{log(N)-wt(bin(i))}
    // take the minimum of X and Z distance

    int info_size = Kz + Kx - N; 
    assert (info_indices.size() == info_size);
    map<int, vector<int>> equiv_class; // expect there to be (at most) 2^{Kz+Kx-N} classes when looking at the list
    map<int, vector<int>> flips;       // weight of noise patterns 
    vector<int> path_info(info_size);  // store the binary representation of desired_class_idx (to which class desired_Z belongs)
    bool is_in_one_class; 

    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(px, 0);
    // weight used in calculating probability of noise, 
    // the SCL-C error is counted using the true weight px/(1-px)
    vector<double> weight = {px/(1.0-px), 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] is set to the true weight
    vector<int> noise_X(N, 0); // Z-type noise, dummy, needed because of the channel add noise interface
    vector<int> noise_Z(N, 0); // X-type noise
    vector<int> noise_Z_diff(N, 0);

    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    //  SCL Frame ER     , SCL-E logical X ER   , SCL-C logical X ER
    int SCL_num_Z_err = 0, SCL_num_Z_err_deg = 0, SCL_num_Z_err_equal_weight = 0;
    int equal_flips_err = 0, SCL_smaller = 0; // SCL frame errors where an MLD will also be wrong
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
        /* In order to simulate the Steane error correction, 
           at each step, a random codeword instead of the all-zero codeword should be generated. 
           This is because SCL has a bias for the decision of the all-zero codeword. 
           At each decision step, SCL first inserts decision bit 0 then decision bit 1. 
           Finally, the front-most codeword is selected. If we use the all-zero codeword, 
           then SCL is naturally biased towards the correct answer and the error rate is lower than usual.
        */
        for (int i = 0; i < Kz; i++) info_bits_Z[i] = d(gen);       // generate random info vector
        encoder_Z->encode(info_bits_Z.data(), desired_Z.data(), 0); // encode to get random codeword
        // generate noise_Z (X-type noise)
        chn_bsc_q->add_noise(noise_Z.data(), noise_X.data(), 0);
        num_flips = count_weight(noise_Z);
        total_flips += num_flips;
        xor_vec(N, desired_Z, noise_Z, noisy_codeword_Z); // noisy_codeword_Z = desired_Z + noise_Z
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-px)/px) : log((1-px)/px); // 0 -> 1.0; 1 -> -1.0
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        if (!verify(N, SCL_denoised_codeword_Z, desired_Z)) { is_SCL_wrong = true; SCL_num_Z_err++; }
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z);

        // for each codewords in the list, determine the equiv class it belongs to by looking at it's u values at the info indices
        // calculate its distance to the noisy codeword and put this distance into that equiv class
        equiv_class.clear(); flips.clear();
        SCL_decoder_Z->partition(info_indices, equiv_class, flips, noisy_codeword_Z, SCL_best_class_idx);
        encoder_Z->light_encode(desired_Z.data()); // find pre-image u
        for (int i = 0; i < info_size; i++)
            path_info[i] = desired_Z[info_indices[i]];  // look at u_{\mathcal{A}_X \bigcap \mathcal{A}_Z}
        desired_class_idx = binary2decimal(path_info, info_size); // to determine which class desired_Z belongs to
        if (SCL_best_class_idx != desired_class_idx) { is_SCL_deg_wrong = true; SCL_num_Z_err_deg++; }
        if (is_SCL_wrong) {
            if (num_flips == SCL_num_flips) equal_flips_err++;
            if (SCL_num_flips < num_flips)  SCL_smaller++;
        }

        // obtain the equiv class with max error probability (SCL-C answer)
        std::fill(max_class_prob.begin(), max_class_prob.end(), 0.0);
        for (int i = 0; i < max_class_idx.size(); i++) max_class_idx[i].clear();
        min_num_flips = (SCL_num_flips > num_flips) ? num_flips : SCL_num_flips; // used for normalization

        for (auto& x : equiv_class) {
            auto& ec = x.second;
            int k = x.first;
            if (ec.size() == 0) continue;
            auto& wt = flips[k];
            sort(wt.begin(), wt.end());
            // print_wt_dist(wt); // sort is included, uncomment this line to see the error coset weight distribution
            for (int i = 0; i < weight.size(); i++) {
                temp_prob = 0.0;
                cal_wt_dist_prob(wt, temp_prob, min_num_flips, weight[i]); // normalize by the min_num_flips, then add probabilities up
                // Normalization is necessary! 
                if (temp_prob > max_class_prob[i]) {
                    max_class_prob[i] = temp_prob;
                    max_class_idx[i].clear();
                    max_class_idx[i].push_back(k);
                } else if (temp_prob >= (max_class_prob[i] - std::numeric_limits<double>::epsilon())) {
                    max_class_idx[i].push_back(k);
                }
            }
        }

        for (int i = 0; i < weight.size(); i++) {
            is_SCL_weighted_deg_wrong = false;
            auto& c = max_class_idx[i];
            int to_compare = c[0];
            // SCL-C and SCL-E should guess the same if the SCL result is among the closest to the noise, in order to avoid fluctuations.
            if (c.size() > 1 && (std::find(c.begin(), c.end(), SCL_best_class_idx) != c.end()))
                to_compare = SCL_best_class_idx;
            if (to_compare != desired_class_idx) {
                SCL_num_Z_err_deg_list[i]++;
                is_SCL_weighted_deg_wrong = true;
            }
            if (!is_SCL_weighted_deg_wrong && is_SCL_deg_wrong && (max_class_idx[i].size() == 1)) 
                degeneracy_helps[i]++;
            if (is_SCL_weighted_deg_wrong && !is_SCL_deg_wrong && (max_class_idx[i].size() == 1)) 
                degeneracy_worse[i]++;
        }

        if (turn_idx % print_interval == (print_interval-1)) {
            cerr << "SCL   #err: " << SCL_num_Z_err << " / " << (turn_idx+1) << endl;
            cerr << "SCL-E #err: " << SCL_num_Z_err_deg << " / " << (turn_idx+1) << endl;
            cerr << "SCL-C #err: ";
            for (int i : SCL_num_Z_err_deg_list) cerr << i << " ";
            cerr << endl << "degeneracy helps: ";
            for (int i : degeneracy_helps) cerr << i << " ";
            cerr << endl << "degeneracy worse: ";
            for (int i : degeneracy_worse) cerr << i << " ";
            cerr << endl;
        }
    }
    cerr << "average #flips: " << (double)total_flips / num_total << endl;
    cerr << "SCL   #err: " << SCL_num_Z_err << ", SCL Frame Error Rate: " << (double)SCL_num_Z_err / num_total << endl;
    cerr << "SCL-E #err: " << SCL_num_Z_err_deg << ", SCL-E Logical X Error Rate: " << (double)SCL_num_Z_err_deg / num_total << endl;
    cerr << "SCL-C #err: " << SCL_num_Z_err_deg_list[0] << ", SCL-C Logical X Error Rate: " << (double)SCL_num_Z_err_deg_list[0] / num_total << endl;
    cerr << "SCL   #err due to equal: " << equal_flips_err << endl;
    cerr << "SCL   #err due to smaller: " << SCL_smaller << endl;
}