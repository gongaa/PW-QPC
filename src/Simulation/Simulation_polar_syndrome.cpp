#include "Simulation/Simulation.hpp"
using namespace std;
void simulation_polar_syndrome(int N, int Kz, int Kx, int list_size, double px, int num_total, CONSTRUCTION con, double beta, int seed, int interval)
{
    // this is an interpolated method of codeword decoding and direct syndrome decoding (next function)
    // it is omitted from the paper, but here is how it works:
    // generate a random noise n each time, and obtain syndrome s for it (using columns \mathcal{A}_Z^c of the encoding circuit E)
    // put the syndrome at \mathcal{A}_Z^c and fill in zeros at \mathcal{A}_Z, obtain \bar{s}
    // find the preimage of this extended syndrome by \tilde{c}=\bar{s}E, which is a noisy codeword compatible with the syndrome s
    // feed this noisy codeword \tilde{c} to the codeword SCL decoder and solve for a codeword \hat{c}
    // \tilde{c}-\hat{c} is the estimated noise and the correction operation to apply
    cerr << "Simulation Polar syndrome decoding N=" << N << ", Kx=" << Kx << ", Kz=" << Kz << ", l=" << list_size << ", px=" << px 
         << ", #samples=" << num_total << ", seed=" << seed << endl;
    vector<int>  desired_Z(N);          // n + \bar{s}E 
    vector<int>  noisy_codeword_Z(N);   // \bar{s}E
    vector<double> llr_noisy_codeword_Z(N);
    vector<int>  SCL_denoised_codeword_Z(N);
    vector<bool> Z_code_frozen_bits(N); // \mathcal{A}_Z^c
    vector<bool> X_stab_info_bits(N);   // C_X^{perp}
    vector<int>  Z_stab_syndromes(N);   // s

    vector<int>  input(N, 0);
    vector<int>  output(N, 0);

    construct_frozen_bits(con, N, Kz, Kx, Z_code_frozen_bits, X_stab_info_bits, beta);
    print_mixing_factor(Z_code_frozen_bits);
    // N-Kx of X_stab_info_bits are 1. N-Kz of Z_code_frozen_bits are 1.
    Encoder_polar* encoder_Z         = new Encoder_polar(Kz, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(Kz, N, list_size, Z_code_frozen_bits);
    
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(px, 0);
    vector<double> weight = {px/(1.0-px), 0.1, 0.3, 0.5, 0.7, 1.0}; 
    vector<int> noise_X(N, 0), noise_Z(N, 0); 
    vector<int> noise_Z_diff(N, 0);
    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    vector<int>  info_indices; // \mathcal{A}_X \bigcap \mathcal{A}_Z
    for (int i = 0; i < N; i++)
        if (!X_stab_info_bits[i] && !Z_code_frozen_bits[i])
            info_indices.push_back(i);

    cerr << "info indices: ";
    for (auto i : info_indices) cerr << i << " ";
    cerr << endl; 

    int info_size = Kz + Kx - N; 
    assert (info_indices.size() == info_size);
    map<int, vector<int>> equiv_class; 
    map<int, vector<int>> flips;
    vector<int> path_info(info_size);
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
        chn_bsc_q->add_noise(noise_Z.data(), noise_X.data(), 0);
        num_flips = count_weight(noise_Z);
        total_flips += num_flips;

        // obtain syndrome
        copy(noise_Z.data(), noise_Z.data() + N, Z_stab_syndromes.data());
        encoder_Z->light_encode(Z_stab_syndromes.data());
        // from syndrome back to a noisy codeword
        for (int i = 0; i < N; i++) { // obtain \bar{s}
            if(!Z_code_frozen_bits[i]) noisy_codeword_Z[i] = 0; 
            else noisy_codeword_Z[i] = Z_stab_syndromes[i];
        }
        encoder_Z->light_encode(noisy_codeword_Z.data()); // \tilde{c}=\bar{s}E
        // SCL decode
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-px)/px) : log((1-px)/px);
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z); // wt(\tilde{c}-\hat{c})
        xor_vec(N, noise_Z, noisy_codeword_Z, desired_Z); // desired_Z = \tilde{c}+n

        // obtain the partition and number of flips
        equiv_class.clear(); flips.clear();
        SCL_decoder_Z->partition(info_indices, equiv_class, flips, noisy_codeword_Z, SCL_best_class_idx);
        encoder_Z->light_encode(desired_Z.data()); // find pre-image u
        for (int i = 0; i < info_size; i++)
            path_info[i] = desired_Z[info_indices[i]]; // u_{\mathcal{A}_X \bigcap \mathcal{A}_Z}
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
            // print_wt_dist(wt); // sort is included, uncomment this line to see the error coset weight distribution
            for (int i = 0; i < weight.size(); i++) {
                temp_prob = 0.0;
                cal_wt_dist_prob(wt, temp_prob, min_num_flips, weight[i]);
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
            // SCL-C and SCL-E should guess the same if the SCL result is among the closest to the noise
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

        if (turn_idx % interval == (interval-1)) {
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
    cerr << "SCL-E #err: " << SCL_num_Z_err_deg << ", SCL-E Logical X Error Rate: " << (double)SCL_num_Z_err_deg / num_total << endl;
    cerr << "SCL-C #err: " << SCL_num_Z_err_deg_list[0] << ", SCL-C Logical X Error Rate: " << (double)SCL_num_Z_err_deg_list[0] / num_total << endl;
    cerr << "SCL-E #err due to equal: " << equal_flips_err << endl;
    cerr << "SCL-E #err due to smaller: " << SCL_smaller << endl;
}

void simulation_polar_syndrome_direct(int N, int Kz, int Kx, int list_size, double px, int num_total, CONSTRUCTION con, double beta, int seed, int interval)
{
    cerr << "Simulation Polar direct syndrome decoding N=" << N << ", Kx=" << Kx << ", Kz=" << Kz << ", l=" << list_size << ", px=" << px 
         << ", #samples=" << num_total << ", seed=" << seed << endl;
    vector<int>  desired_Z(N, 1);            // generated bit-flip noise, always the same as noise_Z
    vector<int>  noisy_codeword_Z(N, 0);     // the noisy codeword is always the all-zero codeword
    vector<double> llr_noisy_codeword_Z(N, log((1-px)/px)); // the llr of the all-zero codeword
    vector<int>  SCL_denoised_codeword_Z(N); // the estimated noise by the syndrome SCL decoder
    vector<bool> Z_code_frozen_bits(N, 0);   // Z code C_Z (X-type codewords) frozen set \mathcal{A}_Z^c
    vector<bool> X_stab_info_bits(N, 0);     // info bits for C_X^{\perp} (X-type stabilizers)
    vector<int>  Z_stab_syndromes(N, 0);     // measurements results of columns \mathcal{A}_Z^c of the encoding circuit E
    vector<int>  Z_code_frozen_values(N, 0); // Z_stab_syndromes will be the pre-agreed frozen values of the syndrome SCL decoder

    vector<int>  input(N, 0);
    vector<int>  output(N, 0);

    construct_frozen_bits(con, N, Kz, Kx, Z_code_frozen_bits, X_stab_info_bits, beta);
    print_mixing_factor(Z_code_frozen_bits);

    Encoder_polar* encoder_Z         = new Encoder_polar(Kz, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(Kz, N, list_size, Z_code_frozen_bits);

    vector<int>  info_indices; // \mathcal{A}_X \bigcap \mathcal{A}_Z
    for (int i = 0; i < N; i++)
        if (!X_stab_info_bits[i] && !Z_code_frozen_bits[i])
            info_indices.push_back(i);

    cerr << "info indices: "; 
    for (auto i : info_indices) cerr << i << " ";
    cerr << endl; 

    int info_size = Kz + Kx - N;
    assert (info_indices.size() == info_size);
    
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(0, px);
    vector<double> weight = {px/(1.0-px), 0.1, 0.3, 0.5, 0.7, 1.0}; 
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_Z_diff(N, 0);

    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    map<int, vector<int>> equiv_class; 
    map<int, vector<int>> flips;
    vector<int> path_info(info_size); 
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
        chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
        num_flips = count_weight(noise_Z);
        total_flips += num_flips;

        // obtain syndrome
        copy(noise_Z.data(), noise_Z.data() + N, Z_stab_syndromes.data());
        encoder_Z->light_encode(Z_stab_syndromes.data());
        // from syndrome back to a noisy codeword
        for (int i = 0; i < N; i++) 
            if(Z_code_frozen_bits[i]) Z_code_frozen_values[i] = Z_stab_syndromes[i];
        SCL_decoder_Z->set_frozen_values(Z_code_frozen_values);
        // SCL decode
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z);

        // obtain the partition and number of flips
        equiv_class.clear(); flips.clear();
        SCL_decoder_Z->partition(info_indices, equiv_class, flips, noisy_codeword_Z, SCL_best_class_idx);

        xor_vec(N, noise_Z, noisy_codeword_Z, desired_Z); // desired_Z = noise_Z
        encoder_Z->light_encode(desired_Z.data()); // find pre-image u
        for (int i = 0; i < info_size; i++)
            path_info[i] = desired_Z[info_indices[i]]; // u_{\mathcal{A}_X \bigcap \mathcal{A}_Z}
        desired_class_idx = binary2decimal(path_info, info_size);
        if (SCL_best_class_idx != desired_class_idx) { is_SCL_deg_wrong = true; SCL_num_Z_err_deg++; }
        if (is_SCL_deg_wrong) {
            if (num_flips == SCL_num_flips) equal_flips_err++;
            if (SCL_num_flips < num_flips)  SCL_smaller++;
        }
        // cerr << "num_flips: " << num_flips << " , SCL_num_flips: " << SCL_num_flips << endl;

        std::fill(max_class_prob.begin(), max_class_prob.end(), 0.0);
        for (int i = 0; i < max_class_idx.size(); i++) max_class_idx[i].clear();
        min_num_flips = (SCL_num_flips > num_flips) ? num_flips : SCL_num_flips;

        for (auto& x : equiv_class) {
            auto& ec = x.second;
            int k = x.first;
            if (ec.size() == 0) continue;

            auto& wt = flips[k];
            sort(wt.begin(), wt.end());
            // print_wt_dist(wt); // sort is included, uncomment this line to see the error coset weight distribution
            for (int i = 0; i < weight.size(); i++) {
                temp_prob = 0.0;
                cal_wt_dist_prob(wt, temp_prob, min_num_flips, weight[i]);
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
            // SCL-C and SCL-E should guess the same if the SCL result is among the closest to the noise
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

        if (turn_idx % interval == (interval-1)) {
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
    cerr << "SCL-E #err: " << SCL_num_Z_err_deg << ", SCL-E Logical X Error Rate: " << (double)SCL_num_Z_err_deg / num_total << endl;
    cerr << "SCL-C #err: " << SCL_num_Z_err_deg_list[0] << ", SCL-C Logical X Error Rate: " << (double)SCL_num_Z_err_deg_list[0] / num_total << endl;
    cerr << "SCL-E #err due to equal: " << equal_flips_err << endl;
    cerr << "SCL-E #err due to smaller: " << SCL_smaller << endl;
}