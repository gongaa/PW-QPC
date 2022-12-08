#include "Simulation/Simulation.hpp"
using namespace std;
void simulation_polar_syndrome(int N, int K, int list_size, double pz, int num_total, CONSTRUCTION con, int exact_t, int seed)
{
    cerr << "Simulation Polar syndrome N=" << N << ", K=" << K << ", l=" << list_size << ", pz=" << pz << ", #samples=" << num_total << endl;
    vector<int>  desired_Z(N, 1);
    vector<int>  noisy_codeword_Z(N);
    vector<double> llr_noisy_codeword_Z(N);
    vector<int>  SCL_denoised_codeword_Z(N);
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> X_code_frozen_bits(N, 0);
    vector<bool> X_stab_frozen_bits(N, 1);
    vector<int>  X_stab_syndromes(N-K, 0);

    vector<vector<int>> X_stab(N-K, vector<int>(N,0));
    vector<int>  input(N, 0);
    vector<int>  output(N, 0);

    construct_frozen_bits(con, N, K, Z_code_frozen_bits);

    for (int i = 0; i < N; i++) X_code_frozen_bits[i] = Z_code_frozen_bits[N-1-i];
    for (int i = 0; i < N; i++) 
        if (Z_code_frozen_bits[i] == 0 && Z_code_frozen_bits[N-1-i] == 1) 
            X_stab_frozen_bits[i] = 0;

    Encoder_polar* encoder_Z         = new Encoder_polar(K, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(K, N, list_size, Z_code_frozen_bits);
    Encoder_polar* X_stab_encoder    = new Encoder_polar(N-K, N, X_stab_frozen_bits);

    int j = 0;
    for (int i = 0; i < N; i++) {
        if (!X_stab_frozen_bits[i]) {
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
    vector<double> weight = {1.0, 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] will be set to the real p 
    weight[0] = pz / (1.0-pz);
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_Z_diff(N, 0);
    vector<vector<int>> Z_list(list_size, vector<int>(N, 0));
    vector<double> pm_Z_list(list_size, 0.0);

    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    vector<vector<int>> equiv_class; // expect there to be 2*K-N classes
    bool is_in_one_class; int largest_class_size; int* largest_class;

    int equal_flips_err = 0, SCL_smaller = 0;
    int SCL_num_Z_err_deg = 0, SCL_num_Z_err_equal_weight = 0;
    int SCL_equal_weight_guess_correct = 0;
    vector<int> SCL_equal_weight_deg_guess_correct(weight.size(), 0);
    vector<int> SCL_num_Z_err_deg_list(weight.size(), 0);
    vector<int> degeneracy_helps(weight.size(), 0);
    vector<int> degeneracy_worse(weight.size(), 0);
    vector<double> max_class_prob(weight.size());
    vector<vector<int>> max_class_idx(weight.size());
    int desired_class_idx, SCL_class_idx, largest_class_idx;
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
        for (int i = 0; i < N-K; i++) {
            X_stab_syndromes[i] = dot_product(N, X_stab[i], noise_Z);
        }
        // from syndrome back to a noisy codeword
        j = 0;
        for (int i = 0; i < N; i++) {
            if (X_code_frozen_bits[i]) {
                noisy_codeword_Z[i] = X_stab_syndromes[j++];
            } else noisy_codeword_Z[i] = 0; // TODO: make it random
        }
        encoder_Z->transpose_encode(noisy_codeword_Z.data());
        // SCL decode
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-pz)/pz) : log((1-pz)/pz); // 0 -> 1.0; 1 -> -1.0
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        SCL_decoder_Z->copy_codeword_list(Z_list, pm_Z_list);
        // post-processing
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z);
        xor_vec(N, noise_Z, noisy_codeword_Z, desired_Z);
        xor_vec(N, SCL_denoised_codeword_Z, desired_Z, SCL_denoised_codeword_Z);
        if (!X_stab_encoder->is_codeword(SCL_denoised_codeword_Z.data())) { is_SCL_deg_wrong = true; SCL_num_Z_err_deg++; }
        if (is_SCL_deg_wrong) {
            if (num_flips == SCL_num_flips) equal_flips_err++;
            if (SCL_num_flips < num_flips)  SCL_smaller++;
        }
        for (auto& dz : Z_list) xor_vec(N, dz, desired_Z, dz);
        // partition into equivalence classes
        equiv_class.clear();
        equiv_class.push_back({0});
        for (int i = 1; i < list_size; i++) {
            is_in_one_class = false;
            for (auto& ec : equiv_class) {
                xor_vec(N, Z_list[ec[0]], Z_list[i], noise_Z_diff);
                if (X_stab_encoder->is_codeword(noise_Z_diff.data())) {
                    ec.push_back(i);
                    is_in_one_class = true;
                    break;
                }
            }
            if (!is_in_one_class)
                equiv_class.push_back({i});
        }
        // cerr << "num_flips: " << num_flips << " , SCL_num_flips: " << SCL_num_flips << endl;
        // cerr << "there are " << equiv_class.size() << " equiv classes" << endl;
        // cerr << "Weight Distribution with added noise:" << endl;
        // determine max class with respect to weighted degeneracy
        std::fill(max_class_prob.begin(), max_class_prob.end(), 0.0);
        for (int i = 0; i < max_class_idx.size(); i++) max_class_idx[i].clear();
        desired_class_idx = -1; // desired class may not be in the list
        SCL_class_idx = -1;     // this is guaranteed in the list
        min_num_flips = (SCL_num_flips > num_flips) ? num_flips : SCL_num_flips;
        largest_class_size = 0; largest_class_idx = 0;

        for (int k = 0; k < equiv_class.size(); k++) {
            auto& ec = equiv_class[k];
            vector<int> wt(ec.size());
            for (int l = 0; l < wt.size(); l++) {
                wt[l] = count_flip(N, Z_list[ec[l]], noise_Z); 
            }
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

            if (X_stab_encoder->is_codeword(Z_list[ec[0]].data())) desired_class_idx = k;

            if (ec.size() > largest_class_size) {
                largest_class_size = ec.size();
                largest_class = ec.data();
                largest_class_idx = k;
            }
        }
        // cerr << "largest class size: " << largest_class_size << ". Weight distribution ";
        // vector<int> wt(largest_class_size);
        // for (int i = 0; i < largest_class_size; i++) 
        //     wt[i] = count_weight(Z_list[largest_class[i]]);
        // print_wt_dist(wt);
        
        // if (largest_class_idx != desired_class_idx) {
        //     if (desired_class_idx >= 0) {
        //         auto& ec = equiv_class[desired_class_idx];
        //         cerr << "stabilizer class size: " << ec.size() << ". Weight distribution ";
        //         vector<int> wt(ec.size());
        //         for (int i = 0; i < ec.size(); i++) 
        //             wt[i] = count_weight(Z_list[ec[i]]);
        //         print_wt_dist(wt);
        //     }
        // }
        // else cerr << "stabilizer class is the largest class" << endl;
       

        // cerr << "max class prob (normalized) : ";
        // for (auto i : max_class_prob) cerr << i << " ";
        // cerr << endl;
        for (int i = 0; i < weight.size(); i++) {
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

        if (turn_idx % 100 == 99) {
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

void simulation_polar_syndrome_fast(int N, int K, int list_size, double pz, int num_total, CONSTRUCTION con, int exact_t, int seed)
{
    cerr << "Simulation Polar syndrome FASTER N=" << N << ", K=" << K << ", l=" << list_size << ", pz=" << pz << ", #samples=" << num_total << endl;
    vector<int>  desired_Z(N, 1);
    vector<int>  noisy_codeword_Z(N);
    vector<double> llr_noisy_codeword_Z(N);
    vector<int>  SCL_denoised_codeword_Z(N);
    vector<bool> Z_code_frozen_bits(N, 0);
    vector<bool> X_code_frozen_bits(N, 0);
    vector<bool> X_stab_frozen_bits(N, 1);
    vector<int>  X_stab_info_indices; // expect to have size 2*K-N
    vector<int>  X_stab_syndromes(N-K, 0);

    vector<vector<int>> X_stab(N-K, vector<int>(N,0));
    vector<int>  input(N, 0);
    vector<int>  output(N, 0);

    construct_frozen_bits(con, N, K, Z_code_frozen_bits);

    for (int i = 0; i < N; i++) X_code_frozen_bits[i] = Z_code_frozen_bits[N-1-i];
    for (int i = 0; i < N; i++) 
        if (Z_code_frozen_bits[i] == 0 && Z_code_frozen_bits[N-1-i] == 1) 
            X_stab_frozen_bits[i] = 0;

    Encoder_polar* encoder_Z         = new Encoder_polar(K, N, Z_code_frozen_bits);
    Decoder_polar_SCL* SCL_decoder_Z = new Decoder_polar_SCL(K, N, list_size, Z_code_frozen_bits);
    Encoder_polar* X_stab_encoder    = new Encoder_polar(N-K, N, X_stab_frozen_bits);

    int j = 0;
    for (int i = 0; i < N; i++) {
        if (!X_stab_frozen_bits[i]) {
            input[i] = 1;
            output = input;
            encoder_Z->light_encode(output.data());
            X_stab[j++] = output;
            input[i] = 0;
        }
    } 

    for (int i = 0; i < N; i++)
        if (X_stab_frozen_bits[i] && !Z_code_frozen_bits[i])
            X_stab_info_indices.push_back(i);

    int info_size = 2*K - N;
    int num_equiv_classes = 1 << info_size;
    assert (X_stab_info_indices.size() == info_size);
    
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0, 0, seed);
    chn_bsc_q->set_prob(0, pz);
    int floor_Z = N * pz, ceil_Z = floor_Z + 1;
    if (exact_t)
        cerr << "Noise #flips in range [" << floor_Z-exact_t << ", " << ceil_Z+exact_t << "]" << endl;
    vector<double> weight = {1.0, 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] will be set to the real p 
    weight[0] = pz / (1.0-pz);
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_Z_diff(N, 0);

    int total_flips = 0, num_flips = 0, SCL_num_flips = 0, min_num_flips = 0;
    double pm_best;

    vector<vector<int>> equiv_class(num_equiv_classes, vector<int>()); // expect there to be 2*K-N classes
    vector<vector<int>> flips(num_equiv_classes, vector<int>());
    vector<int> path_info(info_size);
    bool is_in_one_class; int largest_class_size; int* largest_class;

    int equal_flips_err = 0, SCL_smaller = 0;
    int SCL_num_Z_err_deg = 0, SCL_num_Z_err_equal_weight = 0;
    int SCL_equal_weight_guess_correct = 0;
    vector<int> SCL_equal_weight_deg_guess_correct(weight.size(), 0);
    vector<int> SCL_num_Z_err_deg_list(weight.size(), 0);
    vector<int> degeneracy_helps(weight.size(), 0);
    vector<int> degeneracy_worse(weight.size(), 0);
    vector<double> max_class_prob(weight.size());
    vector<vector<int>> max_class_idx(weight.size());
    int desired_class_idx, largest_class_idx, SCL_best_class_idx;
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
        for (int i = 0; i < N-K; i++) {
            X_stab_syndromes[i] = dot_product(N, X_stab[i], noise_Z);
        }
        // from syndrome back to a noisy codeword
        j = 0;
        for (int i = 0; i < N; i++) {
            if (X_code_frozen_bits[i]) {
                noisy_codeword_Z[i] = X_stab_syndromes[j++];
            } else noisy_codeword_Z[i] = 0; // TODO: make it random
        }
        encoder_Z->transpose_encode(noisy_codeword_Z.data());
        // SCL decode
        for (int i = 0; i < N; i++) llr_noisy_codeword_Z[i] = noisy_codeword_Z[i] ? -log((1-pz)/pz) : log((1-pz)/pz); // 0 -> 1.0; 1 -> -1.0
        pm_best = SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        SCL_num_flips = count_flip(N, SCL_denoised_codeword_Z, noisy_codeword_Z);

        // obtain the partition and number of flips
        for (int i = 0; i < num_equiv_classes; i++) { equiv_class[i].clear(); flips[i].clear(); }
        SCL_decoder_Z->partition(X_stab_info_indices, equiv_class, flips, noisy_codeword_Z, SCL_best_class_idx);

        xor_vec(N, noise_Z, noisy_codeword_Z, desired_Z);
        encoder_Z->light_encode(desired_Z.data()); // find pre-image
        for (int i = 0; i < info_size; i++)
            path_info[i] = desired_Z[X_stab_info_indices[i]];
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
        largest_class_size = 0; largest_class_idx = 0;

        for (int k = 0; k < num_equiv_classes; k++) {
            auto& ec = equiv_class[k];
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

            if (ec.size() > largest_class_size) {
                largest_class_size = ec.size();
                largest_class = ec.data();
                largest_class_idx = k;
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

        if (turn_idx % 100 == 99) {
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