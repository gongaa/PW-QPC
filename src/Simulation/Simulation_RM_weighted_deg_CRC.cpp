#include "Simulation/Simulation.hpp"
using namespace std;
int simulation_RM_CSS_weighted_degeneracy(int m, int rx, int rz, int list_size, int p_min, int p_max, int num_total = 100, bool use_crc = true) {
    // int m = 5, rx = 3, rz = 2, list_size = 16384;
    Encoder_RM_CSS* encoder = new Encoder_RM_CSS(m, rx, rz);
    Encoder_RM* encoder_X = new Encoder_RM(m, rz);
    int K_X = encoder_X->get_K();
    // int crc_size = 8;
    // CRC_polynomial<int>* crc_X = new CRC_polynomial<int>(K_X, "8-WCDMA", crc_size);
    int crc_size = 4;
    CRC_polynomial<int>* crc_X = new CRC_polynomial<int>(K_X, "4-ITU", crc_size);
    vector<int> info_bits_X(K_X, 0);
    vector<int> info_bits_X_crc(K_X + crc_size, 0);
    // if two X-type errors differ by a codeword of RM(m, rz)
    // they will have the same syndrome
    Decoder_RM_SCL* SCL_decoder_X = new Decoder_RM_SCL(m, rz, list_size);
    Decoder_RM_SCL* SCL_decoder_Z = new Decoder_RM_SCL(m, rx, list_size);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", rx="<< rx << ", rz=" << rz
         << ", K=" << K << ", N=" << N << endl;
    cerr << "List size=" << list_size << endl;
    // ideally, the noise_X should be decoded to the all-zero codeword
    // but it also can be any codeword of RM(m, rz)
    vector<int> noise_X(N, 0);
    vector<int> noise_X_diff(N, 0);
    vector<int> desired_X(N, 0);
    vector<int> noise_Z(N, 0);
    vector<int> noisy_codeword_X(N, 0);
    vector<double> llr_noisy_codeword_X(N, 0);
    vector<double> llr_noisy_codeword_Z(N, 0);
    vector<int> SCL_denoised_codeword_X(N, 0);
    vector<vector<int>> X_list(list_size, vector<int>(N, 0));
    vector<double> pm_X_list(list_size, 0.0);
    vector<bool> crc_match_X(list_size, false);
    // vector<int> SCL_denoised_codeword_Z(N, 0);
    // vector<vector<int>> Z_list(list_size, vector<int>(N, 0));
    // vector<double> pm_Z_list(list_size, 0.0);
    vector<vector<int>> equiv_class;
    bool is_in_one_class; int largest_class_size; int *largest_class;
    int SCL_num_X_err = 0, SCL_num_Z_err = 0;
    int SCL_num_X_crc_err = 0, SCL_num_Z_crc_err = 0;
    int SCL_num_X_err_deg = 0, SCL_num_X_err_equal_weight = 0;
    int SCL_equal_weight_guess_correct = 0;
    int all_wrong_crc;
    int num_Z_list_err = 0, SCL_num_Z_err_deg = 0, SCL_num_Z_err_deg_list = 0;
    int num_flips, SCL_num_flips;
    bool is_SCL_wrong = false, is_SCL_deg_wrong = false, is_SCL_weighted_deg_wrong = false;
    // vector<double> weight = {1.0, 0.1, 0.2, 0.3, 0.4, 0.5}; // weight[0] will be set to the real p 
    vector<double> weight = {1.0, 0.1, 0.3, 0.5, 0.7, 1.0}; // weight[0] will be set to the real p 
    vector<int> SCL_equal_weight_deg_guess_correct(weight.size(), 0);
    vector<int> SCL_num_X_err_deg_list(weight.size(), 0);
    vector<int> degeneracy_helps(weight.size(), 0);
    vector<int> degeneracy_worse(weight.size(), 0);
    vector<vector<double>> pow(weight.size(), vector<double>(N+1));
    for (int i = 0; i < weight.size(); i++) {
        pow[i][0] = 1;
        for (int j = 1; j < N+1; j++)
            pow[i][j] = pow[i][j-1] * weight[i];
    }
    vector<double> max_class_prob(weight.size());
    vector<vector<int>> max_class_idx(weight.size());
    vector<double> pxs;
    int desired_class_idx, SCL_class_idx;
    // vector<double> class_prob(list_size+1, 0.0);
    int temp_flips;
    double temp_prob;
    int total_flips;
    double pm_best;
    int best_path;
    for (int i = p_min; i <= p_max; i++) pxs.push_back((double)i/100.0);
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, 0.0, 0.0, 41);
    for (double px : pxs) {
        cerr << "px = " << px << endl;
        chn_bsc_q->set_prob(px, 0.0);
        weight[0] = px / (1.0-px);
        for (int j = 1; j < N+1; j++) pow[0][j] = pow[0][j-1] * weight[0];
        SCL_num_X_err = 0; SCL_num_X_err_deg = 0; SCL_num_X_err_equal_weight = 0;
        SCL_num_X_crc_err = 0;
        SCL_equal_weight_guess_correct = 0; all_wrong_crc = 0;
        std::fill(degeneracy_helps.begin(), degeneracy_helps.end(), 0);
        std::fill(degeneracy_worse.begin(), degeneracy_worse.end(), 0);
        std::fill(SCL_num_X_err_deg_list.begin(), SCL_num_X_err_deg_list.end(), 0);
        std::fill(SCL_equal_weight_deg_guess_correct.begin(), SCL_equal_weight_deg_guess_correct.end(), 0);
        total_flips = 0;

        for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
            // cerr << "******* turn " << turn_idx << endl;
            generate_random(K_X, info_bits_X.data());
            if (use_crc) {
                copy(info_bits_X.begin(), info_bits_X.end(), info_bits_X_crc.begin());
                crc_X->build(info_bits_X_crc.data(), info_bits_X_crc.data(), 0);
            }
            encoder_X->encode(info_bits_X.data(), desired_X.data(), 1);
            chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
            num_flips = count_weight(noise_X);
            total_flips += num_flips;
            xor_vec(N, desired_X.data(), noise_X.data(), noisy_codeword_X.data());
            // cerr << "noise adds " << num_flips << " flips" << endl;
            for (int i = 0; i < N; i++) {
                llr_noisy_codeword_X[i] = noisy_codeword_X[i] ? -log((1-px)/px) : log((1-px)/px); // 0 -> 1.0; 1 -> -1.0
            }
            pm_best = SCL_decoder_X->decode(llr_noisy_codeword_X.data(), SCL_denoised_codeword_X.data(), 0);
            SCL_decoder_X->copy_codeword_list(X_list, pm_X_list);
            if (!verify(N, SCL_denoised_codeword_X.data(), desired_X.data())) SCL_num_X_err++;
            if (use_crc) {
                for (int i = 0; i < list_size; i++) {
                    encoder_X->decode(X_list[i].data(), info_bits_X_crc.data());
                    crc_match_X[i] = crc_X->check(info_bits_X_crc.data(), 0);
                }
                // cerr << "There are " << std::count(crc_match_X.begin(), crc_match_X.end(), true) << " cases that survived CRC" << endl;
                if (std::none_of(crc_match_X.begin(), crc_match_X.end(), [](bool v) { return v; })) {
                    all_wrong_crc++;
                    SCL_num_X_crc_err++;
                    // for (auto& x : X_list) assert(!verify(N, x.data(), desired_X.data()));
                    continue;
                }
                pm_best = std::numeric_limits<double>::max(); best_path = -1;
                for (int i = 0; i < list_size; i++) {
                    if (crc_match_X[i] && (pm_X_list[i] < pm_best)) {
                        best_path = i;
                        pm_best = pm_X_list[i];
                    }
                }
                assert (best_path >= 0);
                std::copy(X_list[best_path].begin(), X_list[best_path].end(), SCL_denoised_codeword_X.begin());
            }
            for (auto& dx : X_list) xor_vec(N, dx.data(), desired_X.data(), dx.data());
            is_SCL_wrong = false;
            if (!verify(N, SCL_denoised_codeword_X.data(), desired_X.data())) {
                SCL_num_X_crc_err++;
                is_SCL_wrong = true;
            }
            xor_vec(N, SCL_denoised_codeword_X.data(), desired_X.data(), SCL_denoised_codeword_X.data());
            SCL_num_flips = count_flip(N, SCL_denoised_codeword_X.data(), noise_X.data());
            if (is_SCL_wrong && (num_flips == SCL_num_flips)) 
                for (int i = 0; i < X_list.size(); i++) {
                    if (count_weight(X_list[i]) == 0 && abs(pm_X_list[i] - pm_best) < std::numeric_limits<double>::epsilon()) 
                        cerr << "Found a case where num_flips=SCL_num_flips=" << num_flips <<  " but path metrics are different" 
                        << ", desired pm=" << pm_X_list[i] << ", SCL best pm=" << pm_best << endl;

                    break;
                }
            if (is_SCL_wrong && (SCL_num_flips != num_flips)) 
                // not fail due to equal weight (even ML cannot decide)
                SCL_num_X_err_equal_weight++;
            if (!is_SCL_wrong && (SCL_num_flips == num_flips)) 
                    SCL_equal_weight_guess_correct++;
            // cerr << "plain SCL is " << (is_SCL_wrong ? "wrong" : "correct") << endl;
            is_SCL_deg_wrong = false;
            if (!encoder->is_X_stabilizer(SCL_denoised_codeword_X.data())) {
                SCL_num_X_err_deg++;
                is_SCL_deg_wrong = true;
            }
            equiv_class.clear();
            equiv_class.push_back({0});
            for (int i = 1; i < list_size; i++) {
                is_in_one_class = false;
                for (auto& ec : equiv_class) {
                    xor_vec(N, X_list[ec[0]].data(), X_list[i].data(), noise_X_diff.data());
                    if (encoder->is_X_stabilizer(noise_X_diff.data())) {
                        ec.push_back(i);
                        is_in_one_class = true;
                        break;
                    }
                }
                if (!is_in_one_class)
                    equiv_class.push_back({i});
            }
            // cerr << "there are " << equiv_class.size() << " equiv classes" << endl;
            std::fill(max_class_prob.begin(), max_class_prob.end(), 0.0);
            for (int i = 0; i < max_class_idx.size(); i++) max_class_idx[i].clear();
            desired_class_idx = -1; // desired class may not be in the list
            SCL_class_idx = -1;     // this is guaranteed in the list
            for (int k = 0; k < equiv_class.size(); k++) {
                auto& ec = equiv_class[k];
                // kick out this class from the competition if no member passed CRC
                // if (use_crc && std::none_of(ec.begin(), ec.end(), [crc_match_X](int i) { return crc_match_X[i]; })) continue;
                for (int i = 0; i < weight.size(); i++) {
                    temp_prob = 0.0;
                    for (int ec_idx : ec) {
                        if (!crc_match_X[ec_idx]) continue;
                        temp_flips = count_flip(N, X_list[ec_idx].data(), noise_X.data());
                        temp_prob += pow[i][temp_flips];
                    }
                    // ******************* This line is causing trouble
                    // if(i == 0) class_prob[k] = temp_prob;
                    // **************************************************
                    if (temp_prob > max_class_prob[i]) {
                        max_class_prob[i] = temp_prob;
                        max_class_idx[i].clear();
                        max_class_idx[i].push_back(k);
                    } else if (temp_prob > (max_class_prob[i] - std::numeric_limits<double>::epsilon())) {
                        max_class_idx[i].push_back(k);
                    }
                }

                // if (encoder->is_X_stabilizer(X_list[ec[0]].data())) desired_class_idx = k;
                // xor_vec(N, X_list[ec[0]].data(), SCL_denoised_codeword_X.data(), noise_X_diff.data());
                // if (encoder->is_X_stabilizer(noise_X_diff.data())) SCL_class_idx = k;
            }
            for (int i = 0; i < weight.size(); i++) {
                is_SCL_weighted_deg_wrong = true;
                // random guessing
                if ((max_class_idx[i].size() > 1) && 
                    (std::find(max_class_idx[i].begin(), max_class_idx[i].end(), SCL_class_idx) != max_class_idx[i].end())) {
                    // the same guess as SCL best
                    is_SCL_weighted_deg_wrong = is_SCL_deg_wrong;
                } else if (encoder->is_X_stabilizer(X_list[equiv_class[max_class_idx[i][0]][0]].data())) {
                    is_SCL_weighted_deg_wrong = false;
                }
                if (is_SCL_weighted_deg_wrong) SCL_num_X_err_deg_list[i]++;
                if (num_flips != SCL_num_flips) {
                    if (!is_SCL_weighted_deg_wrong && is_SCL_wrong) degeneracy_helps[i]++;
                    if (!is_SCL_wrong && is_SCL_weighted_deg_wrong) degeneracy_worse[i]++;
                }
                /*
                if ((i == 0) && (is_SCL_wrong ^ is_SCL_weighted_deg_wrong)) {
                    cerr << "num_flips: " << num_flips << ", SCL_num_flips: " << SCL_num_flips << endl;
                    cerr << "there are " << max_class_idx[0].size() << " classes: ";
                    for (int k : max_class_idx[0]) cerr << k << " ";
                    cerr << ", prob " << max_class_prob[0] << endl;
                    if (desired_class_idx >= 0) cerr << "desired equiv class " << desired_class_idx << ", size " << equiv_class[desired_class_idx].size() << ", prob " << class_prob[desired_class_idx] << endl;
                    cerr << "SCL best equiv class " << SCL_class_idx << ", size " << equiv_class[SCL_class_idx].size() << ", prob " << class_prob[SCL_class_idx] << endl;
                    if (!is_SCL_weighted_deg_wrong && is_SCL_wrong) {
                        if (num_flips == SCL_num_flips)
                            cerr << "weighted degeneracy guess correctly but SCL guess wrongly" << endl;
                        else {
                            cerr << "degeneracy helps" << endl;
                        }
                    }
                    if (!is_SCL_wrong && is_SCL_weighted_deg_wrong) {
                        if (num_flips == SCL_num_flips)
                            cerr << "SCL guess correctly but weighted degeneracy guess wrongly" << endl;
                        else 
                            cerr << "degeneracy makes it worse" << endl;
                    }
                }
                */
                if (!is_SCL_weighted_deg_wrong && (SCL_num_flips == num_flips)) SCL_equal_weight_deg_guess_correct[i]++;
            }
        }
        cerr << "average #flips: " << (double)total_flips / num_total << endl;
        cerr << "In " << num_total << " samples, there are " << all_wrong_crc << " cases where no denoised codeword passes CRC." << endl;
        cerr << "SCL fail not due to equal weight: " << SCL_num_X_err_equal_weight << endl;
        cerr << "SCL #err    : " << SCL_num_X_err << ", SCL Frame Error Rate: " << (double)SCL_num_X_err / num_total << endl;
        cerr << "SCL-CRC #err: " << SCL_num_X_crc_err << ", exists one passes CRC: " << (SCL_num_X_crc_err - all_wrong_crc) << endl;
        cerr << "SCL #err considering degeneracy: " << SCL_num_X_err_deg << endl;
        cerr << "SCL #err weighted degeneracy: ";
        for (int i : SCL_num_X_err_deg_list) cerr << i << " ";
        cerr << ", with weights: ";
        for (double i : weight) cerr << i << " ";
        cerr << endl << "SCL guess correct under equal weight: " << SCL_equal_weight_guess_correct;
        cerr << endl << "SCL with degeneracy guess correct under equal weight: ";
        for (int i : SCL_equal_weight_deg_guess_correct) cerr << i << " ";
        cerr << endl << "degeneracy_helps: ";
        for (int i : degeneracy_helps) cerr << i << " ";
        cerr << endl << "degeneracy worse: ";
        for (int i : degeneracy_worse) cerr << i << " ";
        cerr << endl;
    }
    return 0;
}