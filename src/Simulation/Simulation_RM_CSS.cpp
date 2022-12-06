#include "Simulation/Simulation.hpp"
#define PRINT_FLIPS
int simulation_RM_CSS(int m, int rx, int rz, int list_size) {
    // int m = 5, rx = 3, rz = 2, list_size = 16384;
    Encoder_RM_CSS* encoder = new Encoder_RM_CSS(m, rx, rz);
    // if two X-type errors differ by a codeword of RM(m, rz)
    // they will have the same syndrome
    Decoder_RM_SCL* SCL_decoder_X = new Decoder_RM_SCL(m, rz, list_size);
    Decoder_RM_SCL* SCL_decoder_Z = new Decoder_RM_SCL(m, rx, list_size);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", rx="<< rx << ", rz=" << rz
         << ", K=" << K << ", N=" << N << endl;
    cerr << "List size=" << list_size << endl;
    double px = 0.1, pz = 0;
    cerr << "px=" << px << ", pz=" << pz << endl;
    Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, px, pz, 41);
    // ideally, the noise_X should be decoded to the all-zero codeword
    // but it also can be any codeword of RM(m, rz)
    vector<int> noise_X(N, 0);
    vector<int> noise_X_diff(N, 0);
    vector<int> desired_X(N, 0);
    vector<int> noise_Z(N, 0);
    vector<double> llr_noisy_codeword_X(N, 0);
    vector<double> llr_noisy_codeword_Z(N, 0);
    vector<int> SCL_denoised_codeword_X(N, 0);
    vector<vector<int>> X_list(list_size, vector<int>(N, 0));
    vector<double> pm_X_list(list_size, 0.0);
    // vector<int> SCL_denoised_codeword_Z(N, 0);
    // vector<vector<int>> Z_list(list_size, vector<int>(N, 0));
    // vector<double> pm_Z_list(list_size, 0.0);
    vector<vector<int>> equiv_class;
    bool is_in_one_class; int largest_class_size; int *largest_class;
    int num_total = 100, SCL_num_X_err = 0, SCL_num_Z_err = 0;
    int SCL_num_X_err_deg = 0, SCL_num_X_err_deg_list = 0;
    int num_Z_list_err = 0, SCL_num_Z_err_deg = 0, SCL_num_Z_err_deg_list = 0;
    double pm_best;
    for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
        chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0);
        for (int i = 0; i < N; i++) {
            llr_noisy_codeword_X[i] = noise_X[i] ? -log((1-px)/px) : log((1-px)/px); // 0 -> 1.0; 1 -> -1.0
            // llr_noisy_codeword_Z[i] = noise_Z[i] ? -log((1-pz)/px) : log((1-px)/px); // 0 -> 1.0; 1 -> -1.0
        }
        pm_best = SCL_decoder_X->decode(llr_noisy_codeword_X.data(), SCL_denoised_codeword_X.data(), 0);
        SCL_decoder_X->copy_codeword_list(X_list, pm_X_list);
        // SCL_decoder_Z->decode(llr_noisy_codeword_Z.data(), SCL_denoised_codeword_Z.data(), 0);
        // SCL_decoder_Z->copy_codeword_list(Z_list, pm_Z_list);
// #ifdef VANILLA_LIST_DECODING
        cerr << "******* turn " << turn_idx << endl;
        if (!std::all_of(SCL_denoised_codeword_X.begin(), SCL_denoised_codeword_X.end(), [](int i) { return i==0; })) {
            cerr << "SCL was wrong, the codeword it gave has pm=" << pm_best << " and content" << endl;
            for (int k : SCL_denoised_codeword_X) cerr << k;
            cerr << ", pm=" << pm_best << ", #flips=" << 
            count_flip(N, SCL_denoised_codeword_X, noise_X) << endl;
            SCL_num_X_err++;
            for (int i = 0; i < list_size; i++) {
                if (encoder->is_X_stabilizer(X_list[i].data())) {
                    cerr << "wrong but idx=" << i << " differs by only a stabilizer" << endl;
                    for (int k : X_list[i]) cerr << k;
                    cerr << ", pm=" << pm_X_list[i] << ", #flips=" <<
                    count_flip(N, X_list[i], noise_X) << endl;
                }
            }
        }
        if (!encoder->is_X_stabilizer(SCL_denoised_codeword_X.data()))
            SCL_num_X_err_deg++;
// #elif defined DEGENERACY_LIST_DECODING
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
#ifdef CHECK
        for (int i = 0; i < list_size; i++)
            assert (Encoder_RM::is_codeword(X_list[i].data(), m, rz));
        int num_equiv_class = equiv_class.size();
        for (int i = 0; i < num_equiv_class; i++) {
            for (int j = 0; j < equiv_class[i].size(); j++) 
                cerr << equiv_class[i][j] << " ";
            cerr << endl;
        }

        for (int i = 0; i < num_equiv_class; i++) {
            for (int j = i + 1; j < num_equiv_class; j++) {
                xor_vec(N, X_list[equiv_class[i][0]].data(), X_list[equiv_class[j][0]].data(), noise_X_diff.data());
                assert (!encoder->is_X_stabilizer(noise_X_diff.data()));
                assert (encoder->is_logical_X(noise_X_diff.data()));
            }
            int ec_size = equiv_class[i].size();
            if (ec_size > 1) {
                for (int j = 0; j < ec_size; j++) {
                    for (int k = j + 1; k < ec_size; k++) {
                        xor_vec(N, X_list[equiv_class[i][j]].data(), X_list[equiv_class[i][k]].data(), noise_X_diff.data());
                        // assert (Encoder_RM::is_codeword(noise_X_diff.data(), m, rz));
                        assert (!encoder->is_logical_X(noise_X_diff.data()));
                    }
                }
            }
        }
#endif // CHECK
        largest_class_size = 0;
        cerr << "classes sizes are: ";
        for (auto& ec : equiv_class) {
            if (encoder->is_X_stabilizer(X_list[ec[0]].data())) cerr << "stab_";
            cerr << ec.size() << " ";
            if (ec.size() > largest_class_size) {
                largest_class_size = ec.size();
                largest_class = X_list[ec[0]].data();
            }
        }
        cerr << endl;
        // cerr << "larget equiv class size " << largest_class_size << endl;
        if (!encoder->is_X_stabilizer(largest_class))
            SCL_num_X_err_deg_list++;
        // else
            // cerr << "turn " << turn_idx << " 00....0 is correctly decoded considering degeneracy, largest equiv class size " << largest_class_size << endl;
// endif
    }
    cerr << "SCL_num_err: " << SCL_num_X_err << endl;
    // cerr << "SCL_num_X_list_0_err: " << SCL_num_X_list_0_err << endl;
    cerr << "SCL_num_err_deg: " << SCL_num_X_err_deg << endl;
    cerr << "SCL_num_err_deg_list: " << SCL_num_X_err_deg_list << endl;
    cerr << "SCL Frame Error Rate: " << (double)SCL_num_X_err / num_total << endl;
    return 0;
}

int simulation_RM_degeneracy(int m, int rx, int rz, double px_dummy, double pz)
{
    Encoder_RM_CSS* encoder = new Encoder_RM_CSS(m, rx, rz);
    // if two X-type errors differ by a codeword of RM(m, rz)
    // they will have the same syndrome
    int K = encoder->get_K(), N = encoder->get_N(), N_half = N >> 1;
    cerr << "For m=" << m << ", rx="<< rx << ", rz=" << rz
         << ", K=" << K << ", N=" << N << endl;
    // ideally, the noise_X should be decoded to the all-zero codeword
    // but it also can be any codeword of RM(m, rz)
    vector<int> noise_X(N, 0), noise_Z(N, 0);
    vector<int> noise_X_diff(N, 0);
    vector<int> desired_X(N, 0);
    // vector<int> noise_Z(N, 0);
    vector<vector<int>> equiv_class; // equiv_class[0] is the stabilizers
    vector<vector<int>> codewords;
    generate_all_equiv_classes(m, rx, rz, codewords, equiv_class);
    int num_codewords = codewords.size(), num_equiv_class = equiv_class.size(), ec_size = equiv_class[0].size();
    // assert (encoder->is_X_stabilizer(codewords[equiv_class[0][0]].data())); 
    int num_total = 2000, ML_err = 0, ML_deg_err = 0;
    int ML_err_no_rep = 0, ML_deg_err_no_rep = 0;
    int degeneracy_better = 0;
    int num_flips = N, min_num_flips = N;
    vector<int> min_flip_indices;
    vector<int> min_err_prob_classes;
    vector<vector<int>> min_flip_array;
    double ec_err_prob, max_ec_err_prob;
    vector<int> flips(N + 1, 0);
    // flips[i] stores how many codewords in this equiv class differ
    // by exactly i positions to the noise_X.
    vector<double> flips_err_prob(N + 1, 0);
    bool exists_stab = false;
    vector<double> pxs;
    for (int i = 9; i < 40; i++) pxs.push_back((double)(i + 1) / 100.0);
    for (double px : pxs) {
        cerr << "px = " << px << endl;
        double log_err = N * log(1-px), log_err_diff = log(px) - log(1-px);
        flips_err_prob[0] = exp(log_err);
        for (int i = 1; i < N + 1; i++) {
            log_err += log_err_diff;
            flips_err_prob[i] = exp(log_err);
            // c++ double exponent min is -1022, may underflow
        }
        Channel_BSC_q* chn_bsc_q = new Channel_BSC_q(N, px, pz, 41);
        ML_err = 0; ML_deg_err = 0; ML_err_no_rep = 0; ML_deg_err_no_rep = 0;
        degeneracy_better = 0;
        for (int turn_idx = 0; turn_idx < num_total; turn_idx++) {
            chn_bsc_q->add_noise(noise_X.data(), noise_Z.data(), 0); // only use noise_X for now
            // ML_err is just take the codeword with smallest flip, 
            // see whether it is in the all-zero's equiv class (stabilizers)
            min_num_flips = N; min_flip_indices.clear();
            for (int i = 0; i < num_codewords; i++) {
                num_flips = count_flip(N, noise_X, codewords[i]);
                if (num_flips < min_num_flips) {
                    min_flip_indices.clear();
                    min_flip_indices.push_back(i);
                    min_num_flips = num_flips;
                } else if (num_flips == min_num_flips) {
                    min_flip_indices.push_back(i);
                }
            }
            // cerr << "there are " << min_flip_indices.size() << " codewords that have min flips" << endl;
            // when px=0.1, a lot of times min_flip_indices.size() > 1
            exists_stab = false;
            for (int i : min_flip_indices) {
                if (encoder->is_X_stabilizer(codewords[i].data())) {
                    exists_stab = true;
                    break;
                }
            }
            if (!exists_stab) ML_err++;
            if (!exists_stab || min_flip_indices.size() > 1) ML_err_no_rep++; 

            // ML_deg_err is to add up the error-probability for each equiv class
            // and choose the largest one
            max_ec_err_prob = 0;
            for (int i = 0; i < num_equiv_class; i++) {
                for (int k = 0; k < N + 1; k++) flips[k] = 0;
                for (int idx : equiv_class[i]) {
                    num_flips = count_flip(N, noise_X, codewords[idx]);
                    flips[num_flips]++;
                }
                ec_err_prob = 0.0; 
                for (int k = 0; k < N + 1; k++)
                    ec_err_prob += flips[k] * flips_err_prob[k];
                // cerr << "ec " << i << " has error prob " << ec_err_prob << endl;
                // for (int k = 0; k < N + 1; k++)
                //     cerr << flips[k] << " ";
                // cerr << endl;
                if (ec_err_prob > max_ec_err_prob) {
                    min_err_prob_classes.clear();
                    min_err_prob_classes.push_back(i);
                    max_ec_err_prob = ec_err_prob;
                    #ifdef PRINT_FLIPS
                    min_flip_array.clear();
                    min_flip_array.push_back(flips);
                    #endif // PRINT_FLIPS
                } else if (ec_err_prob > max_ec_err_prob - std::numeric_limits<double>::epsilon()) {
                    min_err_prob_classes.push_back(i);
                    #ifdef PRINT_FLIPS
                    min_flip_array.push_back(flips);
                    #endif
                }
            }
            // cerr << "ec " << min_err_prob_classes[0] << " is the class that has the min error prob" << endl;
            // cerr << "there are " << min_err_prob_classes.size() << " equiv classes that have the same min error prob" << endl;
            if (std::none_of(min_err_prob_classes.begin(), min_err_prob_classes.end(), [](int i) { return i == 0; })) {
                ML_deg_err++; ML_deg_err_no_rep++;
            } else if (min_err_prob_classes.size() > 1) {
                ML_deg_err_no_rep++;
                if (min_err_prob_classes.size() < min_flip_indices.size()) degeneracy_better++;
                if (min_err_prob_classes.size() != min_flip_indices.size()) cerr << "two sizes disagree" << endl;
                #ifdef PRINT_FLIPS
                if (min_num_flips == 4 && min_flip_indices.size() > 2) {
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
                /*
                cerr <<"***********************" << endl;
                cerr << "min_flip_indices size=" << min_flip_indices.size() <<
                ",\t min_err_prob_classes size=" << min_err_prob_classes.size() << endl;
                cerr << "min_num_flips=" << min_num_flips << endl;
                for (auto a : min_flip_array) {
                    for (int f : a) cerr << f << " ";
                    cerr << endl;
                }
                cerr <<"***********************" << endl;
                */
                #endif
            }
            #ifdef PRINT_FLIPS
              else {
                if (min_num_flips >= 4)
                    cerr << "can correctly decode when #flips >= 4, this shouldn't happen???" << endl;
            }
            #endif
        }
        cerr << "ML_err: " << ML_err << ". ML_deg_err: " << ML_deg_err << endl;
        cerr << "ML_err_no_rep: " << ML_err_no_rep << ". ML_deg_err_no_rep: " << ML_deg_err_no_rep << endl;
        cerr << "degeneracy_better: " << degeneracy_better << endl;;
        cerr << "ML FER: " << (double)ML_err / num_total << endl;
        cerr << "ML degeneracy FER: " << (double)ML_deg_err / num_total << endl;
        cerr << "ML no repetition FER: " << (double)ML_err_no_rep / num_total << endl;
        cerr << "ML degeneracy no repetition FER: " << (double)ML_deg_err_no_rep / num_total << endl;
    }
    return 0;
}
