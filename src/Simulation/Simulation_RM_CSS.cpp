#include "Simulation/Simulation.hpp"

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
            count_flip(N, SCL_denoised_codeword_X.data(), noise_X.data()) << endl;
            SCL_num_X_err++;
            for (int i = 0; i < list_size; i++) {
                if (encoder->is_X_stabilizer(X_list[i].data())) {
                    cerr << "wrong but idx=" << i << " differs by only a stabilizer" << endl;
                    for (int k : X_list[i]) cerr << k;
                    cerr << ", pm=" << pm_X_list[i] << ", #flips=" <<
                    count_flip(N, X_list[i].data(), noise_X.data()) << endl;
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


// check SCL decoder behaves symmetrically around every codeword under BSC
void test_RM_SCL_symmetry() {
    int m = 7, r = 3, list_size = 128;
    Encoder* encoder = new Encoder_RM(m, r);
    // Decoder_RM_SC* SC_decoder = new Decoder_RM_SC(m ,r, 1);
    Decoder_RM_SCL* SCL_decoder = new Decoder_RM_SCL(m ,r, list_size);
    int K = encoder->get_K(), N = encoder->get_N();
    cerr << "For m=" << m << ", r="<< r << ", K=" << K << ", N=" << N << endl;
    cerr << "List size=" << list_size << endl;

    double p = 0.01;
    cerr << "p=" << p << endl;
    Channel_BSC* chn_bsc = new Channel_BSC(N, p, 42);
    vector<int> info_bits(K, 1);
    vector<int> codeword(N, 0);
    vector<int> all_zero(N, 0);
    vector<int> noise(N, 0); // all-zero + noise
    vector<int> noisy_codeword(N, 0);
    vector<int> diff_to_noise(N, 0);
    vector<double> llr_noisy_codeword_1(N, 0);
    vector<double> llr_noisy_codeword_2(N, 0);
    vector<int> SCL_denoised_codeword_1(N, 0);
    vector<int> SCL_denoised_codeword_2(N, 0);
    vector<vector<int>> X_list_1(list_size, vector<int>(N, 0));
    vector<double> pm_X_list_1(list_size, 0.0);
    vector<vector<int>> X_list_2(list_size, vector<int>(N, 0));
    vector<double> pm_X_list_2(list_size, 0.0);
    int num_total = 3, SCL_num_err = 0;
    int SCL_num_flips_1 = 0, SCL_num_flips_2 = 0, bsc_flips = 0;
    for (int i = 0; i < num_total; i++) {
        generate_random(K, info_bits.data());
        encoder->encode(info_bits.data(), codeword.data(), 1);
        bsc_flips = chn_bsc->add_noise(all_zero.data(), noise.data(), 0);
        cerr << "BSC #flips=" << bsc_flips << endl;
        xor_vec(N, noise.data(), codeword.data(), noisy_codeword.data());
        for (int i = 0; i < N; i++) {
            llr_noisy_codeword_1[i] = noise[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            llr_noisy_codeword_2[i] = noisy_codeword[i] ? -log((1-p)/p) : log((1-p)/p); // 0 -> 1.0; 1 -> -1.0
            // llr_noisy_codeword[i] = noisy_codeword[i] ? -1.0 : 1.0; // 0 -> 1.0; 1 -> -1.0
        }
        SCL_decoder->decode(llr_noisy_codeword_1.data(), SCL_denoised_codeword_1.data(), 0);
        SCL_decoder->copy_codeword_list(X_list_1, pm_X_list_1);
        SCL_decoder->decode(llr_noisy_codeword_2.data(), SCL_denoised_codeword_2.data(), 0);
        SCL_decoder->copy_codeword_list(X_list_2, pm_X_list_2);
        for (int i = 0; i < list_size; i++) {
            SCL_num_flips_1 = count_flip(N, noise.data(), X_list_1[i].data());
            cerr << "1: pm=" << pm_X_list_1[i] << "\t";
            cerr << "#flips=" << SCL_num_flips_1 << "\t";
            for (int j : X_list_1[i]) cerr << j;
            cerr << endl;
            SCL_num_flips_2 = count_flip(N, noisy_codeword.data(), X_list_2[i].data());
            cerr << "2: pm=" << pm_X_list_2[i] << "\t";
            cerr << "#flips=" << SCL_num_flips_2 << "\t";
            for (int j = 0; j < N; j++) X_list_2[i][j] ^= codeword[j];
            for (int j : X_list_2[i]) cerr << j;
            cerr << endl;
        }
    }
}