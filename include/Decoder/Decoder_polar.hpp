#ifndef DECODER_POLAR_SCL_HPP_
#define DECODER_POLAR_SCL_HPP_
// based on the paper "LLR-Based Successive Canncellation List Decoding of Polar Codes"
#include <vector>
#include <set>
#include "Decoder.hpp"
#include "Tree.hpp"
#include "Util/Util.hpp"
using namespace std;

class Contents_SCL
{
public:
    vector<double> l;     // Log-Likelihood Ratio array
    vector<int> s;       // partial sum array u

    bool is_frozen_bit;
    int max_depth_llrs;

    explicit Contents_SCL(int size) : l(size, 0), s(size, 0), is_frozen_bit(false) {}
    virtual ~Contents_SCL() {}
};

// temporarily single kernel F={{1,0},{1,1}}
class Decoder_polar_SCL : Decoder
{
protected:
    // const double metric_init; // init value of the metrics in the trees
    const int L;             // maximum path number
    std::set<int> active_paths;
    int best_path;

    vector<bool> frozen_bits;
    vector<Tree_metric<Contents_SCL>*> polar_trees;
    vector<vector<Node<Contents_SCL>*>> leaves_array;   

    // vector<vector<function<double(const vector<double> &LLRs, const vector<int> &bits)>>> lambdas;

public:
    Decoder_polar_SCL(const int& K, const int& N, const int& L, const vector<bool>& frozen_bits);
    virtual ~Decoder_polar_SCL();
    virtual int decode(const double *Y_N, int *V_K, const size_t frame_id);
    void decode_SC(const double *Y_N, int *V_K, const size_t frame_id);
    void get_llr_for_frozen_bits(double *Y_N);
    void copy_codeword_list(vector<vector<int>>& c_list, vector<double>& pm_list);
    void partition(vector<int>& info_indices, vector<vector<int>>& par, vector<vector<int>>& flips, vector<int>& noisy_codeword, int& best_path_class_idx);

protected:
    void _load(const double *Y_N);
    void _decode(const size_t frame_id);
    void _store(int *V_K) const;
    void _decode_SC(Node<Contents_SCL>* node_cur);

private:
    void recursive_compute_llr(Node<Contents_SCL>* node_cur, int depth);
    void recursive_propagate_sums(const Node<Contents_SCL>* node_cur);
    void duplicate_path(int path, int leaf_index, vector<vector<Node<Contents_SCL>*>> leaves_array);

    void recursive_duplicate_tree_llr(Node<Contents_SCL>* node_a, Node<Contents_SCL>* node_b);
    void recursive_duplicate_tree_sums(Node<Contents_SCL>* node_a, Node<Contents_SCL>* node_b, Node<Contents_SCL>* node_caller);

protected:
    virtual void select_best_path(const size_t frame_id);
    
    void recursive_allocate_nodes_contents(Node<Contents_SCL>* node_curr, const int vector_size, int& max_depth_llrs);
    void recursive_initialize_frozen_bits(const Node<Contents_SCL>* node_curr, const std::vector<bool>& frozen_bits);
    void recursive_deallocate_nodes_contents(Node<Contents_SCL>* node_curr);

};

#endif /* DECODER_POLAR_SCL_HPP_ */
