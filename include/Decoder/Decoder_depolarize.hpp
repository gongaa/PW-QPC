#ifndef DECODER_DEPOLARIZE_HPP_
#define DECODER_DEPOLARIZE_HPP_
#include <vector>
#include <set>
#include "Decoder.hpp"
#include "Tree.hpp"
#include "Util/Util.hpp"
using namespace std;

class Contents_depolarize
{
public:
    vector<vector<double>> l;     // probability (TODO: make it LLR) array (p_I, p_X, p_Z, p_Y)
    vector<vector<int>> s;       // partial sum array (u_x, u_z)

    bool is_frozen_bit;
    int max_depth_llrs;

    explicit Contents_depolarize(int size) : l(size, vector<double>(4,0)), s(size, vector<int>(2,0)), is_frozen_bit(false) {}
    virtual ~Contents_depolarize() {}
};

class Decoder_depolarize : Decoder
{
protected:
    const int L;             // maximum path number
    std::set<int> active_paths;
    int best_path;
    vector<vector<int>> vs;  // for precoded polar codes

    vector<bool> frozen_bits;
    vector<int> frozen_values; // support arbitrary frozen values (not all-zero)
    vector<Tree_metric<Contents_depolarize>*> polar_trees;
    vector<vector<Node<Contents_depolarize>*>> leaves_array;   

public:
    Decoder_depolarize(const int& K, const int& N, const int& L, const vector<bool>& frozen_bits);
    virtual ~Decoder_depolarize();
    virtual int decode(const vector<vector<double>>& Y_N, vector<vector<int>>& V_K);
    void decode_SC(const vector<vector<double>>& Y_N, vector<vector<int>>& V_K);
    void set_frozen_values(const vector<int>& fv);

protected:
    void _load(const vector<vector<double>>& Y_N);
    void _decode();
    void _store(vector<vector<int>>& V_K) const;
    void _decode_SC(Node<Contents_depolarize>* node_cur);

private:
    void recursive_compute_llr(Node<Contents_depolarize>* node_cur, int depth);
    void recursive_propagate_sums(const Node<Contents_depolarize>* node_cur);
    void duplicate_path(int path, int leaf_index, vector<vector<Node<Contents_depolarize>*>> leaves_array, vector<int>& decisions);

    void recursive_duplicate_tree_llr(Node<Contents_depolarize>* node_a, Node<Contents_depolarize>* node_b);
    void recursive_duplicate_tree_sums(Node<Contents_depolarize>* node_a, Node<Contents_depolarize>* node_b, Node<Contents_depolarize>* node_caller);

protected:
    virtual void select_best_path();
    
    void recursive_allocate_nodes_contents(Node<Contents_depolarize>* node_curr, const int vector_size, int& max_depth_llrs);
    void recursive_initialize_frozen_bits(const Node<Contents_depolarize>* node_curr, const std::vector<bool>& frozen_bits);
    void recursive_deallocate_nodes_contents(Node<Contents_depolarize>* node_curr);

};

#endif /* DECODER_DEPOLARIZE_HPP_ */
