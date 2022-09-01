#ifndef DECODER_POLAR_SCL_HPP_
#define DECODER_POLAR_SCL_HPP_
// based on the paper "LLR-Based Successive Canncellation List Decoding of Polar Codes"
#include <vector>
#include <set>
#include "Decoder.hpp"
using namespace std;

class Contents_SCL
{
public:
    vector<double> l;     // Log-Likelihood Ratio array, size N
    vector<int> s;       // partial sum array u, size N (but may not be fully used). 
    // If root: effective size N. In total, Nlog(N) partial sums.
    // update ruls for u_n^{(i)} := \hat{u}_i:
    // u_{s-1}^{(2i-[i mode 2^{s-1}])} = u_s^{(2i)}\oplus u_s^{(2i+1)}
    // u_{s-1}^{(2^2+2i-[i mod 2^{s-1}])} = u_s^{(2i+1)}
    bool is_frozen_bit;
    int max_depth_llrs;

    explicit Contents_SCL(int size) : l(size), s(size), is_frozen_bit(false), max_depth_llrs(-1) {}
    virtual ~Contents_SCL() {}
};

// temporarily single kernel F={{1,0},{1,1}}
class Decoder_polar_SCL : Decoder
{
protected:
    const double metric_init; // init value of the metrics in the trees
    const int L;             // maximum path number
    vector<uint32_t> stages; // length = log_2(N) = n
    // at each stage s\in {1,...,n}, use f- and f+ to update LLRs, 
    // the decoder uses the rule above to keep track of NlogN partial sums u_s^{i},
    // and update them after decoding each bit \hat{u}_i
    std::set<int> active_paths;

    vector<bool> frozen_bits;
    vector<Tree_metric<Contents_SCL>> polar_trees;
    vector<vector<Node<Contents_SCL>*>> leaves_array;   // possible leaf nodes to split on
    vector<double> LLRs;
    vector<int> bits; 

    // temporary arrays used to update partial sums, 
    // the only usage is in recursive_propagate_sums,
    // to avoid allocating new temp memory every time
    vector<uint32_t> idx;
    vector<int> u;
    // linearized kernel = {1,1,0,1}. TODO: maybe {1,0,1,1}
    vector<int> Ke;

    // refactor lambdas into the Decoder class, since RM needs this too
    // vector<vector<function<double(const vector<double> &LLRs, const vector<int> &bits)>>> lambdas;

public:
    Decoder_polar_SCL(const int& K, const int& N, const int& L, const vector<bool>& frozen_bits);
            // vector<function<double(const vector<double> &LLRS, const vector<int> &bits)>> lambdas);
    virtual ~Decoder_polar_SCL();
    virtual int decode(const double *Y_N, int *V_K, const size_t frame_id);

protected:
    void _load(const double *Y_N);
    void _decode(const size_t frame_id);
    void _store(int *V_K) const;

private:
    void recursive_compute_llr(Node<Contents_SCL>* node_cur, int depth);
    void recursive_propagate_sums(const Node<Contents_SCL>* node_cur);
    void duplicate_path(int path, int leaf_index, vector<vector<Node<Contents_SCL>*>> leaves_array);

    void recursive_duplicate_tree_llr(Node<Contents_SCL>* node_a, Node<Contents_SCL>* node_b);
    void recursive_duplicate_tree_sums(Node<Contents_SCL>* node_a, Node<Contents_SCL>* node_b, Node<Contents_SCL>* node_caller);

protected:
    virtual void select_best_path(const size_t frame_id);
    
    void recursive_allocate_nodes_contents(Node<Contents_SCL>* node_curr, const int vector_size, int &max_depth_llrs);
    void recursive_initialize_frozen_bits(const Node<Contents_SCL>* node_curr, const std::vector<bool>& frozen_bits);
    void recursive_deallocate_nodes_contents(Node<Contents_SCL>* node_curr);

};

#endif /* DECODER_POLAR_SCL_HPP_ */
