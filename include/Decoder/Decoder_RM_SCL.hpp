#ifndef DECODER_RM_SCL_HPP_
#define DECODER_RM_SCL_HPP_
#include <vector>
#include <set>
#include "Decoder.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Tree.hpp"
using namespace std;

enum Status {UNUSED, USED, DETERMINED};
class Contents_RM_SCL
{
public:
    int m, r;
    vector<double> l;     // Log-Likelihood Ratio array, size 2^m 
    vector<int> s;       // estimated codeword array u, size 2^m (but may not be fully used). 
    Status status = UNUSED;

    explicit Contents_RM_SCL(int m, int r, int size) 
    : m(m), r(r), l(size), s(size) {}
    virtual ~Contents_RM_SCL() {}
};

// temporarily single kernel F={{1,0},{1,1}}
class Decoder_RM_SCL : Decoder
{
protected:
    const int m, r;
    const double metric_init; // init value of the metrics in the trees
    const int L;             // maximum path number
    vector<uint32_t> stages; // length = log_2(N) = n
    // at each stage s\in {1,...,n}, use f- and f+ to update LLRs, 
    // the decoder uses the rule above to keep track of NlogN partial sums u_s^{i},
    // and update them after decoding each bit \hat{u}_i
    std::set<int> active_paths;

    vector<bool> frozen_bits;
    vector<RM_Tree_metric<Contents_RM_SCL>> rm_trees;
    // vector<vector<Node<Contents_RM_SCL>*>> leaves_array;   // possible leaf nodes to split on
    vector<Node<Contents_RM_SCL>*> active_nodes;  // it has the same size as the active_paths array
    // and stores for each rm_tree, at which node to split next
    vector<double> LLRs;
    vector<int> bits; 
	// tuples to be sorted. <Path,estimated codeword,metric>
	vector<tuple<Node<Contents_RM_SCL>*, vector<int>, double>> metrics_vec;

    // temporary arrays used to update partial sums, 
    // the only usage is in recursive_propagate_sums,
    // to avoid allocating new temp memory every time
    vector<uint32_t> idx;
    vector<int> u;
    // linearized kernel = {1,1,0,1}. TODO: maybe {1,0,1,1}
    vector<int> Ke;

public:
    Decoder_RM_SCL(const int& K, const int& N, const int& L, const vector<bool>& frozen_bits);
            // vector<function<double(const vector<double> &LLRS, const vector<int> &bits)>> lambdas);
    virtual ~Decoder_RM_SCL();
    virtual int decode(const double *Y_N, int *V_K, const size_t frame_id);

protected:
    void _load(const double *Y_N);
    void _decode(const size_t frame_id);
    // void _decode_llr(const int& m, const int& r); // recursive decoding
    void _decode_llr(Node<Contents_RM_SCL>* node_curr); // recursive decoding
    void _decode_mm(Node<Contents_RM_SCL>* node_curr, bool skip = false);
    void _decode_11(Node<Contents_RM_SCL>* node_curr);
    void _decode_m0(Node<Contents_RM_SCL>* node_curr, bool skip = false);
    void _store(int *V_K) const;

private:
    void recursive_compute_llr(Node<Contents_RM_SCL>* node_cur);
    void recursive_propagate_sums(const Node<Contents_RM_SCL>* node_cur);
    void duplicate_path(int path, int leaf_index, vector<vector<Node<Contents_RM_SCL>*>> leaves_array);

    void recursive_duplicate_tree_llr(Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b);
    void recursive_duplicate_tree_sums(Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b, Node<Contents_RM_SCL>* node_caller);

protected:
    virtual void select_best_path(const size_t frame_id);
    
    void recursive_allocate_nodes_contents(Node<Contents_RM_SCL>* node_curr, const int vector_size, int m, int r);
    void recursive_deallocate_nodes_contents(Node<Contents_RM_SCL>* node_curr);

};

#endif /* DECODER_RM_SCL_HPP_ */