#ifndef DECODER_RM_SCL_HPP_
#define DECODER_RM_SCL_HPP_
#include <vector>
#include <set>
#include "Decoder.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Tree.hpp"
using namespace std;

class Contents_RM_SCL
{
public:
    int m, r;
    vector<double> l;     // Log-Likelihood Ratio array, size 2^m 
    vector<int> s;       // estimated codeword array u, size 2^m 

    explicit Contents_RM_SCL(int m, int r, int size) 
    : m(m), r(r), l(size), s(size) {}
    virtual ~Contents_RM_SCL() {}
};

// temporarily single kernel F={{1,0},{1,1}}
class Decoder_RM_SCL : public Decoder
{
protected:
    const int m, r;
    const int L;             // maximum path number
    // at each stage s\in {1,...,n}, use f- and f+ to update LLRs, 
    // the decoder uses the rule above to keep track of NlogN partial sums u_s^{i},
    // and update them after decoding each bit \hat{u}_i
    std::set<int> active_paths;
    // the active_paths may not store contiguous number, but all path indices are within [0, L)

    vector<RM_Tree_metric<Contents_RM_SCL>*> rm_trees;
    // possible leaf nodes to split on, it has the same size as the active_paths array
    // and stores for each rm_tree, at which node to split next
    vector<Node<Contents_RM_SCL>*> active_nodes;  

	// tuples to be sorted. <node, path idx, estimated codeword, metric>
	vector<tuple<Node<Contents_RM_SCL>*, int, vector<int>, double>> metrics_vec;

public:
    Decoder_RM_SCL(const int& m, const int& r, const int& L);
    virtual ~Decoder_RM_SCL();
    // input: an LLR array. output: a codeword (or a list of codewords)
    // NOTE: output codeword instead of information bits 
    virtual int decode(const double *Y_N, int *V_K, const size_t frame_id);
    void test_copy_until();
    void test_assign_path_idx();

protected:
    void _decode_leaf_llr(Node<Contents_RM_SCL>* node_curr, int i, double& pm); 
    void _decode_mm(Node<Contents_RM_SCL>* node_curr, int i, double& pm); 
    void _decode_11(Node<Contents_RM_SCL>* node_curr, int i, double& pm);
    void _decode_m0(Node<Contents_RM_SCL>* node_curr, int i, double& pm);

private:
    void recursive_compute_llr(int path_idx, Node<Contents_RM_SCL>* node_cur);
    Node<Contents_RM_SCL>* recursive_propagate_sums(const Node<Contents_RM_SCL>* node_cur);
    bool partition_and_copy();
    Node<Contents_RM_SCL>* copy_until(Node<Contents_RM_SCL>* node_stop, Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b);

protected:
    int  select_best_path(const size_t frame_id);
    void recursive_allocate_nodes_contents(Node<Contents_RM_SCL>* node_curr, const int vector_size, int m, int r);
    void recursive_deallocate_nodes_contents(Node<Contents_RM_SCL>* node_curr);

};

#endif /* DECODER_RM_SCL_HPP_ */