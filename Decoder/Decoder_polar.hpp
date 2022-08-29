#ifndef DECODER_POLAR_SCL_HPP_
#define DECODER_POLAR_SCL_HPP_

#include <vector>
#include <set>
#include "Decoder.hpp"
using namespace std;

class Contents_SCL
{
public:
    vector<float> l;     // probability pair array P
    vector<int> s;       // bits array B
    bool is_frozen_bit;
    int max_depth_llrs;

    explicit Contents_SCL(int size) : l(size), s(size), is_frozen_bit(false), max_depth_llrs(-1) {}
    virtual ~Contents_SCL() {}
};

// temporarily single kernel F=((1,1),(0,1))
class Decoder_polar_SCL : Decoder
{
protected:
    const float metric_init; // init value of the metrics in the trees
    const int L;             // maximum path number
    vector<uint32_t> stages; // number of stages = log_2(N)
    std::set<int> active_paths;

    vector<bool> frozen_bits;
    vector<Tree_metric<Contents_SCL>> polar_trees;
    vector<vector<Node<Contents_SCL>*>> leaves_array;

    vector<float> LLRs;
    vector<int> bits; 
    // vector<vector<function<float(const vector<float> &LLRs, const vector<int> &bits)>>> lambdas;
    vector<function<float(const vector<float> &LLRs, const vector<int> &bits)>> lambdas;

public:
    Decoder_polar_SCL(const int& K, const int& N, const int& L, const vector<bool>& frozen_bits);
            // vector<function<float(const vector<float> &LLRS, const vector<int> &bits)>> lambdas);
    virtual ~Decoder_polar_SCL();
    vector<uint32_t> get_stages() const { return stages; }

protected:
    void _load(const float *Y_N);
    void _decode(const size_t frame_id);

private:
    void recursive_compute_llr(Node<Contents_SCL>* node_cur, int depth);
    void recursive_propagate_sums(const Node<Contents_SCL>* node_cur);
    void duplicate_path(int path, int leaf_index, vector<vector<Node<Contents_SCL>*>> leaves_array);

protected:
    virtual void select_best_path(const size_t frame_id);
    
    void recursive_allocate_nodes_contents(Node<Contents_SCL>* node_curr, const int vector_size, int &max_depth_llrs);
    void recursive_initialize_frozen_bits(const Node<Contents_SCL>* node_curr, const std::vector<bool>& frozen_bits);
    void recursive_deallocate_nodes_contents(Node<Contents_SCL>* node_curr);

};

#endif /* DECODER_POLAR_SCL_HPP_ */