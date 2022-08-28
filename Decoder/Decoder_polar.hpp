#ifndef DECODER_POLAR_SCL
#define DECODER_POLAR_SCL

#include <vector>
#include <set>
using namespace std;

template <typename T = float>
class Node
{
    friend Tree_metric<T>;

public:
    Node<T>* father;       // nullptr if is root
    vector<Node<T>*> children;

    T* contents;
    const int depth;    // vertical indexing
    const int lane_id;  // horizontal indexing
    int child_id;

    // implicitly inline
    T* get_c() const { return contents; }
    Node<T>* get_father() const { return father; }
    bool is_root() const { return (father == nullptr); }
    int get_child_id() const { return child_id; }
    int get_depth() const { return depth; }
    vector<Node<T>*> get_children() const { return children; }
    T* get_contents() const { return contents; }
};

class Contents_SCL
{
public:
    vector<int> l;
    vector<float> s;
    bool is_frozen_bit;
    int max_depth_llrs;

    explicit Contents_SCL(int size) : l(size), s(size), is_frozen_bit(false), max_depth_llrs(-1) {}
    virtual ~Contents_SCL() {}
};

template <typename T = float>
class Tree_metric
{
// protected:
    const vector<uint32_t> sequence;
    int depth;
    Node<T>* root;

    float path_metric;

public:
    explicit Tree_metric(const int depth, const int base, float path_metric);
    explicit Tree_metric(vector<uint32_t> &sequence, float path_metric);
    virtual ~Tree_metric();
    float get_path_metric() const { return path_metric; }
    void set_path_metric(float path_metric) { this->path_metric = path_metric; }
    int get_depth() const { return depth; }
    Node<T>* get_root() const { return root; }
};

class Decoder
{
protected:
    const int K, N;
    vector<float> Y_N;
};
// TODO: move the aboce to Decoder.hpp

class Decoder_polar_SCL : Decoder
{
protected:
    const float metric_init; // init value of the metrics in the trees
    const int L;             // maximum path number
    std::set<int> active_paths;

    vector<bool> frozen_bits;
    vector<Tree_metric<Contents_SCL>> polar_trees;
    vector<vector<Node<Contents_SCL>*>> leaves_array;

    vector<float> LLRs;
    vector<int> bits; 
    vector<vector<function<float(const vector<float> &LLRs, const vector<int> &bits)>>> lambdas;

protected:
    void _load(const float *Y_N);
    void _decode(const size_t frame_id);

private:
    void recursive_compute_llr(Node<Contents_SCL>* node_cur, int depth);

protected:
    virtual void select_best_path(const size_t frame_id);

};

#endif /* DECODER_POLAR_SCL */