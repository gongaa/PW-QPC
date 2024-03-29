#ifndef TREE_HPP_
#define TREE_HPP_
#include <vector>
#include <set>
using namespace std;
template <typename T = double> class Tree_metric;

template <typename T = double>
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
    Node(Node<T>*father, vector<Node<T>*> children, T* contents, int depth, int lane_id, int child_id = -1)
    : father(father), children(children), contents(contents), depth(depth), lane_id(lane_id), child_id(child_id) {}
    T* get_c() const { return contents; }
    Node<T>* get_father() const { return father; }
    bool is_root() const { return (father == nullptr); }
    int get_child_id() const { return child_id; }
    int get_depth() const { return depth; }
    vector<Node<T>*> get_children() const { return children; }
    int get_lane_id() const { return lane_id; }
    T* get_contents() const { return contents; }
    void set_contents(T* contents) { this->contents = contents; }
    bool is_leaf() const { return !children.size(); }
    bool is_empty() const { return (contents == nullptr); }
};

template <typename T>
class Tree_metric
{
public:
    int depth;
    Node<T>* root;
    vector<Node<T>*> leaves;
    double path_metric;

    explicit Tree_metric(const int depth, double path_metric);
    virtual ~Tree_metric();
    double get_path_metric() const { return path_metric; }
    void set_path_metric(double path_metric) { this->path_metric = path_metric; }
    int get_depth() const { return depth; }
    Node<T>* get_root() const { return root; }
    vector<Node<T>*> get_leaves() const { return leaves; }

protected:
    void init();
    void recursive_get_leaves(Node<T>* cur_node);
    void create_nodes(Node<T>* cur_node, int cur_depth, vector<int>& lanes);
    void delete_nodes(Node<T>* cur_node);
};

template <typename T>
class RM_Tree_metric 
{
public:
    int m, r;
    Node<T>* root;
    double path_metric;
    explicit RM_Tree_metric(const int m, const int r, double path_metric);
    virtual ~RM_Tree_metric();
    Node<T>* get_root() const { return root; }
    double get_path_metric() const { return path_metric; }
    void set_path_metric(double path_metric) { this->path_metric = path_metric; }

private:
    void RM_init();
    void create_RM_nodes(Node<T>* cur_node, int cur_depth, int m, int r, vector<int>& lanes);
    void delete_RM_nodes(Node<T>* cur_node);
};

// template should have methods' implementation inline
// a workaround to this is
#include "Tree.hxx"

#endif /* Tree.hpp */