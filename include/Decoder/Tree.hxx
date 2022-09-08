#include <vector>
#include <cmath>
#include "Tree.hpp"
using namespace std;

template <typename T>
Tree_metric<T>::Tree_metric(const int depth, const int base, double path_metric)
: sequence(vector<uint32_t>(depth, base)), depth(depth), root(nullptr), path_metric(path_metric)
{
    this->init();
}

template <typename T>
Tree_metric<T>::Tree_metric(vector<uint32_t> &sequence, double path_metric)
: sequence(sequence), depth(sequence.size() + 1), root(nullptr), path_metric(path_metric)
{
    this->init();
}

template <typename T>
void Tree_metric<T>::init()
{
    vector<int> lanes(depth + 1, 0);

    auto cur_depth = 0;
    this->root = new Node<T>(nullptr, vector<Node<T>*>(), nullptr, cur_depth, lanes[cur_depth]);
    this->create_nodes(this->root, cur_depth + 1, lanes, sequence);
    recursive_get_leaves(this->get_root());
}

template <typename T>
void Tree_metric<T>::create_nodes(Node<T>* cur_node, int cur_depth, vector<int> &lanes, const vector<uint32_t> &sequence)
{
    if (cur_depth < this->depth) {
        for (auto c = 0; c < (int)sequence[cur_depth - 1]; c++) {
            auto child = new Node<T>(cur_node, vector<Node<T>*>(), nullptr, cur_depth, lanes[cur_depth]++, c);
            cur_node->children.push_back(child);
            this->create_nodes(child, cur_depth + 1, lanes, sequence);
        }
    }
}

template <typename T>
void Tree_metric<T>::recursive_get_leaves(Node<T>* cur_node) 
{
    if (cur_node->is_leaf())
        leaves.push_back(cur_node);
    else
        for (auto c : cur_node->children)
            recursive_get_leaves(c);
}


template <typename T>
Tree_metric<T>::~Tree_metric()
{
    this->delete_nodes(this->root);
}

template <typename T>
void Tree_metric<T>::delete_nodes(Node<T>* cur_node)
{
    if (cur_node != nullptr) {
        for (auto c : cur_node->children)
            this->delete_nodes(c);
        delete cur_node;
    }
}

template <typename T>
void RM_Tree_metric<T>::RM_init()
{
    vector<int> lanes(depth + 1, 0);

    auto cur_depth = 0;
    this->root = new Node<T>(nullptr, vector<Node<T>*>(), nullptr, cur_depth, lanes[cur_depth]);
    this->create_RM_nodes(this->root, cur_depth + 1, this->m, this->r, lanes, sequence);
    recursive_get_leaves(this->get_root());
}

template <typename T>
void RM_Tree_metric<T>::create_RM_nodes(Node<T>* cur_node, int cur_depth, int m, int r, vector<int> &lanes)
{
    if (m == 0 || r == 0 || m == r) {
        return; // leaves are either the repetition code or the full code
    } else {
        auto left_child = new Node<T>(cur_node, vector<Node<T>*>(), nullptr, cur_depth, lanes[cur_depth]++, 0);
        cur_node->children.push_back(left_child);
        this->create_nodes(left_child, cur_depth + 1, m-1, r-1, lanes);

        auto right_child = new Node<T>(cur_node, vector<Node<T>*>(), nullptr, cur_depth, lanes[cur_depth]++, 1);
        cur_node->children.push_back(right_child);
        this->create_nodes(right_child, cur_depth + 1, m-1, r, lanes);
    }

}

template <typename T>
RM_Tree_metric<T>::RM_Tree_metric(const int m, const int r, double path_metric)
: Tree_metric<T>(m, 2, path_metric), m(m), r(r)
{
    this->RM_init();
}

template <typename T>
RM_Tree_metric<T>::~RM_Tree_metric()
{
    this->delete_nodes(this->root);
}