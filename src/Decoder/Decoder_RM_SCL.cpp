/*
for every LLR bit of (u+v) + u
greedily choose whether it is 0 or 1 (by comparing path metric)

then partition and cut the list
finalize the remaining branch

(0, m)x2, (m, m)x4, (1,1)x4
(m, m) add to path likelihoods.
Optimized branching.
List as a tree.
*/

#include "Decoder/Decoder_RM_SCL.hpp"
#include <iostream>
#include <numeric>
using namespace std;

Decoder_RM_SCL::Decoder_RM_SCL(const int& m, const int& r, const int& L) 
: Decoder(Encoder_RM::calculate_K(m,r), 1 << m), m(m), r(r), L(L), active_nodes(L, nullptr)
{
    for (auto i = 0; i < L; i++) {
        auto new_tree = new RM_Tree_metric<Contents_RM_SCL>(m, r, numeric_limits<double>::min());
        this->recursive_allocate_nodes_contents(new_tree->get_root(), this->N, this->m, this->r);
        this->rm_trees.push_back(new_tree);
    }
}

void Decoder_RM_SCL::recursive_allocate_nodes_contents(Node<Contents_RM_SCL>* node_curr, const int vector_size, int m, int r)
{
	node_curr->set_contents(new Contents_RM_SCL(m, r, vector_size));

    auto children = node_curr->get_children();
	if (!node_curr->is_leaf())
	{
        // assert(children.size() == 2);
		const auto new_vector_size = vector_size / 2;
		this->recursive_allocate_nodes_contents(children[0], new_vector_size, m-1, r-1);
		this->recursive_allocate_nodes_contents(children[1], new_vector_size, m-1, r);
	} 
}

Decoder_RM_SCL::~Decoder_RM_SCL()
{
    for (auto i = 0; i < L; i++)
        this->recursive_deallocate_nodes_contents(this->rm_trees[i]->get_root());
}

void Decoder_RM_SCL::recursive_deallocate_nodes_contents(Node<Contents_RM_SCL>* node_curr)
{
	for (auto c : node_curr->get_children())
		this->recursive_deallocate_nodes_contents(c); 

	delete node_curr->get_contents();
	node_curr->set_contents(nullptr);
}

void Decoder_RM_SCL::select_best_path(const size_t frame_id)
{   // select the best one, not the best L ones.
	int best_path = 0;
	if (active_paths.size() >= 1)
		best_path = *active_paths.begin();

	for (int path : active_paths)
		if(rm_trees[path]->get_path_metric() < rm_trees[best_path]->get_path_metric())
			best_path = path;

	active_paths.clear();
	active_paths.insert(best_path);
}

void Decoder_RM_SCL::recursive_compute_llr(int path_idx, Node<Contents_RM_SCL>* node_curr)
{
    if (node_curr->is_leaf()) {
        active_nodes[path_idx] = node_curr;
        return;
    }
    Contents_RM_SCL* c = node_curr->get_c();
    int m = c->m, N = 1 << m, N_half = N / 2;
    vector<double> llr = c->l;
    vector<Node<Contents_RM_SCL>*> children = node_curr->get_children();
    auto left_child = children[0];
    Decoder::f_plus(llr.data(), llr.data() + N_half, N_half, left_child->get_c()->l.data());
    recursive_compute_llr(path_idx, left_child);
}


int Decoder_RM_SCL::decode(const double *Y_N, int* X_N, const size_t frame_id)
{
    // initialization: only one active path, 
    // and it's active node is RM(m-r, 0), obtained by recursive_compute_llr
    this->active_paths.clear();
    this->active_paths.insert(0);
    for (auto i = 0; i < L; i++) {
        rm_trees[i]->get_root()->get_c()->l = vector<double>(Y_N, Y_N + this->N);
        // this->recursive_compute_llr(rm_trees[i]->get_root());
    }
    this->recursive_compute_llr(0, rm_trees[0]->get_root());
    // the active node of each rm_tree is in the same position, but with different content
    while (true) {
        cerr << "active_nodes size " << active_nodes.size() << endl;
        cerr << "active_path size " << active_paths.size() << endl;
        for (int path_idx : active_paths) {
            cerr << "path_idx=" << path_idx << endl;
            assert(active_nodes[path_idx] != nullptr);
            assert(rm_trees[path_idx] != nullptr);
            _decode_leaf_llr(active_nodes[path_idx], path_idx, rm_trees[path_idx]->path_metric); // split nodes
        }
        // an active node is always a left child, except RM(1,1)
        // change active nodes' status from used to determined
        cerr << "metrics_vec size " << metrics_vec.size() << endl;
        active_paths.clear();
        if (partition_and_copy()) break;
        // cerr << "updated active_nodes size " << active_nodes.size() << endl;
        metrics_vec.clear();
    }

	// this->select_best_path(frame_id);
    auto s = this->rm_trees[0]->get_root()->get_c()->s;
    for (int i : s)
        cerr << i << " ";
    cerr << endl;
    copy(s.begin(), s.end(), X_N);
    return 0;
}


bool Decoder_RM_SCL::partition_and_copy()
{
    // sort path metric in increasing order
    std::sort(metrics_vec.begin(), metrics_vec.end(),
        [](std::tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> x, std::tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> y){
            return std::get<3>(x) < std::get<3>(y);
        });
    // remove worst metrics from list
    if (metrics_vec.size() > L)
        metrics_vec.resize(L);
    vector<int> all_indices(L);
    iota(all_indices.begin(), all_indices.end(), 0);
    set<int> used_path_indices;
    transform(metrics_vec.begin(), metrics_vec.end(), inserter(used_path_indices, used_path_indices.begin()),
            [] (auto& m) { return get<1>(m); });
    vector<int> path_indices_to_reuse;
    set_difference(all_indices.begin(), all_indices.end(), used_path_indices.begin(), used_path_indices.end(), 
            inserter(path_indices_to_reuse, path_indices_to_reuse.begin()));
    vector<bool> used(L, false);
    cerr << "partition_and_copy metrics_vec size " << metrics_vec.size() << endl;
    bool reach_root = false;
    for (auto m : metrics_vec) {
        auto cur_path_node = get<0>(m);
        int cur_path_idx = get<1>(m);
        if (!used[cur_path_idx]) {
            used[cur_path_idx] = true;
            active_paths.insert(cur_path_idx);
            cerr << cur_path_idx << " uses idx " << cur_path_idx << endl;
            // copy s array
            get<0>(m)->get_contents()->s = get<2>(m);
            rm_trees[cur_path_idx]->set_path_metric(get<3>(m));
        } else {
            auto idx_it = path_indices_to_reuse.begin();
            int new_path_idx = *idx_it;
            path_indices_to_reuse.erase(idx_it);
            // assert(!used[idx]);
            cerr << cur_path_idx << " uses new idx " << new_path_idx << endl;
            used[new_path_idx] = true;
            active_paths.insert(new_path_idx);
            rm_trees[new_path_idx]->set_path_metric(get<3>(m));
            assert(rm_trees[cur_path_idx]->get_root() != nullptr);
            assert(rm_trees[new_path_idx]->get_root() != nullptr);
            Node<Contents_RM_SCL>* new_path_node = copy_until(cur_path_node, rm_trees[cur_path_idx]->get_root(), rm_trees[new_path_idx]->get_root());
            assert(new_path_node != nullptr);
            new_path_node->get_contents()->s = get<2>(m);
            cur_path_node = new_path_node;
            cur_path_idx = new_path_idx;
        }
        // TODO: propagate back LLR to obtain bit estimation and return the node to stop
        auto next_node = recursive_propagate_sums(cur_path_node);
        if (next_node == nullptr) {
            reach_root = true;
            continue;
        }
        // TODO: recursive_compute_llr to this node, active_nodes will be updated in this step for each tree
        recursive_compute_llr(cur_path_idx, next_node);
    }
    return reach_root;
}

Node<Contents_RM_SCL>* Decoder_RM_SCL::copy_until(Node<Contents_RM_SCL>* node_stop, Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b)
{   // copy tree 1 (to which node_a and node_stop belong) into copy tree 2 (to which node_b belongs) until node_stop
    // return: whether reached or not
    auto c_a = node_a->get_contents();
    auto c_b = node_b->get_contents();
    // cerr << "in node with m=" << c->m << ", r=" << c->r << endl;
    auto children_a = node_a->get_children();
    auto children_b = node_b->get_children();
    if (node_a->is_root()) {
        assert(node_stop != node_a);
        auto new_node_b = copy_until(node_stop, children_a[0], children_b[0]);
        if (new_node_b != nullptr) return new_node_b;
        else return copy_until(node_stop, children_a[1], children_b[1]);
    } else if (node_a->is_leaf()) {
        if (node_a == node_stop) {
            cerr << "leaf reached m=" << c_b->m << ", r=" << c_b->r << endl;
            return node_b;
        }
        else {
            c_b->l = c_a->l;
            c_b->s = c_a->s;
            return nullptr;
        }
    } else {
        auto new_node_b = copy_until(node_stop, children_a[0], children_b[0]);
        if (new_node_b != nullptr) return new_node_b;
        else  {
            new_node_b = copy_until(node_stop, children_a[1], children_b[1]);
            if (new_node_b != nullptr) return new_node_b;
        }
        // copy self
        if (node_stop == node_a) {
            cerr << "reached m=" << c_b->m << ", r=" << c_b->r << endl;
            return node_b; 
        } else {
            c_b->l = c_a->l;
            c_b->s = c_a->s;
            return nullptr;
        }

    }
}

Node<Contents_RM_SCL>* Decoder_RM_SCL::recursive_propagate_sums(const Node<Contents_RM_SCL>* node_cur)
{   // return the next node to consider
    if (node_cur->is_root()) {
        // reached the root, finish
        return nullptr;
    } else {
        auto p = node_cur->get_father();
        auto p_c = p->get_contents();
        auto children = p->get_children();
        int N_half = 1 << (node_cur->get_contents()->m);
        if (node_cur == children[0]) { // if is the left child
            auto right_child = children[1];
            f_minus(p_c->l.data(), p_c->l.data() + N_half, node_cur->get_contents()->s.data(), N_half, right_child->get_contents()->l.data());
            return right_child;
        } else {                       // if is the right child
            auto left_child = children[0];
            vector<int> left_s = left_child->get_contents()->s;
            vector<int> right_s = node_cur->get_contents()->s;
            assert(N_half == left_s.size());
            copy(right_s.begin(), right_s.end(), p_c->s.begin() + N_half);
            for (int i = 0; i < N_half; i++) {
                p_c->s[i] = left_s[i] ^ right_s[i];
            }
            return recursive_propagate_sums(p); 
        }
    }
}

void Decoder_RM_SCL::_decode_leaf_llr(Node<Contents_RM_SCL>* node_curr, int i, double& pm)
{   // an active node is always a leaf (RM(m,0), RM(m,m), RM(1,1))
    // because otherwise can propogate back
    auto c = node_curr->get_c();
    int m = c->m, r = c->r, N = 1 << m;
    cerr << "at leaf node m=" << m << ", r=" << r << endl;
    if (m == r) {
        if (r == 1) {
            _decode_11(node_curr, i, pm);
        } else {
            _decode_mm(node_curr, i, pm);
        }
    } else if (r == 0) {
        _decode_m0(node_curr, i, pm);
    } 
}

void Decoder_RM_SCL::_decode_m0(Node<Contents_RM_SCL>* node_curr, int i, double& pm)
{   // repetition code, take both 0 and 1
    // calculate path metric for estimate to 0
    // calculate path metric for estimate to 1
    cerr << "in decode m0" << endl;
    auto c = node_curr->get_c();
    int m = c->m, N = 1 << m;
    vector<double> llr = c->l;
    double pm0 = pm, pm1 = pm;
    for (int i = 0; i < N; i++) {
        pm0 = phi(pm0, llr[i], 0);
        pm1 = phi(pm1, llr[i], 1);
    }
    cerr << "after phi" << endl;
    metrics_vec.push_back(make_tuple(node_curr, i, vector(N, 0), pm0));
    metrics_vec.push_back(make_tuple(node_curr, i, vector(N, 1), pm1));
    cerr << "at the end of decode_m0" << endl;
}

void Decoder_RM_SCL::_decode_11(Node<Contents_RM_SCL>* node_curr, int i, double& pm)
{   // push back all four possiblities: 00, 01, 10, 11
    double pm00 = pm, pm01 = pm, pm10 = pm, pm11 = pm;
    auto c = node_curr->get_c();
    int m = c->m, N = 1 << m;
    vector<double> llr = c->l;  // expect size 2
    pm00 = phi(pm00, llr[0], 0); pm00 = phi(pm00, llr[1], 0);
    pm01 = phi(pm01, llr[0], 0); pm01 = phi(pm01, llr[1], 1);
    pm10 = phi(pm10, llr[0], 1); pm10 = phi(pm10, llr[1], 0);
    pm11 = phi(pm11, llr[0], 1); pm11 = phi(pm11, llr[1], 1);
    metrics_vec.push_back(make_tuple(node_curr, i, vector<int>{0,0}, pm00));
    metrics_vec.push_back(make_tuple(node_curr, i, vector<int>{0,1}, pm01));
    metrics_vec.push_back(make_tuple(node_curr, i, vector<int>{1,0}, pm10));
    metrics_vec.push_back(make_tuple(node_curr, i, vector<int>{1,1}, pm11));
}

void Decoder_RM_SCL::_decode_mm(Node<Contents_RM_SCL>* node_curr, int i, double& pm)
{
    auto c = node_curr->get_c();
    int m = c->m, N = 1 << m;
    vector<double> llr = c->l; 
    vector<int> tmp(N);
    double pm_min = pm;
    // temporarily only take 1 instead of 4. TODO: take 4
    for (int i = 0; i < N; i++) {
        if (llr[i] > 0) {
            tmp[i] = 0;
            pm_min = phi(pm_min, llr[i], tmp[i]); // expect \Delta PM = 0
        } else {
            tmp[i] = 1;
            pm_min = phi(pm_min, llr[i], tmp[i]); // expect \Delta PM = 0
        }
    }
    metrics_vec.push_back(make_tuple(node_curr, i, tmp, pm_min));
}

//*************************** TESTS **********************************************
void Decoder_RM_SCL::test_copy_until()
{
    auto r1 = rm_trees[0]->get_root(), r2 = rm_trees[1]->get_root();
    // auto c = r1->get_contents();
    // auto node_stop = (((r1->get_children()[0])->get_children()[0])->get_children()[0])->get_children()[1];
    auto node_stop = r1;
    while (!node_stop->is_leaf())
        node_stop = node_stop->get_children()[0];
    // cerr << "in test_copy_until" << endl;
    auto node_reached = copy_until(node_stop->get_father(), r1, r2);
    // auto c_reached = node_reached->get_contents();
    // cerr << "reaced node m=" << c_reached->m << ", r=" << c_reached->r << endl;

}

void Decoder_RM_SCL::test_assign_path_idx()
{
    metrics_vec.clear();
    metrics_vec.push_back(make_tuple(nullptr, 1, vector<int>(), 0.0));
    metrics_vec.push_back(make_tuple(nullptr, 3, vector<int>(), 0.0));
    metrics_vec.push_back(make_tuple(nullptr, 5, vector<int>(), 0.0));
    metrics_vec.push_back(make_tuple(nullptr, 8, vector<int>(), 0.0));
    metrics_vec.push_back(make_tuple(nullptr, 3, vector<int>(), 0.0));
    metrics_vec.push_back(make_tuple(nullptr, 4, vector<int>(), 0.0));
    metrics_vec.push_back(make_tuple(nullptr, 1, vector<int>(), 0.0));



    vector<int> all_indices(L);
    iota(all_indices.begin(), all_indices.end(), 0);
    set<int> used_path_indices;
    transform(metrics_vec.begin(), metrics_vec.end(), inserter(used_path_indices, used_path_indices.begin()),
            [] (auto& m) { return get<1>(m); });
    vector<int> path_indices_to_reuse;
    set_difference(all_indices.begin(), all_indices.end(), used_path_indices.begin(), used_path_indices.end(), 
            inserter(path_indices_to_reuse, path_indices_to_reuse.begin()));
    for (int i : path_indices_to_reuse)
        cerr << i << " ";
    cerr << endl;
    // path_indices_to_reuse = {0, 2, 6, 7, 9}
    vector<bool> used(L, false);
    for(auto m : metrics_vec) {
        int cur_path_idx = get<1>(m);
        if (!used[cur_path_idx]) {
            used[cur_path_idx] = true;
            cerr << cur_path_idx << " uses idx " << cur_path_idx << endl;
            // copy s array
        } else {
            auto idx_it = path_indices_to_reuse.begin();
            int idx = *idx_it;
            path_indices_to_reuse.erase(idx_it);
            assert(!used[idx]);
            cerr << cur_path_idx << " uses idx " << idx << endl;
            used[idx] = true;
        }
    }

}

//*************************** UNUSED **********************************************
void Decoder_RM_SCL::duplicate_path(int path, Node<Contents_RM_SCL>* leaf_node)
{

}

void Decoder_RM_SCL::recursive_duplicate_tree_llr(Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b)
{   // duplicate starts from a leaf, and propagate up until the root
    node_b->get_c()->l = node_a->get_c()->l;

    if (!node_a->get_father()->is_root())
        this->recursive_duplicate_tree_llr(node_a->get_father(), node_b->get_father());

}

void Decoder_RM_SCL::recursive_duplicate_tree_sums(Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b, Node<Contents_RM_SCL>* node_caller)
{   // again, propagate upwards until the root
    // do not copy if the current node is a leaf
    if (!node_a->is_leaf()) {
        auto c2 = node_b->get_children().begin();
        for (auto c1 : node_a->get_children()) {
            if (c1 != node_caller)
                (*c2)->get_c()->s = c1->get_c()->s;
            ++c2;
        }
    }

    if (!node_a->is_root())
        recursive_duplicate_tree_sums(node_a->get_father(), node_b->get_father(), node_a);

}

/*
void Decoder_RM_SCL::get_pre_in_post_order(Node<Contents_RM_SCL>* node_curr, Node<Contents_RM_SCL>* node_caller)
{   // get all the nodes before the current node in the post-order traversal of the tree
    // post-order: left, right, father 
    reversed_post_order(node_curr);
    auto p = node_curr->get_father();
    if (p != nullptr) { // not the root
        auto children = p->get_children();
        if (node_curr == children[0]) {  
            // if node_curr is a left child
            // node_curr, node_curr->p (do not print), 
        } else if (node_curr == children[1]) {
            // if node_curr is a right child
            // node_curr->p->left_child, node_curr->p (do not print)
        } else {
            // this should not happen, so add an assert to check
            assert(0);
            cerr << "more than two childern" << endl;
        }
    } else {
        // is the root, the last node in the post-order traversal
    }


}

void Decoder_RM_SCL::reversed_post_order(Node<Contents_RM_SCL>* node_curr) 
{
    // copy self, right, left
    auto children = node_curr->get_children();
    if (children.size() == 0) {
        // is a leaf, copy itself and return
        return;
    }
    reversed_post_order(children[1]);
    reversed_post_order(children[0]);
}
*/