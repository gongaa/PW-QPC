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
#include <algorithm>
#include <set>
#include <tuple>
#include <cassert>
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
		const int new_vector_size = vector_size / 2;
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

int Decoder_RM_SCL::select_best_path(const size_t frame_id)
{   // select the best one, not the best L ones.
	int best_path = 0;
	if (active_paths.size() >= 1)
		best_path = *active_paths.begin();

	for (int path : active_paths)
		if (rm_trees[path]->get_path_metric() < rm_trees[best_path]->get_path_metric())
			best_path = path;
    return best_path;
}

int Decoder_RM_SCL::is_codeword_in_list(const int* X_N)
{
    bool is_same = true;
    int num_flips = 0;
    int path_idx_same = -1;
    for (int path : active_paths) {
        is_same = true;
        num_flips = 0;
        const int* X = rm_trees[path]->get_root()->get_contents()->s.data();
        for (int i = 0; i < N; i++) {
            if (X[i] != X_N[i]) {
                is_same = false;
                num_flips++;
                // break;
            }
        }
        // cerr << "path_idx=" << path << ", metric=" << rm_trees[path]->get_path_metric() << ", num_flips=" << num_flips << endl;
        if (is_same) {
            path_idx_same = path;
            // cerr << "codeword matches, path_idx=" << path_idx_same << endl;
            // return path;
        }
    }
    return path_idx_same;
}

void Decoder_RM_SCL::copy_codeword_list(vector<vector<int>>& c_list, vector<double>& pm_list)
{
    int i = 0;
    for (int path : active_paths) {
        const int* X = rm_trees[path]->get_root()->get_contents()->s.data();
        copy(X, X+N, c_list[i].data());
        pm_list[i] = rm_trees[path]->get_path_metric(); 
        ++i;
    }
}

void Decoder_RM_SCL::recursive_compute_llr(int path_idx, Node<Contents_RM_SCL>* node_curr)
{
    if (node_curr->is_leaf()) {
        active_nodes[path_idx] = node_curr;
        return;
    }
    Contents_RM_SCL* c = node_curr->get_c();
    int m = c->m, N = 1 << m, N_half = N / 2;
    vector<double>& llr = c->l;
    vector<Node<Contents_RM_SCL>*> children = node_curr->get_children();
    auto left_child = children[0];
    Decoder::f_plus(llr.data(), llr.data() + N_half, N_half, left_child->get_c()->l.data());
    recursive_compute_llr(path_idx, left_child);
}


double Decoder_RM_SCL::decode(const double *Y_N, int* X_N, const size_t frame_id)
{
    // initialization: only one active path, 
    // and it's active node is RM(m-r, 0), obtained by recursive_compute_llr
    metrics_vec.clear();
    this->active_paths.clear();
    this->active_paths.insert(0);
    for (auto i = 0; i < L; i++) {
        rm_trees[i]->get_root()->get_c()->l = vector<double>(Y_N, Y_N + this->N);
        this->recursive_compute_llr(i, rm_trees[i]->get_root());
        rm_trees[i]->set_path_metric(0.0);
    }
    // this->recursive_compute_llr(0, rm_trees[0]->get_root());
    // the active node of each rm_tree is in the same position, but with different content
    while (true) {
        for (int path_idx : active_paths) {
            _decode_leaf_llr(active_nodes[path_idx], path_idx, rm_trees[path_idx]->path_metric); // split nodes
        }
        active_paths.clear();
        if (partition_and_copy()) break; // update active_paths
        metrics_vec.clear();
    }
    // cerr << "active_path size " << active_paths.size() << endl;
	int best_path = this->select_best_path(frame_id);
    auto s = this->rm_trees[best_path]->get_root()->get_c()->s;
    copy(s.begin(), s.end(), X_N);
    return this->rm_trees[best_path]->get_path_metric();
}

bool Decoder_RM_SCL::partition_and_copy()
{
    // sort path metric in increasing order
    std::sort(metrics_vec.begin(), metrics_vec.end(),
        [](std::tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> x, std::tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> y) {
            return std::get<3>(x) < std::get<3>(y);
        });
    // remove worst metrics from list
    if (metrics_vec.size() > L)
        metrics_vec.resize(L);
    vector<int> all_indices(L); // set to 0,...,L-1
    iota(all_indices.begin(), all_indices.end(), 0);
    set<int> used_path_indices;
    transform(metrics_vec.begin(), metrics_vec.end(), inserter(used_path_indices, used_path_indices.begin()),
            [] (auto& m) { return get<1>(m); });
    vector<int> path_indices_to_reuse; // = all_indices \ used_path_indices
    set_difference(all_indices.begin(), all_indices.end(), used_path_indices.begin(), used_path_indices.end(), 
            inserter(path_indices_to_reuse, path_indices_to_reuse.begin()));
    vector<bool> used(L, false);
    bool reach_root = false;
    for (auto m : metrics_vec) {
        auto cur_path_node = get<0>(m);
        int cur_path_idx = get<1>(m);
        if (!used[cur_path_idx]) {
            used[cur_path_idx] = true;
            active_paths.insert(cur_path_idx);
            cur_path_node->get_contents()->s = get<2>(m);
            rm_trees[cur_path_idx]->set_path_metric(get<3>(m));
        } else {
            auto idx_it = path_indices_to_reuse.begin();
            int new_path_idx = *idx_it;
            path_indices_to_reuse.erase(idx_it);
            used[new_path_idx] = true;
            active_paths.insert(new_path_idx);
            Node<Contents_RM_SCL>* new_path_node = copy_until(cur_path_node, rm_trees[cur_path_idx]->get_root(), rm_trees[new_path_idx]->get_root());
            new_path_node->get_contents()->s = get<2>(m);
            rm_trees[new_path_idx]->set_path_metric(get<3>(m));
            cur_path_node = new_path_node;
            cur_path_idx = new_path_idx;
        }
        // propagate back LLR to obtain bit estimation and return the node to stop
        auto next_node = recursive_propagate_sums(cur_path_node);
        if (next_node == nullptr) 
            reach_root = true;
        else  
            // active_nodes will be updated in this step for each tree
            recursive_compute_llr(cur_path_idx, next_node);
    }
    return reach_root;
}

Node<Contents_RM_SCL>* Decoder_RM_SCL::copy_until(Node<Contents_RM_SCL>* node_stop, Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b)
{   // copy tree 1 (to which node_a and node_stop belong) into copy tree 2 (to which node_b belongs) until node_stop
    // return: whether reached or not
    // instead of copying the whole tree until a node
    // only need to copy the s array for all the leaves before this node, and the left child if in its right subtree
    // and copy the l array along the path that was used to compute the llr for the current code
    auto c_a = node_a->get_contents();
    auto c_b = node_b->get_contents();
    auto children_a = node_a->get_children();
    auto children_b = node_b->get_children();
    if (node_a->is_root()) {
        assert(node_stop != node_a);
        auto new_node_b = copy_until(node_stop, children_a[0], children_b[0]);
        if (new_node_b != nullptr) return new_node_b;
        else {
            children_b[0]->get_contents()->s = children_a[0]->get_contents()->s;
            return copy_until(node_stop, children_a[1], children_b[1]);
        }
    } else if (node_a->is_leaf()) {
        if (node_a == node_stop) return node_b;
        else {
            c_b->s = c_a->s;
            return nullptr;
        }
    } else {
        auto new_node_b = copy_until(node_stop, children_a[0], children_b[0]);
        if (new_node_b != nullptr) { // found in the left subtree
            c_b->l = c_a->l;
            return new_node_b;
        } else  {
            new_node_b = copy_until(node_stop, children_a[1], children_b[1]);
            if (new_node_b != nullptr) { // found in the right subtree
                children_b[0]->get_contents()->s = children_a[0]->get_contents()->s;
                c_b->l = c_a->l;
                return new_node_b;
            }
        }
        if (node_stop == node_a) return node_b; // this line is at the bottom because I need to copy s at the leaves
        else return nullptr; // not itself, not in its subtrees
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
            vector<int> left_s = left_child->get_contents()->s; // v
            vector<int> right_s = node_cur->get_contents()->s;  // u
            copy(right_s.begin(), right_s.end(), p_c->s.begin() + N_half); // copy u to (u+v, u)
            for (int i = 0; i < N_half; i++) { // copy u+v to (u+v, u)
                p_c->s[i] = left_s[i] ^ right_s[i];
            }
            return recursive_propagate_sums(p); 
        }
    }
}

void Decoder_RM_SCL::_decode_leaf_llr(Node<Contents_RM_SCL>* node_curr, int i, const double& pm)
{   // an active node is always a leaf (RM(m,0), RM(m,m), RM(1,1))
    // because otherwise can propogate back
    auto c = node_curr->get_c();
    int m = c->m, r = c->r;
    // cerr << "node m=" << m << ", r=" << r << ", pm=" << pm << endl;
    assert (!isnan(pm));
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

void Decoder_RM_SCL::_decode_m0(Node<Contents_RM_SCL>* node_curr, int i, const double& pm)
{   // repetition code, take both 0 and 1
    // calculate path metric for estimate to 0
    // calculate path metric for estimate to 1
    auto c = node_curr->get_c();
    int m = c->m, N = 1 << m;
    // cerr << "in m=" << m << ", r=0" << endl;
    vector<double> llr = c->l;
    double pm0 = pm, pm1 = pm;
    // double sum_llr = 0;
    for (int j = 0; j < N; j++) {
        // sum_llr += llr[i];
        pm0 = phi(pm0, llr[j], 0);
        pm1 = phi(pm1, llr[j], 1);
    }
    // cerr << "pm0=" << pm0 << ", pm1=" << pm1 << ", decision=" << (sum_llr>0.0 ? 0 : 1) << endl;
    metrics_vec.push_back(make_tuple(node_curr, i, vector(N, 0), pm0));
    metrics_vec.push_back(make_tuple(node_curr, i, vector(N, 1), pm1));
}

void Decoder_RM_SCL::_decode_11(Node<Contents_RM_SCL>* node_curr, int i, const double& pm)
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

void Decoder_RM_SCL::_decode_mm(Node<Contents_RM_SCL>* node_curr, int i, const double& pm)
{
    // if u = 0, flip to 1, PM = PM - log(1+exp(-lambda)) + log(1+exp(lambda)) = PM + lambda
    // if u = 1, flip to 0, PM = PM - lambda
    auto c = node_curr->get_c();
    int m = c->m, N = 1 << m;
    vector<double> llr = c->l; 
    vector<int> tmp(N);
    double pm_min = pm;
    double llr_min = abs(llr[0]), llr_snd_min = abs(llr[0]), llr_third_min = abs(llr[0]);
    int min_idx = 0, snd_min_idx = 0, third_min_idx = 0;
    double l, l_abs;
    // temporarily only take 1 instead of 4. TODO: take 4
    for (int i = 0; i < N; i++) {
        l = llr[i];
        if (l > 0) {
            tmp[i] = 0;
            pm_min = phi(pm_min, l, tmp[i]); // expect \Delta PM = 0
        } else {
            tmp[i] = 1;
            pm_min = phi(pm_min, l, tmp[i]); // expect \Delta PM = 0
        }
        l_abs = abs(l);
        if (l_abs < llr_min) {
            llr_third_min = llr_snd_min; third_min_idx = snd_min_idx;
            llr_snd_min = llr_min; snd_min_idx = min_idx;
            llr_min = l_abs; min_idx = i;
        } else if (l_abs < llr_snd_min) {
            llr_third_min = llr_snd_min; third_min_idx = snd_min_idx;
            llr_snd_min = l_abs; snd_min_idx = i;
        } else if (l_abs < llr_third_min) {
            llr_third_min = l_abs; third_min_idx = i;
        }
    }
    vector<int> snd_best(tmp); double pm_snd_min = pm_min;
    snd_best[min_idx] = !snd_best[min_idx]; pm_snd_min += llr_min;
    vector<int> third_best(tmp); double pm_third_min = pm_min;
    third_best[snd_min_idx] = !third_best[snd_min_idx]; pm_third_min += llr_snd_min;
    vector<int> fourth_best(tmp); double pm_fourth_min = pm_min;
    if (llr_min + llr_snd_min < llr_third_min) {
        fourth_best[min_idx] = !fourth_best[min_idx]; fourth_best[snd_min_idx] = !fourth_best[snd_min_idx];
        pm_fourth_min += llr_min + llr_snd_min;
    } else {
        fourth_best[third_min_idx] = !fourth_best[third_min_idx];
        pm_fourth_min += llr_third_min;
    }
    // assert (pm_min == pm);
    // cerr << "decode_mm " << pm_min << endl;
    metrics_vec.push_back(make_tuple(node_curr, i, tmp, pm_min));
    metrics_vec.push_back(make_tuple(node_curr, i, snd_best, pm_snd_min));
    metrics_vec.push_back(make_tuple(node_curr, i, third_best, pm_third_min));
    metrics_vec.push_back(make_tuple(node_curr, i, fourth_best, pm_fourth_min));
}
