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
using namespace std;

Decoder_RM_SCL::Decoder_RM_SCL(const int& m, const int& r, const int& L) 
: Decoder(Encoder_RM::calculate_K(m,r), 1 << m), m(m), r(r), metric_init(std::numeric_limits<double>::min()), L(L)
{
    for (auto i = 0; i < L; i++) {
        auto new_tree = new RM_Tree_metric<Contents_RM_SCL>(m, r, metric_init);
        this->recursive_allocate_nodes_contents(new_tree->get_root(), this->N, this->m, this->r);
        this->rm_trees.push_back(new_tree);
    }
    // TODO: only one active_path
    // TODO: recursive compute LLR for it
    // TODO: the active_nodes_array is the RM(m-r, 0) node of this active path
    // for (auto i = 0; i < L; i++) 
    //     leaves_array.push_back(this->rm_trees[i].get_leaves());
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

void Decoder_RM_SCL::recursive_compute_llr(Node<Contents_RM_SCL>* node_curr)
{
    if (node_curr->is_leaf()) {
        active_nodes.push_back(node_curr);
        return;
    }
    Contents_RM_SCL* c = node_curr->get_c();
    int m = c->m, N = 1 << m, N_half = N / 2;
    vector<double> llr = c->l;
    vector<Node<Contents_RM_SCL>*> children = node_curr->get_children();
    auto left_child = children[0], right_child = children[1];
    Decoder::f_plus(llr.data(), llr.data() + N_half, N_half, left_child->get_c()->l.data());
    left_child->get_c()->status = USED;
    recursive_compute_llr(left_child);
}


int Decoder_RM_SCL::decode(const double *Y_N, int* X_N, const size_t frame_id)
{
    // initialization: only one active path, 
    // and it's active node is RM(m-r, 0), obtained by recursive_compute_llr
    this->active_paths.clear();
    this->active_paths.insert(0);
    for (auto i = 0; i < L; i++) {
        rm_trees[i]->get_root()->get_c()->l = vector<double>(Y_N, Y_N + this->N);
        this->recursive_compute_llr(rm_trees[i]->get_root());
    }
    // the active node of each rm_tree is in the same position, but with different content
    while (active_nodes.size() != 0) {
        vector<Node<Contents_RM_SCL>*> active_nodes_tmp(active_nodes);
        active_nodes.clear();
        for (int i : active_paths) {
            _decode_leaf_llr(active_nodes_tmp[i], i, rm_trees[i]->path_metric); // split nodes
        }
        // an active node is always a left child, except RM(1,1)
        // update and decode parent's right child
        // TODO: clean up both the active_paths and active_nodes array
        //       through partitioning
        // TODO: copy and delete some rm_tree, 
        // update s for the active nodes,
        // change active nodes' status from used to determined
        metrics_vec.clear();



        partition_copy_delete();

		// propagate sums
		for (int i : active_paths)
			this->recursive_propagate_sums(active_nodes[i]);

        for (int i : active_paths) {
            this->recursive_compute_llr(active_nodes[i]); // for the right child
            // cannot start from the root,
            // because this function always only looks at the left node
        }
    }


	this->select_best_path(frame_id);
    return 0;
}

void Decoder_RM_SCL::partition_copy_delete()
{
    // sort hypothetic metrics
    std::sort(metrics_vec.begin(), metrics_vec.end(),
        [](std::tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> x, std::tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> y){
            return std::get<3>(x) < std::get<3>(y);
        });
    // sorted so that PM is in increasing order
    // can there be duplicate path?
    // search in worst metrics. If a path is found twice, erase it
    for (auto it = metrics_vec.begin() + metrics_vec.size() / 2; it != metrics_vec.end(); ++it)
    {
        auto cur_path_node = get<0>(*it);
        int cur_path_idx = get<1>(*it);

        auto it_double = find_if(it + 1, metrics_vec.end(),
            [cur_path_node, cur_path_idx](tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> x){
                return (get<0>(x) == cur_path_node && get<1>(x) == cur_path_idx);
            });

        if (it_double != metrics_vec.end()) {
            cerr << "found duplicated path" << endl;
            active_paths.erase(std::get<1>(*it));
        }
    }

    // remove worst metrics from list
    if (metrics_vec.size() > L) {
        metrics_vec.resize(L);
    } else {
        // TODO: only do copy
    }

    for (auto it = metrics_vec.begin(); it != metrics_vec.end(); ++it)
    {
        auto cur_path_node = get<0>(*it);
        int cur_path_idx = get<1>(*it);

        auto it_double = find_if(it + 1, metrics_vec.end(),
            [cur_path_node, cur_path_idx](tuple<Node<Contents_RM_SCL>*, int, vector<int>, double> x){
                return (get<0>(x) == cur_path_node && get<1>(x) == cur_path_idx);
            });

        if (it_double != metrics_vec.end())
        {
            // duplicate
            metrics_vec.erase(it_double);
            duplicate_path(get<1>(*it), get<0>(*it));
        }
        else
        {
            // choose
            get<0>(*it)->get_contents()->s = get<2>(*it);
            rm_trees[get<1>(*it)]->set_path_metric(get<3>(*it));
        }
    }
}

bool Decoder_RM_SCL::copy_until(Node<Contents_RM_SCL>* node_stop, Node<Contents_RM_SCL>* node_a, Node<Contents_RM_SCL>* node_b)
{   // copy tree 1 (to which node_a and node_stop belong) into copy tree 2 (to which node_b belongs) until node_stop
    // return: whether reached or not
    auto c = node_a->get_contents();
    cerr << "in node with m=" << c->m << ", r=" << c->r << endl;
    auto children_a = node_a->get_children();
    auto children_b = node_b->get_children();
    if (node_a->is_root()) {
        assert(node_stop != node_a);
        if (copy_until(node_stop, children_a[0], children_b[0])) return true;
        else return copy_until(node_stop, children_a[1], children_b[1]);
    } else if (node_a->is_leaf()) {
        // copy and return;
        if (node_a == node_stop) {
            cerr << "reached m=" << c->m << ", r=" << c->r << endl;
            return true;
        }
        else return false;
    } else {
        if (copy_until(node_stop, children_a[0], children_b[0])) return true;
        else if (copy_until(node_stop, children_a[1], children_b[1])) return true;
        else {
            // copy self
            if (node_stop == node_a) {
                cerr << "reached m=" << c->m << ", r=" << c->r << endl;
                return true; 
            } else return false;
        }

    }
}

void Decoder_RM_SCL::test_copy_until()
{
    auto r1 = rm_trees[0]->get_root(), r2 = rm_trees[1]->get_root();
    cerr << "root children size " << r1->get_children().size() << endl;
    auto c = r1->get_contents();
    cerr << "root m=" << c->m << ", r=" << c->r << endl;
    // auto node_stop = (((r1->get_children()[0])->get_children()[0])->get_children()[0])->get_children()[1];
    auto node_stop = r1;
    while (!node_stop->is_leaf())
        node_stop = node_stop->get_children()[0];
    cerr << "in test_copy_until" << endl;
    copy_until(node_stop->get_father(), r1, r2);
}

void Decoder_RM_SCL::recursive_propagate_sums(const Node<Contents_RM_SCL>* node_cur)
{

}

void Decoder_RM_SCL::_decode_leaf_llr(Node<Contents_RM_SCL>* node_curr, int i, double& pm)
{   // an active node is always a leaf (RM(m,0), RM(m,m), RM(1,1))
    // because otherwise can propogate back
    auto c = node_curr->get_c();
    int m = c->m, r = c->r, N = 1 << m;
    if (m == r) {
        if (active_paths.size() >= L) { // I need a talbe of list sizes at the border nodes 
            _decode_mm(node_curr, i, pm);
        } else if (r == 1) {
            _decode_11(node_curr, i, pm);
        } else {
            _decode_mm(node_curr, i, pm);
        }
    } else if (r == 0) {
        if (active_paths.size() >= L) { // again
            _decode_m0(node_curr, i, pm);
        } else {
            _decode_m0(node_curr, i, pm);
        }
    } 
    node_curr->get_c()->status = DETERMINED;
    auto llr_v = c->l;  // llr for (u+v) + u = v
    auto bits = c->s; // bit estimate for \hat{v}
    auto p = node_curr->get_father();
    auto llr_uv = p->get_c()->l; // llr for u+v, size N
    auto right_child = p->get_children()[1];
    // llr for u = (u+v) + \hat{v}
    f_minus(llr_uv.data(), llr_v.data(), bits.data(), N, right_child->get_c()->l.data());
    active_nodes.push_back(right_child);
    // update parent's right child
    // split nodes: change status from used to determined
}

void Decoder_RM_SCL::_decode_m0(Node<Contents_RM_SCL>* node_curr, int i, double& pm)
{   // repetition code, take both 0 and 1
    // calculate path metric for estimate to 0
    // calculate path metric for estimate to 1
    auto c = node_curr->get_c();
    int m = c->m, N = 1 << m;
    vector<double> llr = c->l;
    double pm0 = pm, pm1 = pm;
    for (int i = 0; i < N; i++) {
        pm0 = phi(pm0, llr[i], 0);
        pm1 = phi(pm1, llr[i], 1);
    }
    metrics_vec.push_back(make_tuple(node_curr, i, vector(N, 0), pm0));
    metrics_vec.push_back(make_tuple(node_curr, i, vector(N, 1), pm1));
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