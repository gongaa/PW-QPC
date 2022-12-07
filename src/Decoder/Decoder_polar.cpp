#include "Decoder/Decoder_polar.hpp"
#include <vector>
#include <set>
#include <iostream>
#include <cmath>
using namespace std;

Decoder_polar_SCL::Decoder_polar_SCL(const int& K, const int& N, const int& L, const vector<bool>& frozen_bits) 
: Decoder(K, N), L(L), frozen_bits(frozen_bits)
{
    this->active_paths.insert(0);
	int n = log2(N);
	int max_depth_llrs = n - 1;
    for (auto i = 0; i < L; i++) {
		auto new_tree = new Tree_metric<Contents_SCL>(n, 0);
		this->polar_trees.push_back(new_tree);
		this->recursive_allocate_nodes_contents(new_tree->get_root(), this->N, max_depth_llrs);
		this->recursive_initialize_frozen_bits(new_tree->get_root(), frozen_bits);
    }
    for (auto i = 0; i < L; i++) 
        leaves_array.push_back(this->polar_trees[i]->get_leaves());
}

void Decoder_polar_SCL::recursive_allocate_nodes_contents(Node<Contents_SCL>* node_curr, const int vector_size, int& max_depth_llrs)
{
	node_curr->set_contents(new Contents_SCL(vector_size));
	if (!node_curr->is_leaf())
	{
		const int new_vector_size = vector_size / 2;
		for (auto c : node_curr->get_children())
			this->recursive_allocate_nodes_contents(c, new_vector_size, max_depth_llrs);
	} else node_curr->get_c()->max_depth_llrs = max_depth_llrs;
	max_depth_llrs = this->polar_trees[0]->get_depth() - node_curr->get_depth();
	// the max_depth_llrs is the (path length -1) from the current leaf
	// to the common ancestor of itself and the previous immediate leaf.
	// max_depth_llrs for even leaf is always 0
}

void Decoder_polar_SCL::recursive_initialize_frozen_bits(const Node<Contents_SCL>* node_curr, const std::vector<bool>& frozen_bits)
{
	if (!node_curr->is_leaf()) 
	{
		for (auto c : node_curr->get_children())
			this->recursive_initialize_frozen_bits(c, frozen_bits); 
	}
	else 
		node_curr->get_contents()->is_frozen_bit = frozen_bits[node_curr->get_lane_id()];
}

Decoder_polar_SCL::~Decoder_polar_SCL()
{
    for (auto i = 0; i < L; i++) 
        this->recursive_deallocate_nodes_contents(this->polar_trees[i]->get_root());
}

void Decoder_polar_SCL::recursive_deallocate_nodes_contents(Node<Contents_SCL>* node_curr)
{
	for (auto c : node_curr->get_children())
		this->recursive_deallocate_nodes_contents(c); 

	delete node_curr->get_contents();
	node_curr->set_contents(nullptr);
}

void Decoder_polar_SCL::copy_codeword_list(vector<vector<int>>& c_list, vector<double>& pm_list)
{
	int i = 0;
	for (int path : this->active_paths) {
		const int* X = polar_trees[i]->get_root()->get_c()->s.data();
		copy(X, X+N, c_list[i].data());
		pm_list[i] = polar_trees[path]->get_path_metric();
		++i;
	}
}

void Decoder_polar_SCL::recursive_compute_llr(Node<Contents_SCL>* node_cur, int depth)
{
	auto node_father = node_cur->get_father();

	if (depth != 0)
		recursive_compute_llr(node_father, --depth);

	if (!node_cur->is_root()) {
		int child_id = node_cur->get_child_id();
		auto p = node_cur->get_father();
		auto p_c = p->get_contents();
		int size = p_c->l.size(), size_half = size / 2;
		if (child_id == 1) {
			// cerr << "update llr for right child at depth " << node_cur->get_depth() << " , lane_id " << node_cur->get_lane_id() << endl;
			auto left_child = p->get_children()[0];
			f_minus(p_c->l.data(), p_c->l.data() + size_half, left_child->get_contents()->s.data(), size_half, node_cur->get_contents()->l.data());
		} else if (child_id == 0) {
			// cerr << "update llr for left child at depth " << node_cur->get_depth() << " , lane_id " << node_cur->get_lane_id() << endl;
			f_plus(p_c->l.data(), p_c->l.data() + size_half, size_half, node_cur->get_contents()->l.data());
		}
	}
}

void Decoder_polar_SCL::partition(vector<int>& info_indices, vector<vector<int>>& par, vector<vector<int>>& flips, vector<int>& noisy_codeword, int& best_path_class_idx)
{
	// partition into equivalence classes according to info bits
	int size = info_indices.size();
	vector<bool> path_info(size);
	auto& best_path_leaves = leaves_array[this->best_path];
	for (int i = 0; i < size; i++)
		path_info[i] = best_path_leaves[info_indices[i]]->get_c()->s[0];
	best_path_class_idx = binary2decimal(path_info, size);

	int class_idx;
	int wt;
	for (int path : active_paths) {
		auto& path_leaves = leaves_array[path];
		for (int i = 0; i < size; i++)
			path_info[i] = path_leaves[info_indices[i]]->get_c()->s[0];
		class_idx = binary2decimal(path_info, size);
		par[class_idx].push_back(path);
		auto& codeword = polar_trees[path]->get_root()->get_c()->s;
		wt = count_flip(N, noisy_codeword, codeword);
		flips[class_idx].push_back(wt);
	}
}

void Decoder_polar_SCL::select_best_path(const size_t frame_id)
{   // select the best one, not the best L ones.
	int best_path = 0;
	if (active_paths.size() >= 1)
		best_path = *active_paths.begin();

	for (int path : active_paths)
		if(polar_trees[path]->get_path_metric() < polar_trees[best_path]->get_path_metric())
			best_path = path;

	this->best_path = best_path;
}

void Decoder_polar_SCL::_load(const double *Y_N)
{
	for (int path = 0; path < this->L; path++) {
		std::copy(Y_N, Y_N + this->N, this->polar_trees[path]->get_root()->get_contents()->l.data());
		// polar_trees[path]->set_path_metric(numeric_limits<double>::min());
		polar_trees[path]->set_path_metric(0);
	}

	// initialization
	active_paths.clear();
	active_paths.insert(0);
}

void Decoder_polar_SCL::_decode(const size_t frame_id)
{
	std::set<int> last_active_paths;
	int cur_path;
	int depth = log2(N);

	// tuples to be sorted. <Path,estimated bit,metric>
	std::vector<std::tuple<int,int,double>> metrics_vec;

	// run through each leaf
	for (auto leaf_index = 0; leaf_index < this->N; leaf_index++)
	{
		// compute LLR for current leaf
		for (auto path : active_paths) 
			this->recursive_compute_llr(leaves_array[path][leaf_index], 
										leaves_array[path][leaf_index]->get_c()->max_depth_llrs);
		// only need to compute llr starting from the common ancestor of itself and the previous leaf

		// if current leaf is a frozen bit
		if (leaves_array[0][leaf_index]->get_c()->is_frozen_bit) {
		    // penalize if the prediction for frozen bit is wrong, frozen value is 0
			// frozen bit should not be set it to CRC of info bits
			// the CRC checksums should be placed at the best channels
			auto min_phi = std::numeric_limits<double>::max();
			for (auto path : active_paths) {
				auto cur_leaf = leaves_array[path][leaf_index];
				cur_leaf->get_c()->s[0] = 0;
				auto phi_cur = phi(polar_trees[path]->get_path_metric(), cur_leaf->get_c()->l[0], 0);
				this->polar_trees[path]->set_path_metric(phi_cur);
				min_phi = std::min<double>(min_phi, phi_cur);
			}

			// normalization
			for (auto path : active_paths)
				this->polar_trees[path]->set_path_metric(this->polar_trees[path]->get_path_metric() - min_phi);

		} else {
			// metrics vec used to store values of hypothetic path metrics
			metrics_vec.clear();
			auto min_phi = std::numeric_limits<double>::max();
			for (auto path : active_paths) {
				auto cur_leaf = leaves_array[path][leaf_index];
				double phi0 = phi(polar_trees[path]->get_path_metric(), cur_leaf->get_c()->l[0], 0);
				double phi1 = phi(polar_trees[path]->get_path_metric(), cur_leaf->get_c()->l[0], 1);
				metrics_vec.push_back(std::make_tuple(path, 0, phi0));
				metrics_vec.push_back(std::make_tuple(path, 1, phi1));

				min_phi = std::min<double>(min_phi, phi0);
				min_phi = std::min<double>(min_phi, phi1);
			}

			
			for (auto& vec : metrics_vec) // normalization
				std::get<2>(vec) -= min_phi;

			if (active_paths.size() <= (unsigned)(L / 2)) {
				last_active_paths = active_paths;
				for (auto path : last_active_paths)
					this->duplicate_path(path, leaf_index, leaves_array);
			} else {
				// sort hypothetic path metrics in increasing order
				std::sort(metrics_vec.begin(), metrics_vec.end(),
					[](std::tuple<int,int,double> x, std::tuple<int,int,double> y){
						return std::get<2>(x) < std::get<2>(y);
					});

				// search in worst metrics. If a path is found twice, erase it
				for (auto it = metrics_vec.begin() + metrics_vec.size() / 2; it != metrics_vec.end(); ++it) {
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it + 1, metrics_vec.end(),
						[cur_path](std::tuple<int,int,double> x){ return std::get<0>(x) == cur_path; });

					if (it_double != metrics_vec.end())
						active_paths.erase(std::get<0>(*it));
				}

				// remove worst metrics from list
				metrics_vec.resize(metrics_vec.size() / 2);

				for (auto it = metrics_vec.begin(); it != metrics_vec.end(); ++it) {
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it + 1, metrics_vec.end(),
						[cur_path](std::tuple<int,int,double> x){ return std::get<0>(x) == cur_path; });

					if (it_double != metrics_vec.end()) { // duplicate
						metrics_vec.erase(it_double);
						duplicate_path(std::get<0>(*it), leaf_index, leaves_array);
					} else { // choose
						leaves_array[std::get<0>(*it)][leaf_index]->get_c()->s[0] = std::get<1>(*it);
						polar_trees[std::get<0>(*it)]->set_path_metric(std::get<2>(*it));
					}
				}
			}
		}
		// right node keeps propagating sums, until itself becomes a left node
		for (auto path : active_paths)
			this->recursive_propagate_sums(leaves_array[path][leaf_index]);
	}

	this->select_best_path(frame_id);
}

void Decoder_polar_SCL::recursive_propagate_sums(const Node<Contents_SCL>* node_cur)
{
	auto children = node_cur->get_children();

	if (children.size() > 0) { // not leaf
		auto n_c = node_cur->get_contents();
		int size = n_c->s.size(), size_half = size / 2;
		auto left_s = children[0]->get_contents()->s;
		auto right_s = children[1]->get_contents()->s;
		copy(right_s.begin(), right_s.end(), n_c->s.begin() + size_half);
		for (int i = 0; i < size_half; i++)
			n_c->s[i] = left_s[i] ^ right_s[i];
	}
	if (!node_cur->is_root() && (node_cur->get_child_id() == 1)) // is the right child
		this->recursive_propagate_sums(node_cur->get_father());
	// else cerr << "recursive propagate sums ends at node at depth " << node_cur->get_depth() << " , lane_id " << node_cur->get_lane_id() << endl;
}

void Decoder_polar_SCL::duplicate_path(int path, int leaf_index, vector<vector<Node<Contents_SCL>*>> leaves_array)
{
	vector<Node<Contents_SCL>*> path_leaves, newpath_leaves;
	int newpath = 0;
	while (active_paths.find(newpath++) != active_paths.end()){};
	newpath--;
	active_paths.insert(newpath);
	path_leaves = leaves_array[path];
	newpath_leaves = leaves_array[newpath];
	for (auto i = 0; i < leaf_index; i++)
		newpath_leaves[i]->get_c()->s = path_leaves[i]->get_c()->s;

	// the cleverer way
	recursive_duplicate_tree_sums(leaves_array[path][leaf_index], leaves_array[newpath][leaf_index], nullptr);
	if (leaf_index < this->N - 1)
		recursive_duplicate_tree_llr(leaves_array[path][leaf_index + 1], leaves_array[newpath][leaf_index + 1]);
	// do not need to copy the whole tree, only copy the necessary part to compute llr for all the future leaf nodes

	leaves_array[newpath][leaf_index]->get_c()->s[0] = 1;
	polar_trees[newpath]->set_path_metric(phi(polar_trees[path]->get_path_metric(),
	                                                     leaves_array[path][leaf_index]->get_c()->l[0], 1));

	leaves_array[path][leaf_index]->get_c()->s[0] = 0;
	polar_trees[path]->set_path_metric(phi(polar_trees[path]->get_path_metric(),
	                                                  leaves_array[path][leaf_index]->get_c()->l[0], 0));
}

void Decoder_polar_SCL::recursive_duplicate_tree_llr(Node<Contents_SCL>* node_a, Node<Contents_SCL>* node_b)
{
	node_b->get_c()->l = node_a->get_c()->l;

	if (!node_a->get_father()->is_root())
		this->recursive_duplicate_tree_llr(node_a->get_father(), node_b->get_father());
}

void Decoder_polar_SCL::recursive_duplicate_tree_sums(Node<Contents_SCL>* node_a, Node<Contents_SCL>* node_b, Node<Contents_SCL>* node_caller)
{
	if (!node_a->is_leaf()) { // if called by its right child
		auto left_child = (node_a->get_children())[0];
		if (left_child != node_caller) {
			node_b->get_children()[0]->get_c()->s = left_child->get_c()->s;
		}
	}
	if (!node_a->is_root())
		this->recursive_duplicate_tree_sums(node_a->get_father(), node_b->get_father(), node_a);
}

void Decoder_polar_SCL::_store(int *V) const
{
    auto *root = this->polar_trees[this->best_path]->get_root();
    std::copy(root->get_c()->s.begin(), root->get_c()->s.begin() + this->N, V);
}

int Decoder_polar_SCL::decode(const double *Y_N, int *V_K, const size_t frame_id)
{
    this->_load(Y_N);
    this->_decode(frame_id);
    this->_store(V_K);
	return 0;
}

void Decoder_polar_SCL::decode_SC(const double *Y_N, int *V_K, const size_t frame_id) 
{
	this->_load(Y_N);
	this->_decode_SC(polar_trees[0]->get_root());
	this->_store(V_K);
}

void Decoder_polar_SCL::_decode_SC(Node<Contents_SCL>* node_cur) 
{
	auto n_c = node_cur->get_c();
	if (node_cur->is_leaf()) {
		// cerr << "node " << node_cur->get_lane_id() << " is frozen bit " << n_c->is_frozen_bit << endl;
		n_c->s[0] = (n_c->is_frozen_bit || n_c->l[0] >= 0) ? 0 : 1;
		return;
	}
	auto left_child = node_cur->get_children()[0];
	auto right_child = node_cur->get_children()[1];
	int size = n_c->l.size(), size_half = size / 2;
	// cerr << "update llr for left child at depth " << node_cur->get_depth() << " , lane_id " << node_cur->get_lane_id() << endl;
	f_plus(n_c->l.data(), n_c->l.data() + size_half, size_half, left_child->get_c()->l.data());
	_decode_SC(left_child);
	// cerr << "update llr for right child at depth " << node_cur->get_depth() << " , lane_id " << node_cur->get_lane_id() << endl;
	f_minus(n_c->l.data(), n_c->l.data() + size_half, left_child->get_c()->s.data(), size_half, right_child->get_c()->l.data());
	_decode_SC(right_child);
	auto left_s = left_child->get_c()->s;
	auto right_s = right_child->get_c()->s;
	copy(right_s.begin(), right_s.end(), n_c->s.begin() + size_half);
	for (int i = 0; i < size_half; i++)
		n_c->s[i] = left_s[i] ^ right_s[i];
}

void Decoder_polar_SCL::get_llr_for_frozen_bits(double *Y_N)
{
	int i = 0;
	for (auto n : this->leaves_array[0]) 
		Y_N[i++] = n->get_c()->l[0];
}