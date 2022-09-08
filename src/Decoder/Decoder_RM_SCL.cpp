/*
for every LLR bit of (u+v) + u
greedily choose whether it is 0 or 1 (by comparing path metric)

then partition and cut the list
finalize the remaining branch

(0, m)x2, (m, m)x4.
(m, m) add to path likelihoods.
Optimized branching.
List as a tree.

what I need
- a table of list sizes at the leaf nodes
  and a counter for it
- slist and plist



*/
#include "Decoder/Decoder_RM_SCL.hpp"
using namespace std;

Decoder_RM_SCL::Decoder_RM_SCL(const int& m, const int& r, const int& L, const vector<bool>& frozen_bits) 
                //   vector<function<double(const vector<double> &LLRS, const vector<int> &bits)>> lambdas)
: Decoder(Encoder_RM::calculate_K(m,r), 1 << m), m(m), r(r), metric_init(std::numeric_limits<double>::min()), L(L), stages((int)std::log(N), 0),
  frozen_bits(frozen_bits), LLRs(2), bits(1), idx(2), u(2), Ke(4)
{
    this->Ke[0] = 1; this->Ke[1] = 1; this->Ke[2] = 0; this->Ke[3] = 1;
    this->active_paths.insert(0);

    for (auto i = 0; i < L; i++) {
        this->rm_trees.push_back(RM_Tree_metric<Contents_RM_SCL>(m, r, metric_init));
        this->recursive_allocate_nodes_contents(this->rm_trees[i].get_root(), this->N, this->m, this->r);
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

	if (!node_curr->is_leaf())
	{
		const auto new_vector_size = vector_size / 2;
		for (auto c : node_curr->get_children())
			this->recursive_allocate_nodes_contents(c, new_vector_size, m-1, r--);
	}
}

Decoder_RM_SCL::~Decoder_RM_SCL()
{
    for (auto i = 0; i < L; i++)
        this->recursive_deallocate_nodes_contents(this->rm_trees[i].get_root());
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
		if(rm_trees[path].get_path_metric() < rm_trees[best_path].get_path_metric())
			best_path = path;

	active_paths.clear();
	active_paths.insert(best_path);
}

void Decoder_RM_SCL::recursive_compute_llr(Node<Contents_RM_SCL>* node_curr)
{
    if (node_curr->is_leaf()) {
        active_nodes.push_back(node_curr);
    }
    Contents_RM_SCL* c = node_curr->get_c();
    int m = c->m, r = c->r, N = 1 << m, N_half = N / 2;
    vector<double> llr = c->l;
    vector<Node<Contents_RM_SCL>*> children = node_curr->get_children();
    auto left_child = children[0], right_child = children[1];
    Decoder::f_plus(llr.data(), llr.data() + N_half, N_half, left_child->get_c()->l.data());
    left_child->get_c()->status = USED;
    recursive_compute_llr(left_child);
}

void Decoder_RM_SCL::_decode(const size_t frame_id)
{
    // while (active_nodes.size() != 0)
    for (int i : active_paths) {
        this->recursive_compute_llr(rm_trees[i].get_root());
    }
    for (auto active_node : active_nodes) {
        this->_decode_llr(active_node);  // split nodes
    }
    // TODO: clean up both the active_paths and active_nodes array
    //       through partitioning
    // TODO: copy and delete some rm_tree, 
    // update s for the active nodes,
    // change active nodes' status from used to determined
    metrics_vec.clear();
}

void Decoder_RM_SCL::_decode_llr(Node<Contents_RM_SCL>* node_curr)
{
    auto c = node_curr->get_c();
    int m = c->m, r = c->r;
    if (m == r) {
        if (active_paths.size() >= L) { // I need a talbe of list sizes at the border nodes 
            _decode_mm(node_curr, true);
        } else if (r == 1) {
            _decode_11(node_curr);
        } else {
            _decode_mm(node_curr, false);
        }
    }

    if (r == 0) {
        if (active_paths.size() >= L) { // again
            _decode_m0(node_curr, true);
        } else {
            _decode_m0(node_curr, false);
        }
    }

    // split nodes: change status from used to determined
}

void Decoder_RM_SCL::_decode_m0(Node<Contents_RM_SCL>* node_curr, bool skip = false)
{   // repetition code, take both 0 and 1
    // calculate path metric for estimate to 0
    // calculate path metric for estimate to 1
    auto c = node_curr->get_c();
    int m = c->m, N = 1 << m;
    vector<double> llr = c->l;
    double pm0 = 0.0, pm1 = 0.0;
    for (int i = 0; i < N; i++) {
        pm0 = phi(pm0, llr[i], 0);
        pm1 = phi(pm1, llr[i], 1);
    }
    metrics_vec.push_back(make_tuple(node_curr, vector(N, 0), pm0));
    metrics_vec.push_back(make_tuple(node_curr, vector(N, 1), pm1));
}

void Decoder_RM_SCL::_decode_11(Node<Contents_RM_SCL>* node_curr)
{

}

void Decoder_RM_SCL::_decode_mm(Node<Contents_RM_SCL>* node_curr, bool skip = false)
{
    int m = node_curr->get_c()->m, N = 1 << m;
}

void Decoder_RM_SCL::_decode(const size_t frame_id)
{
	std::set<int> last_active_paths;
	int cur_path;

	// tuples to be sorted. <Path,estimated bit,metric>
	std::vector<std::tuple<int,int,double>> metrics_vec;

	// run through each leaf
	for (auto leaf_index = 0 ; leaf_index < this->N; leaf_index++)
	{
		// compute LLR for current leaf
		for (auto path : active_paths)
			this->recursive_compute_llr(leaves_array[path][leaf_index],
			                            leaves_array[path][leaf_index]->get_c()->max_depth_llrs);

		// if current leaf is a frozen bit
		if (leaves_array[0][leaf_index]->get_c()->is_frozen_bit)
		{   // penalize if the prediction for frozen bit is wrong, TODO: why defalut frozen value is 0?
			auto min_phi = std::numeric_limits<double>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
				cur_leaf->get_c()->s[0] = 0; // TODO: shouldn't s (u_i) be set to the frozen value? why zero?
				auto phi_cur = phi(polar_trees[path].get_path_metric(), cur_leaf->get_c()->l[0], 0);
				this->polar_trees[path].set_path_metric(phi_cur);
				min_phi = std::min<double>(min_phi, phi_cur);
			}

			// normalization
			for (auto path : active_paths)
				this->polar_trees[path].set_path_metric(this->polar_trees[path].get_path_metric() - min_phi);
		}
		else
		{
			// metrics vec used to store values of hypothetic path metrics
			metrics_vec.clear();
			auto min_phi = std::numeric_limits<double>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
				double phi0 = phi(polar_trees[path].get_path_metric(), cur_leaf->get_c()->l[0], 0);
				double phi1 = phi(polar_trees[path].get_path_metric(), cur_leaf->get_c()->l[0], 1);
				metrics_vec.push_back(std::make_tuple(path, 0, phi0));
				metrics_vec.push_back(std::make_tuple(path, 1, phi1));

				min_phi = std::min<double>(min_phi, phi0);
				min_phi = std::min<double>(min_phi, phi1);
			}

			// normalization
			for (auto vec : metrics_vec)
				std::get<2>(vec) -= min_phi;

			if (active_paths.size() <= (unsigned)(L / 2))
			{
				last_active_paths = active_paths;
				for (auto path : last_active_paths)
					this->duplicate_path(path, leaf_index, leaves_array);
			}
			else
			{
				// sort hypothetic metrics
				std::sort(metrics_vec.begin(), metrics_vec.end(),
					[](std::tuple<int,int,double> x, std::tuple<int,int,double> y){
						return std::get<2>(x) < std::get<2>(y);
					});

				// search in worst metrics. If a path is found twice, erase it
				for (auto it = metrics_vec.begin() + metrics_vec.size() / 2; it != metrics_vec.end(); ++it)
				{
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it + 1, metrics_vec.end(),
						[cur_path](std::tuple<int,int,double> x){
							return std::get<0>(x) == cur_path;
						});

					if (it_double != metrics_vec.end())
						active_paths.erase(std::get<0>(*it));
				}

				// remove worst metrics from list
				metrics_vec.resize(metrics_vec.size() / 2);

				for (auto it = metrics_vec.begin(); it != metrics_vec.end(); ++it)
				{
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it +1, metrics_vec.end(),
						[cur_path](std::tuple<int,int,double> x){
							return std::get<0>(x) == cur_path;
						});

					if (it_double != metrics_vec.end())
					{
						// duplicate
						metrics_vec.erase(it_double);
						duplicate_path(std::get<0>(*it), leaf_index, leaves_array);
					}
					else
					{
						// choose
						leaves_array[std::get<0>(*it)][leaf_index]->get_c()->s[0] = std::get<1>(*it);
						polar_trees[std::get<0>(*it)].set_path_metric(std::get<2>(*it));
					}
				}
			}
		}

		// propagate sums
		for (auto path : active_paths)
			this->recursive_propagate_sums(leaves_array[path][leaf_index]);
	}

	this->select_best_path(frame_id);
}