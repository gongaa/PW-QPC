#include "Decoder_polar.hpp"
#include <vector>
#include <set>
#include <cmath>
using namespace std;

Decoder_polar_SCL::Decoder_polar_SCL(const int& K, const int& N, const int& L, const vector<bool>& frozen_bits) 
                //   vector<function<float(const vector<float> &LLRS, const vector<int> &bits)>> lambdas)
: Decoder(K, N), metric_init(std::numeric_limits<float>::min()), L(L), stages((int)std::log(N), 0),
  frozen_bits(frozen_bits), LLRs(2), bits(1), lambdas(2)
{
    this->active_paths.insert(0);
    vector<uint32_t> sequence(2);
    // TODO: figure out what sequence does
    sequence[0] = 2;
    sequence[1] = 2;
    // ------------

    for (auto i = 0; i < L; i++) {
        this->polar_trees.push_back(Tree_metric<Contents_SCL>(sequence, metric_init));
        int max_depth_llrs = 1;
        this->recursive_allocate_nodes_contents(this->polar_trees[i].get_root(), this->N, max_depth_llrs);
        this->recursive_initialize_frozen_bits(this->polar_trees[i].get_root(), frozen_bits);
    }
    for (auto i = 0; i < L; i++) 
        leaves_array.push_back(this->polar_trees[i].get_leaves());

    this->lambdas[0] = [](const vector<float> &LLRs, const vector<int> &bits) -> float
    {
        auto sign = std::signbit((float)LLRs[0]) ^ std::signbit((float)LLRs[1]);
        auto abs0 = std::abs(LLRs[0]);
        auto abs1 = std::abs(LLRs[1]);
        auto min = std::min(abs0, abs1);
        return sign ? -min : min;
    };

    this->lambdas[1] = [](const vector<float> &LLRs, const vector<int> &bits) -> float
    {
        return ((bits[0] == 0) ? LLRs[0] : -LLRs[0]) + LLRs[1];
    };
}

void Decoder_polar_SCL::recursive_allocate_nodes_contents(Node<Contents_SCL>* node_curr, const int vector_size, int &max_depth_llrs)
{
	node_curr->set_contents(new Contents_SCL(vector_size));

	if (!node_curr->is_leaf())
	{
		const auto new_vector_size = vector_size / 2;
		for (auto c : node_curr->get_children())
			this->recursive_allocate_nodes_contents(c, new_vector_size, max_depth_llrs);
	}
	else
		node_curr->get_c()->max_depth_llrs = max_depth_llrs;

	max_depth_llrs = this->polar_trees[0].get_depth() - node_curr->get_depth() - 1;
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
        this->recursive_deallocate_nodes_contents(this->polar_trees[i].get_root());
}

void Decoder_polar_SCL::recursive_deallocate_nodes_contents(Node<Contents_SCL>* node_curr)
{
	for (auto c : node_curr->get_children())
		this->recursive_deallocate_nodes_contents(c); 

	delete node_curr->get_contents();
	node_curr->set_contents(nullptr);
}

void Decoder_polar_SCL::recursive_compute_llr(Node<Contents_SCL>* node_cur, int depth)
{
	auto node_father = node_cur->get_father();

	if (depth != 0)
		recursive_compute_llr(node_father, --depth);

	if (!node_cur->is_root())
	{
		const auto child = node_cur->get_child_id();
		const auto kern_size = (int)node_father->get_children().size(); // e.g. 2 in my case
		const auto size = (int)node_father->get_c()->l.size();
		const auto n_kernels = size / kern_size;

		for (auto k = 0; k < n_kernels; k++)
		{
			for (auto l = 0; l < kern_size; l++) LLRs[l] = node_father->get_c()->l[l * n_kernels + k];
			for (auto c = 0; c < child;     c++) bits[c] = node_father->get_children()[c]->get_c()->s[k];
			node_cur->get_c()->l[k] = lambdas[child](LLRs, bits);
		}
	}
}

void Decoder_polar_SCL::select_best_path(const size_t frame_id)
{
	int best_path = 0;
	if (active_paths.size() >= 1)
		best_path = *active_paths.begin();

	for (int path : active_paths)
		if(polar_trees[path].get_path_metric() < polar_trees[best_path].get_path_metric())
			best_path = path;

	active_paths.clear();
	active_paths.insert(best_path);
}

inline float phi(const float& mu, const float& lambda, const int& u)
{
	float new_mu;

	if (u == 0 && lambda < 0)
		new_mu = mu - lambda;
	else if (u != 0 && lambda > 0)
		new_mu = mu + lambda;
	else
		new_mu = mu;

	return new_mu;
}

void Decoder_polar_SCL::_load(const float *Y_N)
{
	for (auto path = 0; path < this->L; path++) {
		std::copy(Y_N, Y_N + this->N, this->polar_trees[path].get_root()->get_contents()->l.data());
		polar_trees[path].set_path_metric(metric_init);
	}

	// initialization
	active_paths.clear();
	active_paths.insert(0);
}

void Decoder_polar_SCL::_decode(const size_t frame_id)
{
	std::set<int> last_active_paths;
	int cur_path;

	// tuples to be sorted. <Path,estimated bit,metric>
	std::vector<std::tuple<int,int,float>> metrics_vec;

	// run through each leaf
	for (auto leaf_index = 0 ; leaf_index < this->N; leaf_index++)
	{
		// compute LLR for current leaf
		for (auto path : active_paths)
			this->recursive_compute_llr(leaves_array[path][leaf_index],
			                            leaves_array[path][leaf_index]->get_c()->max_depth_llrs);

		// if current leaf is a frozen bit
		if (leaves_array[0][leaf_index]->get_c()->is_frozen_bit)
		{
			auto min_phi = std::numeric_limits<float>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
				cur_leaf->get_c()->s[0] = 0;
				auto phi_cur = phi(polar_trees[path].get_path_metric(), cur_leaf->get_c()->l[0], 0);
				this->polar_trees[path].set_path_metric(phi_cur);
				min_phi = std::min<float>(min_phi, phi_cur);
			}

			// normalization
			for (auto path : active_paths)
				this->polar_trees[path].set_path_metric(this->polar_trees[path].get_path_metric() - min_phi);
		}
		else
		{
			// metrics vec used to store values of hypothetic path metrics
			metrics_vec.clear();
			auto min_phi = std::numeric_limits<float>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
				float phi0 = phi(polar_trees[path].get_path_metric(), cur_leaf->get_c()->l[0], 0);
				float phi1 = phi(polar_trees[path].get_path_metric(), cur_leaf->get_c()->l[0], 1);
				metrics_vec.push_back(std::make_tuple(path, 0, phi0));
				metrics_vec.push_back(std::make_tuple(path, 1, phi1));

				min_phi = std::min<float>(min_phi, phi0);
				min_phi = std::min<float>(min_phi, phi1);
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
					[](std::tuple<int,int,float> x, std::tuple<int,int,float> y){
						return std::get<2>(x) < std::get<2>(y);
					});

				// search in worst metrics. If a path is found twice, erase it
				for (auto it = metrics_vec.begin() + metrics_vec.size() / 2; it != metrics_vec.end(); ++it)
				{
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it + 1, metrics_vec.end(),
						[cur_path](std::tuple<int,int,float> x){
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
						[cur_path](std::tuple<int,int,float> x){
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
