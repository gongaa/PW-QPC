#include <utility>
#include <algorithm>

#include "Channel.hpp"
#include "Channel_BSC_q.hpp"
#include "Channel_depolarize_q.hpp"

using namespace std;

const string Channel_name = "Channel";
const string Channel_prefix = "chn";

Channel_c::Channel_c(const int N, const int seed = 0)
:N(N)
{
    mt19937.seed(seed);
}

Channel_q::Channel_q(const int N, const int seed = 0) 
:N(N)
{
    this->set_seed(seed);
}

void Channel_q::set_seed(const int seed) 
{
    this->event_generator->set_seed(seed);
}

void Channel_q::_add_noise(const float *CP, int *Y1_N, int *Y2_N, const size_t frame_id)
{
    this->event_generator->generate(Y1_N, Y2_N, this->N);
}

Channel_BSC::Channel_BSC(const int N, float p, const int seed = 0)
:Channel_c(N, seed), p(p)
{

}

void Channel_BSC::_add_noise(const float *CP, const int *X_N, const int *Y_N, const size_t frame_id)
{
    mt19937.randf_cc() <= this->p;
}