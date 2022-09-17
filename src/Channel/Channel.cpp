#include <utility>
#include <algorithm>
#include <iostream>

#include "Channel/Channel.hpp"

using namespace std;

const string Channel_name = "Channel";
const string Channel_prefix = "chn";

Channel_c::Channel_c(const int N, const int seed)
: N(N)
{
    mt19937.seed(seed);
}

Channel_q::Channel_q(const int N, Event_generator_unitary* event_generator, const int seed) 
: event_generator(event_generator), N(N)
{
    this->set_seed(seed);
}

void Channel_q::set_seed(const int seed) 
{
    this->event_generator->set_seed(seed);
}

void Channel_q::add_noise(int *Y1_N, int *Y2_N, const size_t frame_id)
{
    this->event_generator->generate(Y1_N, Y2_N, this->N);
}

Channel_BSC::Channel_BSC(const int N, float p, const int seed)
:Channel_c(N, seed), p(p)
{
    if (p>1 || p<0)
        cerr << "Channel_BSC 0<=p<=1 violated" << endl; 
}

int Channel_BSC::add_noise(const int* X_N, int *Y_N, const size_t frame_id)
{
    // X_N[i] = ((mt19937.randf_cc() >= this->p) != X_N[i]) ? 0 : 1;
    int num_flips = 0;
    for(int i = 0; i < N; i++) {
        if (mt19937.randf_cc() <= this->p) {
            num_flips++;
            Y_N[i] = !X_N[i]; 
        } else Y_N[i] = X_N[i];
    }
    cerr << "num_flips " << num_flips << endl;
    return num_flips;
}