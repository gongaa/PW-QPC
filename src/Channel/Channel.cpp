#include <utility>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "Channel/Channel.hpp"

using namespace std;

const string Channel_name = "Channel";
const string Channel_prefix = "chn";

Channel_c::Channel_c(const int N, const int seed)
: N(N)
{
    mt19937.seed(seed);
}

Channel_q::Channel_q(const int N, Event_generator_unitary* event_generator) 
: event_generator(event_generator), N(N)
{
}

void Channel_q::set_seed(const int seed) 
{
    this->event_generator->set_seed(seed);
}

void Channel_q::add_noise(int *Y1_N, int *Y2_N, const size_t frame_id)
{
    this->event_generator->generate(Y1_N, Y2_N, this->N);
}
Channel_AWGN::Channel_AWGN(const int N, double sigma, double design_sigma, const int seed)
: Channel_c(N, seed), sigma(sigma), design_sigma(design_sigma)
{

}

int Channel_AWGN::add_noise(const int *X_N, double *Y_N, const size_t frame_id)
{
    // bpsk: 0 -> -1, 1 -> +1
    int i = 0;
    double v1, v2, r;
    while (i < N) {
        do {
        v1 = 2.0 * mt19937.randf_cc() - 1.0;
        v2 = 2.0 * mt19937.randf_cc() - 1.0;
        r = v1 * v1 + v2 * v2;
        } while (r >= 1.0);
        r = sqrt((-2.0 * log(r)) / r);
        Y_N[i] = (X_N[i] << 1) - 1 + v1 * r * sigma;
        i++;
        Y_N[i] = (X_N[i] << 1) - 1 + v2 * r * sigma;
        i++;
    }
    // dec_input <-- 2 * c_out / sg^2.
    for (int i = 0; i < N; i++) Y_N[i] = -Y_N[i] * 2 / (design_sigma * design_sigma);
    return 0;
}

Channel_BSC::Channel_BSC(const int N, double p, const int seed)
: Channel_c(N, seed), p(p)
{
    if (p>1 || p<0)
        cerr << "Channel_BSC 0<=p<=1 violated" << endl; 
}

int Channel_BSC::add_noise(const int *X_N, int *Y_N, const size_t frame_id)
{
    // X_N[i] = ((mt19937.randf_cc() >= this->p) != X_N[i]) ? 0 : 1;
    int num_flips = 0;
    for(int i = 0; i < N; i++) {
        if (mt19937.randd_cc() <= this->p) {
            num_flips++;
            Y_N[i] = !X_N[i]; 
        } else Y_N[i] = X_N[i];
    }
    // cerr << "num_flips " << num_flips << endl;
    return num_flips;
}