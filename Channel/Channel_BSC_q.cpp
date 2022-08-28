#include "Channel_BSC_q.hpp"

Channel_BSC_q::Channel_BSC_q(const int N, float px, float pz, const int seed = 0)
:N(N), px(px), pz(pz), seed(seed), type("BSC_q"), 
event_generator(new Event_generator_unitary(seed, px*(1-pz), (1-px)*(1-pz), (1-px)*pz))
{

}