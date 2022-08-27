#include "Channel_depolarize_q.hpp"

Channel_depolarize_q::Channel_depolarize_q(const int N, float p, const int seed = 0)
:N(N), p(p), seed(seed), type("Depolarize_q")
{

}