#include "Channel_depolarize_q.hpp"
#include <iostream>
using namespace std;

Channel_depolarize_q::Channel_depolarize_q(const int N, float p, const int seed = 0)
:N(N), p(p), seed(seed), type("Depolarize_q")
{
    if (3*p > 1)
        cerr << "Depolarizing channel 3*p>1" << endl;
    else
        this->event_generator = new Event_generator_unitary(seed, 2*p, p, 2*p);
}