#ifndef CHANNEL_HPP
#define CHANNEL_HPP

#include <string>
#include <map>

using namespace std;
#include "Algo/Event_generator_unitary.hpp"
#include "PRNG_MT19937.hpp"

class Channel_c // classical, has codeword
{
public:
    string type = "BSC_c";
    int N = 0;
	PRNG_MT19937 mt19937;

    explicit Channel_c(int N, const int seed = 0);
    virtual ~Channel_c() = default;
    // TODO: copy constructor
    // Channel_c* clone() const; 
};

class Channel_q // classical, has codeword
{
public:
    shared_ptr<Event_generator_unitary> event_generator;
    string type = "BSC";
    int N = 0;
    int seed = 0;

    explicit Channel_q(const int N, const int seed = 0);
    virtual ~Channel_q() = default;
    // Channel_q* clone() const;
    virtual void set_seed(const int seed);

protected:
    virtual void _add_noise(const float *CP, int *Y1_N, int *Y2_N, const size_t frame_id);
};

class Channel_BSC : Channel_c
{
protected:
    float p;

public:
    Channel_BSC(const int N, float p, const int seed = 0);
    virtual ~Channel_BSC() = default;
    // virtual Channel_BSC* clone() const;

protected:
    void _add_noise(const float *CP, const int *X_N, const int *Y_N, const size_t frame_id);

};

#endif /* CHANNEL_HPP */