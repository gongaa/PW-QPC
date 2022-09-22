#ifndef CHANNEL_HPP
#define CHANNEL_HPP

#include <string>
#include <map>
#include <iostream>

using namespace std;
#include "Algo/Event_generator_unitary.hpp"
#include "Algo/PRNG_MT19937.hpp"

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

class Channel_BSC : Channel_c
{
protected:
    float p;

public:
    Channel_BSC(const int N, float p, const int seed = 0);
    virtual ~Channel_BSC() = default;
    // virtual Channel_BSC* clone() const;

    // TODO: support a probability vector, this should go hand in hand with SIMD
    // void _add_noise(const float *CP, const int *X_N, const int *Y_N, const size_t frame_id);
    int add_noise(const int *X_N, int *Y_N, const size_t frame_id);

};

class Channel_AWGN : Channel_c
{
protected:
    double sigma;

public: 
    Channel_AWGN(const int N, double sigma, const int seed = 0);
    virtual ~Channel_AWGN() = default;
    int add_noise(const int *X_N, double* Y_N, const size_t frame_id);
};

class Channel_q // classical, has codeword
{
public:
    // shared_ptr<Event_generator_unitary> event_generator;
    Event_generator_unitary* event_generator;
    int N;
    int seed;

    explicit Channel_q(const int N, Event_generator_unitary* event_generator, const int seed = 0);
    virtual ~Channel_q() = default;
    // Channel_q* clone() const;
    virtual void set_seed(const int seed);

    // TODO: support a probability vector, this should go hand in hand with SIMD
    // virtual void _add_noise(const float *CP, int *Y1_N, int *Y2_N, const size_t frame_id);
    virtual void add_noise(int *Y1_N, int *Y2_N, const size_t frame_id);
};

class Channel_depolarize_q : Channel_q
{
protected:
    float p;
    string type = "Depolarize_q";

public:
    Channel_depolarize_q(const int N, float p, const int seed = 0) 
    : Channel_q(N, new Event_generator_unitary(seed, 2*p, p , 2*p), seed), p(p) {
        if (3*p > 1)
            cerr << "Depolarizing channel 3*p>1" << endl;
    }
};

class Channel_BSC_q : Channel_q
{
protected:
    float px, pz;
    string type = "BSC_q";

public:
    Channel_BSC_q(const int N, float px, float pz, const int seed = 0) 
    : Channel_q(N, new Event_generator_unitary(seed, px*(1-pz), (1-px)*(1-pz), (1-px)*pz)), px(px), pz(pz) {}
};

#endif /* CHANNEL_HPP */