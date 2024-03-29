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
    double p;

public:
    Channel_BSC(const int N, double p, const int seed = 0);
    virtual ~Channel_BSC() = default;
    // virtual Channel_BSC* clone() const;

    // TODO: support a probability vector, this should go hand in hand with SIMD
    // void _add_noise(const double *CP, const int *X_N, const int *Y_N, const size_t frame_id);
    int add_noise(const int *X_N, int *Y_N, const size_t frame_id);

};

class Channel_AWGN : Channel_c
{
protected:
    double sigma;
    double design_sigma;

public: 
    Channel_AWGN(const int N, double sigma, double design_sigma, const int seed = 0);
    virtual ~Channel_AWGN() = default;
    int add_noise(const int *X_N, double* Y_N, const size_t frame_id);
};

class Channel_q // classical, has codeword
{
public:
    // shared_ptr<Event_generator_unitary> event_generator;
    Event_generator_unitary* event_generator;
    int N;

    explicit Channel_q(const int N, Event_generator_unitary* event_generator);
    virtual ~Channel_q() = default;
    // Channel_q* clone() const;
    void set_seed(const int seed);

    // TODO: support a probability vector, this should go hand in hand with SIMD
    // always add to the all zero codeword
    void add_noise(int *Y1_N, int *Y2_N, const size_t frame_id);
};

class Channel_depolarize_q : public Channel_q
{
protected:
    double p;
    string type = "Depolarize_q";

public:
    Channel_depolarize_q(const int N, double p, const int seed) 
    : Channel_q(N, new Event_generator_unitary(seed, 2*p/3, p/3 , 2*p/3)), p(p) {
        if (p > 1)
            cerr << "Depolarizing channel p>1" << endl;
    }
    virtual void set_prob(double p) { this->event_generator->set_prob(2*p/3, p/3, 2*p/3); }
};

class Channel_BSC_q : public Channel_q
{
protected:
    double px, pz;
    string type = "BSC_q";

public:
    Channel_BSC_q(const int N, double px, double pz, const int seed) 
    : Channel_q(N, new Event_generator_unitary(seed, px, px*pz, pz)), px(px), pz(pz) {}
    virtual void set_prob(double px, double pz) { this->event_generator->set_prob(px, px*pz, pz); }
};

#endif /* CHANNEL_HPP */