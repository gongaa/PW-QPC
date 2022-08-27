#ifndef CHANNEL_DEPOLARIZE_HPP
#define CHANNEL_DEPOLARIZE_HPP

#include "Channel.hpp"
#include "Algo/Event_generator_unitary.hpp"

class Channel_depolarize_q : Channel_q
{
protected:
    std::shared_ptr<Event_generator_unitary> event_generator;
    float p;

public:
    Channel_depolarize_q(const int N, float p, const int seed = 0);
    virtual ~Channel_depolarize_q() = default;
    virtual Channel_depolarize_q* clone() const;

protected:
    // no implicit codeword, only errors and syndrome (handled by decoder)
    // CSS code, X and Z decoding are always separate
    void _add_noise(const float *CP, const int *Y1_N, const int *Y2_N, const size_t frame_id);
};

#endif 