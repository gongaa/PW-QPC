#ifndef CHANNEL_BSC_HPP
#define CHANNEL_BSC_HPP

#include <memory>

#include "Channel.hpp"

class Channel_BSC_q : Channel_q
{
protected:
    float px, pz;

public:
    Channel_BSC_q(const int N, float px, float pz, const int seed = 0);
    virtual ~Channel_BSC_q() = default;
    virtual Channel_BSC_q* clone() const;

protected:
    // no implicit codeword, only errors and syndrome (handled by decoder)
    void _add_noise(const float *CP, const int *Y1_N, const int *Y2_N, const size_t frame_id);

};

#endif