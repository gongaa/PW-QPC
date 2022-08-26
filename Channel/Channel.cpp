#include <utility>
#include <algorithm>

#include "Channel.hpp"
#include "Channel_BSC.hpp"
#include "Channel_depolarize.hpp"

using namespace std;

const string Channel_name = "Channel";
const string Channel_prefix = "chn";

Channel::Channel(const string &prefix)
{

}

Channel* Channel::build_event() const
{
    if (type == "BEC") return new Channel_BSC(this->N, this->seed);
}