#ifndef EVENT_GENERATOR_UNITARY_HPP
#define EVENT_GENERATOR_UNITARY_HPP

// #include <random>
#include "PRNG_MT19937.hpp"

class Event_generator_unitary
{
protected:
	// std::mt19937 rd_engine; // Mersenne Twister 19937
	PRNG_MT19937 mt19937;
	float px, py, pz;
	float tx, ty, tz;

public:
	explicit Event_generator_unitary(const int seed = 0, float px = 0, float py = 0, float pz = 0);
	virtual ~Event_generator_unitary() = default;
	virtual Event_generator_unitary* clone() const;

	virtual void set_seed(const int seed);
	virtual void generate(int *draw_x, int *draw_z, const unsigned length);
};

#endif //EVENT_GENERATOR_UNITARY_HPP
