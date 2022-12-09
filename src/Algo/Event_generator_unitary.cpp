#include "Algo/Event_generator_unitary.hpp"
#include <iostream>
using namespace std;
enum ERROR_TYPE { NO_ERR = 0, ERR_X = 1, ERR_Y = 2, ERR_Z = 3};

Event_generator_unitary::Event_generator_unitary(const int seed, double px, double py, double pz)
:px(px), py(py), pz(pz), tx(px - py), ty(px), tz(px + pz - py)
{
	if(py > px || py > pz)
		cerr << "Event generator unitary: py>px or py>pz" << endl;
	this->set_seed(seed);
}

Event_generator_unitary* Event_generator_unitary::clone() const
{
	return new Event_generator_unitary(*this);
}


void Event_generator_unitary::set_seed(const int seed)
{
	mt19937.seed(seed);
	cerr << "Event generator random seed: " << seed << endl;;
}


void Event_generator_unitary::generate(int *draw_x, int *draw_z, const unsigned length)
{
	// 0 < p < px - py:        X
	// px - py < p < px:       XZ
	// px < p < px + pz - py:  Z
	// px + pz - py < p < 1:   I
	for (unsigned i = 0; i < length; i++) {
		double temp = mt19937.randd_cc();
		if (temp > tz)
			draw_x[i] = draw_z[i] = 0;
		else if (temp > ty) {
			draw_x[i] = 0;
			draw_z[i] = 1;
		}
		else if (temp > tx) 
			draw_x[i] = draw_z[i] = 1;
		else {
			draw_x[i] = 1;
			draw_z[i] = 0;
		}
	}
}

void Event_generator_unitary::set_prob(double px, double py, double pz)
{
	this->px = px; this->py = py; this->pz = pz;
	this->tx = px-py; this->ty = px; this->tz = px+pz-py;
}
