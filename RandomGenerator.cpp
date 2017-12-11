/*
 * RandomGenerator.cpp
 *
 * TODO: lower and upper boundary checks
 *
 */

#include "RandomGenerator.h"

Random::Random (double lb, double ub)
	: _lb(0.0),
	  _ub(1.0)
{
	_lb = lb;
	_ub = ub;
	std::mt19937::result_type seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	real_rand = std::bind(std::uniform_real_distribution<double>(_lb, _ub), std::mt19937(seed));
}

double Random::operator () ()
{
	double rand = 0.0;
	rand = real_rand();
	return rand;
}
