/*
 * RandomGenerator.cpp
 *
 * TODO:
 * lower and upper boundary checks
 * add exception handling
 */

#include <stdexcept>
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

NormalGenerator::NormalGenerator(double mean, double stddev, double min, double max)
	: seed(std::chrono::system_clock::now().time_since_epoch().count()),
	  generator(seed),
	  _min(min),
	  _max(max),
	  distribution(mean, stddev)
{
	if (stddev < 0)
		throw std::logic_error("Attempt to initialize NormalGenerator with negative standard deviation");
	if (min > max)
		throw std::logic_error("Attempt to initialize NormalGenerator with min > max");
}

double NormalGenerator::operator() () {
	double x = 0.0;
	while (true) {
		x = this->distribution(generator);
		if (x >= this->_min && x <= this->_max)
			return x;
	}
}
