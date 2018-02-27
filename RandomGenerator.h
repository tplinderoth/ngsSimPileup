/*
 * RandomGenerator.h
 */

#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_

#include <chrono>
#include <random>
#include <functional>

class Random {
public:
	// PUBLIC FUNCTIONS
	Random(double lb = 0.0, double ub = 1.0); // initializer
	double operator () ();

	// PUBLIC MEMBERS
	std::function<double()> real_rand;

	// PRIVATE MEMBERS
	double _lb;
	double _ub;
};

class NormalGenerator {
	std::mt19937::result_type seed;
	std::mt19937 generator;
	double _min;
	double _max;
	std::normal_distribution<double> distribution;
public:
	NormalGenerator(double mean, double stddev, double min, double max); // initializer
	double operator()();
};

#endif /* RANDOMGENERATOR_H_ */
