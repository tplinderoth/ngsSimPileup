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
//	double generate(); // spawn a random real number in range (lb,ub)
	double operator () ();

	// PUBLIC MEMBERS
	std::function<double()> real_rand;

	// PRIVATE MEMBERS
	double _lb;
	double _ub;
};



#endif /* RANDOMGENERATOR_H_ */
