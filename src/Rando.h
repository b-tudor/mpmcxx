#pragma once
#include <random>
class Rando
{
public:
	static void seed( unsigned int seed ) { mt.seed(seed); }

	static double rand()        { return uniform_distribution(mt); }
	static double rand_normal() { return  normal_distribution(mt); }
	
private:
	static std::mt19937 mt;
	static std::normal_distribution<double>       normal_distribution;
	static std::uniform_real_distribution<double> uniform_distribution;
};
