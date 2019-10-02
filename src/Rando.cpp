#include "Rando.h"
#include <chrono>
#include <time.h>

std::normal_distribution<double>       Rando::normal_distribution(0.0, 1.0);
std::uniform_real_distribution<double> Rando::uniform_distribution( 0.0, 1.0 );
std::mt19937 Rando::mt((unsigned int) std::chrono::high_resolution_clock::now().time_since_epoch().count() );
