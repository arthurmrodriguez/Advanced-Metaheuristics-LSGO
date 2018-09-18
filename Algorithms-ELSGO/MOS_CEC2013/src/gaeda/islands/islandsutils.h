#ifndef UTIL_H_
#define UTIL_H_

#include "../GAPopulation.h"
#include "../genomes/GAGenome.h"
#include <vector>

using namespace std;

vector<GAGenome*>* convertPopToVector(GAPopulation* pop);



#endif /*UTIL_H_*/
