#ifndef FITNESSSTATVALUE_H_
#define FITNESSSTATVALUE_H_

#include "CollectionLogStat.h"

class FitnessLogStat : public CollectionLogStat
{
public:
	FitnessLogStat(const Algorithm& alg) : CollectionLogStat("fitness",alg) {}
	~FitnessLogStat(){}

	long double indValue(GAGenome& ind){return ind.fitness();}
};

#endif /*FITNESSSTATVALUE_H_*/
