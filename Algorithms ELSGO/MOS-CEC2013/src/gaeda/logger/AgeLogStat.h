#ifndef AGESTATVALUE_H_
#define AGESTATVALUE_H_

#include "CollectionLogStat.h"

class AgeLogStat : public CollectionLogStat {
public:
	AgeLogStat(const Algorithm& alg) : CollectionLogStat("age",alg) {}
	~AgeLogStat(){}
	
	double indValue(GAGenome& ind){return ind.age();}
};

#endif /*AGESTATVALUE_H_*/
