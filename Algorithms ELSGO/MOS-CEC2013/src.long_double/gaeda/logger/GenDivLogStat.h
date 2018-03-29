#ifndef GENDIVSTATVALUE_H_
#define GENDIVSTATVALUE_H_

#include "CollectionLogStat.h"

class GenDivLogStat : public CollectionLogStat {
public:
	GenDivLogStat(const Algorithm& alg);
	~GenDivLogStat();
	
	void   computeValues(long double& max, long double& min, long double& avg, long double& dev);
	long double indValue(GAGenome& ind);
};

#endif /*GENDIVSTATVALUE_H_*/
