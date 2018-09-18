#ifndef GENDIVSTATVALUE_H_
#define GENDIVSTATVALUE_H_

#include "CollectionLogStat.h"

class GenDivLogStat : public CollectionLogStat {
public:
	GenDivLogStat(const Algorithm& alg);
	~GenDivLogStat();
	
	void   computeValues(double& max, double& min, double& avg, double& dev);
	double indValue(GAGenome& ind);
};

#endif /*GENDIVSTATVALUE_H_*/
