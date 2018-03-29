#ifndef FITINCSTATVALUE_H_
#define FITINCSTATVALUE_H_

#include "SingleLogStat.h"
#include "../GAStatistics.h"

class FitIncLogStat : public SingleLogStat {
public:
	FitIncLogStat(const Algorithm& alg) : SingleLogStat("fit_inc",alg) {}
	
	~FitIncLogStat(){}

	long double computeValue(){ return alg_.statistics().fitInc();	}
};

#endif /*FITINCSTATVALUE_H_*/
