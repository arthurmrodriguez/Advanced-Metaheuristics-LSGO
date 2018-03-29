#ifndef GREFENSTETTELOGSTAT_H_
#define GREFENSTETTELOGSTAT_H_

#include "SingleLogStat.h"

class GrefenstetteLogStat : public SingleLogStat {
public:
	GrefenstetteLogStat (const Algorithm& alg) : SingleLogStat("gref_bias",alg){}
	~GrefenstetteLogStat(                    ) {}

	long double computeValue() { return alg_.population().grefenstetteBias(); } 
};

#endif /*GREFENSTETTELOGSTAT_H_*/
