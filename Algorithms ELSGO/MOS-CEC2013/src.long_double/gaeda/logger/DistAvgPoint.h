#ifndef DISTAVGPOINT_H_
#define DISTAVGPOINT_H_

#include "SingleLogStat.h"

/*
 * Measure explained in article http://portal.acm.org/citation.cfm?id=669270
*/
class DistAvgPoint : public SingleLogStat {
public:
	DistAvgPoint(const Algorithm& alg) : SingleLogStat("distAvgPoint",alg){}
	virtual ~DistAvgPoint() {}
	
	long double computeValue() { return alg_.population().distAvgPoint(); }
};

#endif /*DISTAVGPOINT_H_*/
