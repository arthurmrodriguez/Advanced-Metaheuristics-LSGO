#ifndef NATIVEPRCNTSTATVALUE_H_
#define NATIVEPRCNTSTATVALUE_H_

#include "SingleLogStat.h"

class NativePrcntLogStat : public SingleLogStat {
  int rank_;
public:
	NativePrcntLogStat(const Algorithm& alg,int rank) : SingleLogStat("native_prcnt",alg), rank_(rank) {}
	
	~NativePrcntLogStat() {}
	
	long double computeValue() { return alg_.population().nativePrcnt(rank_); }
};

#endif /*NATIVEPRCNTSTATVALUE_H_*/
